"""foamgci.report — end-to-end V&V report assembly.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence


from ._version import __version__
from .reader import read_fieldminmax
from .stats import window_stats, WindowStats
from .gci import gci_over_hierarchy, GCIResult


def _fmt_kpss_p(p: float) -> str:
    """KPSS p-values are interpolated on the published critical values and
    clamped to [0.01, 0.10]; render the clamped ends honestly."""
    if p >= 0.10:
        return ">=0.100"
    if p <= 0.01:
        return "<=0.010"
    return f"{p:.3f}"


# ---------------------------------------------------------------------------
# Analytical reference: Rayleigh–Pitot stagnation pressure ratio
# ---------------------------------------------------------------------------

def rayleigh_pitot(M1: float, gamma: float = 1.4) -> float:
    if M1 <= 1.0:
        raise ValueError("Rayleigh-Pitot requires supersonic upstream (M1 > 1).")
    g = float(gamma)
    M2 = float(M1) ** 2
    term1 = ((g + 1.0) ** 2 * M2 / (4.0 * g * M2 - 2.0 * (g - 1.0))) ** (g / (g - 1.0))
    term2 = (1.0 - g + 2.0 * g * M2) / (g + 1.0)
    return float(term1 * term2)


# ---------------------------------------------------------------------------
# Container
# ---------------------------------------------------------------------------

@dataclass
class GridCase:
    """One OpenFOAM case in a refinement hierarchy."""

    label: str
    path: Path
    h: float
    n_cells: int | None = None


@dataclass
class ReportTable:
    cases: list[GridCase]
    field: str
    quantity: str
    window: tuple[float, float]
    stats: list[WindowStats] = field(default_factory=list)
    gcis: list[GCIResult] = field(default_factory=list)
    reference_value: float | None = None
    reference_label: str = ""

    # ------------------------------------------------------------------
    # Pretty printers
    # ------------------------------------------------------------------

    def as_text(self) -> str:
        """Human-readable plain-text report."""
        lines: list[str] = []
        lines.append("=" * 72)
        lines.append(f"foamgci V&V report — field {self.field!r}, "
                     f"quantity {self.quantity!r}, "
                     f"window [{self.window[0]:g}, {self.window[1]:g}]")
        lines.append("=" * 72)
        # Per-grid table
        lines.append("")
        lines.append("Per-grid time-averaged statistics:")
        lines.append("")
        hdr = (f"  {'label':<12}"
               f"{'N_cells':>10}"
               f"{'h':>10}"
               f"{'N':>7}"
               f"{'mean':>12}"
               f"{'std':>10}"
               f"{'τ_int':>8}"
               f"{'SEM':>10}"
               f"{'N_eff':>8}"
               f"{'KPSS_p':>9}")
        lines.append(hdr)
        lines.append("  " + "-" * (len(hdr) - 2))
        for c, s in zip(self.cases, self.stats):
            lines.append(
                f"  {c.label:<12}"
                f"{(c.n_cells if c.n_cells else 0):>10d}"
                f"{c.h:>10.4g}"
                f"{s.n:>7d}"
                f"{s.mean:>12.6g}"
                f"{s.std:>10.4g}"
                f"{s.tau_int:>8.2f}"
                f"{s.sem:>10.4g}"
                f"{s.n_eff:>8.1f}"
                f"{_fmt_kpss_p(s.kpss_p):>9}"
            )

        # GCI block
        lines.append("")
        lines.append("Roache GCI on consecutive triplets (coarse → medium → fine):")
        lines.append("")
        for g in self.gcis:
            lines.append(
                f"  triplet ({g.label_coarse}, {g.label_medium}, {g.label_fine}):"
            )
            lines.append(f"      regime                   = {g.regime}")
            lines.append(f"      apparent order p̂        = {g.p_apparent:.3f}")
            lines.append(f"      Richardson φ_exact       = {g.phi_exact:.6g}")
            lines.append(f"      GCI_fine_21              = {g.gci_fine_21_pct:.4f} %")
            lines.append(f"      GCI_medium_32            = {g.gci_medium_32_pct:.4f} %")
            lines.append(f"      asymptotic ratio (≈1)    = {g.asymptotic_ratio:.3f}")
            if g.regime != "monotonic":
                lines.append(f"      note                     = {g.note}")
            if g.regime == "oscillatory" and g.u_oscillatory_pct == g.u_oscillatory_pct:
                lines.append(
                    f"      U_oscillatory (Celik)    = "
                    f"{g.u_oscillatory_pct:.4f} %  (half solution span)")

        # Reference cross-check
        if self.reference_value is not None and self.gcis:
            phi_ext = self.gcis[-1].phi_exact
            rel = 100.0 * abs(phi_ext - self.reference_value) / abs(self.reference_value)
            lines.append("")
            lines.append(f"Analytical reference ({self.reference_label}):")
            lines.append(f"      reference value         = {self.reference_value:.6g}")
            lines.append(f"      Richardson extrapolation = {phi_ext:.6g}  "
                         f"(error = {rel:.4f} %)")
            lines.append(f"      finest-grid mean         = "
                         f"{self.stats[-1].mean:.6g}  "
                         f"(error = "
                         f"{100*abs(self.stats[-1].mean - self.reference_value)/abs(self.reference_value):.4f} %)")
            g_last = self.gcis[-1]
            band = g_last.gci_fine_21_pct
            if band == band:  # not NaN
                covered = rel <= band
                lines.append(
                    f"      |phi_ext - ref| vs GCI_21 = {rel:.4f} % vs "
                    f"{band:.4f} %  -> "
                    + ("reference COVERED by GCI band"
                       if covered else
                       "reference OUTSIDE GCI band: GCI likely understates "
                       "total uncertainty, or the reference is not directly "
                       "comparable to this QoI"))
        lines.append("=" * 72)
        return "\n".join(lines)

    def as_latex(self) -> str:
        """LaTeX ``tabular`` block suitable for direct insertion into a paper."""
        rows = []
        for c, s in zip(self.cases, self.stats):
            ncells = f"{c.n_cells:,}" if c.n_cells else "—"
            kp_tex = _fmt_kpss_p(s.kpss_p)
            kp_tex = kp_tex.replace(">=", "$\\geq$").replace("<=", "$\\leq$")
            rows.append(
                f"{c.label} & {ncells} & {c.h:.4g} & "
                f"{s.n:d} & {s.mean:.4f} & {s.std:.4f} & "
                f"{s.tau_int:.2f} & {s.sem:.4f} & "
                f"{kp_tex}"
                r" \\"
            )
        body = "\n        ".join(rows)
        # GCI table
        gci_rows = []
        for g in self.gcis:
            gci_rows.append(
                f"({g.label_coarse},{g.label_medium},{g.label_fine}) & "
                f"{g.p_apparent:.3f} & {g.phi_exact:.4f} & "
                f"{g.gci_fine_21_pct:.4f} & {g.gci_medium_32_pct:.4f} & "
                f"{g.asymptotic_ratio:.3f}"
                r" \\"
            )
        gci_body = "\n        ".join(gci_rows)

        ref_line = ""
        if self.reference_value is not None and self.gcis:
            phi_ext = self.gcis[-1].phi_exact
            rel = 100.0 * abs(phi_ext - self.reference_value) / abs(self.reference_value)
            ref_line = (
                f"\n% Analytical reference ({self.reference_label}): "
                f"phi_ref = {self.reference_value:.4f}; "
                f"|phi_exact - phi_ref|/phi_ref = {rel:.4f} %"
            )

        latex = rf"""% foamgci v{__version__} — auto-generated table
% field = {self.field!r}, quantity = {self.quantity!r},
% window = [{self.window[0]:g}, {self.window[1]:g}]{ref_line}
\begin{{table}}[t]
\centering
\caption{{Time-averaged extrema and Roache GCI on a four-grid hierarchy.}}
\label{{tab:gci}}
\begin{{tabular}}{{lrrrrrrrr}}
    \toprule
    Grid & $N_{{\text{{cells}}}}$ & $h$ & $N_t$ & $\langle\phi\rangle$ &
    $\sigma$ & $\hat\tau_{{\mathrm{{int}}}}$ & SEM & KPSS-$p$ \\
    \midrule
        {body}
    \bottomrule
\end{{tabular}}
\vspace{{0.5em}}
\begin{{tabular}}{{lrrrrr}}
    \toprule
    Triplet & $\hat p$ & $\phi_{{\mathrm{{exact}}}}$ &
    GCI$_{{21}}^{{\mathrm{{fine}}}}$ (\%) &
    GCI$_{{32}}^{{\mathrm{{medium}}}}$ (\%) &
    $R_{{\mathrm{{asym}}}}$ \\
    \midrule
        {gci_body}
    \bottomrule
\end{{tabular}}
\end{{table}}
"""
        return latex


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def full_report(
    cases: Sequence[GridCase],
    field: str = "p",
    quantity: str = "max",
    window: tuple[float, float] = (3.0, 10.0),
    reference_value: float | None = None,
    reference_label: str = "",
    kpss_regression: str = "c",
) -> ReportTable:
    if quantity not in ("max", "min"):
        raise ValueError("quantity must be 'max' or 'min'.")
    if len(cases) < 2:
        raise ValueError("Need at least 2 grids for a refinement study; "
                         "3+ for GCI.")
    if any(cases[i].h <= cases[i + 1].h for i in range(len(cases) - 1)):
        raise ValueError(
            "Cases must be coarse-to-fine (h strictly decreasing). "
            "Got h's: " + ", ".join(f"{c.h:g}" for c in cases)
        )

    stats_list: list[WindowStats] = []
    means: list[float] = []
    hs: list[float] = []
    labels: list[str] = []
    for c in cases:
        data = read_fieldminmax(c.path, field=field)
        series = data.max if quantity == "max" else data.min
        ws = window_stats(data.time, series, window[0], window[1],
                          kpss_regression=kpss_regression)
        stats_list.append(ws)
        means.append(ws.mean)
        hs.append(c.h)
        labels.append(c.label)

    gcis: list[GCIResult] = []
    if len(cases) >= 3:
        # gci_over_hierarchy expects coarse-to-fine (h strictly decreasing) — same order as `cases`.
        gcis = gci_over_hierarchy(means, hs, labels)

    return ReportTable(
        cases=list(cases),
        field=field,
        quantity=quantity,
        window=tuple(window),
        stats=stats_list,
        gcis=gcis,
        reference_value=reference_value,
        reference_label=reference_label,
    )
