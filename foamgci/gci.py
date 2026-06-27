from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np


# ---------------------------------------------------------------------------
# Apparent order
# ---------------------------------------------------------------------------

def apparent_order(
    phi1: float, phi2: float, phi3: float,
    r21: float, r32: float,
    p_min: float = 0.1, p_max: float = 6.0, tol: float = 1e-8,
    max_iter: int = 100,
) -> float:
    
    if r21 <= 1.0 or r32 <= 1.0:
        raise ValueError("Refinement ratios must be > 1.")
    eps32 = phi3 - phi2
    eps21 = phi2 - phi1
    if eps21 == 0.0:
        raise ValueError("phi2 - phi1 = 0; refinement is exact (or stalled).")
    if eps32 == 0.0:
        raise ValueError("phi3 - phi2 = 0; coarsening is exact (or stalled).")

    s = float(np.sign(eps32 / eps21))
    p = 2.0  # benign starting guess
    for _ in range(max_iter):
        try:
            num = r21 ** p - s
            den = r32 ** p - s
            if num == 0.0 or den == 0.0:
                # Avoid log(0); fall back to constant-r limit.
                q = 0.0
            else:
                q = np.log(num / den)
        except (OverflowError, FloatingPointError):
            q = 0.0
        p_new = (1.0 / np.log(r21)) * np.abs(np.log(abs(eps32 / eps21)) + q)
        p_new = float(np.clip(p_new, p_min, p_max))
        if abs(p_new - p) < tol:
            return p_new
        p = p_new
    return p


# ---------------------------------------------------------------------------
# Richardson + GCI
# ---------------------------------------------------------------------------

def richardson_extrapolation(
    phi1: float, phi2: float, r21: float, p: float
) -> float:
    """Exact-grid value via Richardson extrapolation."""
    denom = r21 ** p - 1.0
    if denom == 0.0:
        return float(phi1)
    return float(phi1 + (phi1 - phi2) / denom)


@dataclass
class GCIResult:
    """Result of Roache GCI on a triplet (coarse → medium → fine)."""

    label_coarse: str
    label_medium: str
    label_fine: str
    h_coarse: float
    h_medium: float
    h_fine: float
    phi_coarse: float
    phi_medium: float
    phi_fine: float
    r21: float
    r32: float
    p_apparent: float
    phi_exact: float
    gci_fine_21_pct: float      # GCI_21 on the fine grid (%)
    gci_medium_32_pct: float    # GCI_32 on the medium grid (%)
    asymptotic_ratio: float     # ≈ 1 for asymptotic range (valid for r21 ≠ r32)
    Fs: float                   # safety factor (1.25 for ≥ 3 grids)
    e21_relative: float         # |φ1-φ2|/|φ1|
    e32_relative: float         # |φ2-φ3|/|φ2|
    regime: str = "monotonic"   # monotonic | oscillatory | divergent | degenerate
    note: str = ""
    u_oscillatory_pct: float = float("nan")  # Celik (2008) oscillatory-convergence
                                # uncertainty 0.5*(phi_max - phi_min)/|phi_fine|, in %

    def __repr__(self) -> str: 
        return (
            f"GCIResult(triplet={self.label_coarse}/"
            f"{self.label_medium}/{self.label_fine}, "
            f"p̂={self.p_apparent:.3f}, "
            f"φ_exact={self.phi_exact:.6g}, "
            f"GCI_fine={self.gci_fine_21_pct:.3f}%, "
            f"R_asym={self.asymptotic_ratio:.3f})"
        )


def roache_gci(
    phi_coarse: float, phi_medium: float, phi_fine: float,
    h_coarse: float, h_medium: float, h_fine: float,
    label_coarse: str = "coarse",
    label_medium: str = "medium",
    label_fine: str = "fine",
    Fs: float = 1.25,
) -> GCIResult:
    """Roache GCI on a single refinement triplet.

    Parameters
    ----------
    phi_coarse, phi_medium, phi_fine : float
        Quantity of interest on each grid.
    h_coarse, h_medium, h_fine : float
        Representative cell size on each grid (units irrelevant — only
        ratios enter).
    label_* : str
        Used only for human-readable output.
    Fs : float
        Safety factor. Celik et al. (2008) recommend ``1.25`` when three
        or more grids are available, ``3.0`` otherwise.

    Returns
    -------
    GCIResult
    """
    if not (h_coarse > h_medium > h_fine > 0):
        raise ValueError(
            "Cell sizes must satisfy h_coarse > h_medium > h_fine > 0, "
            f"got h_c={h_coarse}, h_m={h_medium}, h_f={h_fine}."
        )
    r21 = h_medium / h_fine
    r32 = h_coarse / h_medium

    d21 = phi_fine - phi_medium
    d32 = phi_medium - phi_coarse
    eps21 = abs(d21)
    eps32 = abs(d32)
    e21_rel = eps21 / max(abs(phi_fine), 1e-300)
    e32_rel = eps32 / max(abs(phi_medium), 1e-300)

    if eps21 < 1e-14 or eps32 < 1e-14:
        regime, note = "degenerate", "vanishing residual; observed order undefined"
    else:
        R_conv = d21 / d32
        if R_conv < 0.0:
            regime = "oscillatory"
            note = f"R={R_conv:.3f} < 0: sign change between successive differences"
        elif R_conv > 1.0:
            regime = "divergent"
            note = f"R={R_conv:.3f} > 1: residual grows under refinement (pre-asymptotic)"
        else:
            regime, note = "monotonic", "monotonic convergence; GCI valid"

    if regime != "monotonic":
        nan = float("nan")
        # Celik et al. (2008): for oscillatory convergence, estimate the
        # uncertainty as half the span of the solutions on the triplet.
        u_osc = nan
        if regime == "oscillatory":
            span = max(phi_coarse, phi_medium, phi_fine) - min(
                phi_coarse, phi_medium, phi_fine)
            u_osc = 100.0 * 0.5 * span / max(abs(phi_fine), 1e-300)
        return GCIResult(
            label_coarse=label_coarse, label_medium=label_medium,
            label_fine=label_fine,
            h_coarse=h_coarse, h_medium=h_medium, h_fine=h_fine,
            phi_coarse=phi_coarse, phi_medium=phi_medium, phi_fine=phi_fine,
            r21=r21, r32=r32,
            p_apparent=nan, phi_exact=nan,
            gci_fine_21_pct=nan, gci_medium_32_pct=nan, asymptotic_ratio=nan,
            Fs=Fs, e21_relative=e21_rel, e32_relative=e32_rel,
            regime=regime, note=note, u_oscillatory_pct=u_osc,
        )

    p = apparent_order(phi_fine, phi_medium, phi_coarse, r21, r32)
    phi_exact = richardson_extrapolation(phi_fine, phi_medium, r21, p)
    gci21 = Fs * e21_rel / (r21 ** p - 1.0)
    gci32 = Fs * e32_rel / (r32 ** p - 1.0)
    # Asymptotic-range diagnostic, exact for non-constant refinement
    # ratios: in the asymptotic range eps32 = r21^p (r32^p-1)/(r21^p-1) eps21,
    # so R -> 1. For r21 == r32 this reduces to r^p * eps21/eps32.
    R = (r21 ** p) * eps21 * (r32 ** p - 1.0) / (
        max(eps32, 1e-300) * (r21 ** p - 1.0))

    return GCIResult(
        label_coarse=label_coarse,
        label_medium=label_medium,
        label_fine=label_fine,
        h_coarse=h_coarse, h_medium=h_medium, h_fine=h_fine,
        phi_coarse=phi_coarse, phi_medium=phi_medium, phi_fine=phi_fine,
        r21=r21, r32=r32,
        p_apparent=p,
        phi_exact=phi_exact,
        gci_fine_21_pct=100.0 * gci21,
        gci_medium_32_pct=100.0 * gci32,
        asymptotic_ratio=R,
        Fs=Fs,
        e21_relative=e21_rel,
        e32_relative=e32_rel,
        regime="monotonic",
        note=note,
    )


def gci_over_hierarchy(
    phis: Sequence[float],
    hs: Sequence[float],
    labels: Sequence[str] | None = None,
    Fs: float = 1.25,
) -> list[GCIResult]:
    
    if len(phis) != len(hs):
        raise ValueError("phis and hs must have the same length.")
    n = len(phis)
    if n < 3:
        raise ValueError("Need at least 3 grids.")
    if labels is None:
        labels = [f"grid{i}" for i in range(n)]
    if any(hs[i] <= hs[i + 1] for i in range(n - 1)):
        raise ValueError("hs must be strictly decreasing (coarse → fine).")

    out: list[GCIResult] = []
    for i in range(n - 2):
        out.append(
            roache_gci(
                phi_coarse=phis[i], phi_medium=phis[i + 1], phi_fine=phis[i + 2],
                h_coarse=hs[i], h_medium=hs[i + 1], h_fine=hs[i + 2],
                label_coarse=labels[i], label_medium=labels[i + 1],
                label_fine=labels[i + 2], Fs=Fs,
            )
        )
    return out
