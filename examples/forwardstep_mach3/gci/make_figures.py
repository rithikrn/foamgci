"""
make_figures.py
---------------
Generate grid-convergence figures for the manuscript.

Reads gci_summary.json produced by analyze.py; if that file is
missing the script tells you to run analyze.py first.

Figures generated:
  - fig_grid_convergence.*       : primary pressure QoI, p_max
  - fig_rho_grid_convergence.*   : secondary density QoI, rho_max

Emits both vector PDF and 600-dpi PNG for every figure.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Iterable

sys.path.insert(0, str(Path(__file__).parent))

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


HERE = Path(__file__).parent
SUMMARY = HERE / "gci_summary.json"
OUTDIR = HERE / "figures"

EXTS: tuple[str, ...] = ("pdf", "png")
PNG_DPI = 600

mpl.rcParams.update({
    # --- typography ---
    "font.family": "serif",
    "font.serif": [
        "STIX Two Text",
        "STIXGeneral",
        "Times New Roman",
        "DejaVu Serif",
    ],
    "mathtext.fontset": "stix",
    "font.size": 9.0,
    "axes.labelsize": 9.0,
    "axes.titlesize": 9.0,
    "legend.fontsize": 8.0,
    "xtick.labelsize": 8.0,
    "ytick.labelsize": 8.0,

    # --- line art ---
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.width": 0.4,
    "ytick.minor.width": 0.4,
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "xtick.minor.size": 1.7,
    "ytick.minor.size": 1.7,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "legend.frameon": False,
    "legend.handletextpad": 0.4,
    "legend.handlelength": 1.8,

    # --- output ---
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "savefig.transparent": False,
    "savefig.dpi": PNG_DPI,
    "figure.dpi": 150,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def save_all(fig: plt.Figure, stem: Path, exts: Iterable[str] = EXTS) -> None:
    """Save `fig` to `<stem>.<ext>` for every extension in `exts`."""
    for ext in exts:
        out = stem.with_suffix(f".{ext}")
        fig.savefig(out, dpi=PNG_DPI if ext == "png" else None)
        print(f"  wrote {out}")


def load_summary() -> dict:
    if not SUMMARY.is_file():
        raise SystemExit(
            f"\n  ERROR: {SUMMARY} not found.\n"
            "  Run `python3 analyze.py` first to produce it.\n"
        )
    return json.loads(SUMMARY.read_text())


def _is_number(x: object) -> bool:
    if x is None:
        return False
    try:
        v = float(x)
    except (TypeError, ValueError):
        return False
    return math.isfinite(v)


def _qoi_cases(summary: dict, qoi_key: str) -> list[dict]:
    """Return cases for one QoI, sorted coarse -> fine."""
    try:
        cases = list(summary["qoi_results"][qoi_key]["cases"])
    except KeyError as exc:
        raise SystemExit(
            f"\n  ERROR: qoi_results['{qoi_key}'] not found in {SUMMARY}.\n"
            "  Run the updated analyze.py first.\n"
        ) from exc

    return sorted(cases, key=lambda c: -float(c["dx"]))


def _qoi_report(summary: dict, qoi_key: str) -> dict:
    try:
        return summary["qoi_results"][qoi_key]
    except KeyError as exc:
        raise SystemExit(
            f"\n  ERROR: qoi_results['{qoi_key}'] not found in {SUMMARY}.\n"
            "  Run the updated analyze.py first.\n"
        ) from exc


def _plot_convergence(
    summary: dict,
    qoi_key: str,
    stem: Path,
    y_column: str,
    sem_column: str,
    ylabel: str,
    err_ylabel: str,
    connect_from: int,
    legend_loc_left: str = "lower right",
    legend_loc_right: str = "lower right",
) -> None:
    """Generic two-panel convergence plot for one scalar QoI."""
    report = _qoi_report(summary, qoi_key)
    cases = _qoi_cases(summary, qoi_key)

    h_vals = np.array([float(c["dx"]) for c in cases])
    phi_vals = np.array([float(c[y_column]) for c in cases])
    sem_vals = np.array([float(c[sem_column]) for c in cases])
    labels = [str(c["label"]) for c in cases]
    markers = ["o", "s", "D", "^"]

    phi_ext = report.get("phi_star")
    if not _is_number(phi_ext):
        raise SystemExit(f"ERROR: no finite phi_star for {qoi_key}")

    phi_ext = float(phi_ext)
    err = np.abs(phi_vals - phi_ext)

    tB = report.get("triplet_B_MFXF", {})
    p_obs = tB.get("p_obs")
    has_p = tB.get("regime") == "monotonic" and _is_number(p_obs)
    p_obs = float(p_obs) if has_p else float("nan")

    fig, axes = plt.subplots(1, 2, figsize=(6.85, 2.95))

    # ------------------------------------------------------------------
    # Left: QoI vs h
    # ------------------------------------------------------------------
    ax = axes[0]

    ax.axhline(
        phi_ext,
        color="0.55",
        lw=0.7,
        linestyle=(0, (4, 2)),
        label=fr"Richardson $\phi_{{\mathrm{{ext}}}}={phi_ext:.4f}$",
    )

    if 0 <= connect_from < len(h_vals) - 1:
        ax.plot(
            h_vals[connect_from:],
            phi_vals[connect_from:],
            color="C0",
            lw=0.6,
            zorder=2,
        )

    for i in range(len(cases)):
        ax.errorbar(
            h_vals[i],
            phi_vals[i],
            yerr=sem_vals[i],
            marker=markers[i],
            markersize=4.5,
            markerfacecolor="white",
            markeredgewidth=0.8,
            color="k",
            linestyle="none",
            capsize=2,
            elinewidth=0.6,
            zorder=5,
        )
        ax.annotate(
            labels[i],
            xy=(h_vals[i], phi_vals[i]),
            xytext=(6, -2),
            textcoords="offset points",
            fontsize=7.5,
        )

    ax.set_xscale("log")
    ax.set_xlabel(r"grid spacing $h$")
    ax.set_ylabel(ylabel)

    pad = max(0.02 * abs(phi_vals.mean()), 3.0 * sem_vals.max())
    ymin = min(phi_vals.min(), phi_ext) - pad
    ymax = max(phi_vals.max(), phi_ext) + pad
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(h_vals.min() * 0.7, h_vals.max() * 1.7)

    ax.legend(loc=legend_loc_left, borderaxespad=0.4)
    ax.xaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, subs=(1.0,)))
    ax.xaxis.set_minor_locator(
        mpl.ticker.LogLocator(base=10.0, subs=(0.2, 0.4, 0.6, 0.8))
    )
    ax.text(0.025, 0.96, "(a)", transform=ax.transAxes, fontsize=9, va="top")

    # ------------------------------------------------------------------
    # Right: error vs h
    # ------------------------------------------------------------------
    ax = axes[1]

    positive_err = err[err > 0.0]
    if positive_err.size == 0:
        raise SystemExit(f"ERROR: all errors are zero for {qoi_key}")

    h_ref = np.array([h_vals.min() * 0.7, h_vals.max() * 1.4])
    anchor_h = h_vals[-1]
    anchor_err = err[-1] if err[-1] > 0 else positive_err.min()

    for p_slope, ls, lbl in (
        (1.0, (0, (1, 1.5)), r"$\mathcal{O}(h)$"),
        (2.0, (0, (4, 2)), r"$\mathcal{O}(h^{2})$"),
    ):
        ref = anchor_err * (h_ref / anchor_h) ** p_slope
        ax.plot(h_ref, ref, color="0.55", lw=0.7, linestyle=ls, label=lbl)

    if has_p:
        # Anchor the observed line at the medium grid so it follows the
        # deepest M--F--XF triplet rather than the coarse-grid outlier.
        h_obs = np.array([h_vals[1], h_vals[-1] * 0.9])
        anchor_h2 = h_vals[1]
        anchor_err2 = err[1] if err[1] > 0 else positive_err.min()
        obs = anchor_err2 * (h_obs / anchor_h2) ** p_obs
        ax.plot(
            h_obs,
            obs,
            color="C3",
            lw=0.9,
            label=fr"observed $p={p_obs:.2f}$",
        )

    for i in range(len(cases)):
        if err[i] <= 0:
            continue
        ax.plot(
            h_vals[i],
            err[i],
            marker=markers[i],
            color="k",
            markersize=4.5,
            markerfacecolor="white",
            markeredgewidth=0.8,
            linestyle="none",
            zorder=5,
        )
        ax.annotate(
            labels[i],
            xy=(h_vals[i], err[i]),
            xytext=(5, -3 if i % 2 == 0 else 6),
            textcoords="offset points",
            fontsize=7.5,
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"grid spacing $h$")
    ax.set_ylabel(err_ylabel)
    ax.set_xlim(h_ref[0], h_ref[1])

    ax.legend(loc=legend_loc_right, borderaxespad=0.4)
    ax.xaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, subs=(1.0,)))
    ax.xaxis.set_minor_locator(
        mpl.ticker.LogLocator(base=10.0, subs=(0.2, 0.4, 0.6, 0.8))
    )
    ax.text(0.025, 0.96, "(b)", transform=ax.transAxes, fontsize=9, va="top")

    fig.subplots_adjust(wspace=0.32)
    save_all(fig, stem)
    plt.close(fig)


def fig_pressure(summary: dict) -> None:
    """Primary p_max convergence figure."""
    _plot_convergence(
        summary=summary,
        qoi_key="p_max",
        stem=OUTDIR / "fig_grid_convergence",
        y_column="mean",
        sem_column="sem",
        ylabel=r"$\langle p_{\max}\rangle_{\,t\in[6,10]}$",
        err_ylabel=(
            r"$\,\left|\,\langle p_{\max}\rangle - "
            r"\phi_{\mathrm{ext}}\,\right|$"
        ),
        # Pressure C--M--F is pre-asymptotic, so visually connect only M--F--XF.
        connect_from=1,
        legend_loc_left="lower right",
        legend_loc_right="lower right",
    )


def fig_density(summary: dict) -> None:
    """Secondary rho_max convergence figure."""
    _plot_convergence(
        summary=summary,
        qoi_key="rho_max",
        stem=OUTDIR / "fig_rho_grid_convergence",
        y_column="mean",
        sem_column="sem",
        ylabel=r"$\langle \rho_{\max}\rangle_{\,t\in[6,10]}$",
        err_ylabel=(
            r"$\,\left|\,\langle \rho_{\max}\rangle - "
            r"\phi_{\mathrm{ext}}\,\right|$"
        ),
        # Density is monotonic on both triplets, so connect all grids.
        connect_from=0,
        legend_loc_left="upper right",
        legend_loc_right="lower right",
    )


def main() -> None:
    OUTDIR.mkdir(exist_ok=True)
    summary = load_summary()

    fig_pressure(summary)
    fig_density(summary)


if __name__ == "__main__":
    main()
