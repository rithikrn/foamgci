"""foamgci.plot — optional convergence plot.

Importing this module requires matplotlib. The rest of foamgci does
not, by design: plotting is a downstream concern and we keep core
verification pure NumPy / standard-library to maximise portability.
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    from .report import ReportTable


def plot_convergence(
    report: "ReportTable",
    out_path: str | Path | None = None,
    dpi: int = 150,
    show: bool = False,
):
    """Two-panel convergence figure: ⟨φ⟩ vs h, and GCI vs h.

    Left panel
        Time-averaged mean at the bow-shock stagnation point vs cell
        size, with autocorrelation-corrected SEM error bars and a
        horizontal line at the analytical reference if available.
    Right panel
        Roache GCI on the fine grid of each consecutive triplet, on a
        log-log scale, with a slope-:math:`\\hat p` reference line.

    Returns the matplotlib Figure object so callers can post-process.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "foamgci.plot requires matplotlib; install with "
            "`pip install foamgci[plot]` or `pip install matplotlib`."
        ) from exc

    hs = np.array([c.h for c in report.cases])
    means = np.array([s.mean for s in report.stats])
    sems = np.array([s.sem for s in report.stats])

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.0, 3.6))

    # Left ---------------------------------------------------------------
    axL.errorbar(hs, means, yerr=sems, fmt="o-", color="C0",
                 capsize=3, label=r"$\langle\phi\rangle\pm\mathrm{SEM}_\tau$")
    if report.reference_value is not None:
        axL.axhline(report.reference_value, color="k", lw=1.0, ls="--",
                    label=report.reference_label)
    if report.gcis:
        phi_ext = report.gcis[-1].phi_exact
        axL.axhline(phi_ext, color="C3", lw=1.0, ls=":",
                    label=r"Richardson $\phi_{\mathrm{exact}}$")
    axL.set_xscale("log")
    axL.set_xlabel(r"cell size $h$")
    axL.set_ylabel(rf"$\langle {report.field}_{{\rm {report.quantity}}}\rangle$")
    axL.grid(alpha=0.3)
    axL.legend(fontsize=8, loc="best")
    axL.invert_xaxis()  # finer grids to the right

    # Right --------------------------------------------------------------
    if report.gcis:
        gci_h = np.array([g.h_fine for g in report.gcis])
        gci_v = np.array([g.gci_fine_21_pct for g in report.gcis])
        axR.loglog(gci_h, gci_v, "s-", color="C2",
                   label=r"GCI$_{21}^{\mathrm{fine}}$")
        # Reference slope p̂ line
        if len(gci_h) >= 1:
            p_hat = report.gcis[-1].p_apparent
            h_ref = gci_h[0]
            g_ref = gci_v[0]
            xs = np.linspace(gci_h.min(), gci_h.max(), 50)
            ys = g_ref * (xs / h_ref) ** p_hat
            axR.loglog(xs, ys, "k--", lw=1.0,
                       label=rf"slope $\hat p={p_hat:.2f}$")
        axR.set_xlabel(r"cell size $h$")
        axR.set_ylabel("GCI (%)")
        axR.grid(which="both", alpha=0.3)
        axR.legend(fontsize=8, loc="best")
        axR.invert_xaxis()

    fig.tight_layout()
    if out_path is not None:
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    if show:  # pragma: no cover
        plt.show()
    return fig
