"""
make_aux_figures.py
-------------------
Auxiliary figure for the wedge case:
  fig_domain.pdf : 15-degree ramp geometry, boundary conditions, and the
                   analytical oblique shock (angle beta) overlaid.

Note: the forward-step example also produced a peak-location convergence
panel. That panel tracked the (x, y) of a pointwise pressure maximum. This
case's QoI is an area-averaged wall pressure with no pointwise location, so
the peak-location panel does not apply and is intentionally omitted.

Reads gci_summary.json only to label beta / p2/p1 if present; the figure is
produced even without it.
"""
from __future__ import annotations

import json
import sys
import subprocess
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrowPatch

from data import DOMAIN_X, DOMAIN_Y, THETA_DEG, ramp_y
from oblique_shock import oblique_shock
from data import M_INFLOW, GAMMA


HERE = Path(__file__).parent
SUMMARY = HERE / "gci_summary.json"
OUTDIR = HERE / "figures"


mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIX Two Text", "STIXGeneral", "Times New Roman",
                   "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 9.0,
    "axes.labelsize": 9.0,
    "axes.titlesize": 9.0,
    "legend.fontsize": 8.0,
    "xtick.labelsize": 8.0,
    "ytick.labelsize": 8.0,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "legend.frameon": False,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42,
})


def fig_domain(out: Path, ref) -> None:
    x0, x1 = DOMAIN_X
    y0, y1 = DOMAIN_Y
    xt = 0.0                       # ramp tip
    beta = np.radians(ref.beta_deg)

    fig, ax = plt.subplots(figsize=(5.6, 2.7))

    # Fluid boundary outline.
    outline = [
        (x0, y0), (xt, y0), (x1, ramp_y(x1)), (x1, y1), (x0, y1), (x0, y0),
    ]
    xs, ys = zip(*outline)
    ax.plot(xs, ys, color="k", lw=0.9)

    # Solid wedge body (below the ramp), shaded.
    body = [(xt, y0), (x1, y0), (x1, ramp_y(x1))]
    ax.add_patch(Polygon(body, closed=True, facecolor="0.85",
                         edgecolor="k", lw=0.0, hatch="//", zorder=0))
    ax.plot([xt, x1], [y0, ramp_y(x1)], color="k", lw=1.1)   # ramp surface

    # Inflow arrows.
    for y in np.linspace(0.15 * y1, 0.85 * y1, 5):
        ax.add_patch(FancyArrowPatch(
            (x0 - 0.05, y), (x0 - 0.005, y),
            arrowstyle="-|>", mutation_scale=8, lw=0.7, color="0.2",
            shrinkA=0, shrinkB=0))
    ax.text(x0 - 0.05, y1 * 1.04, fr"inflow: $M={ref.M1:.0f}$, $p=T=1$",
            fontsize=8.5, ha="left")
    ax.text(x1 + 0.004, 0.5 * y1, "outflow", fontsize=8.5, ha="left",
            rotation=-90, va="center")
    ax.text(0.5 * (x0 + x1), y1 * 1.02, "top (symmetry)", fontsize=8.0,
            ha="center")
    ax.text(0.5 * x0, -0.06 * y1, "upstream wall (symmetry)", fontsize=8.0,
            ha="center", va="top")
    ax.text(0.55 * x1, 0.30 * ramp_y(x1), r"15$^\circ$ ramp (slip)",
            fontsize=8.0, ha="center", va="center", color="0.25",
            rotation=THETA_DEG)

    # Analytical oblique shock from the ramp tip.
    x_sh = np.array([xt, x1])
    y_sh = (x_sh - xt) * np.tan(beta)
    y_sh = np.minimum(y_sh, y1)
    ax.plot(x_sh, y_sh, color="C3", lw=1.0, linestyle=(0, (5, 2)))
    ax.text(0.62 * x1, 0.62 * x1 * np.tan(beta),
            fr"oblique shock, $\beta={ref.beta_deg:.1f}^\circ$",
            fontsize=8.0, color="C3", rotation=ref.beta_deg, va="bottom")
    ax.text(0.40 * x1, 0.78 * y1,
            fr"$p_2/p_1={ref.p2_p1:.3f}$, $M_2={ref.M2:.2f}$",
            fontsize=8.0, ha="left")

    ax.set_xlim(x0 - 0.07, x1 + 0.06)
    ax.set_ylim(-0.12 * y1, y1 * 1.18)
    ax.set_aspect("equal")
    ax.axis("off")

    fig.savefig(out)
    plt.close(fig)
    print(f"  wrote {out}")


def main() -> None:
    OUTDIR.mkdir(exist_ok=True)
    ref = oblique_shock(M_INFLOW, THETA_DEG, GAMMA)
    fig_domain(OUTDIR / "fig_domain.pdf", ref)

    for f in OUTDIR.glob("*.pdf"):
        try:
            subprocess.run(["pdftoppm", "-r", "200", "-png",
                            str(f), str(f.with_suffix(""))],
                           check=True, capture_output=True)
        except (FileNotFoundError, subprocess.CalledProcessError):
            pass


if __name__ == "__main__":
    main()
