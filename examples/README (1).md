"""
make_contour_figure.py
----------------------
Produce a 4-panel pressure-contour figure of the steady wedge field,
one panel per grid, with the solid wedge body shaded and the analytical
oblique-shock line (angle beta) overlaid.

Reads snapshots/snap_<label>_t0.200.npz produced by extract_snapshot.py.

The captured shock sharpens from Coarse to Extra-fine and should line up
with the analytical beta from oblique_shock.py.
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
from matplotlib.patches import Polygon

from data import (DOMAIN_X, DOMAIN_Y, THETA_DEG, ramp_y, M_INFLOW, GAMMA)
from oblique_shock import oblique_shock


HERE = Path(__file__).parent
SNAPSHOTS = HERE / "snapshots"
OUTDIR = HERE / "figures"
SUMMARY = HERE / "gci_summary.json"


mpl.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIX Two Text", "STIXGeneral", "Times New Roman",
                   "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 9.0,
    "axes.labelsize": 9.0,
    "legend.fontsize": 8.0,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42,
})


# Pressure contour band: free-stream p1 = 1 up to a little above p2.
N_LEVELS = 25
PANEL_ORDER = ["Coarse", "Medium", "Fine", "Extra-fine"]


def find_snapshot(label: str) -> Path | None:
    key = label.lower().replace("-", "")
    matches = sorted(SNAPSHOTS.glob(f"snap_{key}_t*.npz"))
    return matches[-1] if matches else None


def draw_panel(ax, data, label: str, n_cells: int, levels, ref) -> None:
    x, y, p = data["x"], data["y"], data["p"]
    x0, x1 = DOMAIN_X
    y0, y1 = DOMAIN_Y

    ax.contour(x, y, p, levels=levels, colors="k", linewidths=0.35)

    # Solid wedge body.
    body = [(0.0, y0), (x1, y0), (x1, ramp_y(x1))]
    ax.add_patch(Polygon(body, closed=True, facecolor="0.88",
                         edgecolor="k", lw=0.4, zorder=2))

    # Analytical shock line.
    xs = np.array([0.0, x1])
    ys = np.minimum(xs * np.tan(np.radians(ref.beta_deg)), y1)
    ax.plot(xs, ys, color="C3", lw=0.7, linestyle=(0, (5, 2)), zorder=3)

    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_xticks([-0.1, 0.0, 0.1, 0.2, 0.3])
    ax.set_yticks([0.0, 0.075, 0.15])

    ax.text(0.02, 0.92, f"{label}  ($N={n_cells:,}$)",
            transform=ax.transAxes, fontsize=8.5, va="top")
    pmax = np.nanmax(p)
    ax.text(0.98, 0.92, fr"$p_{{\max}}={pmax:.2f}$",
            transform=ax.transAxes, fontsize=7.0, va="top", ha="right")


def main() -> None:
    OUTDIR.mkdir(exist_ok=True)
    ref = oblique_shock(M_INFLOW, THETA_DEG, GAMMA)
    levels = np.linspace(1.0, 1.05 * ref.p2_p1, N_LEVELS)

    n_cells = {}
    if SUMMARY.is_file():
        for c in json.loads(SUMMARY.read_text())["cases"]:
            n_cells[c["label"]] = c["n_cells"]

    paths = {label: find_snapshot(label) for label in PANEL_ORDER}
    missing = [k for k, v in paths.items() if v is None]
    if missing:
        print(f"  Missing snapshot files for: {missing}")
        print(f"  Run: python3 extract_snapshot.py "
              f"(or for one grid: --grid {missing[0]})")
        if all(v is None for v in paths.values()):
            return

    fig, axes = plt.subplots(4, 1, figsize=(5.6, 6.2),
                             sharex=True, sharey=True)
    for ax, label in zip(axes.flat, PANEL_ORDER):
        p = paths[label]
        if p is None:
            ax.text(0.5, 0.5, f"{label}\n(snapshot missing)", ha="center",
                    va="center", transform=ax.transAxes, fontsize=9)
            ax.set_xticks([]); ax.set_yticks([])
            continue
        with np.load(p, allow_pickle=True) as d:
            draw_panel(ax, d, label, n_cells.get(label, -1), levels, ref)

    axes[-1].set_xlabel(r"$x$")
    for ax in axes:
        ax.set_ylabel(r"$y$")

    fig.suptitle(
        fr"Steady pressure contours: {N_LEVELS} levels in "
        fr"$p \in [1,\,{1.05 * ref.p2_p1:.2f}]$; "
        fr"red dashed = analytical shock $\beta={ref.beta_deg:.1f}^\circ$",
        fontsize=8.5, y=0.995)
    fig.subplots_adjust(hspace=0.18, top=0.95, left=0.10, right=0.99,
                        bottom=0.07)

    out = OUTDIR / "fig_pressure_contours.pdf"
    fig.savefig(out)
    plt.close(fig)
    print(f"  wrote {out}")

    try:
        subprocess.run(["pdftoppm", "-r", "200", "-png",
                        str(out), str(out.with_suffix(""))],
                       check=True, capture_output=True)
    except (FileNotFoundError, subprocess.CalledProcessError):
        pass


if __name__ == "__main__":
    main()
