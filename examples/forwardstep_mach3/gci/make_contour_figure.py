"""
make_contour_figure.py
----------------------
Produce a 4-panel Greenshields-style density contour figure at t = 4.

Reads snapshots/snap_<label>_t4.000.npz produced by extract_snapshot.py.

Matches the style of Greenshields et al. (2010), Fig. 6:
    30 contour lines uniformly distributed over 0.2568 <= rho <= 6.067,
    full-domain view (0 <= x <= 3, 0 <= y <= 1), step body shaded.
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
from matplotlib.patches import Rectangle


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


# Greenshields 2010 Fig. 6: 30 levels over [0.2568, 6.067].
RHO_MIN, RHO_MAX = 0.2568, 6.067
N_LEVELS = 30
LEVELS = np.linspace(RHO_MIN, RHO_MAX, N_LEVELS)

# Order of panels (top-left -> bottom-right).
PANEL_ORDER = ["Coarse", "Medium", "Fine", "Extra-fine"]


def find_snapshot(label: str) -> Path | None:
    key = label.lower().replace("-", "")
    matches = sorted(SNAPSHOTS.glob(f"snap_{key}_t*.npz"))
    if not matches:
        return None
    return matches[-1]


def draw_panel(ax, data: np.lib.npyio.NpzFile, label: str,
               n_cells: int) -> None:
    x = data["x"]
    y = data["y"]
    rho = data["rho"]

    cs = ax.contour(x, y, rho, levels=LEVELS,
                    colors="k", linewidths=0.35)
    # Step body (slip wall) shaded.
    ax.add_patch(Rectangle((0.6, 0.0), 2.4, 0.2,
                           facecolor="0.88", edgecolor="k", lw=0.4,
                           zorder=2))

    ax.set_xlim(0.0, 3.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_aspect("equal")
    ax.set_xticks([0, 1, 2, 3])
    ax.set_yticks([0, 0.5, 1])

    rho_actual = np.nanmax(rho)
    rho_min_actual = np.nanmin(rho)
    # Annotate panel.
    ax.text(0.02, 0.92,
            f"{label}  ($N={n_cells:,}$)",
            transform=ax.transAxes, fontsize=8.5, va="top")
    ax.text(0.98, 0.92,
            (fr"$\rho \in [{rho_min_actual:.3f},\,{rho_actual:.3f}]$"),
            transform=ax.transAxes, fontsize=7.0,
            va="top", ha="right")


def main() -> None:
    OUTDIR.mkdir(exist_ok=True)

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

    fig, axes = plt.subplots(2, 2, figsize=(7.0, 3.0),
                             sharex=True, sharey=True)
    for ax, label in zip(axes.flat, PANEL_ORDER):
        p = paths[label]
        if p is None:
            ax.text(0.5, 0.5, f"{label}\n(snapshot missing)",
                    ha="center", va="center",
                    transform=ax.transAxes, fontsize=9)
            ax.set_xticks([]); ax.set_yticks([])
            continue
        with np.load(p, allow_pickle=True) as d:
            draw_panel(ax, d, label, n_cells.get(label, -1))

    for ax in axes[1, :]:
        ax.set_xlabel(r"$x$")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$y$")

    fig.suptitle(
        fr"Density contours at $t = 4$: {N_LEVELS} equispaced levels in "
        fr"$\rho \in [{RHO_MIN},\,{RHO_MAX}]$ "
        "(cf. Greenshields et al. 2010, Fig. 6)",
        fontsize=8.5, y=1.0,
    )
    fig.subplots_adjust(wspace=0.05, hspace=0.25, top=0.92,
                        left=0.07, right=0.99, bottom=0.13)

    out = OUTDIR / "fig_density_contours.pdf"
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
