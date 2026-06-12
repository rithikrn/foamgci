"""
make_aux_figures.py
-------------------
Auxiliary JFM-style figures:
  fig_domain.pdf       : forward-step geometry + boundary conditions
  fig_peak_location.pdf: convergence of (x,y) of peak-pressure on step face
Reads gci_summary.json for the peak-location data.
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
from matplotlib.patches import Rectangle, FancyArrowPatch


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
    "legend.handletextpad": 0.4,
    "legend.handlelength": 1.5,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42,
})


def fig_domain(out: Path) -> None:
    fig, ax = plt.subplots(figsize=(5.5, 2.4))

    fluid = [
        (0.0, 1.0), (3.0, 1.0), (3.0, 0.2),
        (0.6, 0.2), (0.6, 0.0), (0.0, 0.0), (0.0, 1.0),
    ]
    xs, ys = zip(*fluid)
    ax.plot(xs, ys, color="k", lw=0.9)
    ax.add_patch(Rectangle((0.6, 0.0), 2.4, 0.2,
                           facecolor="0.85", edgecolor="k", lw=0.0,
                           hatch="//", zorder=0))
    ax.plot([0.6, 0.6, 3.0], [0.2, 0.2, 0.2], color="k", lw=0.9)
    ax.plot([0.6, 0.6], [0.0, 0.2], color="k", lw=0.9)

    for y in np.linspace(0.10, 0.90, 5):
        ax.add_patch(FancyArrowPatch(
            (-0.2, y), (-0.02, y),
            arrowstyle="-|>", mutation_scale=8,
            lw=0.7, color="0.2", shrinkA=0, shrinkB=0,
        ))
    ax.text(-0.2, 1.06, r"inflow:  $M=3$, $\rho=p=1$",
            fontsize=8.5, ha="left")
    ax.text(3.05, 0.5, "outflow", fontsize=8.5, ha="left",
            rotation=-90, va="center")
    ax.text(1.5, 1.05, "top wall (slip)", fontsize=8.5, ha="center")
    ax.text(0.3, -0.07, "bottom wall (slip)", fontsize=8.5,
            ha="center", va="top")
    ax.text(1.8, 0.10, "step (slip)", fontsize=8.5,
            ha="center", va="center", color="0.25")
    ax.plot(0.6, 0.0, "o", markersize=3, markerfacecolor="C3",
            markeredgecolor="C3")
    ax.annotate("peak-$p$ location\n(post-bow-shock\nstagnation)",
                xy=(0.6, 0.0), xytext=(0.05, 0.45),
                fontsize=7.5, ha="left",
                arrowprops=dict(arrowstyle="->", color="C3", lw=0.6))

    ax.annotate("", xy=(0.0, -0.18), xytext=(3.0, -0.18),
                arrowprops=dict(arrowstyle="<->", color="0.3", lw=0.5))
    ax.text(1.5, -0.22, r"$L_x = 3$", fontsize=8, ha="center", va="top")
    ax.annotate("", xy=(3.2, 0.0), xytext=(3.2, 1.0),
                arrowprops=dict(arrowstyle="<->", color="0.3", lw=0.5))
    ax.text(3.27, 0.5, r"$L_y = 1$", fontsize=8, ha="left", va="center")
    ax.annotate("", xy=(0.0, 1.18), xytext=(0.6, 1.18),
                arrowprops=dict(arrowstyle="<->", color="0.3", lw=0.5))
    ax.text(0.3, 1.22, r"$0.6$", fontsize=8, ha="center", va="bottom")
    ax.annotate("", xy=(3.5, 0.0), xytext=(3.5, 0.2),
                arrowprops=dict(arrowstyle="<->", color="0.3", lw=0.5))
    ax.text(3.57, 0.10, r"$0.2$", fontsize=8, ha="left", va="center")

    ax.set_xlim(-0.45, 4.0)
    ax.set_ylim(-0.35, 1.45)
    ax.set_aspect("equal")
    ax.axis("off")

    fig.savefig(out)
    plt.close(fig)
    print(f"  wrote {out}")


def fig_peak_location(out: Path, summary: dict) -> None:
    cases = sorted(summary["cases"], key=lambda c: -c["dx"])
    h      = np.array([c["dx"]         for c in cases])
    x_peak = np.array([c["peak_loc_x"] for c in cases])
    y_peak = np.array([c["peak_loc_y"] for c in cases])
    labels = [c["label"] for c in cases]
    markers = ["o", "s", "D", "^"]

    fig, axes = plt.subplots(1, 2, figsize=(6.85, 2.85))

    # (a) x of peak p -> step face x = 0.6
    ax = axes[0]
    ax.axhline(0.6, color="0.55", lw=0.7, linestyle=(0, (4, 2)),
               label=r"step face $x = 0.6$")
    for i in range(len(cases)):
        ax.plot(h[i], x_peak[i], marker=markers[i], color="k",
                markersize=4.5, markerfacecolor="white",
                markeredgewidth=0.8, linestyle="none", zorder=5)
        ax.annotate(labels[i], xy=(h[i], x_peak[i]),
                    xytext=(6, -2), textcoords="offset points",
                    fontsize=7.5)
    ax.set_xscale("log")
    ax.set_xlabel(r"grid spacing $h$")
    ax.set_ylabel(r"$x$ of peak $p$")
    ax.set_xlim(h.min() * 0.7, h.max() * 1.7)
    ymin = min(x_peak.min(), 0.6) - 0.005
    ymax = max(x_peak.max(), 0.6) + 0.005
    ax.set_ylim(ymin, ymax)
    ax.legend(loc="lower right", borderaxespad=0.4)
    ax.text(0.025, 0.96, "(a)", transform=ax.transAxes,
            fontsize=9, va="top")

    # (b) y of peak p ~ h/2 (cell centre of first row)
    ax = axes[1]
    h_ref = np.array([h.min() * 0.7, h.max() * 1.4])
    ax.plot(h_ref, h_ref / 2.0, color="0.55", lw=0.7,
            linestyle=(0, (4, 2)), label=r"$y = h/2$  (first cell)")
    for i in range(len(cases)):
        ax.plot(h[i], y_peak[i], marker=markers[i], color="k",
                markersize=4.5, markerfacecolor="white",
                markeredgewidth=0.8, linestyle="none", zorder=5)
        ax.annotate(labels[i], xy=(h[i], y_peak[i]),
                    xytext=(6, -4), textcoords="offset points",
                    fontsize=7.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"grid spacing $h$")
    ax.set_ylabel(r"$y$ of peak $p$")
    ax.set_xlim(h_ref[0], h_ref[1])
    ax.legend(loc="lower right", borderaxespad=0.4)
    ax.text(0.025, 0.96, "(b)", transform=ax.transAxes,
            fontsize=9, va="top")

    fig.subplots_adjust(wspace=0.32)
    fig.savefig(out)
    plt.close(fig)
    print(f"  wrote {out}")


def main() -> None:
    OUTDIR.mkdir(exist_ok=True)
    summary = (json.loads(SUMMARY.read_text())
               if SUMMARY.is_file() else None)
    fig_domain(OUTDIR / "fig_domain.pdf")
    if summary is None:
        print("  skipped fig_peak_location (run analyze.py first)")
        return
    fig_peak_location(OUTDIR / "fig_peak_location.pdf", summary)

    for f in OUTDIR.glob("*.pdf"):
        try:
            subprocess.run(["pdftoppm", "-r", "200", "-png",
                            str(f), str(f.with_suffix(""))],
                           check=True, capture_output=True)
        except (FileNotFoundError, subprocess.CalledProcessError):
            pass


if __name__ == "__main__":
    main()
