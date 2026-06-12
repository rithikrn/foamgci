"""
make_figures.py
---------------
Generate the primary grid-convergence figure for the manuscript,
in a style consistent with the Journal of Fluid Mechanics.

Reads gci_summary.json produced by analyze.py; if that file is
missing the script tells you to run analyze.py first.

Emits both vector PDF (preferred for journal line art) and 600-dpi
PNG (preview / Markdown embedding) for every figure.
"""
from __future__ import annotations

import json
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

# Output formats. PDF stays first so it's the canonical/preferred copy.
EXTS: tuple[str, ...] = ("pdf", "png")
PNG_DPI = 600  # JFM/CUP raster spec for line art with text

mpl.rcParams.update({
    # --- typography ---
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

    # --- hairline rules / journal-grade line widths ---
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
    "figure.dpi": 150,            # screen preview only

    # Embedded TrueType in PDF (Type 42); harmless for PNG.
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def save_all(fig: plt.Figure, stem: Path, exts: Iterable[str] = EXTS) -> None:
    """Save `fig` to `<stem>.<ext>` for every extension in `exts`."""
    for ext in exts:
        out = stem.with_suffix(f".{ext}")
        # dpi only matters for raster; vector formats ignore it.
        fig.savefig(out, dpi=PNG_DPI if ext == "png" else None)
        print(f"  wrote {out}")


def load_summary() -> dict:
    if not SUMMARY.is_file():
        raise SystemExit(
            f"\n  ERROR: {SUMMARY} not found.\n"
            "  Run `python3 analyze.py` first to produce it.\n"
        )
    return json.loads(SUMMARY.read_text())


def fig_combined(summary: dict, stem: Path) -> None:
    cases = summary["cases"]
    cases.sort(key=lambda c: -c["dx"])  # coarse -> fine
    h_vals     = np.array([c["dx"]         for c in cases])
    phi_vals   = np.array([c["p_max_mean"] for c in cases])
    sigma_vals = np.array([c["p_max_std"]  for c in cases])
    labels     = [c["label"] for c in cases]
    markers    = ["o", "s", "D", "^"]

    phi_ext = summary["phi_star"]
    err = np.abs(phi_vals - phi_ext)

    tB = summary["triplet_B_MFXF"]
    p_obs = tB.get("p_obs", float("nan"))
    has_p = tB.get("regime") == "monotonic" and p_obs == p_obs

    fig, axes = plt.subplots(1, 2, figsize=(6.85, 2.95))

    # ----- left: phi vs h -----
    ax = axes[0]
    ax.axhline(phi_ext, color="0.55", lw=0.7, linestyle=(0, (4, 2)),
               label=fr"Richardson $\phi_{{\mathrm{{ext}}}}={phi_ext:.4f}$")
    ax.plot(h_vals[1:], phi_vals[1:], color="C0", lw=0.6, zorder=2)
    for i in range(len(cases)):
        ax.errorbar(h_vals[i], phi_vals[i], yerr=sigma_vals[i],
                    marker=markers[i],
                    markersize=4.5, markerfacecolor="white",
                    markeredgewidth=0.8, color="k", linestyle="none",
                    capsize=2, elinewidth=0.6, zorder=5)
        ax.annotate(labels[i], xy=(h_vals[i], phi_vals[i]),
                    xytext=(6, -2), textcoords="offset points",
                    fontsize=7.5)
    ax.set_xscale("log")
    ax.set_xlabel(r"grid spacing $h$")
    ax.set_ylabel(r"$\langle p_{\max}\rangle_{\,t\in[3,10]}$")
    ymin = min(phi_vals.min(), phi_ext) - 0.02
    ymax = max(phi_vals.max(), phi_ext) + 0.02
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(h_vals.min() * 0.7, h_vals.max() * 1.7)
    ax.legend(loc="lower right", borderaxespad=0.4)
    ax.xaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, subs=(1.0,)))
    ax.xaxis.set_minor_locator(mpl.ticker.LogLocator(
        base=10.0, subs=(0.2, 0.4, 0.6, 0.8)))
    ax.text(0.025, 0.96, "(a)", transform=ax.transAxes,
            fontsize=9, va="top")

    # ----- right: log-log error -----
    ax = axes[1]
    h_ref = np.array([h_vals.min() * 0.7, h_vals.max() * 1.4])
    anchor_h = h_vals[-1]
    anchor_err = err[-1] if err[-1] > 0 else err[err > 0].min()
    for p_slope, ls, lbl in (
        (1.0, (0, (1, 1.5)), r"$\mathcal{O}(h)$"),
        (2.0, (0, (4, 2)),  r"$\mathcal{O}(h^{2})$"),
    ):
        ref = anchor_err * (h_ref / anchor_h) ** p_slope
        ax.plot(h_ref, ref, color="0.55", lw=0.7, linestyle=ls, label=lbl)

    if has_p:
        h_obs = np.array([h_vals[1], h_vals[-1] * 0.9])
        anchor_h2 = h_vals[1]
        anchor_err2 = err[1] if err[1] > 0 else err[err > 0].min()
        obs = anchor_err2 * (h_obs / anchor_h2) ** p_obs
        ax.plot(h_obs, obs, color="C3", lw=0.9,
                label=fr"observed $p={p_obs:.2f}$")

    for i in range(len(cases)):
        if err[i] <= 0:
            continue
        ax.plot(h_vals[i], err[i], marker=markers[i], color="k",
                markersize=4.5, markerfacecolor="white",
                markeredgewidth=0.8, linestyle="none", zorder=5)
        ax.annotate(labels[i], xy=(h_vals[i], err[i]),
                    xytext=(5, -3 if i % 2 == 0 else 6),
                    textcoords="offset points", fontsize=7.5)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"grid spacing $h$")
    ax.set_ylabel(r"$\,\left|\,\langle p_{\max}\rangle - "
                  r"\phi_{\mathrm{ext}}\,\right|$")
    ax.set_xlim(h_ref[0], h_ref[1])
    ax.legend(loc="lower right", borderaxespad=0.4)
    ax.xaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, subs=(1.0,)))
    ax.xaxis.set_minor_locator(mpl.ticker.LogLocator(
        base=10.0, subs=(0.2, 0.4, 0.6, 0.8)))
    ax.text(0.025, 0.96, "(b)", transform=ax.transAxes,
            fontsize=9, va="top")

    fig.subplots_adjust(wspace=0.32)
    save_all(fig, stem)
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(exist_ok=True)
    summary = load_summary()
    # Pass a stem (no suffix); save_all() appends .pdf and .png.
    fig_combined(summary, OUTDIR / "fig_grid_convergence")


if __name__ == "__main__":
    main()