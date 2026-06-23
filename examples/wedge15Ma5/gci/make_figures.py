"""
make_figures.py
---------------
Grid-convergence figure for the wedge case.

Reads gci_summary.json produced by analyze.py; if that file is missing the
script tells you to run analyze.py first.

Figure generated:
  - fig_grid_convergence.*  : wall-pressure ratio p_wall/p_inf vs h, with the
                              exact analytical p2/p1 and the Richardson
                              extrapolate; plus an error-vs-h panel measured
                              against the analytical reference.

Emits both vector PDF and 600-dpi PNG.
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
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "savefig.dpi": PNG_DPI,
    "figure.dpi": 150,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def save_all(fig: plt.Figure, stem: Path, exts: Iterable[str] = EXTS) -> None:
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
        return math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def fig_convergence(summary: dict) -> None:
    cases = sorted(summary["cases"], key=lambda c: -float(c["dx"]))
    h = np.array([float(c["dx"]) for c in cases])
    phi = np.array([float(c["mean"]) for c in cases])
    sem = np.array([float(c["sem"]) for c in cases])
    labels = [str(c["label"]) for c in cases]
    markers = ["o", "s", "D", "^"]

    p2p1 = float(summary["reference"]["p2_p1"])
    phi_star = summary.get("phi_star")
    phi_star = float(phi_star) if _is_number(phi_star) else float("nan")

    tB = summary.get("triplet_B_MFXF", {})
    p_obs = tB.get("p_obs")
    has_p = tB.get("regime") == "monotonic" and _is_number(p_obs)
    p_obs = float(p_obs) if has_p else float("nan")

    # Error measured against the EXACT analytical reference.
    err = np.abs(phi - p2p1)

    fig, axes = plt.subplots(1, 2, figsize=(6.85, 2.95))

    # (a) QoI vs h
    ax = axes[0]
    ax.axhline(p2p1, color="C3", lw=0.8, linestyle=(0, (4, 2)),
               label=fr"analytical $p_2/p_1={p2p1:.4f}$")
    if _is_number(phi_star):
        ax.axhline(phi_star, color="0.55", lw=0.7, linestyle=(0, (1, 1.5)),
                   label=fr"Richardson $\phi_{{\mathrm{{ext}}}}={phi_star:.4f}$")
    ax.plot(h, phi, color="C0", lw=0.6, zorder=2)
    for i in range(len(cases)):
        ax.errorbar(h[i], phi[i], yerr=sem[i], marker=markers[i],
                    markersize=4.5, markerfacecolor="white",
                    markeredgewidth=0.8, color="k", linestyle="none",
                    capsize=2, elinewidth=0.6, zorder=5)
        ax.annotate(labels[i], xy=(h[i], phi[i]), xytext=(6, -2),
                    textcoords="offset points", fontsize=7.5)
    ax.set_xscale("log")
    ax.set_xlabel(r"grid spacing $h$")
    ax.set_ylabel(r"$\langle p_{\mathrm{wall}}\rangle / p_\infty$")
    pad = max(0.02 * abs(phi.mean()), 3.0 * sem.max())
    ax.set_ylim(min(phi.min(), p2p1) - pad, max(phi.max(), p2p1) + pad)
    ax.set_xlim(h.min() * 0.7, h.max() * 1.7)
    ax.legend(loc="best", borderaxespad=0.4)
    ax.xaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, subs=(1.0,)))
    ax.xaxis.set_minor_locator(
        mpl.ticker.LogLocator(base=10.0, subs=(0.2, 0.4, 0.6, 0.8)))
    ax.text(0.025, 0.96, "(a)", transform=ax.transAxes, fontsize=9, va="top")

    # (b) error vs h, against the analytical reference
    ax = axes[1]
    positive = err[err > 0.0]
    if positive.size == 0:
        raise SystemExit("ERROR: all errors are zero (exact match?)")
    h_ref = np.array([h.min() * 0.7, h.max() * 1.4])
    anchor_h, anchor_err = h[-1], (err[-1] if err[-1] > 0 else positive.min())
    for p_slope, ls, lbl in ((1.0, (0, (1, 1.5)), r"$\mathcal{O}(h)$"),
                             (2.0, (0, (4, 2)), r"$\mathcal{O}(h^{2})$")):
        ax.plot(h_ref, anchor_err * (h_ref / anchor_h) ** p_slope,
                color="0.55", lw=0.7, linestyle=ls, label=lbl)
    if has_p:
        h_obs = np.array([h[1], h[-1] * 0.9])
        a_err = err[1] if err[1] > 0 else positive.min()
        ax.plot(h_obs, a_err * (h_obs / h[1]) ** p_obs, color="C3", lw=0.9,
                label=fr"$p_{{\mathrm{{obs}}}}={p_obs:.2f}$")
    for i in range(len(cases)):
        if err[i] <= 0:
            continue
        ax.plot(h[i], err[i], marker=markers[i], color="k", markersize=4.5,
                markerfacecolor="white", markeredgewidth=0.8,
                linestyle="none", zorder=5)
        ax.annotate(labels[i], xy=(h[i], err[i]),
                    xytext=(5, -3 if i % 2 == 0 else 6),
                    textcoords="offset points", fontsize=7.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"grid spacing $h$")
    ax.set_ylabel(r"$\left|\,\langle p_{\mathrm{wall}}\rangle/p_\infty - p_2/p_1\,\right|$")
    ax.set_xlim(h_ref[0], h_ref[1])
    ax.legend(loc="lower right", borderaxespad=0.4)
    ax.xaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, subs=(1.0,)))
    ax.xaxis.set_minor_locator(
        mpl.ticker.LogLocator(base=10.0, subs=(0.2, 0.4, 0.6, 0.8)))
    ax.text(0.025, 0.96, "(b)", transform=ax.transAxes, fontsize=9, va="top")

    fig.subplots_adjust(wspace=0.32)
    save_all(fig, OUTDIR / "fig_grid_convergence")
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(exist_ok=True)
    fig_convergence(load_summary())


if __name__ == "__main__":
    main()
