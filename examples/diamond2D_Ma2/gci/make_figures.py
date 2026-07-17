"""Figures for the diamond example, read from gci_summary.json.

Produces, in figures/:
  fig_convergence.pdf   smooth QoIs (Cd, p_front, p_rear) vs h with references
  fig_coupling.pdf      Cd from forceCoeffs vs Cd from facet pressures
  fig_meshlock.pdf      max(p) peak location vs h (mesh-locking, slope 1)

Skips quietly if matplotlib or the summary is missing.
"""
from __future__ import annotations

import json
from pathlib import Path

GCI = Path(__file__).parent
FIG = GCI / "figures"


def _load():
    p = GCI / "gci_summary.json"
    if not p.exists():
        print("no gci_summary.json; run analyze.py first")
        return None
    return json.loads(p.read_text())


def main() -> int:
    out = _load()
    if out is None:
        return 0
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception as e:  # noqa: BLE001
        print("matplotlib unavailable; skipping figures:", e)
        return 0
    plt.rcParams.update({"font.family": "serif", "font.size": 10,
                         "mathtext.fontset": "cm", "figure.dpi": 150})
    FIG.mkdir(exist_ok=True)

    # ---- smooth-QoI convergence ----------------------------------------
    sq = out["smooth_qois"]
    fig, axs = plt.subplots(1, len(sq), figsize=(4.2 * len(sq), 3.4))
    if len(sq) == 1:
        axs = [axs]
    labels = {"Cd": r"$C_d$", "p_front": r"$p_{\rm front}/p_\infty$",
              "p_rear": r"$p_{\rm rear}/p_\infty$"}
    for ax, (name, blk) in zip(axs, sq.items()):
        h = np.array([c["h"] for c in blk["cases"]])
        m = np.array([c["mean"] for c in blk["cases"]])
        ax.axhline(blk["reference"], color="0.4", ls="--", lw=1.0,
                   label="shock-expansion")
        ax.plot(h, m, "o-", color="#1f4e79", mfc="white", lw=1.4, label="CFD")
        if blk.get("deepest"):
            ax.plot(0.0, blk["deepest"]["phi_ext"], "s", color="#c0392b",
                    label="Richardson")
        ax.set_xlabel(r"representative $h$")
        ax.set_title(labels.get(name, name), fontsize=10)
        ax.set_xlim(left=-0.05 * h.max())
        ax.legend(frameon=False, fontsize=7.8)
    fig.tight_layout()
    fig.savefig(FIG / "fig_convergence.pdf")
    plt.close(fig)

    # ---- coupling: Cd two ways -----------------------------------------
    cp = out.get("coupling")
    if cp and cp.get("per_grid"):
        rows = [r for r in cp["per_grid"] if r["cd_force"] is not None]
        if rows:
            lab = [r["label"] for r in rows]
            cf = [r["cd_force"] for r in rows]
            cpr = [r["cd_pressure"] for r in rows]
            x = np.arange(len(lab))
            fig, ax = plt.subplots(figsize=(5.0, 3.4))
            ax.plot(x, cf, "o-", color="#1f4e79", mfc="white", label="from forceCoeffs")
            ax.plot(x, cpr, "s--", color="#c0392b", mfc="white",
                    label="from facet pressures")
            ax.axhline(out["reference"]["Cd_shock_expansion"], color="0.4", ls=":",
                       lw=1.0, label="shock-expansion")
            ax.set_xticks(x)
            ax.set_xticklabels(lab, rotation=20, ha="right", fontsize=8)
            ax.set_ylabel(r"$C_d$ (full airfoil)")
            ax.set_title("Coupled cross-check: two routes to drag", fontsize=10)
            ax.legend(frameon=False, fontsize=8)
            fig.tight_layout()
            fig.savefig(FIG / "fig_coupling.pdf")
            plt.close(fig)

    # ---- mesh-lock: max(p) location ------------------------------------
    dp = out.get("diagnostic_pmax")
    if dp and dp.get("cases"):
        cases = dp["cases"]
        # location wander is per-grid; reconstruct loc_x from across_grid note is
        # not available, so plot value and annotate the across-grid verdict.
        h = np.array([c["h"] for c in cases])
        m = np.array([c["mean"] for c in cases])
        fig, ax = plt.subplots(figsize=(5.0, 3.4))
        ax.axhline(dp["reference_p2"], color="0.4", ls="--", lw=1.0,
                   label=r"post-shock $p_2/p_\infty$")
        ax.plot(h, m, "o-", color="#6b6b6b", mfc="white", lw=1.4,
                label=r"$\max(p)$")
        ax.set_xlabel(r"representative $h$")
        ax.set_ylabel(r"$\max(p)/p_\infty$")
        verdict = (dp.get("across_grid") or {}).get("verdict", "")
        ax.set_title(f"max(p): {dp['pct_above_post_shock']:.0f}% above $p_2$ "
                     f"({verdict})", fontsize=10)
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(FIG / "fig_meshlock.pdf")
        plt.close(fig)

    # ---- entropy volume integral (smooth VOLUME QoI; convergence-only) --
    sv = (out.get("volume_qoi") or {}).get("S_vol") or {}
    if sv.get("cases") and sv.get("status") != "absent":
        h = np.array([c["h"] for c in sv["cases"]])
        m = np.array([c["mean"] for c in sv["cases"]])
        fig, ax = plt.subplots(figsize=(5.0, 3.4))
        ax.plot(h, m, "o-", color="#1f4e79", mfc="white", lw=1.4,
                label=r"$S_{\rm vol}=\int \Delta s\,dV$")
        deep = sv.get("deepest")
        if deep and deep.get("phi_ext") is not None:
            ax.plot(0.0, deep["phi_ext"], "s", color="#c0392b",
                    label="Richardson")
        ax.set_xlabel(r"representative $h$")
        ax.set_ylabel(r"entropy integral $S_{\rm vol}$")
        ax.set_xlim(left=-0.05 * h.max())
        p_obs = None if not deep else deep.get("p_obs")
        ttl = "Entropy volume integral (Oswatitsch wave-drag proxy)"
        if p_obs is not None:
            ttl += fr"; $p_{{\rm obs}}={p_obs:.2f}$"
        ax.set_title(ttl, fontsize=9.5)
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(FIG / "fig_entropy.pdf")
        plt.close(fig)

    print("figures written to", FIG)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
