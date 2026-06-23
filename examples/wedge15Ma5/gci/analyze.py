from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

from data import (  # noqa: E402
    GRIDS, T_STAT, M_INFLOW, THETA_DEG, GAMMA, P_INF,
)
from oblique_shock import oblique_shock  # noqa: E402
from foamgci.reader import (  # noqa: E402
    read_surface_field_value,
    read_fieldminmax,
)
from foamgci.stats import window_stats  # noqa: E402
from foamgci.gci import gci_over_hierarchy  # noqa: E402

# order_consistency + extremum_location_drift are library functions from
# foamgci v3.5+. Inline compact equivalents so this file is swap-and-run on
# an older (3.4.x) install too; the library versions are used when present.
try:
    from foamgci.gci import order_consistency, extremum_location_drift  # noqa
except ImportError:
    from types import SimpleNamespace as _NS

    def order_consistency(results, rel_tol=0.15, p_bounds=(0.1, 6.0)):
        mono = [r for r in results if r.regime == "monotonic"]
        labs = [f"{r.label_coarse}-{r.label_medium}-{r.label_fine}" for r in mono]
        orders = [float(r.p_apparent) for r in mono]
        lo, hi = p_bounds
        clip = any(o <= lo + 1e-9 or o >= hi - 1e-9 for o in orders)
        if len(mono) < 2:
            return _NS(orders=orders, triplet_labels=labs, n_monotonic=len(mono),
                       p_spread=float("nan"), p_rel_spread=float("nan"),
                       rel_tol=rel_tol, any_order_clipped=clip, consistent=None,
                       verdict="indeterminate",
                       note=f"need >=2 monotonic triplets; orders={[round(o,3) for o in orders]}")
        spread = max(orders) - min(orders)
        rel = spread / (sum(abs(o) for o in orders) / len(orders))
        ok = bool(rel <= rel_tol)
        return _NS(orders=orders, triplet_labels=labs, n_monotonic=len(mono),
                   p_spread=spread, p_rel_spread=rel, rel_tol=rel_tol,
                   any_order_clipped=clip, consistent=ok,
                   verdict="consistent" if ok else "inconsistent",
                   note=f"orders {[round(o,3) for o in orders]}; spread "
                        f"{100*rel:.1f}% vs tol {100*rel_tol:.0f}% "
                        f"(necessary, not sufficient)")

    def extremum_location_drift(locations, hs, labels=None, tol_cells=2.0):
        L = [np.asarray(p, float) for p in locations]; h = np.asarray(hs, float)
        n = len(L); labels = list(labels) if labels else [f"g{i}" for i in range(n)]
        if n < 2 or len(h) != n:
            return _NS(cell_displacements=[], physical_displacements=[],
                       pair_labels=[], n_grids=n, finest_cell_disp=float("nan"),
                       max_cell_disp=float("nan"), shrinking=False,
                       tol_cells=tol_cells, mesh_locked=None,
                       verdict="indeterminate", note="need >=2 grids")
        dim = min(p.size for p in L); phys = []; cells = []; pl = []
        for k in range(n - 1):
            d = float(np.linalg.norm(L[k+1][:dim] - L[k][:dim])); hk = float(h[k+1])
            phys.append(d); cells.append(d/hk if hk > 0 else float("nan"))
            pl.append(f"{labels[k]}->{labels[k+1]}")
        finest = cells[-1]
        shrink = all(cells[i+1] <= cells[i]*1.05+1e-12 for i in range(len(cells)-1))
        locked = bool(finest > tol_cells)
        note = (f"mesh-locked: two finest grids {finest:.1f} cells apart "
                f"(> {tol_cells:g}); per-pair {[round(c,1) for c in cells]} -> "
                f"location scales with h, tracks a mesh feature" if locked else
                f"converging: two finest grids agree to {finest:.1f} cells "
                f"(<= {tol_cells:g}); per-pair {[round(c,1) for c in cells]}")
        return _NS(cell_displacements=cells, physical_displacements=phys,
                   pair_labels=pl, n_grids=n, finest_cell_disp=finest,
                   max_cell_disp=max(cells), shrinking=shrink, tol_cells=tol_cells,
                   mesh_locked=locked,
                   verdict="mesh_locked" if locked else "converging", note=note)


# ----------------------------------------------------------------------------
# This case reads TWO OpenFOAM outputs per grid and reports both QoIs:
#
#   PRIMARY  (p_wall_ratio): area-averaged ramp-surface static pressure from
#            the surface-region area-average (`surfaceRegion` in OpenFOAM-4.x,
#            `surfaceFieldValue` in v5.0+) -- a different output type than the
#            forward step. It is an integrated functional with a definite
#            continuum value and is robust to per-cell grid noise, so it is a
#            better-posed Richardson/GCI target than a pointwise extremum. Its
#            observed order is still set by the shock-capturing scheme (vanLeer
#            + Euler), so analyze.py measures p_obs rather than assuming it.
#            Reference-anchored against the exact oblique-shock p2/p1. This QoI
#            carries the verdict and owns the summary keys the figures consume.
#
#   SECONDARY (p_max): global max(p) from `fieldMinMax` (the SAME output type
#            the forward step used; here for cross-case consistency and as a
#            multi-file demonstration). The four-grid run shows this is NOT a
#            valid verification target: max(p) tracks a near-tip mesh feature,
#            not the post-shock plateau. Its VALUE (~6.8) sits ~42% ABOVE the
#            analytical post-shock p2, and its LOCATION marches to the sharp
#            ramp tip at a roughly constant ~13-cell offset (position scales
#            with h). It is steady within each window -- so the within-window
#            wander check reports "localized" (false comfort) -- yet the
#            across-grid location-drift check flags it "mesh_locked", and the
#            analytical reference shows the 42% gap. That is the real lesson:
#            a clean-looking monotone GCI on a pointwise extremum can be junk;
#            an independent reference AND across-grid location drift are what
#            expose it. Reported as a diagnostic / cautionary contrast with the
#            well-posed surface integral, not as a convergent QoI.
# ----------------------------------------------------------------------------
QOI_NAME = "p_wall_ratio"
QOI_DESCRIPTION = "area-averaged ramp-surface static pressure / p_inf"

# Max in-window extremum wander (in cell widths) tolerated before the pointwise
# QoI is flagged not-localized. Same threshold as the forward-step example.
LOC_WANDER_CELLS_MAX = 5.0


def _clean_json(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {k: _clean_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean_json(v) for v in obj]
    if isinstance(obj, np.generic):
        return _clean_json(obj.item())
    if isinstance(obj, float):
        if math.isnan(obj) or math.isinf(obj):
            return None
        return obj
    return obj


def _fmt(v: Any, spec: str = ".4f") -> str:
    if v is None:
        return "--"
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    if math.isnan(f) or math.isinf(f):
        return "--"
    return format(f, spec)


# ----------------------------------------------------------------------------
# Per-grid statistics for each source
# ----------------------------------------------------------------------------

def _case_stats_wall(grid) -> dict:
    """Primary QoI: time-averaged area-averaged wall pressure (surfaceFieldValue)."""
    data = read_surface_field_value(grid.sfv_path, column="p")
    ws = window_stats(data.time, data.value, T_STAT[0], T_STAT[1])
    return {
        "label": grid.label,
        "n_cells": grid.n_cells,
        "dx": grid.dx,
        "column": data.column,
        "n_samples_stat": ws.n,
        "mean": ws.mean / P_INF,
        "std": ws.std / P_INF,
        "sem": ws.sem / P_INF,
        "tau_int": ws.tau_int,
        "n_eff": ws.n_eff,
        "kpss_stat": ws.kpss_stat,
        "kpss_p": ws.kpss_p,
        "kpss_stationary": bool(ws.kpss_stationary_5pct),
    }


def _case_stats_pmax(grid) -> dict:
    """Secondary QoI: time-averaged global max(p) and its localization
    (fieldMinMax). Mirrors the forward-step pointwise-extremum diagnostics."""
    data = read_fieldminmax(grid.fmm_path, field="p")
    ws = window_stats(data.time, data.max, T_STAT[0], T_STAT[1])

    # In-window extremum wander, in cell widths (5th-95th pct spread of loc_max).
    loc_wander_cells = float("nan")
    loc_x = loc_y = float("nan")
    if data.loc_max is not None:
        m = (data.time >= T_STAT[0]) & (data.time <= T_STAT[1])
        if m.any():
            locw = data.loc_max[m]
            loc_x, loc_y = (float(v) for v in np.median(locw, axis=0)[:2])
            span_x = float(np.percentile(locw[:, 0], 95) - np.percentile(locw[:, 0], 5))
            span_y = float(np.percentile(locw[:, 1], 95) - np.percentile(locw[:, 1], 5))
            loc_wander_cells = float(np.hypot(span_x, span_y) / grid.dx)

    return {
        "label": grid.label,
        "n_cells": grid.n_cells,
        "dx": grid.dx,
        "field": data.field,
        "n_samples_stat": ws.n,
        "mean": ws.mean / P_INF,
        "std": ws.std / P_INF,
        "sem": ws.sem / P_INF,
        "tau_int": ws.tau_int,
        "n_eff": ws.n_eff,
        "kpss_stat": ws.kpss_stat,
        "kpss_p": ws.kpss_p,
        "kpss_stationary": bool(ws.kpss_stationary_5pct),
        "loc_x": loc_x,
        "loc_y": loc_y,
        "loc_wander_cells": loc_wander_cells,
    }


def _triplet_dict(g, label: str) -> dict:
    return {
        "label": label,
        "phi1": g.phi_fine, "phi2": g.phi_medium, "phi3": g.phi_coarse,
        "h1": g.h_fine, "h2": g.h_medium, "h3": g.h_coarse,
        "r": g.r21,
        "regime": g.regime,
        "p_obs": g.p_apparent,
        "phi_ext": g.phi_exact,
        "gci_fine_pct": g.gci_fine_21_pct,
        "gci_med_pct": g.gci_medium_32_pct,
        "asymptotic_ratio": g.asymptotic_ratio,
        "note": g.note,
    }


def _oc_dict(gcis) -> dict:
    oc = order_consistency(gcis)
    return {
        "orders": oc.orders,
        "triplet_labels": oc.triplet_labels,
        "n_monotonic": oc.n_monotonic,
        "p_spread": oc.p_spread,
        "p_rel_spread": oc.p_rel_spread,
        "rel_tol": oc.rel_tol,
        "any_order_clipped": oc.any_order_clipped,
        "consistent": oc.consistent,
        "verdict": oc.verdict,
        "note": oc.note,
    }


def _gci_blocks(cases: list[dict]) -> tuple[dict, dict, list]:
    """Run GCI on the C-M-F and M-F-XF triplets for a set of per-grid means."""
    order = ["Coarse", "Medium", "Fine", "Extra-fine"]
    by = {c["label"]: c for c in cases}
    means = [by[k]["mean"] for k in order]
    hs = [by[k]["dx"] for k in order]
    gcis = gci_over_hierarchy(means, hs, order)
    tA = _triplet_dict(gcis[0], "A : C-M-F")
    tB = _triplet_dict(gcis[1], "B : M-F-XF")
    return tA, tB, gcis


def _pick_phi_star(gcis, by, tA, tB):
    """Richardson phi_star from the deepest monotonic triplet."""
    if gcis[1].regime == "monotonic":
        return gcis[1].phi_exact, tB["label"], tB["gci_fine_pct"], gcis[1].phi_fine
    if gcis[0].regime == "monotonic":
        return gcis[0].phi_exact, tA["label"], tA["gci_fine_pct"], gcis[0].phi_fine
    xf = by["Extra-fine"]["mean"]
    return xf, "Extra-fine (no monotonic triplet)", float("nan"), xf


def _error_table(cases, ref_value):
    return [
        {
            "label": c["label"], "dx": c["dx"], "mean": c["mean"],
            "err": abs(c["mean"] - ref_value),
            "rel_pct": 100.0 * abs(c["mean"] - ref_value) / abs(ref_value),
        }
        for c in cases
    ]


def _make_diagnostics_figure(out: dict) -> None:
    """Three-panel results figure: primary convergence, wall-pressure error,
    and the max(p) mesh-locking plot (location distance-to-tip vs h). Saved to
    figures/fig_diagnostics.{pdf,png}. Skips quietly if matplotlib is absent."""
    try:
        import matplotlib as mpl
        mpl.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:  # noqa: BLE001
        print(f"  (diagnostics figure skipped: {e})")
        return

    ref = out["reference"]["p2_p1"]
    pc = out["cases"]
    dx = np.array([c["dx"] for c in pc]); m = np.array([c["mean"] for c in pc])
    lab = [c["label"] for c in pc]
    sec = out["secondary_qoi_fieldminmax"]
    sc = [c for c in sec["cases"]
          if isinstance(c.get("loc_x"), (int, float))
          and not math.isnan(c["loc_x"])]
    phi_ext = out.get("phi_star")

    mpl.rcParams.update({"font.family": "serif", "font.size": 9,
                         "axes.linewidth": 0.6, "xtick.direction": "in",
                         "ytick.direction": "in"})
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.2))

    ax[0].axhline(ref, ls="--", c="crimson", lw=1.2,
                  label=fr"analytical $p_2/p_1$={ref:.4f}")
    if phi_ext is not None:
        ax[0].axhline(phi_ext, ls=":", c="0.4", lw=1.1,
                      label=f"Richardson {phi_ext:.4f}")
    ax[0].plot(dx, m, "o-", c="k", ms=6, lw=1)
    for l, x, y in zip(lab, dx, m):
        ax[0].annotate(l, (x, y), fontsize=6.5, xytext=(0, 6),
                       textcoords="offset points", ha="center")
    ax[0].set_xscale("log"); ax[0].set_xlabel("grid spacing $h$")
    ax[0].set_ylabel(r"$\langle p_{\rm wall}\rangle/p_\infty$")
    ax[0].set_title("(a) primary QoI -> $p_2/p_1$", fontsize=8.5)
    ax[0].legend(fontsize=6.5, loc="lower right")

    err = np.abs(m - ref)
    xx = np.array([dx.min(), dx.max()])
    ax[1].loglog(dx, err, "o-", c="k", ms=6, lw=1)
    ax[1].loglog(xx, err[-1]*(xx/dx[-1])**1, "--", c="0.5", lw=1,
                 label=r"$\mathcal{O}(h)$")
    ax[1].loglog(xx, err[-1]*(xx/dx[-1])**2, "-.", c="0.7", lw=1,
                 label=r"$\mathcal{O}(h^2)$")
    ax[1].set_xlabel("grid spacing $h$")
    ax[1].set_ylabel(r"$|\langle p_{\rm wall}\rangle-p_2/p_1|$")
    ax[1].set_title("(b) wall-pressure error", fontsize=8.5)
    ax[1].legend(fontsize=6.5)

    if len(sc) >= 2:
        dx2 = np.array([c["dx"] for c in sc])
        dist = np.array([math.hypot(c["loc_x"], c["loc_y"]) for c in sc])
        ax[2].loglog(dx2, dist, "s-", c="darkblue", ms=6, lw=1,
                     label="max(p) dist. to tip")
        cc = dist[-1] / dx2[-1]
        xx2 = np.array([dx2.min(), dx2.max()])
        ax[2].loglog(xx2, cc*xx2, "--", c="crimson", lw=1.2,
                     label=fr"$\propto h$ ({cc:.0f} cells)")
        for c, x, y in zip(sc, dx2, dist):
            ax[2].annotate(f"{y/x:.0f}c", (x, y), fontsize=6,
                           xytext=(4, -2), textcoords="offset points")
        ld = sec["localization"].get("across_grid_drift", {})
        verdict = ld.get("verdict", "")
        ax[2].set_title(f"(c) max(p) location vs $h$\n[{verdict}]", fontsize=8.5)
        ax[2].legend(fontsize=6.5, loc="upper left")
    else:
        ax[2].text(0.5, 0.5, "max(p) location\nunavailable", ha="center",
                   va="center", transform=ax[2].transAxes)
    ax[2].set_xlabel("grid spacing $h$")
    ax[2].set_ylabel("dist. of max(p) to ramp tip")

    fig.tight_layout()
    figdir = Path(__file__).parent / "figures"
    figdir.mkdir(exist_ok=True)
    fig.savefig(figdir / "fig_diagnostics.pdf", bbox_inches="tight")
    fig.savefig(figdir / "fig_diagnostics.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {figdir / 'fig_diagnostics.pdf'}")


def main() -> int:
    ref = oblique_shock(M_INFLOW, THETA_DEG, GAMMA)
    print(f"Case: 15-deg wedge, M={M_INFLOW}, gamma={GAMMA} (inviscid)")
    print(f"Analytical oblique shock: beta={ref.beta_deg:.4f} deg, "
          f"p2/p1={ref.p2_p1:.5f} (reference for BOTH QoIs)")
    print(f"Stationary window t in {T_STAT}\n")

    # --- input guards: BOTH files must exist for every grid -----------------
    missing = []
    for grid in GRIDS:
        if not grid.sfv_path.is_file():
            missing.append(str(grid.sfv_path))
        if not grid.fmm_path.is_file():
            missing.append(str(grid.fmm_path))
    if missing:
        print("ERROR: missing input file(s):")
        for p in missing:
            print(f"  - {p}")
        print("Populate gci/data/ from your OpenFOAM runs first "
              "(see gci/data/README.md).")
        return 2

    # ======================= PRIMARY QoI: wall pressure =====================
    print("[PRIMARY] surfaceFieldValue: area-averaged wall pressure on `obstacle`")
    wall_cases = []
    for grid in GRIDS:
        s = _case_stats_wall(grid)
        wall_cases.append(s)
        print(
            f"  {s['label']:11s} N={s['n_samples_stat']:4d} "
            f"p_wall={s['mean']:.6f} sigma={s['std']:.4g} "
            f"tau={s['tau_int']:.2f} SEM={s['sem']:.4g} "
            f"KPSS_p={s['kpss_p']:.2f} stationary={s['kpss_stationary']}"
        )

    by_wall = {c["label"]: c for c in wall_cases}
    tA, tB, gcis = _gci_blocks(wall_cases)
    oc_wall = _oc_dict(gcis)
    for t in (tA, tB):
        print(
            f"\n  Triplet {t['label']}: regime={t['regime']} "
            f"p_obs={_fmt(t['p_obs'])} phi_ext={_fmt(t['phi_ext'], '.5f')} "
            f"GCI_fine={_fmt(t['gci_fine_pct'])}% "
            f"R_fit={_fmt(t['asymptotic_ratio'], '.3f')}"
        )
        print(f"    {t['note']}")

    if oc_wall["orders"]:
        _pairs = ", ".join(f"{lbl}: p={p:.3f}" for lbl, p
                           in zip(oc_wall["triplet_labels"], oc_wall["orders"]))
    else:
        _pairs = "none"
    print(f"\n  Order consistency [{oc_wall['verdict']}] (necessary, not "
          f"sufficient, for asymptotic range): {_pairs}")
    print(f"    {oc_wall['note']}")

    phi_star, src, gci_band, phi_fine_sel = _pick_phi_star(gcis, by_wall, tA, tB)
    err_vs_ref = _error_table(wall_cases, ref.p2_p1)
    ext_err_pct = 100.0 * abs(phi_star - ref.p2_p1) / ref.p2_p1
    # GCI_fine is the Roache band around the selected triplet's FINE-GRID
    # solution, so the coverage test compares the fine-grid value (not the
    # extrapolate) against the reference, normalised the same way as the band.
    fine_err_pct = 100.0 * abs(phi_fine_sel - ref.p2_p1) / abs(phi_fine_sel)
    covered = (not math.isnan(gci_band)) and (fine_err_pct <= gci_band)

    print(f"\n  Richardson phi_star  = {phi_star:.5f} ({src})")
    print(f"  Analytical p2/p1     = {ref.p2_p1:.5f}")
    print(f"  |phi_extrap - ref|   = {ext_err_pct:.4f} %  (extrapolate accuracy)")
    print(f"  |phi_fine - ref|     = {fine_err_pct:.4f} %   vs  GCI_fine = "
          f"{_fmt(gci_band)} %  ->  reference "
          f"{'COVERED' if covered else 'NOT covered'} by fine-grid GCI band")
    print(f"  finest-grid error    = {err_vs_ref[-1]['rel_pct']:.4f} %")
    
    # ============= SECONDARY QoI: pointwise max(p) (fieldMinMax) =============
    print("\n[SECONDARY] fieldMinMax: global max(p) (diagnostic; same output "
          "type as the forward step)")
    pmax_cases = []
    for grid in GRIDS:
        s = _case_stats_pmax(grid)
        pmax_cases.append(s)
        print(
            f"  {s['label']:11s} N={s['n_samples_stat']:4d} "
            f"max_p={s['mean']:.6f} tau={s['tau_int']:.2f} "
            f"KPSS_p={s['kpss_p']:.2f} stationary={s['kpss_stationary']} "
            f"wander={_fmt(s['loc_wander_cells'], '.1f')} cells"
        )

    by_pmax = {c["label"]: c for c in pmax_cases}
    pA, pB, pgcis = _gci_blocks(pmax_cases)
    oc_pmax = _oc_dict(pgcis)
    pphi_star, psrc, pgci_band, pphi_fine_sel = _pick_phi_star(pgcis, by_pmax, pA, pB)
    pmax_err = _error_table(pmax_cases, ref.p2_p1)

    wanders = [c["loc_wander_cells"] for c in pmax_cases
               if isinstance(c["loc_wander_cells"], float)
               and not math.isnan(c["loc_wander_cells"])]
    med_wander = float(np.median(wanders)) if wanders else float("nan")
    localized = (not math.isnan(med_wander)) and (med_wander <= LOC_WANDER_CELLS_MAX)

    # Across-grid location drift: does the extremum converge to a fixed
    # physical point, or march with the mesh?  The within-window wander above
    # is steady on each grid yet cannot see a location that moves across grids
    # -- so it can (and here does) report "localized" for a mesh-singular
    # extremum.  This second check is what actually establishes well-posedness.
    drift_cases = sorted(
        (c for c in pmax_cases
         if isinstance(c["loc_x"], float) and not math.isnan(c["loc_x"])),
        key=lambda c: -c["dx"])  # coarse -> fine
    drift = extremum_location_drift(
        [(c["loc_x"], c["loc_y"]) for c in drift_cases],
        [c["dx"] for c in drift_cases],
        labels=[c["label"] for c in drift_cases],
    )

    # A pointwise extremum is a valid continuum QoI only if it is BOTH steady
    # within each window AND converges to a fixed physical location across
    # grids.  mesh_locked overrides a "localized" within-window verdict.
    well_posed_location = localized and (drift.mesh_locked is False)
    if drift.mesh_locked:
        loc_note = (
            "max(p) is steady within each window but its location is NOT a "
            "convergent physical point: it tracks a mesh feature near the ramp "
            "tip at a roughly constant cell offset (location scales with h). "
            "max(p) value sits well above the analytical post-shock pressure, "
            "so it is not a valid verification QoI against the oblique-shock "
            "reference -- the monotone GCI on it notwithstanding.")
    elif not localized:
        loc_note = ("max(p) location wanders within the window; the extremum "
                    "is not stably resolved on this grid.")
    else:
        loc_note = ("max(p) location is stable within each window and "
                    "converges across grids to a fixed point.")
    for t in (pA, pB):
        print(
            f"\n  Triplet {t['label']}: regime={t['regime']} "
            f"p_obs={_fmt(t['p_obs'])} GCI_fine={_fmt(t['gci_fine_pct'])}%"
        )
    print(f"\n  max(p) phi_star = {pphi_star:.5f} ({psrc}); "
          f"finest-grid |max_p - p2/p1| = {pmax_err[-1]['rel_pct']:.4f} %")
    print(f"  median in-window wander = {_fmt(med_wander, '.1f')} cells "
          f"(threshold {LOC_WANDER_CELLS_MAX:.0f}); within-window localized={localized}")
    print(f"  across-grid location drift [{drift.verdict}]: finest-pair "
          f"displacement {_fmt(drift.finest_cell_disp, '.1f')} cells "
          f"(tol {drift.tol_cells:g}); well_posed_location={well_posed_location}")
    print(f"  -> {loc_note}")

    # ============================= JSON summary =============================
    out = {
        "case": "wedge15Ma5",
        # --- top-level keys below describe the PRIMARY QoI and are kept stable
        #     for the figure scripts (make_figures.py reads these directly) ---
        "qoi": {"name": QOI_NAME, "description": QOI_DESCRIPTION,
                "source": "surfaceFieldValue (areaAverage(p) on patch obstacle)"},
        "stationary_window": list(T_STAT),
        "reference": {
            "kind": "analytical_oblique_shock",
            "M_inflow": M_INFLOW, "theta_deg": THETA_DEG, "gamma": GAMMA,
            "beta_deg": ref.beta_deg, "p2_p1": ref.p2_p1,
            "rho2_rho1": ref.rho2_rho1, "T2_T1": ref.T2_T1, "M2": ref.M2,
        },
        "cases": wall_cases,
        "triplet_A_CMF": tA,
        "triplet_B_MFXF": tB,
        "order_consistency": oc_wall,
        "phi_star": phi_star,
        "phi_star_source": src,
        "phi_star_rel_err_pct": ext_err_pct,
        "fine_grid_rel_err_pct": fine_err_pct,
        "reference_covered_by_gci": bool(covered),
        "error_table_vs_reference": err_vs_ref,
        # --- secondary QoI: pointwise max(p) from fieldMinMax (diagnostic) ---
        "secondary_qoi_fieldminmax": {
            "qoi": {
                "name": "p_max",
                "description": "global max(p) / p_inf",
                "source": "fieldMinMax (max(p), magnitude mode)",
            },
            "cases": pmax_cases,
            "triplet_A_CMF": pA,
            "triplet_B_MFXF": pB,
            "order_consistency": oc_pmax,
            "phi_star": pphi_star,
            "phi_star_source": psrc,
            "phi_star_rel_err_pct": 100.0 * abs(pphi_star - ref.p2_p1) / ref.p2_p1,
            "error_table_vs_reference": pmax_err,
            "localization": {
                "median_wander_cells": med_wander,
                "wander_threshold_cells": LOC_WANDER_CELLS_MAX,
                "within_window_localized": bool(localized),
                "localized": bool(localized),  # kept for backward compat
                "across_grid_drift": {
                    "verdict": drift.verdict,
                    "mesh_locked": drift.mesh_locked,
                    "finest_cell_disp": drift.finest_cell_disp,
                    "cell_displacements": drift.cell_displacements,
                    "pair_labels": drift.pair_labels,
                    "tol_cells": drift.tol_cells,
                    "shrinking": drift.shrinking,
                    "note": drift.note,
                },
                "well_posed_location": bool(well_posed_location),
                "note": loc_note,
            },
        },
        "inputs_per_grid": {
            "primary": "surfaceRegion.dat (area-averaged wall pressure; "
                       "surfaceFieldValue.dat on OpenFOAM v5.0+)",
            "secondary": "fieldMinMax.dat (global max p/rho)",
        },
    }
    out_path = Path(__file__).parent / "gci_summary.json"
    out_path.write_text(json.dumps(_clean_json(out), indent=2, allow_nan=False))
    print(f"\nWrote {out_path}")
    _make_diagnostics_figure(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
