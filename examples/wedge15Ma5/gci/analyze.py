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
#            multi-file demonstration). In VALUE it also approaches p2, but the
#            post-shock field is a near-uniform plateau, so the extremum
#            LOCATION is degenerate and the localization (wander) check is
#            expected to flag it. Reported as diagnostic only -- the contrast
#            between a well-posed surface integral and a degenerate pointwise
#            extremum that target the *same* physical pressure is the point.
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
        "std": ws.std,
        "sem": ws.sem,
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
    for t in (tA, tB):
        print(
            f"\n  Triplet {t['label']}: regime={t['regime']} "
            f"p_obs={_fmt(t['p_obs'])} phi_ext={_fmt(t['phi_ext'], '.5f')} "
            f"GCI_fine={_fmt(t['gci_fine_pct'])}% "
            f"R_asym={_fmt(t['asymptotic_ratio'], '.3f')}"
        )
        print(f"    {t['note']}")

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
    pphi_star, psrc, pgci_band = _pick_phi_star(pgcis, by_pmax, pA, pB)
    pmax_err = _error_table(pmax_cases, ref.p2_p1)

    wanders = [c["loc_wander_cells"] for c in pmax_cases
               if isinstance(c["loc_wander_cells"], float)
               and not math.isnan(c["loc_wander_cells"])]
    med_wander = float(np.median(wanders)) if wanders else float("nan")
    localized = (not math.isnan(med_wander)) and (med_wander <= LOC_WANDER_CELLS_MAX)
    loc_note = (
        "max(p) sits on the uniform post-shock plateau: value approaches p2 but "
        "its location is degenerate"
        if not localized else
        "max(p) location is stable within the window"
    )
    for t in (pA, pB):
        print(
            f"\n  Triplet {t['label']}: regime={t['regime']} "
            f"p_obs={_fmt(t['p_obs'])} GCI_fine={_fmt(t['gci_fine_pct'])}%"
        )
    print(f"\n  max(p) phi_star = {pphi_star:.5f} ({psrc}); "
          f"finest-grid |max_p - p2/p1| = {pmax_err[-1]['rel_pct']:.4f} %")
    print(f"  median in-window wander = {_fmt(med_wander, '.1f')} cells "
          f"(threshold {LOC_WANDER_CELLS_MAX:.0f}); localized={localized}")
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
            "phi_star": pphi_star,
            "phi_star_source": psrc,
            "phi_star_rel_err_pct": 100.0 * abs(pphi_star - ref.p2_p1) / ref.p2_p1,
            "error_table_vs_reference": pmax_err,
            "localization": {
                "median_wander_cells": med_wander,
                "wander_threshold_cells": LOC_WANDER_CELLS_MAX,
                "localized": bool(localized),
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
