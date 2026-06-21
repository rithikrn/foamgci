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
from foamgci.reader import read_surface_field_value  # noqa: E402
from foamgci.stats import window_stats  # noqa: E402
from foamgci.gci import gci_over_hierarchy  # noqa: E402


# Single, integrated QoI for this case: the area-averaged ramp-surface static
# pressure (= post-shock pressure p2 for inviscid flow), normalised by the
# free-stream static pressure. Unlike a fieldMinMax extremum this is an
# integrated functional, for which Richardson/GCI is formally better founded.
QOI_NAME = "p_wall_ratio"
QOI_DESCRIPTION = "area-averaged ramp-surface static pressure / p_inf"


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


def _case_stats(grid) -> dict:
    """Per-grid time-averaged wall-pressure statistics."""
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


def main() -> int:
    ref = oblique_shock(M_INFLOW, THETA_DEG, GAMMA)
    print(f"Case: 15-deg wedge, M={M_INFLOW}, gamma={GAMMA} (inviscid)")
    print(f"Analytical oblique shock: beta={ref.beta_deg:.4f} deg, "
          f"p2/p1={ref.p2_p1:.5f} (reference QoI)")
    print(f"Stationary window t in {T_STAT}\n")

    for grid in GRIDS:
        if not grid.sfv_path.is_file():
            print(f"ERROR: missing input file {grid.sfv_path}")
            print("Populate gci/data/ from your OpenFOAM runs first "
                  "(see gci/data/README.md).")
            return 2

    cases = []
    for grid in GRIDS:
        s = _case_stats(grid)
        cases.append(s)
        print(
            f"  {s['label']:11s} N={s['n_samples_stat']:4d} "
            f"p_wall={s['mean']:.6f} sigma={s['std']:.4g} "
            f"tau={s['tau_int']:.2f} SEM={s['sem']:.4g} "
            f"KPSS_p={s['kpss_p']:.2f} stationary={s['kpss_stationary']}"
        )

    order = ["Coarse", "Medium", "Fine", "Extra-fine"]
    by = {c["label"]: c for c in cases}
    means = [by[k]["mean"] for k in order]
    hs = [by[k]["dx"] for k in order]

    gcis = gci_over_hierarchy(means, hs, order)
    tA = _triplet_dict(gcis[0], "A : C-M-F")
    tB = _triplet_dict(gcis[1], "B : M-F-XF")

    for t in (tA, tB):
        print(
            f"\n  Triplet {t['label']}: regime={t['regime']} "
            f"p_obs={_fmt(t['p_obs'])} phi_ext={_fmt(t['phi_ext'], '.5f')} "
            f"GCI_fine={_fmt(t['gci_fine_pct'])}% "
            f"R_asym={_fmt(t['asymptotic_ratio'], '.3f')}"
        )
        print(f"    {t['note']}")

    # Pick phi_star from the deepest monotonic triplet.
    if gcis[1].regime == "monotonic":
        phi_star, src, gci_band = gcis[1].phi_exact, tB["label"], tB["gci_fine_pct"]
    elif gcis[0].regime == "monotonic":
        phi_star, src, gci_band = gcis[0].phi_exact, tA["label"], tA["gci_fine_pct"]
    else:
        phi_star, src, gci_band = by["Extra-fine"]["mean"], "Extra-fine (no monotonic triplet)", float("nan")

    err_vs_ref = [
        {
            "label": c["label"], "dx": c["dx"], "mean": c["mean"],
            "err": abs(c["mean"] - ref.p2_p1),
            "rel_pct": 100.0 * abs(c["mean"] - ref.p2_p1) / ref.p2_p1,
        }
        for c in cases
    ]

    ext_err_pct = 100.0 * abs(phi_star - ref.p2_p1) / ref.p2_p1
    covered = (
        not math.isnan(gci_band)
        and ext_err_pct <= gci_band
    )

    print(f"\n  Richardson phi_star  = {phi_star:.5f} ({src})")
    print(f"  Analytical p2/p1     = {ref.p2_p1:.5f}")
    print(f"  |phi_star - ref|     = {ext_err_pct:.4f} %   vs  GCI_fine = "
          f"{_fmt(gci_band)} %  ->  reference "
          f"{'COVERED' if covered else 'NOT covered'} by GCI band")
    print(f"  finest-grid error    = {err_vs_ref[-1]['rel_pct']:.4f} %")

    out = {
        "case": "wedge15Ma5",
        "qoi": {"name": QOI_NAME, "description": QOI_DESCRIPTION,
                "source": "surfaceFieldValue (areaAverage(p) on patch obstacle)"},
        "stationary_window": list(T_STAT),
        "reference": {
            "kind": "analytical_oblique_shock",
            "M_inflow": M_INFLOW, "theta_deg": THETA_DEG, "gamma": GAMMA,
            "beta_deg": ref.beta_deg, "p2_p1": ref.p2_p1,
            "rho2_rho1": ref.rho2_rho1, "T2_T1": ref.T2_T1, "M2": ref.M2,
        },
        "cases": cases,
        "triplet_A_CMF": tA,
        "triplet_B_MFXF": tB,
        "phi_star": phi_star,
        "phi_star_source": src,
        "phi_star_rel_err_pct": ext_err_pct,
        "reference_covered_by_gci": bool(covered),
        "error_table_vs_reference": err_vs_ref,
    }
    out_path = Path(__file__).parent / "gci_summary.json"
    out_path.write_text(json.dumps(_clean_json(out), indent=2, allow_nan=False))
    print(f"\nWrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
