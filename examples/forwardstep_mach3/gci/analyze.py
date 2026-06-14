from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any, Iterable

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

from data import GRIDS, GREENSHIELDS_2010, T_STAT  # noqa: E402
from foamgci.reader import read_fieldminmax  # noqa: E402
from foamgci.stats import window_stats  # noqa: E402
from foamgci.gci import gci_over_hierarchy  # noqa: E402
from foamgci.report import rayleigh_pitot  # noqa: E402


# Formal QoIs analyzed from the same fieldMinMax.dat files.
# Keep p_max first because existing figure scripts treat it as the primary QoI.
#
# U_mag_max is intentionally excluded from the formal table for this case:
# its extremum location changes substantially across grids and the GCI
# extrapolation is not physically useful. Velocity can still be added later
# as a diagnostic-only QoI if needed.
# Maximum in-window extremum wander (in cell widths) tolerated before a QoI
# is flagged as not pointwise-localized. A stagnation-point pressure peak
# stays put (~O(1) cell); a density max that trades between the triple point
# and the stagnation foot wanders over tens-to-hundreds of cells.
LOC_WANDER_CELLS_MAX = 5.0


QOIS = [
    {
        "key": "p_max",
        "field_candidates": ("p",),
        "quantity": "max",
        "description": "maximum pressure",
        "reference": "rayleigh_pitot",
        "primary": True,
    },
    {
        "key": "rho_max",
        "field_candidates": ("rho",),
        "quantity": "max",
        "description": "maximum density",
        "reference": None,
        "primary": False,
    },
]


def _clean_json(obj: Any) -> Any:
    """Convert NaN/Inf and NumPy scalars into strict JSON-safe values."""
    if isinstance(obj, dict):
        return {k: _clean_json(v) for k, v in obj.items()}

    if isinstance(obj, list):
        return [_clean_json(v) for v in obj]

    if isinstance(obj, tuple):
        return [_clean_json(v) for v in obj]

    if isinstance(obj, np.generic):
        return _clean_json(obj.item())

    if isinstance(obj, float):
        if math.isnan(obj) or math.isinf(obj):
            return None
        return obj

    return obj


def _fmt(value: Any, spec: str = ".4f") -> str:
    """Format floats cleanly, printing -- for undefined values."""
    if value is None:
        return "--"

    try:
        v = float(value)
    except (TypeError, ValueError):
        return str(value)

    if math.isnan(v) or math.isinf(v):
        return "--"

    return format(v, spec)


def _read_first_available(path: Path, field_candidates: Iterable[str]):
    """Read the first matching field name from a fieldMinMax.dat file."""
    last_error: Exception | None = None

    for field_name in field_candidates:
        try:
            return read_fieldminmax(path, field=field_name), field_name
        except ValueError as exc:
            last_error = exc

    names = ", ".join(repr(f) for f in field_candidates)
    raise ValueError(f"None of the requested fields were found: {names}") from last_error


def _qoi_case_stats(grid, qoi: dict, window) -> dict:
    """Return per-grid statistics for one QoI."""
    path = grid.fieldminmax_path
    data, actual_field = _read_first_available(path, qoi["field_candidates"])

    quantity = qoi["quantity"]
    if quantity == "max":
        series = data.max
        loc = data.loc_max
    elif quantity == "min":
        series = data.min
        loc = data.loc_min
    else:
        raise ValueError(f"Unsupported quantity: {quantity!r}")

    ws = window_stats(data.time, series, window[0], window[1])

    # Median extremum location over the statistical window.
    # Median is used because pointwise extrema can hop by one or two cells.
    m = (data.time >= window[0]) & (data.time <= window[1])
    if loc is not None and m.any():
        locw = loc[m]
        lx, ly = (float(v) for v in np.median(locw, axis=0)[:2])
        # Within-window extremum wander: the robust (5th-95th percentile)
        # spread of the extremum LOCATION inside the averaging window,
        # expressed in cell widths. A localized pointwise QoI (e.g. a
        # stagnation-point pressure peak) keeps its extremum within ~1 cell
        # over the whole window, so loc_wander_cells ~ O(1). A QoI whose
        # extremum migrates between competing flow features (e.g. a density
        # maximum that trades between the triple-point region and the
        # stagnation foot) wanders over tens-to-hundreds of cells, which
        # means its windowed time-average is not a single-feature scalar and
        # the GCI on it is comparing different physical quantities across
        # grids. This is independent of, and complementary to, the KPSS
        # stationarity test on the VALUE series.
        span_x = float(np.percentile(locw[:, 0], 95) - np.percentile(locw[:, 0], 5))
        span_y = float(np.percentile(locw[:, 1], 95) - np.percentile(locw[:, 1], 5))
        loc_wander_cells = float(np.hypot(span_x, span_y) / grid.dx)
    else:
        lx = ly = float("nan")
        loc_wander_cells = float("nan")

    return {
        "label": grid.label,
        "n_cells": grid.n_cells,
        "dx": grid.dx,
        "field": actual_field,
        "quantity": quantity,
        "description": qoi["description"],
        "n_samples_stat": ws.n,
        "mean": ws.mean,
        "std": ws.std,
        "sem": ws.sem,
        "tau_int": ws.tau_int,
        "n_eff": ws.n_eff,
        "kpss_stat": ws.kpss_stat,
        "kpss_p": ws.kpss_p,
        "kpss_stationary": bool(ws.kpss_stationary_5pct),
        "loc_x": lx,
        "loc_y": ly,
        "loc_wander_cells": loc_wander_cells,
    }


def _triplet_dict(g, label: str) -> dict:
    """Map a foamgci GCIResult onto a JSON-friendly schema."""
    eps21 = g.phi_fine - g.phi_medium
    eps32 = g.phi_medium - g.phi_coarse

    if eps32 == 0.0:
        ratio = float("nan")
    else:
        ratio = eps21 / eps32

    return {
        "label": label,
        "phi1": g.phi_fine,
        "phi2": g.phi_medium,
        "phi3": g.phi_coarse,
        "h1": g.h_fine,
        "h2": g.h_medium,
        "h3": g.h_coarse,
        "r": g.r21,
        "eps21": eps21,
        "eps32": eps32,
        "R": ratio,
        "regime": g.regime,
        "p_obs": g.p_apparent,
        "phi_ext": g.phi_exact,
        "gci_fine_pct": g.gci_fine_21_pct,
        "gci_med_pct": g.gci_medium_32_pct,
        "asymptotic_ratio": g.asymptotic_ratio,
        "note": g.note,
    }


def _qoi_quality_flags(report: dict) -> list[str]:
    """Return human-readable caveats for one QoI report."""
    flags: list[str] = []

    cases = report["cases"]
    if not all(c["kpss_stationary"] for c in cases):
        failed = [c["label"] for c in cases if not c["kpss_stationary"]]
        flags.append("KPSS stationarity rejected on: " + ", ".join(failed))

    # Extremum-location localization. If the extremum wanders over many
    # cells within the averaging window, the windowed mean is not a single
    # physical feature and triplet-wise GCI on it compares different maxima
    # across grids. Threshold is deliberately generous: a localized pointwise
    # QoI stays within a handful of cells.
    wanders = [
        c["loc_wander_cells"]
        for c in cases
        if isinstance(c["loc_wander_cells"], float)
        and not math.isnan(c["loc_wander_cells"])
    ]
    if wanders:
        med_wander = float(np.median(wanders))
        report["loc_wander_cells_median"] = med_wander
        if med_wander > LOC_WANDER_CELLS_MAX:
            flags.append(
                f"extremum is not localized: median in-window wander "
                f"{med_wander:.0f} cells (> {LOC_WANDER_CELLS_MAX:.0f}); the "
                f"max relocates between distinct flow features, so its "
                f"time-average is not a single-feature scalar"
            )

    tB = report["triplet_B_MFXF"]
    gci = tB.get("gci_fine_pct")
    p_obs = tB.get("p_obs")

    if gci is not None and not (isinstance(gci, float) and math.isnan(gci)):
        if gci > 10.0:
            flags.append("deepest-triplet GCI exceeds 10%")

    if p_obs is not None and not (isinstance(p_obs, float) and math.isnan(p_obs)):
        if p_obs <= 0.1:
            flags.append("observed order is at or below the solver lower bound")

    return flags


def _qoi_report(qoi: dict, cases: list[dict]) -> dict:
    """Compute GCI block and reference comparison for one QoI."""
    order = ["Coarse", "Medium", "Fine", "Extra-fine"]
    by = {c["label"]: c for c in cases}

    missing = [label for label in order if label not in by]
    if missing:
        raise ValueError(f"Missing grid labels in QoI cases: {missing}")

    means = [by[k]["mean"] for k in order]
    hs = [by[k]["dx"] for k in order]

    gcis = gci_over_hierarchy(means, hs, order)
    triplet_A = _triplet_dict(gcis[0], "A : C-M-F")
    triplet_B = _triplet_dict(gcis[1], "B : M-F-XF")

    if gcis[1].regime == "monotonic":
        phi_star, src = gcis[1].phi_exact, triplet_B["label"]
    elif gcis[0].regime == "monotonic":
        phi_star, src = gcis[0].phi_exact, triplet_A["label"]
    else:
        phi_star, src = by["Extra-fine"]["mean"], "Extra-fine (no monotonic triplet)"

    error_table_vs_phi_star = [
        {
            "label": c["label"],
            "dx": c["dx"],
            "mean": c["mean"],
            "err": abs(c["mean"] - phi_star),
            "rel_pct": 100.0 * abs(c["mean"] - phi_star) / abs(phi_star),
        }
        for c in cases
    ]

    reference_value = None
    reference_label = None
    error_table_vs_reference = None

    if qoi.get("reference") == "rayleigh_pitot":
        reference_value = rayleigh_pitot(
            GREENSHIELDS_2010["M_inflow"],
            GREENSHIELDS_2010["gamma"],
        )
        reference_label = "Rayleigh-Pitot p02"
        error_table_vs_reference = [
            {
                "label": c["label"],
                "dx": c["dx"],
                "mean": c["mean"],
                "err": abs(c["mean"] - reference_value),
                "rel_pct": 100.0 * abs(c["mean"] - reference_value) / abs(reference_value),
            }
            for c in cases
        ]

    report = {
        "key": qoi["key"],
        "description": qoi["description"],
        "field_used": cases[0]["field"],
        "quantity": qoi["quantity"],
        "primary": bool(qoi.get("primary", False)),
        "cases": cases,
        "triplet_A_CMF": triplet_A,
        "triplet_B_MFXF": triplet_B,
        "phi_star": phi_star,
        "phi_star_source": src,
        "error_table_vs_phi_star": error_table_vs_phi_star,
        "reference_label": reference_label,
        "reference_value": reference_value,
        "error_table_vs_reference": error_table_vs_reference,
    }

    flags = _qoi_quality_flags(report)
    report["quality_flags"] = flags
    report["diagnostic_only"] = bool(flags and not report["primary"])

    return report


def _legacy_cases(qoi_results: dict[str, dict]) -> list[dict]:
    """Build backward-compatible top-level cases for existing figure scripts.

    The old schema is pressure-centered but also included rho summary columns.
    U_mag columns were intentionally removed after U_mag_max was excluded from
    the formal forward-step QoI set.
    """
    p_cases = qoi_results["p_max"]["cases"]
    rho_cases = qoi_results["rho_max"]["cases"]

    old_cases = []
    for cp, crho in zip(p_cases, rho_cases):
        old_cases.append(
            {
                "label": cp["label"],
                "n_cells": cp["n_cells"],
                "dx": cp["dx"],
                "n_samples_stat": cp["n_samples_stat"],
                "p_max_mean": cp["mean"],
                "p_max_std": cp["std"],
                "p_max_sem": cp["sem"],
                "p_max_tau_int": cp["tau_int"],
                "p_max_n_eff": cp["n_eff"],
                "kpss_stat": cp["kpss_stat"],
                "kpss_p": cp["kpss_p"],
                "kpss_stationary": cp["kpss_stationary"],
                "rho_max_mean": crho["mean"],
                "rho_max_std": crho["std"],
                "rho_max_sem": crho["sem"],
                "rho_max_tau_int": crho["tau_int"],
                "rho_max_n_eff": crho["n_eff"],
                "rho_kpss_stat": crho["kpss_stat"],
                "rho_kpss_p": crho["kpss_p"],
                "rho_kpss_stationary": crho["kpss_stationary"],
                "peak_loc_x": cp["loc_x"],
                "peak_loc_y": cp["loc_y"],
                "p_max_loc_wander_cells": cp["loc_wander_cells"],
                "rho_max_loc_wander_cells": crho["loc_wander_cells"],
            }
        )

    return old_cases


def main() -> int:
    print(f"Stationary window t in {T_STAT}\n")

    for grid in GRIDS:
        if not grid.fieldminmax_path.is_file():
            print(f"ERROR: missing input file {grid.fieldminmax_path}")
            print("Populate gci/data/ from your OpenFOAM runs first.")
            print("See gci/data/README.md (one `cp` per grid).")
            return 2

    qoi_results: dict[str, dict] = {}

    for qoi in QOIS:
        print(f"\nQoI: {qoi['key']} ({qoi['description']})")
        qoi_cases = []

        for grid in GRIDS:
            s = _qoi_case_stats(grid, qoi, T_STAT)
            qoi_cases.append(s)

            print(
                f"  {s['label']:11s} "
                f"N={s['n_samples_stat']:4d} "
                f"mean={s['mean']:.6f} "
                f"sigma={s['std']:.4g} "
                f"tau={s['tau_int']:.2f} "
                f"SEM={s['sem']:.4g} "
                f"N_eff={s['n_eff']:.1f} "
                f"KPSS_p={s['kpss_p']:.2f} "
                f"stationary={s['kpss_stationary']}"
            )

        report = _qoi_report(qoi, qoi_cases)
        qoi_results[qoi["key"]] = report

        tA = report["triplet_A_CMF"]
        tB = report["triplet_B_MFXF"]

        for t in (tA, tB):
            print(
                f"  Triplet {t['label']}: "
                f"regime={t['regime']} "
                f"R={_fmt(t['R'], '+.4f')} "
                f"p_obs={_fmt(t['p_obs'], '.4f')} "
                f"phi_ext={_fmt(t['phi_ext'], '.5f')} "
                f"GCI_fine={_fmt(t['gci_fine_pct'], '.4f')}%"
            )
            print(f"    {t['note']}")

        if report["quality_flags"]:
            print("  Caveats:")
            for flag in report["quality_flags"]:
                print(f"    - {flag}")

    pmax = qoi_results["p_max"]

    out = {
        "stationary_window": list(T_STAT),
        "cases": _legacy_cases(qoi_results),
        "triplet_A_CMF": pmax["triplet_A_CMF"],
        "triplet_B_MFXF": pmax["triplet_B_MFXF"],
        "phi_star": pmax["phi_star"],
        "phi_star_source": pmax["phi_star_source"],
        "rayleigh_pitot_p02": pmax["reference_value"],
        "error_table": pmax["error_table_vs_phi_star"],
        "error_table_vs_rayleigh_pitot": pmax["error_table_vs_reference"],
        "qoi_results": qoi_results,
        "included_qois": [q["key"] for q in QOIS],
        "excluded_qois": {
            "U_mag_max": (
                "Excluded from the formal forward-step GCI table because the "
                "pointwise velocity maximum changes location across grids and "
                "did not produce a physically useful extrapolation."
            )
        },
        "greenshields_2010": GREENSHIELDS_2010,
    }

    out_path = Path(__file__).parent / "gci_summary.json"
    out_path.write_text(json.dumps(_clean_json(out), indent=2, allow_nan=False))

    print(f"\nPrimary p_max phi_star = {pmax['phi_star']:.5f} ({pmax['phi_star_source']})")
    print(f"Rayleigh-Pitot p02 = {pmax['reference_value']:.5f}")
    print(f"Included QoIs: {', '.join(out['included_qois'])}")
    print(f"Wrote {out_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
