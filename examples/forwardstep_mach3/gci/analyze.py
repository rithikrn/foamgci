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

# order_consistency + extremum_location_drift are library functions from
# foamgci v3.5+. Inline compact equivalents so this file is swap-and-run on an
# older (3.4.x) install too; the library versions are used when present.
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

    # Across-grid drift: a max that is steady within each window can still
    # march with the mesh (location scales with h). mesh_locked means its
    # value is not a continuum quantity, however clean the GCI looks.
    drift = report.get("across_grid_location_drift", {})
    if drift.get("mesh_locked"):
        flags.append(
            f"extremum location is mesh-locked across grids: the two finest "
            f"grids place it {drift.get('finest_cell_disp', float('nan')):.0f} "
            f"cells apart (> {drift.get('tol_cells', 2)}), so it tracks a mesh "
            f"feature rather than a fixed point; GCI on its value is not a "
            f"continuum-error estimate"
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
    oc = order_consistency(gcis)
    order_consistency_block = {
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

    # Across-grid extremum-location drift: does the max converge to a fixed
    # physical point, or march with the mesh? Complements the within-window
    # wander flag (which is steady on each grid yet blind to cross-grid drift).
    order_cf = ["Coarse", "Medium", "Fine", "Extra-fine"]
    drift_cases = [by[k] for k in order_cf
                   if isinstance(by[k]["loc_x"], float)
                   and not math.isnan(by[k]["loc_x"])]
    drift = extremum_location_drift(
        [(c["loc_x"], c["loc_y"]) for c in drift_cases],
        [c["dx"] for c in drift_cases],
        labels=[c["label"] for c in drift_cases],
    )
    across_grid_drift = {
        "verdict": drift.verdict,
        "mesh_locked": drift.mesh_locked,
        "finest_cell_disp": drift.finest_cell_disp,
        "cell_displacements": drift.cell_displacements,
        "pair_labels": drift.pair_labels,
        "tol_cells": drift.tol_cells,
        "shrinking": drift.shrinking,
        "note": drift.note,
    }

    report = {
        "key": qoi["key"],
        "description": qoi["description"],
        "field_used": cases[0]["field"],
        "quantity": qoi["quantity"],
        "primary": bool(qoi.get("primary", False)),
        "cases": cases,
        "triplet_A_CMF": triplet_A,
        "triplet_B_MFXF": triplet_B,
        "order_consistency": order_consistency_block,
        "across_grid_location_drift": across_grid_drift,
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


def _make_diagnostics_figure(out: dict) -> None:
    """Results figure for the forward step: primary p_max convergence vs the
    Rayleigh-Pitot reference, error vs h, and the extremum-location drift that
    shows whether max(p) settles onto a fixed point or marches with the mesh.
    Saved to figures/fig_diagnostics.{pdf,png}; skips if matplotlib is absent."""
    try:
        import matplotlib as mpl
        mpl.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:  # noqa: BLE001
        print(f"  (diagnostics figure skipped: {e})")
        return

    qr = out.get("qoi_results", {})
    p = qr.get("p_max")
    if not p:
        return
    ref = p.get("reference_value")
    pc = p["cases"]
    dx = np.array([c["dx"] for c in pc]); m = np.array([c["mean"] for c in pc])
    lab = [c["label"] for c in pc]
    phi_ext = p.get("phi_star")

    mpl.rcParams.update({"font.family": "serif", "font.size": 9,
                         "axes.linewidth": 0.6, "xtick.direction": "in",
                         "ytick.direction": "in"})
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.2))

    if ref is not None:
        ax[0].axhline(ref, ls="--", c="crimson", lw=1.2,
                      label=f"Rayleigh-Pitot {ref:.3f}")
    if phi_ext is not None:
        ax[0].axhline(phi_ext, ls=":", c="0.4", lw=1.1,
                      label=f"Richardson {phi_ext:.3f}")
    ax[0].plot(dx, m, "o-", c="k", ms=6, lw=1)
    for l, x, y in zip(lab, dx, m):
        ax[0].annotate(l, (x, y), fontsize=6.5, xytext=(0, 6),
                       textcoords="offset points", ha="center")
    ax[0].set_xscale("log"); ax[0].set_xlabel("grid spacing $h$")
    ax[0].set_ylabel(r"$\max(p)/p_\infty$")
    ax[0].set_title("(a) primary $p_{\\max}$ convergence", fontsize=8.5)
    ax[0].legend(fontsize=6.5, loc="best")

    xx = np.array([dx.min(), dx.max()])
    if ref is not None:
        err = np.abs(m - ref)
        ax[1].loglog(dx, err, "o-", c="k", ms=6, lw=1)
        ax[1].loglog(xx, err[-1]*(xx/dx[-1])**1, "--", c="0.5", lw=1,
                     label=r"$\mathcal{O}(h)$")
        ax[1].loglog(xx, err[-1]*(xx/dx[-1])**2, "-.", c="0.7", lw=1,
                     label=r"$\mathcal{O}(h^2)$")
        ax[1].set_ylabel(r"$|\max(p)-p_{02}|$")
        ax[1].legend(fontsize=6.5)
    ax[1].set_xlabel("grid spacing $h$")
    ax[1].set_title("(b) $p_{\\max}$ error vs reference", fontsize=8.5)

    sc = [c for c in pc if isinstance(c.get("loc_x"), (int, float))
          and not math.isnan(c["loc_x"])]
    if len(sc) >= 2:
        dx2 = np.array([c["dx"] for c in sc])
        # distance between consecutive locations, in cells of the finer grid
        xs = np.array([c["loc_x"] for c in sc]); ys = np.array([c["loc_y"] for c in sc])
        disp_cells = [math.hypot(xs[k+1]-xs[k], ys[k+1]-ys[k]) / dx2[k+1]
                      for k in range(len(sc)-1)]
        pairs = [f"{sc[k]['label'][:1]}->{sc[k+1]['label'][:1]}"
                 for k in range(len(sc)-1)]
        ax[2].axhline(2.0, ls="--", c="crimson", lw=1.0,
                      label="tol = 2 cells")
        ax[2].plot(range(len(disp_cells)), disp_cells, "s-", c="darkblue",
                   ms=6, lw=1)
        ax[2].set_xticks(range(len(pairs))); ax[2].set_xticklabels(pairs, fontsize=7)
        ax[2].set_yscale("log")
        ld = p.get("across_grid_location_drift", {})
        ax[2].set_title(f"(c) $p_{{\\max}}$ location drift\n[{ld.get('verdict','')}]",
                        fontsize=8.5)
        ax[2].set_ylabel("inter-grid displacement (cells)")
        ax[2].legend(fontsize=6.5)
    else:
        ax[2].text(0.5, 0.5, "location drift\nunavailable", ha="center",
                   va="center", transform=ax[2].transAxes)

    fig.tight_layout()
    figdir = Path(__file__).parent / "figures"
    figdir.mkdir(exist_ok=True)
    fig.savefig(figdir / "fig_diagnostics.pdf", bbox_inches="tight")
    fig.savefig(figdir / "fig_diagnostics.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {figdir / 'fig_diagnostics.pdf'}")


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

        ocb = report.get("order_consistency")
        if ocb is not None:
            if ocb["orders"]:
                pairs = ", ".join(f"{lbl}: p={p:.3f}" for lbl, p
                                  in zip(ocb["triplet_labels"], ocb["orders"]))
            else:
                pairs = "none"
            print(f"  Order consistency [{ocb['verdict']}] (necessary, not "
                  f"sufficient, for asymptotic range): {pairs}")

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
    _make_diagnostics_figure(out)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
