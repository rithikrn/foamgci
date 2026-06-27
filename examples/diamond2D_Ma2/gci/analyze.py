from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from types import SimpleNamespace as _NS
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

from data import (  # noqa: E402
    GRIDS, T_STAT, M_INFLOW, EPS_DEG, GAMMA, P_INF, T_OVER_C, CD_HALF_TO_FULL,
)
from shock_expansion import diamond_reference  # noqa: E402
from foamgci.reader import read_surface_field_value, read_fieldminmax  # noqa: E402
from foamgci.stats import window_stats  # noqa: E402
from foamgci.gci import gci_over_hierarchy  # noqa: E402

# order_consistency, extremum_location_drift, least_squares_gci are library
# functions in newer foamgci; inline compact equivalents so this driver runs
# on the current install too. The library versions are used when present.
try:
    from foamgci.gci import order_consistency  # noqa
except ImportError:
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
                        f"{100*rel:.1f}% vs tol {100*rel_tol:.0f}% (necessary, not sufficient)")

try:
    from foamgci.gci import extremum_location_drift  # noqa
except ImportError:
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

try:
    from foamgci.gci import least_squares_gci  # noqa
except ImportError:
    def least_squares_gci(phis, hs, p_band=(0.5, 2.1)):
        """Eca-Hoekstra least-squares GCI over all grids (>=4).

        Fit phi_i = phi0 + alpha*h_i^p; for each fixed p this is linear in
        (phi0, alpha), so grid-search p and keep the fit with the smallest
        residual standard deviation. Safety factor 1.25 in [0.5, 2.1], else
        3.0 (Eca & Hoekstra, JCP 2014). Returns phi_ext, p, sigma, u_pct.
        """
        phi = np.asarray(phis, float); h = np.asarray(hs, float)
        n = phi.size
        if n < 4:
            return _NS(phi_ext=float("nan"), p=float("nan"), sigma=float("nan"),
                       fs=float("nan"), u=float("nan"), u_pct=float("nan"),
                       well_behaved=None, reliable=None, n=n,
                       note="need >=4 grids for a least-squares fit")
        p_grid = np.unique(np.concatenate([np.arange(0.5, 3.001, 0.05), [1.0, 2.0]]))
        best = None
        for p in p_grid:
            A = np.vstack([np.ones(n), h**p]).T
            coef, *_ = np.linalg.lstsq(A, phi, rcond=None)
            resid = phi - A @ coef
            sigma = float(np.sqrt(np.sum(resid**2) / (n - 2)))
            if best is None or sigma < best[0]:
                best = (sigma, float(p), float(coef[0]), float(coef[1]))
        sigma, p, phi0, alpha = best
        well = bool(p_band[0] <= p <= p_band[1])
        fs = 1.25 if well else 3.0
        u = fs * sigma
        u_pct = float(100.0 * u / abs(phi0)) if phi0 else float("nan")
        reliable = bool(well and sigma <= 0.02 * abs(phi0))
        return _NS(phi_ext=phi0, p=p, sigma=sigma, alpha=alpha, fs=fs, u=u,
                   u_pct=u_pct, well_behaved=well, reliable=reliable, n=n,
                   note=(f"p={p:.2f}, Fs={fs}, sigma={sigma:.3e}, "
                         f"u={u_pct:.3f}% of phi_ext "
                         f"({'reliable' if reliable else 'treat as indicative'})"))


# ----------------------------------------------------------------------------
# QoI hierarchy (knocks off the single-scalar, surface-only, and smoothness
# limitations in LIMITATIONS.md):
#   Cd        - drag coefficient from forceCoeffs (surface FORCE integral),
#               doubled for the half-domain. Smooth; anchored to shock-expansion.
#   p_front   - area-averaged front-facet pressure (surface average). Anchored.
#   p_rear    - area-averaged rear-facet pressure (surface average). Anchored.
#   Cd_press  - COUPLED cross-check: (p_front - p_rear)*(t/c)/(0.5*g*M^2). Must
#               agree with the independent forceCoeffs Cd, which validates the
#               force normalisation and the pressure integration together.
#   S_vol     - domain entropy integral (VOLUME integral). Smooth; converges to
#               a grid-independent value; physically tied to wave drag.
#   p_max     - global max(p) (pointwise extremum). DIAGNOSTIC: expected to be
#               mesh-locked at the leading-edge tip, the non-smooth contrast.
# ----------------------------------------------------------------------------
LOC_WANDER_CELLS_MAX = 5.0
ORDER = [g.label for g in GRIDS]
TRIPLETS = [("Coarse", "Medium", "Fine"),
            ("Medium", "Fine", "Extra-fine"),
            ("Fine", "Extra-fine", "Ultra-fine")]


def _clean_json(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {k: _clean_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean_json(v) for v in obj]
    if isinstance(obj, np.generic):
        return _clean_json(obj.item())
    if isinstance(obj, float):
        return None if (math.isnan(obj) or math.isinf(obj)) else obj
    return obj


def read_force_coeffs(path: Path, column: str = "Cd"):
    """Read an OpenFOAM forceCoeffs .dat; return (time, values) for one column.

    The 4.x header line is ``# Time  Cm  Cd  Cl  Cl(f)  Cl(r)``; we locate the
    requested column by name. (Inline here; a library reader belongs in
    foamgci.reader alongside read_surface_field_value.)
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"No such file: {path}")
    header = None
    rows = []
    with path.open() as fh:
        for raw in fh:
            s = raw.strip()
            if not s:
                continue
            if s.startswith("#"):
                body = s.lstrip("#").strip()
                if body.split()[:1] == ["Time"]:
                    header = body.split()
                continue
            rows.append([float(x) for x in s.split()])
    if header is None:
        raise ValueError(f"no '# Time ...' header in {path}")
    if column not in header:
        raise ValueError(f"column {column!r} not in {header}")
    j = header.index(column)
    arr = np.asarray(rows, float)
    return arr[:, 0], arr[:, j]


def _scalar_stats(time, value, label, h, n_cells) -> dict:
    ws = window_stats(time, value, T_STAT[0], T_STAT[1])
    return dict(label=label, h=float(h), n_cells=int(n_cells),
                mean=float(ws.mean), sem=float(ws.sem), std=float(ws.std),
                tau=float(ws.tau_int), n_window=int(ws.n),
                kpss_stationary=bool(ws.kpss_stationary_5pct))


def _gci_triplets(cases) -> tuple[list, list[dict]]:
    by = {c["label"]: c for c in cases}
    means = [by[k]["mean"] for k in ORDER]
    hs = [by[k]["h"] for k in ORDER]
    gcis = gci_over_hierarchy(means, hs, ORDER)
    out = []
    for g in gcis:
        out.append(dict(
            triplet=f'{g.label_coarse}-{g.label_medium}-{g.label_fine}',
            regime=g.regime,
            p_obs=None if g.regime != "monotonic" else float(g.p_apparent),
            phi_ext=None if g.regime != "monotonic" else float(g.phi_exact),
            gci_fine_pct=None if g.regime != "monotonic" else float(g.gci_fine_21_pct),
            R_fit=g.asymptotic_ratio,
        ))
    return gcis, out


def _lsq_block(cases) -> dict:
    phis = [c["mean"] for c in cases]
    hs = [c["h"] for c in cases]
    r = least_squares_gci(phis, hs)
    return dict(phi_ext=r.phi_ext, p=r.p, sigma=r.sigma, fs=r.fs,
                u_pct=r.u_pct, well_behaved=r.well_behaved,
                reliable=r.reliable, n_grids=r.n, note=r.note)


def _deepest_monotonic(gcis):
    mono = [g for g in gcis if g.regime == "monotonic"]
    return mono[-1] if mono else None


def _err_pct(value, ref):
    return None if (value is None or ref == 0) else float(100.0 * abs(value - ref) / abs(ref))


def main() -> int:
    ref = diamond_reference(M_INFLOW, EPS_DEG, GAMMA)

    # --- read every QoI on every grid -----------------------------------
    missing = []
    cd_cases, pf_cases, pr_cases, sv_cases, pmax_cases = [], [], [], [], []
    pmax_loc = {}
    for g in GRIDS:
        # drag coefficient (force integral), doubled for the half-domain
        try:
            t, cd = read_force_coeffs(g.force_path, "Cd")
            cd = CD_HALF_TO_FULL * cd
            cd_cases.append(_scalar_stats(t, cd, g.label, g.h, g.n_cells))
        except FileNotFoundError:
            missing.append(str(g.force_path))
        # facet pressures (surface averages)
        for path, bucket in ((g.pfront_path, pf_cases), (g.prear_path, pr_cases)):
            try:
                d = read_surface_field_value(path, column="p")
                bucket.append(_scalar_stats(d.time, d.value, g.label, g.h, g.n_cells))
            except FileNotFoundError:
                missing.append(str(path))
        # entropy volume integral (volFieldValue volIntegrate, same .dat format)
        try:
            d = read_surface_field_value(g.svol_path, column="deltaS")
            sv_cases.append(_scalar_stats(d.time, d.value, g.label, g.h, g.n_cells))
        except FileNotFoundError:
            pass  # optional QoI (coded FO); case still complete without it
        except Exception:
            pass
        # pointwise max(p) + location (diagnostic)
        try:
            d = read_fieldminmax(g.pmax_path, field="p")
            ws = window_stats(d.time, d.max, T_STAT[0], T_STAT[1])
            c = dict(label=g.label, h=float(g.h), n_cells=int(g.n_cells),
                     mean=float(ws.mean), sem=float(ws.sem), tau=float(ws.tau_int),
                     kpss_stationary=bool(ws.kpss_stationary_5pct))
            m = (d.time >= T_STAT[0]) & (d.time <= T_STAT[1])
            if d.loc_max is not None and np.any(m):
                loc = d.loc_max[m]
                c["loc_mean"] = [float(np.mean(loc[:, k])) for k in range(loc.shape[1])]
                xspread = float(np.percentile(loc[:, 0], 95) - np.percentile(loc[:, 0], 5))
                c["loc_wander_cells"] = xspread / g.h
                pmax_loc[g.label] = np.array(c["loc_mean"])
            pmax_cases.append(c)
        except FileNotFoundError:
            missing.append(str(g.pmax_path))

    if len(cd_cases) < 3 and len(pf_cases) < 3:
        print("ERROR: need at least three grids of force/pressure data.")
        for m in missing:
            print("  missing:", m)
        return 1

    out: dict[str, Any] = {
        "case": "diamond2D_Ma2",
        "M_inflow": M_INFLOW, "eps_deg": EPS_DEG, "gamma": GAMMA,
        "window": list(T_STAT),
        "reference": {
            "beta_deg": ref.beta_deg, "p_front": ref.p_front, "p_rear": ref.p_rear,
            "Cd_shock_expansion": ref.Cd_shock_expansion, "Cd_linear": ref.Cd_linear,
            "M2": ref.M2, "M3": ref.M3, "t_over_c": ref.t_over_c,
        },
        "smooth_qois": {}, "coupling": {}, "volume_qoi": {}, "diagnostic_pmax": {},
    }

    # --- smooth, reference-anchored QoIs --------------------------------
    smooth = [("Cd", cd_cases, ref.Cd_shock_expansion),
              ("p_front", pf_cases, ref.p_front),
              ("p_rear", pr_cases, ref.p_rear)]
    for name, cases, refval in smooth:
        if len(cases) < 3:
            continue
        gcis, tri = _gci_triplets(cases)
        deep = _deepest_monotonic(gcis)
        oc = order_consistency(gcis)
        lsq = _lsq_block(cases) if len(cases) >= 4 else None
        fine = cases[-1]
        out["smooth_qois"][name] = dict(
            reference=refval,
            cases=[dict(label=c["label"], h=c["h"], n_cells=c["n_cells"],
                        mean=c["mean"], sem=c["sem"], tau=c["tau"],
                        kpss_stationary=c["kpss_stationary"]) for c in cases],
            triplets=tri,
            deepest=None if deep is None else dict(
                triplet=f'{deep.label_coarse}-{deep.label_medium}-{deep.label_fine}',
                p_obs=float(deep.p_apparent), phi_ext=float(deep.phi_exact),
                gci_fine_pct=float(deep.gci_fine_21_pct)),
            order_consistency=dict(verdict=oc.verdict, note=oc.note),
            least_squares=lsq,
            finest_mean=fine["mean"],
            finest_err_pct=_err_pct(fine["mean"], refval),
            phi_ext_err_pct=(None if deep is None
                             else _err_pct(float(deep.phi_exact), refval)),
        )

    # --- coupled QoI: Cd from forces vs Cd from facet pressures ----------
    if len(pf_cases) >= 3 and len(pr_cases) >= 3:
        pf_by = {c["label"]: c for c in pf_cases}
        pr_by = {c["label"]: c for c in pr_cases}
        cd_by = {c["label"]: c for c in cd_cases}
        q = 0.5 * GAMMA * M_INFLOW ** 2
        rows = []
        for lab in ORDER:
            if lab in pf_by and lab in pr_by:
                cd_press = (pf_by[lab]["mean"] - pr_by[lab]["mean"]) * T_OVER_C / q
                cd_force = cd_by[lab]["mean"] if lab in cd_by else None
                rel = (None if cd_force in (None, 0)
                       else float(100.0 * abs(cd_force - cd_press) / abs(cd_press)))
                rows.append(dict(label=lab, cd_pressure=float(cd_press),
                                 cd_force=cd_force, rel_diff_pct=rel))
        out["coupling"] = dict(
            definition="Cd_pressure = (p_front - p_rear)*(t/c)/(0.5*gamma*M^2); "
                       "Cd_force = 2 * forceCoeffs Cd (half-domain). The two "
                       "independent routes to drag should agree.",
            per_grid=rows,
            finest_rel_diff_pct=rows[-1]["rel_diff_pct"] if rows else None,
        )

    # --- volume QoI: entropy integral (convergence only; no closed form) -
    if len(sv_cases) >= 3:
        gcis, tri = _gci_triplets(sv_cases)
        deep = _deepest_monotonic(gcis)
        lsq = _lsq_block(sv_cases) if len(sv_cases) >= 4 else None
        out["volume_qoi"]["S_vol"] = dict(
            description="domain integral of specific-entropy rise deltaS "
                        "(volFieldValue volIntegrate); smooth volume functional, "
                        "physically tied to wave drag via Oswatitsch.",
            cases=[dict(label=c["label"], h=c["h"], mean=c["mean"],
                        sem=c["sem"], tau=c["tau"]) for c in sv_cases],
            triplets=tri, least_squares=lsq,
            deepest=None if deep is None else dict(
                p_obs=float(deep.p_apparent), phi_ext=float(deep.phi_exact),
                gci_fine_pct=float(deep.gci_fine_21_pct)),
        )
    else:
        out["volume_qoi"]["S_vol"] = {
            "status": "absent",
            "note": "no sVol .dat found; the coded entropy FO may not have run. "
                    "The other QoIs are unaffected (see README fallback)."}

    # --- diagnostic: pointwise max(p), expected mesh-locked at the tip ---
    if len(pmax_cases) >= 3:
        gcis, tri = _gci_triplets(pmax_cases)
        wanders = [c.get("loc_wander_cells") for c in pmax_cases
                   if c.get("loc_wander_cells") is not None]
        med_wander = float(np.median(wanders)) if wanders else float("nan")
        localized = (not math.isnan(med_wander)) and (med_wander <= LOC_WANDER_CELLS_MAX)
        drift = None
        if len(pmax_loc) >= 2:
            labs = [l for l in ORDER if l in pmax_loc]
            drift = extremum_location_drift([pmax_loc[l] for l in labs],
                                            [g.h for g in GRIDS if g.label in labs],
                                            labels=labs)
        fine = pmax_cases[-1]
        out["diagnostic_pmax"] = dict(
            reference_p2=ref.p_front,
            cases=[dict(label=c["label"], h=c["h"], mean=c["mean"],
                        loc_wander_cells=c.get("loc_wander_cells")) for c in pmax_cases],
            triplets=tri,
            finest_value=fine["mean"],
            pct_above_post_shock=float(100.0 * (fine["mean"] - ref.p_front) / ref.p_front),
            within_window_localized=localized,
            across_grid=None if drift is None else dict(
                verdict=drift.verdict, finest_cell_disp=drift.finest_cell_disp,
                note=drift.note),
            lesson="A clean monotone GCI on a pointwise extremum can be junk; "
                   "the analytical reference and across-grid location drift expose it.",
        )

    out_path = Path(__file__).parent / "gci_summary.json"
    out_path.write_text(json.dumps(_clean_json(out), indent=2, allow_nan=False))

    # --- console summary -------------------------------------------------
    print(f"diamond2D_Ma2  M={M_INFLOW}  eps={EPS_DEG} deg")
    print(f"  reference: beta={ref.beta_deg:.2f} deg, Cd={ref.Cd_shock_expansion:.5f}, "
          f"p_front={ref.p_front:.4f}, p_rear={ref.p_rear:.4f}")
    for name, blk in out["smooth_qois"].items():
        d = blk.get("deepest")
        fe = blk.get("finest_err_pct")
        if d:
            print(f"  {name:8s}: finest={blk['finest_mean']:.5f} "
                  f"(err {fe:.3f}% vs ref), p_obs={d['p_obs']:.3f}, "
                  f"GCI_fine={d['gci_fine_pct']:.3f}%")
        ls = blk.get("least_squares")
        if ls:
            print(f"           least-squares (5 grids): {ls['note']}")
    if out["coupling"]:
        print(f"  coupling Cd_force vs Cd_pressure (finest): "
              f"{out['coupling']['finest_rel_diff_pct']}% apart")
    dp = out["diagnostic_pmax"]
    if dp:
        ag = dp.get("across_grid")
        print(f"  max(p) diagnostic: {dp['pct_above_post_shock']:.1f}% above post-shock; "
              f"{'' if ag is None else ag['verdict']}")
    print(f"  wrote {out_path}")
    if missing:
        print("  NOTE missing inputs:")
        for m in missing:
            print("   ", m)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
