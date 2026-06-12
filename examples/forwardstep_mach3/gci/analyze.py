from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from data import GRIDS, GREENSHIELDS_2010, T_STAT  # noqa: E402

from foamgci.reader import read_fieldminmax            # noqa: E402
from foamgci.stats import window_stats                 # noqa: E402
from foamgci.gci import gci_over_hierarchy             # noqa: E402
from foamgci.report import rayleigh_pitot              # noqa: E402

U_FIELD_NAMES = ("mag(U)", "U")


def _mean_std(path: Path, field: str, window) -> tuple[float, float]:
    """Trapezoidal time-mean and sample std of a field's max over window."""
    d = read_fieldminmax(path, field=field)
    ws = window_stats(d.time, d.max, window[0], window[1])
    return ws.mean, ws.std


def summarise(grid, window) -> dict:
    path = grid.fieldminmax_path
    dp = read_fieldminmax(path, field="p")
    ws = window_stats(dp.time, dp.max, window[0], window[1])

    rho_mean, rho_std = _mean_std(path, "rho", window)
    u_mean = u_std = float("nan")
    for name in U_FIELD_NAMES:
        try:
            u_mean, u_std = _mean_std(path, name, window)
            break
        except (ValueError, RuntimeError):
            continue

    # Median peak-pressure location over the window (robust to ties).
    m = (dp.time >= window[0]) & (dp.time <= window[1])
    if dp.loc_max is not None and m.any():
        px, py = (float(v) for v in np.median(dp.loc_max[m], axis=0)[:2])
    else:
        px = py = float("nan")

    return dict(
        label=grid.label,
        n_cells=grid.n_cells, dx=grid.dx,
        n_samples_stat=ws.n,
        p_max_mean=ws.mean, p_max_std=ws.std,
        p_max_sem=ws.sem, p_max_tau_int=ws.tau_int,
        p_max_n_eff=ws.n_eff,
        kpss_stat=ws.kpss_stat, kpss_p=ws.kpss_p,
        kpss_stationary=bool(ws.kpss_stationary_5pct),
        rho_max_mean=rho_mean, rho_max_std=rho_std,
        U_mag_max_mean=u_mean, U_mag_max_std=u_std,
        peak_loc_x=px, peak_loc_y=py,
    )


def triplet_dict(g, label: str) -> dict:
    """Map a foamgci GCIResult onto the JSON schema the figures expect."""
    return dict(
        label=label,
        phi1=g.phi_fine, phi2=g.phi_medium, phi3=g.phi_coarse,
        h1=g.h_fine, h2=g.h_medium, h3=g.h_coarse,
        r=g.r21,
        eps21=g.phi_fine - g.phi_medium,
        eps32=g.phi_medium - g.phi_coarse,
        R=(g.phi_fine - g.phi_medium) / (g.phi_medium - g.phi_coarse),
        regime=g.regime,
        p_obs=g.p_apparent,
        phi_ext=g.phi_exact,
        gci_fine_pct=g.gci_fine_21_pct,
        gci_med_pct=g.gci_medium_32_pct,
        asymptotic_ratio=g.asymptotic_ratio,
        note=g.note,
    )


def main() -> int:
    print(f"Stationary window t in {T_STAT}\n")
    cases = []
    for grid in GRIDS:
        if not grid.fieldminmax_path.is_file():
            print(f"  ERROR: missing input file {grid.fieldminmax_path}")
            print("         Populate gci/data/ from your OpenFOAM runs first.")
            print("         See gci/data/README.md (one `cp` per grid).")
            return 2
        s = summarise(grid, T_STAT)
        cases.append(s)
        print(f"  {s['label']:11s} N={s['n_samples_stat']:4d}  "
              f"<p_max>={s['p_max_mean']:.5f}  sigma={s['p_max_std']:.4f}  "
              f"tau={s['p_max_tau_int']:.2f}  SEM={s['p_max_sem']:.5f}  "
              f"N_eff={s['p_max_n_eff']:.0f}  KPSS_p={s['kpss_p']:.2f}")

    by = {c["label"]: c for c in cases}
    order = ["Coarse", "Medium", "Fine", "Extra-fine"]
    means = [by[k]["p_max_mean"] for k in order]
    hs = [by[k]["dx"] for k in order]

    gcis = gci_over_hierarchy(means, hs, order)   # [(C,M,F), (M,F,XF)]
    tA = triplet_dict(gcis[0], "A : C-M-F")
    tB = triplet_dict(gcis[1], "B : M-F-XF")
    for t in (tA, tB):
        print(f"\nTriplet {t['label']}: regime={t['regime']}  "
              f"R={t['R']:+.4f}  p_obs={t['p_obs']:.4f}  "
              f"phi_ext={t['phi_ext']:.5f}  GCI_fine={t['gci_fine_pct']:.4f}%")
        print(f"  {t['note']}")

    # Pick phi_star from the deepest monotonic triplet, else extra-fine.
    if gcis[1].regime == "monotonic":
        phi_star, src = gcis[1].phi_exact, tB["label"]
    elif gcis[0].regime == "monotonic":
        phi_star, src = gcis[0].phi_exact, tA["label"]
    else:
        phi_star, src = by["Extra-fine"]["p_max_mean"], "Extra-fine (no monotonic triplet)"

    p02 = rayleigh_pitot(GREENSHIELDS_2010["M_inflow"], GREENSHIELDS_2010["gamma"])
    err_table = [
        dict(label=c["label"], dx=c["dx"], p_max=c["p_max_mean"],
             err=abs(c["p_max_mean"] - phi_star),
             rel_pct=100.0 * abs(c["p_max_mean"] - phi_star) / abs(phi_star))
        for c in cases
    ]

    out = dict(
        stationary_window=list(T_STAT),
        cases=cases,
        triplet_A_CMF=tA,
        triplet_B_MFXF=tB,
        phi_star=phi_star,
        phi_star_source=src,
        rayleigh_pitot_p02=p02,           # inviscid normal-shock lower-bound check
        error_table=err_table,
        greenshields_2010=GREENSHIELDS_2010,
    )
    out_path = Path(__file__).parent / "gci_summary.json"
    out_path.write_text(json.dumps(out, indent=2, default=str))
    print(f"\nphi_star = {phi_star:.5f} ({src});  "
          f"Rayleigh-Pitot p02 = {p02:.4f} (lower-bound check)")
    print(f"Wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
