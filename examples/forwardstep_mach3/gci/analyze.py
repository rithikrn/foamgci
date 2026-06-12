"""
analyze.py
----------
End-to-end GCI analysis on the four-grid forward-step study.

For each grid specified in data.py:
  1. parses postProcessing/fieldMinMax1/0/fieldMinMax.dat
  2. time-averages p_max, rho_max, |U|_max over the stationary window
  3. records the late-time location of the peak-pressure cell

Then computes Roache GCI on TWO overlapping triplets:
  A : {Coarse, Medium, Fine}
  B : {Medium, Fine, Extra-fine}
The shift in observed order p between A and B indicates whether the
asymptotic range has been reached (Celik et al. 2008 sec. 2.3).

Writes:
  gci_summary.json  -- structured results consumed by figure scripts.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from data import GRIDS, COARSE, MEDIUM, FINE, EXTRA_FINE, \
                 GREENSHIELDS_2010, T_STAT
from parse_fieldminmax import summarise_case, StationaryStats
from gci import gci_triplet


def banner(title: str) -> None:
    print("=" * 72)
    print(title)
    print("=" * 72)


def main() -> int:
    # ------------------------------------------------------------------
    # Step 1: parse each case's fieldMinMax.dat
    # ------------------------------------------------------------------
    banner(f"Parsing fieldMinMax over stationary window t in {T_STAT}")
    cases: list[StationaryStats] = []
    for g in GRIDS:
        path = g.fieldminmax_path
        print(f"  {g.label:11s} <- {path}")
        if not path.is_file():
            print(f"    ERROR: file missing; check {g.case_dir}")
            return 2
        stats = summarise_case(
            case_label=g.label, folder=g.folder,
            n_cells=g.n_cells, dx=g.dx,
            fieldminmax_path=path,
            window=T_STAT,
        )
        if stats.n_samples_stat < 10:
            print(f"    WARN: only {stats.n_samples_stat} samples in "
                  f"{T_STAT} -- run may not have reached stationary "
                  "regime; numbers below are NOT trustworthy.")
        cases.append(stats)

    # ------------------------------------------------------------------
    # Step 2: per-grid table
    # ------------------------------------------------------------------
    print()
    banner("Per-grid stationary statistics")
    print(f"{'Grid':12s}  {'cells':>8s}  {'dx':>9s}  "
          f"{'<p_max>':>9s}  {'sigma_p':>8s}  "
          f"{'<rho_max>':>10s}  {'<|U|_max>':>10s}")
    for s in cases:
        print(f"{s.label:12s}  {s.n_cells:>8d}  {s.dx:>9.5f}  "
              f"{s.p_max_mean:>9.4f}  {s.p_max_std:>8.4f}  "
              f"{s.rho_max_mean:>10.4f}  {s.U_mag_max_mean:>10.4f}")

    print()
    print(f"{'Grid':12s}  {'peak (x, y)':>18s}")
    for s in cases:
        print(f"{s.label:12s}  ({s.peak_loc_x:.5f}, {s.peak_loc_y:.5f})")

    # Convenient handles by label.
    by = {s.label: s for s in cases}

    # ------------------------------------------------------------------
    # Step 3: Roache GCI on both triplets
    # ------------------------------------------------------------------
    print()
    banner("Roache GCI -- primary metric <p_max>_stat")

    tA = gci_triplet(
        by["Fine"].p_max_mean,
        by["Medium"].p_max_mean,
        by["Coarse"].p_max_mean,
        by["Fine"].dx, by["Medium"].dx, by["Coarse"].dx,
        label="A : C-M-F",
    )
    tB = gci_triplet(
        by["Extra-fine"].p_max_mean,
        by["Fine"].p_max_mean,
        by["Medium"].p_max_mean,
        by["Extra-fine"].dx, by["Fine"].dx, by["Medium"].dx,
        label="B : M-F-XF",
    )

    for t in (tA, tB):
        print(f"\nTriplet {t.label}  ({t.regime})")
        print(f"  phi_1={t.phi1:.5f}  phi_2={t.phi2:.5f}  "
              f"phi_3={t.phi3:.5f}")
        print(f"  eps21={t.eps21:+.6f}  eps32={t.eps32:+.6f}  "
              f"R={t.R:+.4f}")
        if t.regime == "monotonic":
            print(f"  p_obs={t.p_obs:.3f}  phi_ext={t.phi_ext:.5f}")
            print(f"  GCI_12={t.gci_fine_pct:.4f}%  "
                  f"GCI_23={t.gci_med_pct:.4f}%  AR={t.asymptotic_ratio:.4f}")
        print(f"  {t.note}")

    # ------------------------------------------------------------------
    # Step 4: pick best phi_ext for error reporting
    # ------------------------------------------------------------------
    if tB.regime == "monotonic":
        phi_star = tB.phi_ext
        phi_star_source = tB.label
    elif tA.regime == "monotonic":
        phi_star = tA.phi_ext
        phi_star_source = tA.label
    else:
        # Fall back to extra-fine value.
        phi_star = by["Extra-fine"].p_max_mean
        phi_star_source = "Extra-fine (no monotonic triplet available)"

    print()
    banner(f"Relative error vs phi_ext = {phi_star:.5f}  "
           f"(source: {phi_star_source})")
    print(f"{'Grid':12s}  {'dx':>9s}  {'<p_max>':>9s}  "
          f"{'|f-f*|':>10s}  {'rel %':>8s}")
    err_table = []
    for s in cases:
        err = abs(s.p_max_mean - phi_star)
        rel = 100.0 * err / abs(phi_star) if phi_star else float("nan")
        print(f"{s.label:12s}  {s.dx:>9.5f}  {s.p_max_mean:>9.4f}  "
              f"{err:>10.5f}  {rel:>8.4f}")
        err_table.append(dict(label=s.label, dx=s.dx,
                              p_max=s.p_max_mean, err=err, rel_pct=rel))

    # ------------------------------------------------------------------
    # Step 5: dump JSON for downstream figure scripts
    # ------------------------------------------------------------------
    out = {
        "stationary_window": list(T_STAT),
        "cases": [
            dict(
                label=s.label, folder=s.folder,
                n_cells=s.n_cells, dx=s.dx,
                n_samples_stat=s.n_samples_stat,
                p_max_mean=s.p_max_mean, p_max_std=s.p_max_std,
                rho_max_mean=s.rho_max_mean, rho_max_std=s.rho_max_std,
                U_mag_max_mean=s.U_mag_max_mean,
                U_mag_max_std=s.U_mag_max_std,
                peak_loc_x=s.peak_loc_x, peak_loc_y=s.peak_loc_y,
            )
            for s in cases
        ],
        "triplet_A_CMF": tA.as_dict(),
        "triplet_B_MFXF": tB.as_dict(),
        "phi_star": phi_star,
        "phi_star_source": phi_star_source,
        "error_table": err_table,
        "greenshields_2010": GREENSHIELDS_2010,
        "mesh_face_counts": {
            g.label: dict(inlet=g.faces_inlet, outlet=g.faces_outlet,
                          bottom=g.faces_bottom, top=g.faces_top,
                          obstacle=g.faces_obstacle)
            for g in GRIDS
        },
    }
    out_path = Path(__file__).parent / "gci_summary.json"
    out_path.write_text(json.dumps(out, indent=2, default=str))
    print(f"\nWrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
