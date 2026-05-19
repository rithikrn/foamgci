"""verify_abstract.py — regenerate the AIAA-abstract Table 1 from real data.

This script is the bridge between the OpenFOAM Mach-3 forward-facing
step runs and the abstract numbers. Given the four cases on disk, it:

    1. Reads each ``fieldMinMax.dat`` (binary OpenFOAM output).
    2. Computes the time-averaged maximum pressure on the stationary
       window ``[t0, t1]`` for each grid.
    3. Estimates Geyer's τ_int and the autocorrelation-corrected SEM.
    4. Runs the KPSS stationarity test on each window.
    5. Applies Roache GCI to every consecutive triplet.
    6. Compares the Richardson-extrapolated value to the analytical
       Rayleigh-Pitot p_02/p_1 for M=3, γ=1.4.
    7. Prints the LaTeX block to be pasted directly into the abstract.

Run from the repository root::

    python examples/forwardstep_mach3/verify_abstract.py \\
        --coarse     /path/to/coarse/postProcessing/fieldMinMax/0/fieldMinMax.dat \\
        --medium     /path/to/medium/postProcessing/fieldMinMax/0/fieldMinMax.dat \\
        --fine       /path/to/fine/postProcessing/fieldMinMax/0/fieldMinMax.dat \\
        --extra-fine /path/to/extra-fine/postProcessing/fieldMinMax/0/fieldMinMax.dat \\
        --window 3 10

Grid sizes are hard-coded below to match the four runs documented in
the abstract. Cell counts are taken from `log.checkMesh`:

    coarse     :   4,032 cells   h = 0.025
    medium     :  16,128 cells   h = 0.0125
    fine       :  64,512 cells   h = 0.00625
    extra-fine : 258,048 cells   h = 0.003125

Refinement ratio r = 2 between consecutive grids (cell count grows
×4 in 2D, ×8 in 3D — the FFS is 2D so the cell ratio is exactly 4).
"""
from __future__ import annotations

import argparse
from pathlib import Path

from foamgci import (
    GridCase,
    full_report,
    rayleigh_pitot,
)
from foamgci.plot import plot_convergence


# --------------------------------------------------------------------------
# Fixed grid properties of the four runs (from log.checkMesh).
# These are intrinsic to the cases; do not change unless re-meshing.
# --------------------------------------------------------------------------
GRID_SPECS = [
    # (label,        h,        n_cells)
    ("coarse",       0.025,     4_032),
    ("medium",       0.0125,   16_128),
    ("fine",         0.00625,  64_512),
    ("extra-fine",   0.003125, 258_048),
]


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Regenerate the abstract Table 1 from real OpenFOAM runs.",
    )
    p.add_argument("--coarse",     required=True, type=Path,
                   help="path to coarse-grid fieldMinMax.dat")
    p.add_argument("--medium",     required=True, type=Path,
                   help="path to medium-grid fieldMinMax.dat")
    p.add_argument("--fine",       required=True, type=Path,
                   help="path to fine-grid fieldMinMax.dat")
    p.add_argument("--extra-fine", dest="extra_fine", required=True, type=Path,
                   help="path to extra-fine-grid fieldMinMax.dat")
    p.add_argument("--field", default="p",
                   help="field name to read from fieldMinMax.dat (default: p)")
    p.add_argument("--window", nargs=2, type=float, default=[3.0, 10.0],
                   metavar=("T0", "T1"),
                   help="stationary window for time-averaging (default: 3 10)")
    p.add_argument("--mach", type=float, default=3.0,
                   help="freestream Mach number for Rayleigh-Pitot reference")
    p.add_argument("--gamma", type=float, default=1.4,
                   help="ratio of specific heats (default: 1.4)")
    p.add_argument("--out", type=Path, default=Path("out"),
                   help="directory for output files (default: ./out)")
    p.add_argument("--no-plot", action="store_true",
                   help="skip the matplotlib convergence figure")
    return p


def main() -> int:
    args = _build_parser().parse_args()

    paths = [args.coarse, args.medium, args.fine, args.extra_fine]
    for p in paths:
        if not p.exists():
            print(f"ERROR: {p} does not exist.")
            return 2

    cases = [
        GridCase(label=label, path=path, h=h, n_cells=n)
        for path, (label, h, n) in zip(paths, GRID_SPECS)
    ]

    p_ref = rayleigh_pitot(args.mach, args.gamma)
    print(
        f"\nAnalytical Rayleigh-Pitot p_02/p_1 (M={args.mach}, γ={args.gamma}): "
        f"{p_ref:.6f}\n"
    )

    rep = full_report(
        cases=cases,
        field=args.field,
        quantity="max",
        window=tuple(args.window),
        reference_value=p_ref,
        reference_label=f"Rayleigh-Pitot p_02/p_1 (M={args.mach})",
    )

    args.out.mkdir(parents=True, exist_ok=True)
    text = rep.as_text()
    print(text)
    (args.out / "report.txt").write_text(text + "\n")

    latex = rep.as_latex()
    (args.out / "table1.tex").write_text(latex)
    print(f"\n  → wrote {args.out/'report.txt'}")
    print(f"  → wrote {args.out/'table1.tex'} (paste into the abstract)")

    if not args.no_plot:
        try:
            plot_convergence(rep, out_path=args.out / "fig_convergence.pdf")
            print(f"  → wrote {args.out/'fig_convergence.pdf'}")
        except ImportError as exc:
            print(f"  (skipping plot: {exc})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
