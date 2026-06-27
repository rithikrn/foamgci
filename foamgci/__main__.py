"""foamgci CLI.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence

from ._version import __version__
from .report import full_report, rayleigh_pitot, GridCase


def _parse_case(spec: str) -> GridCase:
    parts = spec.split(":")
    if len(parts) not in (3, 4):
        raise argparse.ArgumentTypeError(
            f"--case must be label:path:h[:n_cells]; got {spec!r}"
        )
    label = parts[0]
    path = Path(parts[1])
    try:
        h = float(parts[2])
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"--case spacing must be a number; got {parts[2]!r}"
        ) from exc
    n_cells: int | None = None
    if len(parts) == 4:
        try:
            n_cells = int(parts[3])
        except ValueError as exc:
            raise argparse.ArgumentTypeError(
                f"--case n_cells must be an integer; got {parts[3]!r}"
            ) from exc
    return GridCase(label=label, path=path, h=h, n_cells=n_cells)


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="foamgci",
        description="Verification utilities for unsteady OpenFOAM CFD: "
                    "Roache GCI + autocorrelation-corrected SEM + KPSS.",
    )
    p.add_argument("--version", action="version", version=f"foamgci {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)
    pr = sub.add_parser("report", help="generate a V&V report from fieldMinMax.dat")
    pr.add_argument(
        "--case", action="append", required=True, type=_parse_case,
        help="label:path:h[:n_cells]; pass once per grid, coarse-to-fine",
    )
    pr.add_argument("--field", default="p",
                    help="field name in fieldMinMax.dat (default: p)")
    pr.add_argument("--quantity", choices=["max", "min"], default="max",
                    help="which spatial extremum to track (default: max)")
    pr.add_argument("--window", nargs=2, type=float, default=[3.0, 10.0],
                    metavar=("T0", "T1"),
                    help="stationary window for time-averaging "
                         "(default: 3 10)")
    pr.add_argument("--reference", default=None,
                    help="reference value: a number, or 'rayleigh-pitot'")
    pr.add_argument("--mach", type=float, default=3.0,
                    help="Mach number for rayleigh-pitot reference (default: 3)")
    pr.add_argument("--gamma", type=float, default=1.4,
                    help="γ for rayleigh-pitot reference (default: 1.4)")
    pr.add_argument("--kpss", choices=["c", "ct"], default="c",
                    help="KPSS regression: c=level, ct=trend (default: c)")
    pr.add_argument("--text", type=Path, default=None,
                    help="write plain-text report to this path")
    pr.add_argument("--latex", type=Path, default=None,
                    help="write LaTeX table to this path")
    pr.add_argument("--version", action="version", version=f"foamgci {__version__}")
    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    if args.cmd != "report":
        print(f"Unknown command: {args.cmd}", file=sys.stderr)
        return 2

    ref_value: float | None = None
    ref_label = ""
    if args.reference is not None:
        if args.reference == "rayleigh-pitot":
            ref_value = rayleigh_pitot(args.mach, gamma=args.gamma)
            ref_label = f"Rayleigh-Pitot p_02/p_1 at M={args.mach}, γ={args.gamma}"
        else:
            try:
                ref_value = float(args.reference)
                ref_label = f"user-supplied reference = {ref_value:g}"
            except ValueError:
                print(f"Bad --reference value: {args.reference!r}", file=sys.stderr)
                return 2

    report = full_report(
        cases=args.case,
        field=args.field,
        quantity=args.quantity,
        window=tuple(args.window),
        reference_value=ref_value,
        reference_label=ref_label,
        kpss_regression=args.kpss,
    )

    text = report.as_text()
    print(text)
    if args.text is not None:
        args.text.parent.mkdir(parents=True, exist_ok=True)
        args.text.write_text(text + "\n")
    if args.latex is not None:
        args.latex.parent.mkdir(parents=True, exist_ok=True)
        args.latex.write_text(report.as_latex())
        print(f"\nLaTeX table written to {args.latex}", file=sys.stderr)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
