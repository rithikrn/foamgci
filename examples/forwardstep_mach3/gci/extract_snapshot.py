"""
extract_snapshot.py
-------------------
Sample an OpenFOAM case at a target time (default t=4) onto a uniform
Cartesian grid in the (x, y) plane and save the resampled fields as
a compressed NumPy archive (NPZ).

Usage:
    python3 extract_snapshot.py             # extract all four cases
    python3 extract_snapshot.py --grid Fine # one specific case
    python3 extract_snapshot.py --time 4.0 --nx 600 --ny 200

Output goes to ./snapshots/snap_<label>_t<time>.npz with arrays
    x, y, rho, p, t, (optionally Ux, Uy, T)
Cells inside the step body are masked with NaN.

Requirements: pyvista >= 0.43 (handles both OpenFOAMReader and
POpenFOAMReader transparently via pv.get_reader).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import numpy as np

try:
    import pyvista as pv
except ImportError as e:
    raise SystemExit(
        "PyVista is required for snapshot extraction.  Install with:\n"
        "    pip install --user pyvista\n"
        f"Original error: {e}"
    )

from data import GRIDS, T_BENCHMARK


# Greenshields 2010 used the full domain in their Fig. 6.
DOMAIN_X = (0.0, 3.0)
DOMAIN_Y = (0.0, 1.0)
STEP_X = (0.6, 3.0)
STEP_Y = (0.0, 0.2)


def is_decomposed(case_dir: Path) -> bool:
    return (case_dir / "processor0").is_dir()


def get_internal_mesh(foam_path: Path):
    """Open a .foam stub and return the configured reader.

    Note: POpenFOAMReader exposes 'internalMesh' as an item in
    patch_array_names, so we MUST NOT blanket-disable patches here
    -- doing so silently drops the internal volume.  Disable only
    named boundary patches if you need to.
    """
    reader = pv.get_reader(str(foam_path))
    # Make sure cell/point data are enabled (some readers default off).
    if hasattr(reader, "enable_all_cell_arrays"):
        reader.enable_all_cell_arrays()
    if hasattr(reader, "enable_all_point_arrays"):
        reader.enable_all_point_arrays()
    return reader


def closest_time(reader, t_target: float) -> float:
    times = np.asarray(reader.time_values, dtype=float)
    if times.size == 0:
        raise RuntimeError("Reader reports no time directories.  "
                           "Did the case actually run?")
    idx = int(np.argmin(np.abs(times - t_target)))
    return float(times[idx])


def extract_one(foam_path: Path, label: str, t_target: float,
                nx: int, ny: int, out_dir: Path) -> Path:
    if not foam_path.is_file():
        raise FileNotFoundError(f".foam stub not found: {foam_path}")
    print(f"  {label:11s} <- {foam_path}")

    reader = get_internal_mesh(foam_path)
    t_actual = closest_time(reader, t_target)
    reader.set_active_time_value(t_actual)
    mesh = reader.read()

    # The reader yields a MultiBlock; the internal volume is 'internalMesh'.
    if isinstance(mesh, pv.MultiBlock):
        if "internalMesh" in mesh.keys():
            internal = mesh["internalMesh"]
        else:
            # Take the first non-empty UnstructuredGrid block.
            internal = next(
                b for b in mesh if isinstance(b, pv.UnstructuredGrid)
            )
    else:
        internal = mesh
    if internal is None or internal.n_points == 0:
        raise RuntimeError(f"{label}: internalMesh has zero points")

    # Pick z slice at the midplane of the extruded thickness.
    zmid = 0.5 * (internal.bounds[4] + internal.bounds[5])

    xi = np.linspace(DOMAIN_X[0], DOMAIN_X[1], nx)
    yi = np.linspace(DOMAIN_Y[0], DOMAIN_Y[1], ny)
    X, Y = np.meshgrid(xi, yi)
    pts = np.column_stack([X.ravel(), Y.ravel(),
                           np.full(X.size, zmid)])
    probe = pv.PolyData(pts).sample(internal,
                                    pass_cell_data=False,
                                    pass_point_data=True)

    have = set(probe.array_names)
    fields: dict[str, np.ndarray] = {}
    for k in ("rho", "p", "T"):
        if k in have:
            fields[k] = np.asarray(probe[k]).reshape(ny, nx)
    if "U" in have:
        U = np.asarray(probe["U"]).reshape(ny, nx, 3)
        fields["Ux"] = U[..., 0]
        fields["Uy"] = U[..., 1]
        fields["Umag"] = np.linalg.norm(U, axis=-1)

    # Mask cells inside the step body.
    inside_step = ((X >= STEP_X[0]) & (X <= STEP_X[1]) &
                   (Y >= STEP_Y[0]) & (Y <= STEP_Y[1]))
    valid = ~inside_step
    for k, v in fields.items():
        v[~valid] = np.nan

    out_dir.mkdir(exist_ok=True)
    out = out_dir / f"snap_{label.lower().replace('-', '')}_t{t_actual:.3f}.npz"
    np.savez_compressed(
        out,
        x=xi, y=yi, t=float(t_actual),
        label=label, foam=str(foam_path),
        valid=valid.astype(np.uint8),
        **fields,
    )
    print(f"    wrote {out}  (t_actual = {t_actual:.4f}, "
          f"nx x ny = {nx} x {ny}, fields = {list(fields)})")
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", choices=[g.label for g in GRIDS],
                        default=None,
                        help="Extract only one grid (default: all)")
    parser.add_argument("--time", type=float, default=T_BENCHMARK)
    parser.add_argument("--nx", type=int, default=600)
    parser.add_argument("--ny", type=int, default=200)
    parser.add_argument("--out", type=Path,
                        default=Path(__file__).parent / "snapshots")
    args = parser.parse_args()

    grids = [g for g in GRIDS
             if args.grid is None or g.label == args.grid]
    print(f"Extracting t={args.time} snapshots onto "
          f"{args.nx}x{args.ny} uniform grid")
    print(f"Output dir: {args.out}")

    for g in grids:
        try:
            extract_one(g.foam_path, g.label, args.time,
                        args.nx, args.ny, args.out)
        except (FileNotFoundError, RuntimeError) as e:
            print(f"  ! {g.label}: {e}")
            continue
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
