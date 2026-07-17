"""
extract_snapshot.py
-------------------
Sample an OpenFOAM diamond-airfoil case at a target time (default t = 10, the
steady field) onto a uniform Cartesian grid in the (x, y) plane and save the
resampled fields as a compressed NumPy archive (NPZ), for the density/pressure
contour figure.  Mirrors the wedge15Ma5 / forwardstep_mach3 extractors; only
the domain and the solid-body mask (the airfoil triangle, not a step) differ.

Usage:
    python3 extract_snapshot.py                 # all five grids
    python3 extract_snapshot.py --grid Fine     # one grid
    python3 extract_snapshot.py --time 10 --nx 700 --ny 240
    python3 extract_snapshot.py --entropy       # also write data/<prefix>_entropy.dat

Snapshot output: ./snapshots/snap_<label>_t<time>.npz with arrays
    x, y, t, label, valid, p, rho, (optionally Ux, Uy, Umag, T)
Cells inside the solid airfoil body (below the upper surface) are NaN-masked.

--entropy additionally computes the domain entropy volume integral
    S_vol = integral of deltaS dV,  deltaS = Cv*ln(p/pInf) - Cp*ln(rho/rhoInf)
on the NATIVE cells (no interpolation), the smooth volume QoI analyze.py reads
as <prefix>_entropy.dat.  This replaces the version-fragile coded function
object that was removed from system/controlDict.

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
    ) from e

from data import (
    GRIDS, T_BENCHMARK, DOMAIN_X, DOMAIN_Y, DATA_DIR, in_airfoil_body,
)

# Thermodynamic constants (match 0/ and the removed entropy function object).
CV, CP, P_INF, RHO_INF = 1.78571, 2.5, 1.0, 1.4


def get_internal_mesh(foam_path: Path):
    """Open a .foam stub and return the configured reader.

    POpenFOAMReader exposes 'internalMesh' as an item in patch_array_names,
    so we MUST NOT blanket-disable patches here.  Enable all arrays.
    """
    reader = pv.get_reader(str(foam_path))
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


def _internal_block(mesh):
    if isinstance(mesh, pv.MultiBlock):
        if "internalMesh" in mesh.keys():
            return mesh["internalMesh"]
        return next(b for b in mesh if isinstance(b, pv.UnstructuredGrid))
    return mesh


def _cell_scalar(internal, name: str) -> np.ndarray:
    """Return a per-cell scalar, converting point data if necessary."""
    if name in internal.cell_data:
        return np.asarray(internal.cell_data[name], dtype=float)
    if name in internal.point_data:
        conv = internal.point_data_to_cell_data()
        return np.asarray(conv.cell_data[name], dtype=float)
    raise KeyError(f"field {name!r} not found in cell or point data")


def entropy_integral(internal) -> float:
    """Domain integral of specific-entropy rise on native cells (no interp).

    S_vol = sum_c [Cv*ln(p_c/pInf) - Cp*ln(rho_c/rhoInf)] * V_c
    Matches the OpenFOAM volFieldValue volIntegrate the coded FO used to feed.
    """
    sized = internal.compute_cell_sizes(length=False, area=False, volume=True)
    vol = np.asarray(sized.cell_data["Volume"], dtype=float)
    p = _cell_scalar(internal, "p")
    rho = _cell_scalar(internal, "rho")
    dS = CV * np.log(p / P_INF) - CP * np.log(rho / RHO_INF)
    return float(np.sum(dS * vol))


def extract_one(foam_path: Path, label: str, prefix: str, t_target: float,
                nx: int, ny: int, out_dir: Path, want_entropy: bool) -> Path:
    if not foam_path.is_file():
        raise FileNotFoundError(f".foam stub not found: {foam_path}")
    print(f"  {label:11s} <- {foam_path}")

    reader = get_internal_mesh(foam_path)
    t_actual = closest_time(reader, t_target)
    reader.set_active_time_value(t_actual)
    mesh = reader.read()

    internal = _internal_block(mesh)
    if internal is None or internal.n_points == 0:
        raise RuntimeError(f"{label}: internalMesh has zero points")

    if want_entropy:
        s_vol = entropy_integral(internal)
        DATA_DIR.mkdir(exist_ok=True)
        ent_path = DATA_DIR / f"{prefix}_entropy.dat"
        ent_path.write_text(
            "# Time  volIntegrate(deltaS)\n"
            f"{t_actual:.8g}  {s_vol:.10g}\n"
        )
        print(f"    entropy S_vol = {s_vol:.6g}  ->  {ent_path.name}")

    zmid = 0.5 * (internal.bounds[4] + internal.bounds[5])
    xi = np.linspace(DOMAIN_X[0], DOMAIN_X[1], nx)
    yi = np.linspace(DOMAIN_Y[0], DOMAIN_Y[1], ny)
    X, Y = np.meshgrid(xi, yi)
    pts = np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, zmid)])
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

    # Mask the solid airfoil body (below the upper surface, LE_X..TE_X).
    valid = ~in_airfoil_body(X, Y)
    for v in fields.values():
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
                        default=None, help="Extract only one grid (default: all)")
    parser.add_argument("--time", type=float, default=T_BENCHMARK)
    parser.add_argument("--nx", type=int, default=700)
    parser.add_argument("--ny", type=int, default=240)
    parser.add_argument("--entropy", action="store_true",
                        help="Also write data/<prefix>_entropy.dat (S_vol QoI)")
    parser.add_argument("--out", type=Path,
                        default=Path(__file__).parent / "snapshots")
    args = parser.parse_args()

    grids = [g for g in GRIDS if args.grid is None or g.label == args.grid]
    print(f"Extracting t={args.time} snapshots onto "
          f"{args.nx}x{args.ny} uniform grid")
    print(f"Output dir: {args.out}")

    for g in grids:
        try:
            extract_one(g.foam_path, g.label, g.prefix, args.time,
                        args.nx, args.ny, args.out, args.entropy)
        except (FileNotFoundError, RuntimeError, KeyError) as e:
            print(f"  ! {g.label}: {e}")
            continue
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
