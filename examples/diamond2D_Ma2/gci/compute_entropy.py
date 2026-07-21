"""
compute_entropy.py
------------------
Compute the specific-entropy volume integral for each grid from the native
(cell-centred) t=10 field, and write one scalar per grid to
    data/<prefix>_entropy.dat
which analyze.py ingests as the smooth VOLUME QoI (S_vol).

    deltaS = Cv*ln(p/pInf) - Cp*ln(rho/rhoInf)          (specific-entropy rise)
    S_vol  = sum_cells( deltaS_cell * cellVolume )       (half-domain integral)

Steady flow => one converged field (t=10) per grid is all the GCI needs; there
is NO in-solver function object and NO runtime compilation. Uses native cell
data + cell volumes directly (NOT the uniform-grid probe from
extract_snapshot.py, which would interpolate and mask the body and corrupt a
volume integral).

Usage:
    python3 compute_entropy.py                # all grids
    python3 compute_entropy.py --grid Fine    # one grid
    python3 compute_entropy.py --time 10

Requirements: pyvista >= 0.43. Reads the same <prefix>_grid/case.foam run
directories (Bundle B) as extract_snapshot.py.
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
        "PyVista is required for the entropy integral.  Install with:\n"
        "    pip install --user pyvista\n"
        f"Original error: {e}"
    ) from e

from data import GRIDS, T_BENCHMARK, GAMMA, P_INF, DATA_DIR

# Non-dimensional thermodynamics (a = sqrt(gamma R T) = 1 at T=1):
R = 1.0 / GAMMA                 # 0.71429
CV = R / (GAMMA - 1.0)          # 1.78571
CP = GAMMA * CV                 # 2.5
RHO_INF = P_INF / R             # 1.4  (at T_inf = 1)


def get_internal_mesh(foam_path: Path):
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


def _cell_field(internal, name: str) -> np.ndarray:
    """Return a cell-centred array for `name`, converting from point data
    if the reader only exposed it there."""
    if name in internal.cell_data:
        return np.asarray(internal.cell_data[name])
    if name in internal.point_data:
        conv = internal.point_data_to_cell_data()
        return np.asarray(conv.cell_data[name])
    raise KeyError(f"field {name!r} not found in cell or point data "
                   f"(have cell={list(internal.cell_data)}, "
                   f"point={list(internal.point_data)})")


def entropy_integral(foam_path: Path, t_target: float):
    if not foam_path.is_file():
        raise FileNotFoundError(f".foam stub not found: {foam_path}")
    reader = get_internal_mesh(foam_path)
    t_actual = closest_time(reader, t_target)
    reader.set_active_time_value(t_actual)
    mesh = reader.read()

    if isinstance(mesh, pv.MultiBlock):
        if "internalMesh" in mesh.keys():
            internal = mesh["internalMesh"]
        else:
            internal = next(b for b in mesh
                            if isinstance(b, pv.UnstructuredGrid))
    else:
        internal = mesh
    if internal is None or internal.n_cells == 0:
        raise RuntimeError(f"internalMesh has zero cells: {foam_path}")

    sized = internal.compute_cell_sizes(length=False, area=False, volume=True)
    vol = np.asarray(sized.cell_data["Volume"], dtype=float)
    p = _cell_field(internal, "p").astype(float)
    rho = _cell_field(internal, "rho").astype(float)

    dS = CV * np.log(p / P_INF) - CP * np.log(rho / RHO_INF)
    s_vol = float(np.sum(dS * vol))
    return s_vol, t_actual, int(internal.n_cells)


def write_dat(prefix: str, s_vol: float, t_actual: float, n_cells: int) -> Path:
    DATA_DIR.mkdir(exist_ok=True)
    out = DATA_DIR / f"{prefix}_entropy.dat"
    out.write_text(
        "# specific-entropy volume integral (half-domain), native-cell "
        "volume-weighted\n"
        "# deltaS = Cv*ln(p/pInf) - Cp*ln(rho/rhoInf); "
        f"Cv={CV:.6g} Cp={CP:.6g} pInf={P_INF:.6g} rhoInf={RHO_INF:.6g}\n"
        f"# t={t_actual:.6g}  n_cells={n_cells}\n"
        f"{s_vol:.10g}\n",
        encoding="utf-8",
    )
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--grid", choices=[g.label for g in GRIDS], default=None)
    ap.add_argument("--time", type=float, default=T_BENCHMARK)
    args = ap.parse_args()

    grids = [g for g in GRIDS if args.grid is None or g.label == args.grid]
    print(f"Entropy integral at t={args.time} "
          f"(Cv={CV:.5g}, Cp={CP:.5g}, pInf={P_INF:g}, rhoInf={RHO_INF:g})")
    for g in grids:
        try:
            s_vol, t_actual, n = entropy_integral(g.foam_path, args.time)
        except (FileNotFoundError, RuntimeError, KeyError) as e:
            print(f"  ! {g.label}: {e}")
            continue
        out = write_dat(g.prefix, s_vol, t_actual, n)
        print(f"  {g.label:11s} S_vol={s_vol:.6g}  (t={t_actual:.4g}, "
              f"N={n})  -> {out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())