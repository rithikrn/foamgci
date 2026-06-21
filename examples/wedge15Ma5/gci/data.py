from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

# gci/ directory (this file lives in gci/) and the example root above it.
GCI_DIR = Path(__file__).resolve().parent
DATA_DIR = GCI_DIR / "data"
EXAMPLE_ROOT = GCI_DIR.parent

# Steady-regime time window used for time-averaging the wall-pressure signal.
# If KPSS rejects stationarity on this window for any grid, widen the run or
# narrow the window here and re-run analyze.py.
T_STAT = (0.15, 0.20)

# Benchmark time for snapshots/contours. The wedge is steady by endTime, so
# the snapshot is just the final field.
T_BENCHMARK = 0.20

# ----------------------------------------------------------------------
# Free-stream conditions (from 0/U, 0/p, 0/T and thermophysicalProperties):
#   U = (5,0,0), p = 1, T = 1, gamma = 1.4  =>  a = sqrt(gamma*R*T) = 1,
#   so M_inflow = 5.  Ramp (flow-deflection) angle = 15 deg.  Inviscid.
# ----------------------------------------------------------------------
M_INFLOW = 5.0
THETA_DEG = 15.0
GAMMA = 1.4
P_INF = 1.0   # free-stream static pressure (normalised); QoI is p_wall / P_INF

# ----------------------------------------------------------------------
# Geometry (from system/blockMeshDict), used by the snapshot + contour
# figures. The fluid sits above the ramp; the wedge body below the ramp
# line y = x*tan(theta) for x > 0 is solid and is masked out.
# ----------------------------------------------------------------------
DOMAIN_X = (-0.15242, 0.3048)
DOMAIN_Y = (0.0, 0.1524)
RAMP_TIP_X = 0.0


def ramp_y(x):
    """Height of the ramp surface at station x (0 upstream of the tip)."""
    x = np.asarray(x, dtype=float)
    return np.where(x > RAMP_TIP_X, (x - RAMP_TIP_X) * np.tan(np.radians(THETA_DEG)), 0.0)


def in_wedge_body(X, Y):
    """True where (X, Y) lies inside the solid wedge (below the ramp)."""
    return (X > RAMP_TIP_X) & (Y < ramp_y(X))


@dataclass(frozen=True)
class GridSpec:
    label: str          # human-readable, used in tables/figures
    dat: str            # filename in gci/data/
    case_dirname: str   # run directory (sibling of gci/) for snapshot extraction
    n_cells: int        # total cells (block1 + block2)
    dx: float           # uniform square-cell spacing h = 0.1524 / Ny

    @property
    def sfv_path(self) -> Path:
        """surfaceFieldValue.dat input read by analyze.py."""
        return DATA_DIR / self.dat

    @property
    def foam_path(self) -> Path:
        """`.foam` stub inside this grid's run directory, read by pyvista
        in extract_snapshot.py. The run directory is a sibling of gci/."""
        return EXAMPLE_ROOT / self.case_dirname / "case.foam"


# ----------------------------------------------------------------------
# Four grids.  Refinement ratio r = 2 between successive levels
# (cells x4 in 2D).  Medium == the upstream OpenFOAM-4.x tutorial mesh.
# h = domain_height / Ny = 0.1524 / Ny; cells are square (dx = dy).
# ----------------------------------------------------------------------
COARSE     = GridSpec("Coarse",     "coarse.dat",     "coarse_grid",     1_200, 0.0076200)
MEDIUM     = GridSpec("Medium",     "medium.dat",     "medium_grid",     4_800, 0.0038100)
FINE       = GridSpec("Fine",       "fine.dat",       "fine_grid",      19_200, 0.0019050)
EXTRA_FINE = GridSpec("Extra-fine", "extrafine.dat",  "extrafine_grid", 76_800, 0.0009525)

GRIDS = [COARSE, MEDIUM, FINE, EXTRA_FINE]
