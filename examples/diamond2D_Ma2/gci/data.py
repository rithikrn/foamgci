from __future__ import annotations

from dataclasses import dataclass
from math import sqrt, tan, radians
from pathlib import Path

# gci/ directory (this file lives in gci/) and the example root above it.
GCI_DIR = Path(__file__).resolve().parent
DATA_DIR = GCI_DIR / "data"
EXAMPLE_ROOT = GCI_DIR.parent

# Steady-regime time window for averaging. The M=2 diamond reaches a steady
# attached-shock field; this window sits on the established plateau. If KPSS
# rejects stationarity on any grid, widen the run (endTime) or narrow here.
T_STAT = (6.0, 10.0)
T_BENCHMARK = 10.0

# ----------------------------------------------------------------------
# Free-stream (from 0/U, 0/p, 0/T and thermophysicalProperties):
#   U = (2,0,0), p = 1, T = 1, gamma = 1.4  =>  a = sqrt(gamma*R*T) = 1,
#   so M_inflow = 2.  Symmetric diamond, half-angle 10 deg, chord 1, inviscid.
# Half-domain (upper half, symmetry at y=0): airfoil coefficients are the
# upper-surface contribution and are DOUBLED in post-processing.
# ----------------------------------------------------------------------
M_INFLOW = 2.0
EPS_DEG = 10.0          # diamond half-angle
GAMMA = 1.4
P_INF = 1.0
CHORD = 1.0
SPAN = 0.1             # z-extent in blockMeshDict (-0.05..0.05)
T_OVER_C = tan(radians(EPS_DEG))   # diamond thickness ratio = tan(eps)

# Geometry for figures (from system/blockMeshDict).
DOMAIN_X = (-1.0, 2.5)
DOMAIN_Y = (0.0, 1.2)
PEAK_Y = 0.5 * tan(radians(EPS_DEG))   # 0.0881635 upper-surface peak

# Fluid area (half-domain) = box minus the upper-half airfoil triangle.
# Used for the volume-averaged representative cell size h = sqrt(area / N).
_BOX_AREA = (DOMAIN_X[1] - DOMAIN_X[0]) * (DOMAIN_Y[1] - DOMAIN_Y[0])
_AIRFOIL_UPPER_AREA = 0.5 * CHORD * PEAK_Y
FLUID_AREA = _BOX_AREA - _AIRFOIL_UPPER_AREA


@dataclass(frozen=True)
class GridSpec:
    label: str          # human-readable
    prefix: str         # data/<prefix>_<qoi>.dat
    n_cells: int        # total cells (sum of the four blocks)

    @property
    def h(self) -> float:
        """Volume-averaged representative cell size, h = sqrt(area / N)."""
        return sqrt(FLUID_AREA / self.n_cells)

    # one .dat per QoI per grid (see system/controlDict function objects)
    @property
    def force_path(self) -> Path:
        return DATA_DIR / f"{self.prefix}_forceCoeffs.dat"

    @property
    def pfront_path(self) -> Path:
        return DATA_DIR / f"{self.prefix}_pFront.dat"

    @property
    def prear_path(self) -> Path:
        return DATA_DIR / f"{self.prefix}_pRear.dat"

    @property
    def pmax_path(self) -> Path:
        return DATA_DIR / f"{self.prefix}_fieldMinMax.dat"

    @property
    def svol_path(self) -> Path:
        return DATA_DIR / f"{self.prefix}_sVol.dat"


# ----------------------------------------------------------------------
# Five grids, systematic refinement ratio r = 2 (cells x4 in 2D), so the
# least-squares (Eca-Hoekstra) fit has enough levels to be meaningful while
# the deepest triplet still drives the Roache GCI. Cell counts match the
# four (nx ny 1) tuples documented in system/blockMeshDict.
# ----------------------------------------------------------------------
GRIDS = [
    GridSpec("Coarse",     "coarse",    1_120),
    GridSpec("Medium",     "medium",    4_480),
    GridSpec("Fine",       "fine",     17_920),
    GridSpec("Extra-fine", "extrafine", 71_680),
    GridSpec("Ultra-fine", "ultrafine", 286_720),
]

# Half-domain doubling: the modelled upper-surface drag is half the total.
CD_HALF_TO_FULL = 2.0
