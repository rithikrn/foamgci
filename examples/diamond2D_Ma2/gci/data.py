from __future__ import annotations

from dataclasses import dataclass
from math import sqrt, tan, radians
from pathlib import Path

import numpy as np

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


def upper_surface(x):
    """Upper-surface height of the diamond at station x; 0 outside the chord.
    Front facet (0->0.5) and rear facet (0.5->1) both have slope tan(eps)."""
    x = np.asarray(x, dtype=float)
    t = np.tan(np.radians(EPS_DEG))
    return np.where((x > 0.0) & (x < 1.0), t * np.minimum(x, 1.0 - x), 0.0)


def in_airfoil_body(X, Y):
    """True where (X, Y) lies inside the solid upper-half diamond body
    (below the upper surface, chordwise 0<x<1). Masked out of snapshots."""
    return (X > 0.0) & (X < 1.0) & (Y < upper_surface(X))

# Upper-surface geometry, used by extract_snapshot.py to mask the solid body.
LE_X, TE_X = 0.0, CHORD            # leading / trailing edge x
_HALF_CHORD = 0.5 * CHORD          # peak sits at x = LE_X + half-chord


def airfoil_y(x):
    """Upper-surface height of the diamond for LE_X <= x <= TE_X, else 0."""
    import numpy as np
    x = np.asarray(x, dtype=float)
    inside = (x >= LE_X) & (x <= TE_X)
    return np.where(inside, PEAK_Y * (1.0 - np.abs(x - _HALF_CHORD) / _HALF_CHORD), 0.0)


def in_airfoil_body(X, Y):
    """True where (X, Y) lies inside the solid diamond (below the upper surface)."""
    return (X >= LE_X) & (X <= TE_X) & (Y < airfoil_y(X))

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
        return DATA_DIR / f"{self.prefix}_coefficient.dat"

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

    @property
    def entropy_path(self) -> Path:
        """Native-cell entropy volume integral written by compute_entropy.py
        (one scalar per grid; generated locally from the t=10 field)."""
        return DATA_DIR / f"{self.prefix}_entropy.dat"

    @property
    def foam_path(self) -> Path:
        """`.foam` stub in this grid's run directory (sibling of gci/), read by
        pyvista in extract_snapshot.py. Run dir is <prefix>_grid/."""
        return EXAMPLE_ROOT / f"{self.prefix}_grid" / "case.foam"

    @property
    def foam_path(self) -> Path:
        """`.foam` stub inside this grid's run directory (sibling of the
        example root), read by pyvista in extract_snapshot.py. Matches the
        HPC run-dir naming <prefix>_grid used in submit.sh."""
        return EXAMPLE_ROOT / f"{self.prefix}_grid" / "case.foam"


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
