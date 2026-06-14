from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

# gci/ directory (this file lives in gci/).
GCI_DIR = Path(__file__).resolve().parent
DATA_DIR = GCI_DIR / "data"

# Stationary-regime time window used everywhere for time averaging.
# If KPSS rejects stationarity on this window for any grid, narrow it
# here and re-run analyze.py.
T_STAT = (3.0, 10.0)

# Benchmark comparison time (Greenshields 2010 Fig. 6, WC84 Fig. 4(d)).
T_BENCHMARK = 4.0


@dataclass(frozen=True)
class GridSpec:
    label: str          # human-readable, used in tables/figures
    dat: str            # filename in gci/data/
    n_cells: int        # total cells (from log.checkMesh)
    dx: float           # uniform Cartesian spacing h

    @property
    def fieldminmax_path(self) -> Path:
        return DATA_DIR / self.dat


# ----------------------------------------------------------------------
# The four grids.  Refinement ratio r = 2 between successive levels
# (cells x 4 in 2D).  Medium == Greenshields 2010 "Mesh 1";
# Fine == Greenshields "Mesh 2".
# ----------------------------------------------------------------------
COARSE     = GridSpec("Coarse",     "coarse.dat",     4_032,   0.025)
MEDIUM     = GridSpec("Medium",     "medium.dat",    16_128,   0.0125)
FINE       = GridSpec("Fine",       "fine.dat",      64_512,   0.00625)
EXTRA_FINE = GridSpec("Extra-fine", "extrafine.dat", 258_048,  0.003125)

GRIDS = [COARSE, MEDIUM, FINE, EXTRA_FINE]


# Greenshields et al. (2010) IJNMF 63, 1-21, Section 4.2 -- reference values
GREENSHIELDS_2010 = {
    "rho_contour_min":   0.2568,
    "rho_contour_max":   6.067,
    "n_contours":        30,
    "comparison_time":   4.0,
    "M_inflow":          3.0,
    "gamma":             1.4,
    "Mesh1_cells":       16_128,    # our Medium
    "Mesh2_cells":       64_512,    # our Fine
}
