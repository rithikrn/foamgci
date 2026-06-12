"""
data.py
-------
Grid metadata and case paths for the four-grid forward-step study.

Single source of truth: change paths or grid attributes here and the
rest of the pipeline follows.  No physical measurements are hardcoded;
those come from parsing each case's postProcessing output.

Directory layout expected (relative to the GCI script directory):

    ../coarse_grid/coarsecase.foam
    ../medium_grid/mediumcase.foam
    ../fine_grid/finecase.foam
    ../extrafine_grid/extrafinecase.foam

with each case folder containing the standard OpenFOAM tree plus
    postProcessing/fieldMinMax1/0/fieldMinMax.dat
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


# Root of the four case folders, relative to this file (gci/).
# All four cases live as sister directories of `gci/`.
PARENT = Path(__file__).resolve().parent.parent

# Common postProcessing relative path for fieldMinMax output.
# Override per-grid in GridSpec below if a case writes to a different
# functionObject name (e.g. fieldMinMax2/0/fieldMinMax.dat).
FIELDMINMAX_REL = Path("postProcessing/fieldRange/0/fieldMinMax.dat")

# Stationary-regime time window used everywhere for time averaging.
T_STAT = (3.0, 10.0)

# Benchmark comparison time (Greenshields 2010 Fig. 6, WC84 Fig. 4(d)).
T_BENCHMARK = 4.0


@dataclass(frozen=True)
class GridSpec:
    label: str
    folder: str            # e.g. "coarse_grid"
    foam_file: str         # e.g. "coarsecase.foam"
    n_cells: int
    dx: float              # uniform Cartesian spacing
    # Boundary face counts from checkMesh (mesh statistics for table 1).
    faces_inlet: int
    faces_outlet: int
    faces_bottom: int
    faces_top: int
    faces_obstacle: int
    # Optional override if a particular case used a non-default
    # fieldMinMax function-object name; None -> use FIELDMINMAX_REL.
    fieldminmax_rel: Path | None = None

    @property
    def case_dir(self) -> Path:
        return PARENT / self.folder

    @property
    def foam_path(self) -> Path:
        return self.case_dir / self.foam_file

    @property
    def fieldminmax_path(self) -> Path:
        rel = self.fieldminmax_rel or FIELDMINMAX_REL
        return self.case_dir / rel


# ----------------------------------------------------------------------
# The four grids.  Refinement ratio r = 2 between successive levels
# (cells x 4 in 2D).  Medium == Greenshields 2010 "Mesh 1";
# Fine == Greenshields "Mesh 2".
# ----------------------------------------------------------------------
COARSE = GridSpec(
    label="Coarse",
    folder="coarse_grid",
    foam_file="coarsecase.foam",
    n_cells=4_032,
    dx=0.025,
    faces_inlet=40, faces_outlet=32,
    faces_bottom=24, faces_top=120, faces_obstacle=104,
)

MEDIUM = GridSpec(
    label="Medium",
    folder="medium_grid",
    foam_file="mediumcase.foam",
    n_cells=16_128,
    dx=0.0125,
    faces_inlet=80, faces_outlet=64,
    faces_bottom=48, faces_top=240, faces_obstacle=208,
)

FINE = GridSpec(
    label="Fine",
    folder="fine_grid",
    foam_file="finecase.foam",
    n_cells=64_512,
    dx=0.00625,
    faces_inlet=160, faces_outlet=128,
    faces_bottom=96, faces_top=480, faces_obstacle=416,
)

EXTRA_FINE = GridSpec(
    label="Extra-fine",
    folder="extrafine_grid",
    foam_file="extrafinecase.foam",
    n_cells=258_048,
    dx=0.003125,
    faces_inlet=320, faces_outlet=256,
    faces_bottom=192, faces_top=960, faces_obstacle=832,
)

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
