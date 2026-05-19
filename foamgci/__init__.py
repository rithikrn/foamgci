"""foamgci — verification utilities for unsteady OpenFOAM CFD.

Public API
----------
Reader:
    read_fieldminmax      — parse OpenFOAM ``fieldMinMax.dat``
Statistics:
    window_stats          — time-averaged mean / std on a stationary window
    tau_int_geyer         — Geyer initial-positive-sequence τ_int
    sem_autocorr_corrected— autocorrelation-aware standard error of the mean
    kpss_test             — KPSS stationarity test (level + trend variants)
Grid convergence:
    roache_gci            — Roache GCI from a refinement triplet
    richardson_extrapolation — exact-grid extrapolated value
    apparent_order        — observed order of accuracy p̂
Reporting:
    full_report           — end-to-end: read files → window stats → GCI → table
    rayleigh_pitot        — analytical M-3-shock stagnation pressure ratio

References
----------
- Roache, P.J. (1994), "Perspective: A method for uniform reporting of grid
  refinement studies", J. Fluids Eng., 116, 405-413.
- Celik, I.B. et al. (2008), "Procedure for estimation and reporting of
  uncertainty due to discretization in CFD applications", J. Fluids Eng., 130.
- Geyer, C.J. (1992), "Practical Markov chain Monte Carlo", Stat. Sci., 7, 473.
- Kwiatkowski, D. et al. (1992), "Testing the null hypothesis of stationarity",
  J. Econometrics, 54, 159-178.
"""
from .reader import read_fieldminmax, FieldMinMaxData
from .stats import (
    window_stats,
    tau_int_geyer,
    sem_autocorr_corrected,
    kpss_test,
)
from .gci import (
    roache_gci,
    richardson_extrapolation,
    apparent_order,
    GCIResult,
)
from .report import full_report, rayleigh_pitot, ReportTable

__version__ = "0.2.0"

__all__ = [
    "read_fieldminmax",
    "FieldMinMaxData",
    "window_stats",
    "tau_int_geyer",
    "sem_autocorr_corrected",
    "kpss_test",
    "roache_gci",
    "richardson_extrapolation",
    "apparent_order",
    "GCIResult",
    "full_report",
    "rayleigh_pitot",
    "ReportTable",
    "__version__",
]
