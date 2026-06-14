from .reader import read_fieldminmax, read_timeseries, FieldMinMaxData
from .stats import (
    window_stats,
    tau_int_geyer,
    sem_autocorr_corrected,
    kpss_test,
)
from .gci import (
    roache_gci,
    gci_over_hierarchy,
    richardson_extrapolation,
    apparent_order,
    GCIResult,
)
from .report import full_report, rayleigh_pitot, ReportTable, GridCase

from ._version import __version__

__all__ = [
    "read_fieldminmax",
    "FieldMinMaxData",
    "window_stats",
    "tau_int_geyer",
    "sem_autocorr_corrected",
    "kpss_test",
    "roache_gci",
    "gci_over_hierarchy",
    "richardson_extrapolation",
    "apparent_order",
    "GCIResult",
    "full_report",
    "rayleigh_pitot",
    "ReportTable",
    "GridCase",
    "read_timeseries",
    "__version__",
]
