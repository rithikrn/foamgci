"""
parse_fieldminmax.py
--------------------
Parse OpenFOAM `fieldMinMax` function-object output and compute
time-averaged extrema in the statistically stationary regime.

Expected file format (one row per timestep per field):

    # Time   field  min   location(min)  processor  max   location(max)  processor

The `location()` columns are parenthesised 3-vectors "(x y z)".

The parser is robust to:
  * extra header comment lines
  * extra columns (some OF versions emit "min(magSqr)" etc.)
  * vector fields written with `mode magnitude` (column 'mag(U)') --
    in this case the min/max columns are still scalar magnitudes.

For vector fields written with `mode component`, three triplets appear
per row; that variant is NOT handled here (rarely used in practice).
Switch the case to `mode magnitude` for U.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


# A parenthesised 3-vector of floats like "(0.59 0.012 0.005)".
# This deliberately does NOT match "(U)" or "(magSqr)" appearing in
# field names -- only well-formed location triples.
_FLOAT = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
_LOC = re.compile(
    rf"\(\s*({_FLOAT})\s+({_FLOAT})\s+({_FLOAT})\s*\)"
)
_LOC_PLACE = re.compile(
    rf"\(\s*{_FLOAT}\s+{_FLOAT}\s+{_FLOAT}\s*\)"
)


@dataclass
class FieldMinMax:
    """All time-series data for one tracked field from fieldMinMax.dat."""
    field: str
    t:        np.ndarray    # (N,)
    v_min:    np.ndarray    # (N,)
    v_max:    np.ndarray    # (N,)
    loc_min:  np.ndarray    # (N, 3)
    loc_max:  np.ndarray    # (N, 3)


def parse_file(path: Path | str) -> dict[str, FieldMinMax]:
    """Read a fieldMinMax.dat and return {field_name: FieldMinMax}."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"fieldMinMax.dat not found: {path}")

    buckets: dict[str, dict[str, list]] = {}

    with open(path) as fh:
        for raw in fh:
            s = raw.strip()
            if not s or s.startswith("#"):
                continue
            # Find well-formed location triples (won't match "(U)").
            loc_matches = _LOC.findall(s)
            if len(loc_matches) < 2:
                continue
            # Replace location triples with a placeholder so the rest of
            # the line tokenises cleanly on whitespace.
            s_flat = _LOC_PLACE.sub("LOC", s)
            tok = s_flat.split()
            # Canonical layout: [t, field, min, LOC, proc, max, LOC, proc].
            if len(tok) < 8:
                continue
            try:
                t = float(tok[0])
                field = tok[1]
                v_min = float(tok[2])
                v_max = float(tok[5])
            except (ValueError, IndexError):
                continue
            lmin = np.asarray(loc_matches[0], dtype=float)
            lmax = np.asarray(loc_matches[1], dtype=float)
            b = buckets.setdefault(field, dict(t=[], vmin=[], vmax=[],
                                               lmin=[], lmax=[]))
            b["t"].append(t)
            b["vmin"].append(v_min)
            b["vmax"].append(v_max)
            b["lmin"].append(lmin)
            b["lmax"].append(lmax)

    out: dict[str, FieldMinMax] = {}
    for fname, b in buckets.items():
        out[fname] = FieldMinMax(
            field=fname,
            t=np.asarray(b["t"]),
            v_min=np.asarray(b["vmin"]),
            v_max=np.asarray(b["vmax"]),
            loc_min=np.vstack(b["lmin"]),
            loc_max=np.vstack(b["lmax"]),
        )
    return out


def time_average(t: np.ndarray, y: np.ndarray,
                 window: tuple[float, float]) -> tuple[float, float, int]:
    """Trapezoidal mean and sample std of `y(t)` over `window=(t1,t2)`.

    Returns (mean, std, n_samples).  If fewer than 3 samples fall in
    the window, returns (NaN, NaN, n) so the caller can fall back.
    """
    t1, t2 = window
    mask = (t >= t1) & (t <= t2)
    n = int(mask.sum())
    if n < 3:
        return float("nan"), float("nan"), n
    tm = t[mask]
    ym = y[mask]
    mean = float(np.trapezoid(ym, tm) / (tm[-1] - tm[0]))
    std = float(np.std(ym, ddof=1))
    return mean, std, n


def late_time_location(loc: np.ndarray, t: np.ndarray,
                       window: tuple[float, float]) -> tuple[float, float, float]:
    """Median (x, y, z) of a tracked location over the stationary window.
    Median is robust to occasional sample-flip artifacts when two cells
    tie for the extremum."""
    mask = (t >= window[0]) & (t <= window[1])
    if not mask.any():
        return float("nan"), float("nan"), float("nan")
    med = np.median(loc[mask], axis=0)
    return float(med[0]), float(med[1]), float(med[2])


# ---------------------------------------------------------------------
# High-level helper used by analyze.py.
# ---------------------------------------------------------------------
@dataclass
class StationaryStats:
    """All summary numbers needed for the GCI study, per case."""
    label: str
    folder: str
    n_cells: int
    dx: float
    n_samples_stat: int
    # Primary GCI metric: time-mean of post-bow-shock stagnation pressure.
    p_max_mean: float
    p_max_std: float
    # Reported for the manuscript (not used for GCI).
    rho_max_mean: float
    rho_max_std: float
    U_mag_max_mean: float
    U_mag_max_std: float
    # Late-time stagnation-point location.
    peak_loc_x: float
    peak_loc_y: float


def summarise_case(case_label: str, folder: str, n_cells: int, dx: float,
                   fieldminmax_path: Path,
                   window: tuple[float, float],
                   U_field_names: Iterable[str] = ("mag(U)", "U"),
                   ) -> StationaryStats:
    """Parse one case and return its StationaryStats."""
    fields = parse_file(fieldminmax_path)

    if "p" not in fields:
        raise RuntimeError(f"{case_label}: field 'p' missing from "
                           f"{fieldminmax_path}")
    if "rho" not in fields:
        raise RuntimeError(f"{case_label}: field 'rho' missing from "
                           f"{fieldminmax_path}")

    fp = fields["p"]
    fr = fields["rho"]

    p_mean, p_std, n = time_average(fp.t, fp.v_max, window)
    r_mean, r_std, _ = time_average(fr.t, fr.v_max, window)

    # Velocity magnitude -- try common field names.
    u_mean, u_std = float("nan"), float("nan")
    for name in U_field_names:
        if name in fields:
            u_mean, u_std, _ = time_average(fields[name].t,
                                            fields[name].v_max, window)
            break

    px, py, _ = late_time_location(fp.loc_max, fp.t, window)

    return StationaryStats(
        label=case_label, folder=folder,
        n_cells=n_cells, dx=dx,
        n_samples_stat=n,
        p_max_mean=p_mean, p_max_std=p_std,
        rho_max_mean=r_mean, rho_max_std=r_std,
        U_mag_max_mean=u_mean, U_mag_max_std=u_std,
        peak_loc_x=px, peak_loc_y=py,
    )
