"""foamgci.reader — parse OpenFOAM function-object output.

This module reads the two OpenFOAM function-object outputs that ``foamgci``
treats as first-class QoI sources, plus a generic fallback:

* :func:`read_fieldminmax` — ``fieldMinMax.dat`` (pointwise field extrema and
  their locations); used by the forward-step example.
* :func:`read_surface_field_value` — ``surfaceFieldValue.dat`` (integrated
  surface functionals such as ``areaAverage(p)``); used by the wedge example.
* :func:`read_timeseries` — a generic ``(time, value[, x, y, z])`` reader for
  solver-independent scalar histories.

A single OpenFOAM case can write several of these at once, and one analysis
driver may read more than one per grid (the wedge example reads both a
``surfaceFieldValue.dat`` and a ``fieldMinMax.dat`` from every run).

The ``fieldMinMax`` function object in OpenFOAM writes one row per write
timestep. The format has evolved across versions; this parser handles both
of the two common dialects:

  1. **Combined dialect** (one file holds all requested fields, OpenFOAM ≥ v5)::

         # Time         field   min   location(min)         processor   max   location(max)         processor
         0.0001         p       0.99  (0.012 0.5 0)         0           1.01  (0.123 0.6 0)         0
         0.0001         U       0.0   (0.0 0.5 0)           0           3.0   (1.5  0.5 0)          0
         ...

  2. **Per-field dialect** (one file per field, older versions)::

         # Time         min   location(min)            processor   max   location(max)            processor
         0.0001         0.99  (0.012 0.5 0)            0           1.01  (0.123 0.6 0)            0

The parser returns a :class:`FieldMinMaxData` instance: column-wise NumPy
arrays for *time*, *min*, *max*, and (when available) the locations of the
extrema. Comment lines beginning with ``#`` are skipped.

When the file does not declare which field a row belongs to, the user must
supply the field name explicitly (``field=``).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union

import numpy as np


@dataclass
class FieldMinMaxData:
    """Parsed contents of a single field's ``fieldMinMax.dat`` time-series.

    Attributes
    ----------
    field : str
        Field name (e.g., ``'p'``, ``'U'``, ``'mag(U)'``, ``'rho'``).
    time : np.ndarray, shape (N,)
        Time values at which the extrema were recorded.
    min : np.ndarray, shape (N,)
        Field minimum at each time.
    max : np.ndarray, shape (N,)
        Field maximum at each time.
    loc_min : np.ndarray | None, shape (N, 3)
        Spatial location of the minimum, if present in the file.
    loc_max : np.ndarray | None, shape (N, 3)
        Spatial location of the maximum, if present in the file.
    source : Path
        Source file path (for traceability in reports).
    """

    field: str
    time: np.ndarray
    min: np.ndarray
    max: np.ndarray
    loc_min: Optional[np.ndarray] = None
    loc_max: Optional[np.ndarray] = None
    source: Optional[Path] = None

    def __len__(self) -> int:
        return int(self.time.size)

    def restrict(self, t0: float, t1: float) -> "FieldMinMaxData":
        """Return a new ``FieldMinMaxData`` restricted to ``t0 ≤ t ≤ t1``."""
        m = (self.time >= t0) & (self.time <= t1)
        return FieldMinMaxData(
            field=self.field,
            time=self.time[m],
            min=self.min[m],
            max=self.max[m],
            loc_min=self.loc_min[m] if self.loc_min is not None else None,
            loc_max=self.loc_max[m] if self.loc_max is not None else None,
            source=self.source,
        )


# ---------------------------------------------------------------------------
# Tokenization
# ---------------------------------------------------------------------------

def _tokenize(line: str) -> list[str]:
    """Tokenize an OpenFOAM data line, treating ``(...)`` as a single token.

    This is what makes the format awkward to read with :func:`numpy.loadtxt`:
    locations are written as ``(x y z)`` with internal whitespace.
    """
    out: list[str] = []
    buf: list[str] = []
    depth = 0
    for ch in line.rstrip():
        if ch == "(":
            depth += 1
            buf.append(ch)
        elif ch == ")":
            depth -= 1
            buf.append(ch)
            if depth == 0:
                out.append("".join(buf))
                buf = []
        elif ch.isspace() and depth == 0:
            if buf:
                out.append("".join(buf))
                buf = []
        else:
            buf.append(ch)
    if buf:
        out.append("".join(buf))
    return out


def _parse_loc(token: str) -> np.ndarray:
    """Parse a ``(x y z)`` location token into a length-3 array."""
    inner = token.strip().lstrip("(").rstrip(")")
    parts = inner.split()
    if len(parts) != 3:
        raise ValueError(f"Bad location token: {token!r}")
    return np.array([float(p) for p in parts], dtype=float)


# ---------------------------------------------------------------------------
# Public reader
# ---------------------------------------------------------------------------

def read_fieldminmax(
    path: Union[str, Path],
    field: Optional[str] = None,
) -> FieldMinMaxData:
    """Read an OpenFOAM ``fieldMinMax.dat`` file.

    Parameters
    ----------
    path : str | Path
        Path to the ``fieldMinMax.dat`` file (any OpenFOAM dialect).
    field : str, optional
        Required for the *combined dialect* with multiple fields per file:
        select which field's rows to extract. For the per-field dialect, the
        caller must pass the field name explicitly (the file itself does not
        record it).

    Returns
    -------
    FieldMinMaxData

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    ValueError
        If the file is empty, has an unrecognised header, or ``field`` is
        ambiguous / missing.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"No such file: {path}")

    rows: list[list[str]] = []
    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            rows.append(_tokenize(line))
    if not rows:
        raise ValueError(f"No data rows in {path}")

    # Dialect detection ---------------------------------------------------
    # Combined dialect: token[1] is non-numeric (field name).
    # Per-field dialect: token[1] is numeric (min value).
    first = rows[0]
    try:
        float(first[1])
        per_field_dialect = True
    except ValueError:
        per_field_dialect = False

    times: list[float] = []
    mins: list[float] = []
    maxs: list[float] = []
    lmins: list[np.ndarray] = []
    lmaxs: list[np.ndarray] = []
    have_locations = False

    if per_field_dialect:
        if field is None:
            raise ValueError(
                f"{path} appears to be per-field dialect (no field column); "
                "pass field= explicitly."
            )
        # Layout: time min [loc_min proc] max [loc_max proc]
        for r in rows:
            t = float(r[0])
            vmin = float(r[1])
            idx = 2
            lmin = lmax = None
            if idx < len(r) and r[idx].startswith("("):
                lmin = _parse_loc(r[idx])
                idx += 1
                # skip processor index if present
                try:
                    int(r[idx])
                    idx += 1
                except (ValueError, IndexError):
                    pass
            vmax = float(r[idx])
            idx += 1
            if idx < len(r) and r[idx].startswith("("):
                lmax = _parse_loc(r[idx])
            times.append(t)
            mins.append(vmin)
            maxs.append(vmax)
            if lmin is not None and lmax is not None:
                lmins.append(lmin)
                lmaxs.append(lmax)
                have_locations = True
    else:
        # Combined dialect: time field min [loc_min proc] max [loc_max proc]
        # Filter on the requested field.
        if field is None:
            seen = {r[1] for r in rows}
            if len(seen) == 1:
                field = next(iter(seen))
            else:
                raise ValueError(
                    f"{path} contains multiple fields {sorted(seen)}; "
                    "pass field= to select one."
                )
        for r in rows:
            if r[1] != field:
                continue
            t = float(r[0])
            vmin = float(r[2])
            idx = 3
            lmin = lmax = None
            if idx < len(r) and r[idx].startswith("("):
                lmin = _parse_loc(r[idx])
                idx += 1
                try:
                    int(r[idx])
                    idx += 1
                except (ValueError, IndexError):
                    pass
            vmax = float(r[idx])
            idx += 1
            if idx < len(r) and r[idx].startswith("("):
                lmax = _parse_loc(r[idx])
            times.append(t)
            mins.append(vmin)
            maxs.append(vmax)
            if lmin is not None and lmax is not None:
                lmins.append(lmin)
                lmaxs.append(lmax)
                have_locations = True

    if not times:
        raise ValueError(
            f"No rows for field={field!r} found in {path}."
        )

    return FieldMinMaxData(
        field=field,
        time=np.array(times, dtype=float),
        min=np.array(mins, dtype=float),
        max=np.array(maxs, dtype=float),
        loc_min=np.array(lmins) if have_locations else None,
        loc_max=np.array(lmaxs) if have_locations else None,
        source=path,
    )
# ---------------------------------------------------------------------------
# surfaceFieldValue reader (integrated surface functionals)
# ---------------------------------------------------------------------------

@dataclass
class SurfaceFieldValueData:
    """Parsed contents of a single column of a ``surfaceFieldValue.dat`` file.

    The ``surfaceFieldValue`` (a.k.a. ``fieldValues``/``surfaceRegion``)
    function object writes one row per write timestep and one column per
    requested ``operation(field)`` (e.g. ``areaAverage(p)``,
    ``areaIntegrate(p)``). This dataclass holds the time axis and a single
    selected value column.

    Attributes
    ----------
    column : str
        Full column label as it appears in the file header
        (e.g. ``'areaAverage(p)'``).
    field : str
        The field token inside the operation parentheses (e.g. ``'p'``), or
        the raw column label when no parentheses are present.
    time : np.ndarray, shape (N,)
        Sample times.
    value : np.ndarray, shape (N,)
        The selected functional value at each time.
    source : Path | None
        Source file path, for traceability in reports.
    """

    column: str
    field: str
    time: np.ndarray
    value: np.ndarray
    source: Optional[Path] = None

    def __len__(self) -> int:
        return int(self.time.size)

    def restrict(self, t0: float, t1: float) -> "SurfaceFieldValueData":
        """Return a copy restricted to ``t0 ≤ t ≤ t1``."""
        m = (self.time >= t0) & (self.time <= t1)
        return SurfaceFieldValueData(
            column=self.column,
            field=self.field,
            time=self.time[m],
            value=self.value[m],
            source=self.source,
        )


def _field_token(label: str) -> str:
    """Extract the field name from an operation label.

    ``'areaAverage(p)' -> 'p'``; ``'areaIntegrate(mag(U))' -> 'mag(U)'``;
    a label with no parentheses is returned unchanged.
    """
    i = label.find("(")
    if i == -1:
        return label
    j = label.rfind(")")
    if j <= i:
        return label
    return label[i + 1 : j]


def read_surface_field_value(
    path: Union[str, Path],
    column: Optional[Union[str, int]] = None,
) -> SurfaceFieldValueData:
    """Read an OpenFOAM ``surfaceFieldValue.dat`` file.

    The file is self-describing: a header comment line beginning ``# Time``
    lists the column labels, e.g.::

        # Region type : patch obstacle
        # Faces  : 160
        # Area   : 0.0789...
        # Time             areaAverage(p)
        0.002              4.7791
        0.004              4.7805
        ...

    Parameters
    ----------
    path : str | Path
        Path to the ``surfaceFieldValue.dat`` file.
    column : str | int, optional
        Which value column to return.

        * ``None`` (default): require exactly one non-time column and return
          it; raise if the file has several columns.
        * ``str``: match by field token — ``column='p'`` selects the column
          whose label is ``...(p)`` (e.g. ``areaAverage(p)``). Falls back to a
          case-sensitive substring match on the full label, so passing the
          full label (``'areaAverage(p)'``) also works.
        * ``int``: a 0-based index into the **value** columns (column 0 is the
          first quantity after ``Time``).

    Returns
    -------
    SurfaceFieldValueData

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    ValueError
        If the file is empty, has no parseable ``# Time`` header, or
        ``column`` is ambiguous / not found.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"No such file: {path}")

    header_labels: Optional[list[str]] = None
    data_rows: list[list[str]] = []
    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                # The column header is the last comment line that starts a
                # token sequence beginning with 'Time'.
                body = line.lstrip("#").strip()
                toks = _tokenize(body)
                if toks and toks[0].lower() == "time":
                    header_labels = toks
                continue
            data_rows.append(_tokenize(line))

    if not data_rows:
        raise ValueError(f"No data rows in {path}")

    n_value_cols = len(data_rows[0]) - 1
    if n_value_cols < 1:
        raise ValueError(f"{path}: rows have no value column after Time.")

    # Value-column labels (drop the leading 'Time' from the header if present).
    if header_labels is not None and len(header_labels) >= 1:
        value_labels = header_labels[1:]
    else:
        value_labels = []
    # Pad / fall back to positional labels if the header is missing or short.
    if len(value_labels) != n_value_cols:
        value_labels = [
            value_labels[i] if i < len(value_labels) else f"col{i + 1}"
            for i in range(n_value_cols)
        ]

    # Resolve the requested column to a 0-based value-column index.
    if column is None:
        if n_value_cols != 1:
            raise ValueError(
                f"{path} has {n_value_cols} value columns {value_labels}; "
                "pass column= to select one (by field token, label, or index)."
            )
        sel = 0
    elif isinstance(column, int):
        if not (0 <= column < n_value_cols):
            raise ValueError(
                f"column index {column} out of range for {n_value_cols} "
                f"value columns {value_labels} in {path}."
            )
        sel = column
    else:
        # String: match by field token first, then by substring of the label.
        tokens = [_field_token(lbl) for lbl in value_labels]
        matches = [i for i, tok in enumerate(tokens) if tok == column]
        if not matches:
            matches = [i for i, lbl in enumerate(value_labels) if column in lbl]
        if not matches:
            raise ValueError(
                f"No column matching {column!r} in {path}; "
                f"available labels: {value_labels}."
            )
        if len(matches) > 1:
            raise ValueError(
                f"Ambiguous column {column!r} in {path}: matches "
                f"{[value_labels[i] for i in matches]}. Pass the full label."
            )
        sel = matches[0]

    times: list[float] = []
    values: list[float] = []
    for r in data_rows:
        if len(r) <= sel + 1:
            raise ValueError(
                f"{path}: a data row has fewer columns than the header "
                f"({len(r)} tokens, need >= {sel + 2}). Row: {r}"
            )
        times.append(float(r[0]))
        values.append(float(r[sel + 1]))

    label = value_labels[sel]
    return SurfaceFieldValueData(
        column=label,
        field=_field_token(label),
        time=np.array(times, dtype=float),
        value=np.array(values, dtype=float),
        source=path,
    )


def read_timeseries(
    path,
    *,
    field,
    time_col=0,
    value_col=1,
    loc_cols=None,        # e.g. (2, 3, 4) for x, y, z; or None
    delimiter=None,       # None => any whitespace
    comments="#",
):
    
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"No such file: {path}")
    t, v, loc = [], [], []
    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith(comments):
                continue
            r = line.split(delimiter)
            t.append(float(r[time_col]))
            v.append(float(r[value_col]))
            if loc_cols is not None:
                loc.append([float(r[c]) for c in loc_cols])
    if not t:
        raise ValueError(f"No data rows in {path}")
    t = np.asarray(t, float)
    v = np.asarray(v, float)
    locs = np.asarray(loc, float) if loc_cols is not None else None
    return FieldMinMaxData(
        field=field, time=t, min=v, max=v,
        loc_min=locs, loc_max=locs, source=path,
    )
