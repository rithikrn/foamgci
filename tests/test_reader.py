"""Tests for foamgci.reader."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from foamgci.reader import read_fieldminmax


COMBINED_HEADER = (
    "# Field minMax\n"
    "# Time         field   min   location(min)         processor   max   "
    "location(max)         processor\n"
)

PER_FIELD_HEADER = (
    "# Field minMax for p\n"
    "# Time         min   location(min)         processor   max   "
    "location(max)         processor\n"
)


def _write(p: Path, body: str) -> Path:
    p.write_text(body)
    return p


# ---- Combined dialect ------------------------------------------------------

def test_combined_dialect_single_field(tmp_path: Path) -> None:
    body = COMBINED_HEADER + (
        "0.0001  p  0.99  (0.012 0.5 0)  0  1.01  (0.123 0.6 0)  0\n"
        "0.0002  p  0.98  (0.011 0.5 0)  0  1.03  (0.124 0.6 0)  0\n"
        "0.0003  p  0.97  (0.010 0.5 0)  0  1.05  (0.125 0.6 0)  0\n"
    )
    f = _write(tmp_path / "fieldMinMax.dat", body)
    d = read_fieldminmax(f, field="p")
    assert d.field == "p"
    np.testing.assert_allclose(d.time, [0.0001, 0.0002, 0.0003])
    np.testing.assert_allclose(d.max, [1.01, 1.03, 1.05])
    np.testing.assert_allclose(d.min, [0.99, 0.98, 0.97])
    assert d.loc_max is not None
    np.testing.assert_allclose(d.loc_max[-1], [0.125, 0.6, 0.0])


def test_combined_dialect_multi_field_selects_one(tmp_path: Path) -> None:
    body = COMBINED_HEADER + (
        "0.0001  p  0.99  (0 0 0)  0  1.01  (0 0 0)  0\n"
        "0.0001  U  0.00  (0 0 0)  0  3.00  (1 1 0)  0\n"
        "0.0002  p  0.98  (0 0 0)  0  1.03  (0 0 0)  0\n"
        "0.0002  U  0.01  (0 0 0)  0  3.01  (1 1 0)  0\n"
    )
    f = _write(tmp_path / "fieldMinMax.dat", body)
    dp = read_fieldminmax(f, field="p")
    du = read_fieldminmax(f, field="U")
    assert len(dp) == 2
    assert len(du) == 2
    np.testing.assert_allclose(dp.max, [1.01, 1.03])
    np.testing.assert_allclose(du.max, [3.00, 3.01])


def test_combined_dialect_multi_field_without_field_raises(tmp_path: Path) -> None:
    body = COMBINED_HEADER + (
        "0.0001  p  0.99  (0 0 0)  0  1.01  (0 0 0)  0\n"
        "0.0001  U  0.00  (0 0 0)  0  3.00  (1 1 0)  0\n"
    )
    f = _write(tmp_path / "fieldMinMax.dat", body)
    with pytest.raises(ValueError, match="multiple fields"):
        read_fieldminmax(f)


# ---- Per-field dialect -----------------------------------------------------

def test_per_field_dialect(tmp_path: Path) -> None:
    body = PER_FIELD_HEADER + (
        "0.0001  0.99  (0.012 0.5 0)  0  1.01  (0.123 0.6 0)  0\n"
        "0.0002  0.98  (0.011 0.5 0)  0  1.03  (0.124 0.6 0)  0\n"
    )
    f = _write(tmp_path / "fieldMinMax.dat", body)
    d = read_fieldminmax(f, field="p")
    assert d.field == "p"
    np.testing.assert_allclose(d.max, [1.01, 1.03])


def test_per_field_dialect_without_field_raises(tmp_path: Path) -> None:
    body = PER_FIELD_HEADER + "0.0001  0.99  (0 0 0)  0  1.01  (0 0 0)  0\n"
    f = _write(tmp_path / "fieldMinMax.dat", body)
    with pytest.raises(ValueError, match="per-field"):
        read_fieldminmax(f)


# ---- Edge cases ------------------------------------------------------------

def test_missing_file(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        read_fieldminmax(tmp_path / "nope.dat", field="p")


def test_no_locations(tmp_path: Path) -> None:
    body = (
        "# Time   field  min   max\n"
        "0.0001  p  0.99  1.01\n"
        "0.0002  p  0.98  1.03\n"
    )
    f = _write(tmp_path / "fieldMinMax.dat", body)
    d = read_fieldminmax(f, field="p")
    assert d.loc_max is None
    assert d.loc_min is None


def test_restrict(tmp_path: Path) -> None:
    body = COMBINED_HEADER + "".join(
        f"{t:.4f}  p  0.99  (0 0 0)  0  {1.0 + 0.01*i:.4f}  (0 0 0)  0\n"
        for i, t in enumerate(np.linspace(0.0, 1.0, 11))
    )
    f = _write(tmp_path / "fieldMinMax.dat", body)
    d = read_fieldminmax(f, field="p")
    d2 = d.restrict(0.3, 0.7)
    assert d2.time.min() >= 0.3
    assert d2.time.max() <= 0.7
    assert len(d2) == 5
