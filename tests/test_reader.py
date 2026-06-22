"""Tests for foamgci.reader."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from foamgci.reader import read_fieldminmax, read_timeseries


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
    
def test_read_timeseries_roundtrip(tmp_path):
    p = tmp_path / "series.csv"
    p.write_text("# t value x y z\n0.0 1.0 0.5 0.1 0\n0.1 2.0 0.5 0.1 0\n")
    d = read_timeseries(p, field="cd", time_col=0, value_col=1, loc_cols=(2, 3, 4))
    assert len(d) == 2
    assert d.max[1] == 2.0 and d.min[1] == 2.0
    assert d.loc_max.shape == (2, 3)


# ---- surfaceFieldValue reader ----------------------------------------------

from foamgci.reader import read_surface_field_value  # noqa: E402


SFV_SINGLE = (
    "# Region type : patch obstacle\n"
    "# Faces  : 160\n"
    "# Area   : 7.890000e-02\n"
    "# Time             areaAverage(p)\n"
    "0.002              4.77910\n"
    "0.004              4.78050\n"
    "0.006              4.78081\n"
    "0.008              4.78083\n"
)

SFV_MULTI_DISTINCT = (
    "# Time          areaAverage(p)   areaAverage(rho)\n"
    "0.002           4.77910          2.75300\n"
    "0.004           4.78050          2.75340\n"
)

SFV_MULTI_SAME_FIELD = (
    "# Time          areaAverage(p)   areaIntegrate(p)\n"
    "0.002           4.77910          0.37700\n"
    "0.004           4.78050          0.37711\n"
)


def test_sfv_single_column_auto(tmp_path: Path) -> None:
    f = _write(tmp_path / "surfaceFieldValue.dat", SFV_SINGLE)
    d = read_surface_field_value(f)
    assert d.column == "areaAverage(p)"
    assert d.field == "p"
    assert len(d) == 4
    assert d.value[-1] == pytest.approx(4.78083)
    assert d.time[0] == pytest.approx(0.002)


def test_sfv_token_select(tmp_path: Path) -> None:
    f = _write(tmp_path / "surfaceFieldValue.dat", SFV_MULTI_DISTINCT)
    assert read_surface_field_value(f, column="p").column == "areaAverage(p)"
    rho = read_surface_field_value(f, column="rho")
    assert rho.field == "rho"
    assert rho.value[1] == pytest.approx(2.75340)


def test_sfv_token_ambiguous_raises(tmp_path: Path) -> None:
    f = _write(tmp_path / "surfaceFieldValue.dat", SFV_MULTI_SAME_FIELD)
    with pytest.raises(ValueError, match="[Aa]mbiguous"):
        read_surface_field_value(f, column="p")


def test_sfv_full_label_and_index_select(tmp_path: Path) -> None:
    f = _write(tmp_path / "surfaceFieldValue.dat", SFV_MULTI_SAME_FIELD)
    assert read_surface_field_value(f, column="areaIntegrate(p)").column == "areaIntegrate(p)"
    assert read_surface_field_value(f, column=0).column == "areaAverage(p)"
    assert read_surface_field_value(f, column=1).column == "areaIntegrate(p)"


def test_sfv_multi_column_without_selector_raises(tmp_path: Path) -> None:
    f = _write(tmp_path / "surfaceFieldValue.dat", SFV_MULTI_DISTINCT)
    with pytest.raises(ValueError, match="value columns"):
        read_surface_field_value(f)


def test_sfv_missing_file(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        read_surface_field_value(tmp_path / "nope.dat")


def test_sfv_index_out_of_range_raises(tmp_path: Path) -> None:
    f = _write(tmp_path / "surfaceFieldValue.dat", SFV_SINGLE)
    with pytest.raises(ValueError, match="out of range"):
        read_surface_field_value(f, column=5)


def test_sfv_no_header_single_column(tmp_path: Path) -> None:
    # Header missing entirely: a lone value column must still parse.
    body = "0.002  4.77910\n0.004  4.78050\n"
    f = _write(tmp_path / "surfaceFieldValue.dat", body)
    d = read_surface_field_value(f)
    assert len(d) == 2
    assert d.value[0] == pytest.approx(4.77910)


def test_sfv_restrict(tmp_path: Path) -> None:
    f = _write(tmp_path / "surfaceFieldValue.dat", SFV_SINGLE)
    d = read_surface_field_value(f).restrict(0.004, 0.008)
    assert len(d) == 3
    assert d.time.min() >= 0.004
