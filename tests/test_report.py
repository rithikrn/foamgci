"""Tests for foamgci.report — end-to-end V&V pipeline."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from foamgci.report import full_report, rayleigh_pitot, GridCase


# ---- Rayleigh-Pitot ------------------------------------------------------

def test_rayleigh_pitot_M3_gamma14() -> None:
    """Canonical FFS reference: M=3, γ=1.4 → p_02/p_1 ≈ 12.0610."""
    assert math.isclose(rayleigh_pitot(3.0, 1.4), 12.0610, abs_tol=1e-3)


def test_rayleigh_pitot_M2_gamma14() -> None:
    # Anderson 'Modern Compressible Flow' Table A.2: M=2 → p_02/p_1 = 5.6404
    assert math.isclose(rayleigh_pitot(2.0, 1.4), 5.6404, abs_tol=1e-3)


def test_rayleigh_pitot_rejects_subsonic() -> None:
    with pytest.raises(ValueError, match="supersonic"):
        rayleigh_pitot(0.8)


# ---- End-to-end with synthetic 4-grid hierarchy --------------------------

def _make_synthetic_fieldminmax(
    path: Path,
    field: str,
    mean: float,
    sigma: float,
    n: int = 500,
    t_end: float = 10.0,
    seed: int = 0,
) -> None:
    rng = np.random.default_rng(seed)
    t = np.linspace(0.0, t_end, n)
    series = mean + sigma * rng.standard_normal(n)
    body = ["# Time   field   min   location(min)   processor   max   location(max)   processor"]
    for ti, si in zip(t, series):
        body.append(
            f"{ti:.6f}  {field}  {si - 1.0:.6f}  (0 0 0)  0  {si:.6f}  (0.6 0 0)  0"
        )
    path.write_text("\n".join(body) + "\n")


def test_full_report_end_to_end(tmp_path: Path) -> None:
    """Synthetic 4-grid study converging to 12.06 as h → 0, with r=2."""
    targets = {
        "coarse":     (11.99, 0.025),
        "medium":     (12.02, 0.0125),
        "fine":       (12.05, 0.00625),
        "extra-fine": (12.058, 0.003125),
    }
    cases = []
    for i, (label, (mean, h)) in enumerate(targets.items()):
        f = tmp_path / f"{label}_fieldMinMax.dat"
        # Deterministic per-grid seed. NB: do NOT use hash(label) — Python
        # salts str hashing per process (PYTHONHASHSEED), which would make
        # this test non-reproducible across runs.
        _make_synthetic_fieldminmax(f, "p", mean=mean, sigma=0.01, seed=100 + i)
        cases.append(GridCase(label=label, path=f, h=h))

    rep = full_report(
        cases=cases,
        field="p",
        quantity="max",
        window=(3.0, 10.0),
        reference_value=rayleigh_pitot(3.0, 1.4),
        reference_label="Rayleigh-Pitot M=3",
    )

    assert len(rep.stats) == 4
    assert len(rep.gcis) == 2  # two overlapping triplets

    # Per-grid statistics look right
    for ws, (label, (mean, _)) in zip(rep.stats, targets.items()):
        assert abs(ws.mean - mean) < 0.005, f"{label}: {ws.mean} vs {mean}"
        assert ws.std > 0.0
        assert ws.tau_int >= 1.0
        # The series is iid Gaussian, so τ_int should be close to 1 and the
        # KPSS statistic finite & non-negative. We deliberately do NOT assert
        # "KPSS never rejects": a 5%-level test false-rejects ~5% of the time
        # even under the null, so that assertion would be inherently flaky.
        assert ws.tau_int < 2.0, f"{label}: tau_int={ws.tau_int} too high for iid"
        assert ws.kpss_stat >= 0.0
        assert 0.01 <= ws.kpss_p <= 0.10

    # GCI sanity
    for g in rep.gcis:
        assert 0.0 < g.gci_fine_21_pct < 50.0
        assert g.p_apparent > 0.0

    # Pretty-printing must not raise
    txt = rep.as_text()
    assert "foamgci V&V report" in txt
    tex = rep.as_latex()
    assert r"\begin{table}" in tex
    assert r"\bottomrule" in tex


def test_full_report_rejects_wrong_ordering(tmp_path: Path) -> None:
    f1 = tmp_path / "a.dat"
    f2 = tmp_path / "b.dat"
    _make_synthetic_fieldminmax(f1, "p", 12.0, 0.01, seed=1)
    _make_synthetic_fieldminmax(f2, "p", 12.0, 0.01, seed=2)
    with pytest.raises(ValueError, match="coarse-to-fine"):
        full_report(
            cases=[GridCase("fine", f1, 0.001), GridCase("coarse", f2, 0.01)],
        )
