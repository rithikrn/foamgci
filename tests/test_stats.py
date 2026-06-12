"""Tests for foamgci.stats."""
from __future__ import annotations

import numpy as np
import pytest

from foamgci.stats import (
    tau_int_geyer,
    sem_autocorr_corrected,
    kpss_test,
    window_stats,
)


def _ar1(n: int, phi: float, sigma: float = 1.0, rng=None) -> np.ndarray:
    """Simulate an AR(1) process x_t = phi * x_{t-1} + sigma*eps_t."""
    rng = rng or np.random.default_rng(0)
    eps = rng.standard_normal(n) * sigma
    x = np.empty(n)
    x[0] = eps[0] / np.sqrt(1.0 - phi * phi)
    for t in range(1, n):
        x[t] = phi * x[t - 1] + eps[t]
    return x


# ---- τ_int ---------------------------------------------------------------

def test_tau_int_white_noise() -> None:
    rng = np.random.default_rng(42)
    x = rng.standard_normal(20_000)
    tau = tau_int_geyer(x)
    assert 0.9 < tau < 1.3


def test_tau_int_ar1_matches_analytical() -> None:
    for phi in (0.5, 0.8):
        tau_true = (1.0 + phi) / (1.0 - phi)      # 3.0 at 0.5, 9.0 at 0.8
        x = _ar1(200_000, phi, rng=np.random.default_rng(0))
        tau = tau_int_geyer(x)
        assert tau >= 1.0
        assert abs(tau - tau_true) / tau_true < 0.10, (phi, tau, tau_true)


def test_tau_int_short_series_raises() -> None:
    with pytest.raises(ValueError, match="N ≥ 4"):
        tau_int_geyer(np.array([1.0, 2.0]))


# ---- SEM correction ------------------------------------------------------

def test_sem_correction_inflates_under_correlation() -> None:
    rng = np.random.default_rng(1)
    x = rng.standard_normal(5000)
    sem_naive, n_eff_naive = sem_autocorr_corrected(x, tau=1.0)
    sigma = np.std(x, ddof=1)
    np.testing.assert_allclose(sem_naive, sigma / np.sqrt(x.size), rtol=1e-12)
    np.testing.assert_allclose(n_eff_naive, x.size, rtol=1e-12)

    sem_inflated, n_eff_inflated = sem_autocorr_corrected(x, tau=4.0)
    assert sem_inflated > sem_naive
    assert n_eff_inflated < n_eff_naive  
    np.testing.assert_allclose(sem_inflated / sem_naive, 2.0, rtol=1e-12)

def test_sem_correction_inflates_under_correlation() -> None:
    rng = np.random.default_rng(1)
    x = rng.standard_normal(5000)
    sem_naive, n_eff_naive = sem_autocorr_corrected(x, tau=1.0)
    sigma = np.std(x, ddof=1)
    np.testing.assert_allclose(sem_naive, sigma / np.sqrt(x.size), rtol=1e-12)
    np.testing.assert_allclose(n_eff_naive, x.size, rtol=1e-12)

    sem_inflated, n_eff_inflated = sem_autocorr_corrected(x, tau=4.0)
    assert sem_inflated > sem_naive
    assert n_eff_inflated < n_eff_naive
    np.testing.assert_allclose(sem_inflated / sem_naive, 2.0, rtol=1e-12)


# ---- KPSS ----------------------------------------------------------------

def test_kpss_stationary_white_noise() -> None:
    rng = np.random.default_rng(0)
    x = rng.standard_normal(2000)
    out = kpss_test(x, regression="c")
    assert out["stationary_5pct"] is True
    assert out["p_value"] >= 0.05

    rejects = 0
    for s in range(50):
        rng_s = np.random.default_rng(1000 + s)
        x_s = rng_s.standard_normal(2000)
        if not kpss_test(x_s, regression="c")["stationary_5pct"]:
            rejects += 1
    assert rejects / 50 < 0.20


def test_kpss_rejects_random_walk() -> None:
    """A random walk is *not* level-stationary; KPSS must reject."""
    rng = np.random.default_rng(11)
    rw = np.cumsum(rng.standard_normal(2000))
    out = kpss_test(rw, regression="c")
    assert out["stationary_5pct"] is False
    assert out["p_value"] < 0.05


def test_kpss_trend_stationary_under_ct() -> None:
    """A deterministic linear trend + noise IS trend-stationary,
    so the trend-variant (ct) should NOT reject — even though the
    level variant (c) would."""
    rng = np.random.default_rng(2)
    n = 2000
    x = 0.01 * np.arange(n) + rng.standard_normal(n)
    out_c = kpss_test(x, regression="c")
    out_ct = kpss_test(x, regression="ct")
    assert out_c["stationary_5pct"] is False
    assert out_ct["stationary_5pct"] is True


# ---- window_stats aggregator --------------------------------------------

def test_window_stats_end_to_end() -> None:
    rng = np.random.default_rng(0)
    n = 1000
    t = np.linspace(0, 10, n)
    x = 5.0 + 0.1 * rng.standard_normal(n)
    ws = window_stats(t, x, 3.0, 10.0)
    assert ws.n > 600  # most of the series
    assert 4.9 < ws.mean < 5.1
    assert ws.std > 0.0
    assert ws.tau_int >= 0.5
    assert ws.sem > 0.0
    assert ws.kpss_stationary_5pct is True


def test_window_stats_too_few_samples_raises() -> None:
    t = np.linspace(0, 10, 1000)
    x = np.zeros_like(t)
    with pytest.raises(ValueError, match="contains only"):
        window_stats(t, x, 9.999, 10.0)
