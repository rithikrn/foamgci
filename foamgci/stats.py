"""foamgci.stats — autocorrelation-aware time-series statistics for V&V.

This module implements the **statistical** half of an unsteady-CFD
verification report:

  * :func:`tau_int_geyer` — Geyer's "initial positive sequence" estimator
    of the integrated autocorrelation time τ_int.
  * :func:`sem_autocorr_corrected` — standard error of the mean, corrected
    for serial correlation via N_eff = N / τ_int.
  * :func:`kpss_test` — KPSS test of the null hypothesis that the window
    is (level- or trend-) stationary.
  * :func:`window_stats` — convenience aggregator returning all of the
    above for a stationary window of a single time series. Resamples a
    non-uniformly spaced series onto a uniform time grid first, so that
    lag-k autocorrelation corresponds to a fixed physical time interval.

Convention
----------
The integrated autocorrelation time used throughout is the MCMC /
turbulence-V&V convention

    τ_int = 1 + 2 Σ_{k≥1} ρ_k,

so that iid samples give τ_int = 1 and an AR(1) process with parameter
φ gives τ_int = (1+φ)/(1-φ). The variance of the sample mean is
Var(x̄) = (σ²/N) τ_int, the effective sample size is N_eff = N/τ_int,
and the standard error of the mean is σ √(τ_int/N). The naïve σ/√N
therefore understates the sampling uncertainty by a factor √τ_int —
typically 2–3× for mildly under-resolved shock-dominated unsteady
flows. This is the quantity that, when properly reported, prevents
discretization increments smaller than the sampling noise from being
labelled "grid-converged".
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


# ---------------------------------------------------------------------------
# Geyer integrated autocorrelation time
# ---------------------------------------------------------------------------

def _autocov(x: np.ndarray, max_lag: int) -> np.ndarray:
    """Biased sample autocovariance from lag 0 to ``max_lag`` (inclusive)."""
    n = x.size
    xc = x - x.mean()
    # FFT-based autocorrelation is O(N log N).
    m = 1 << (2 * n - 1).bit_length()  # next power of two ≥ 2n-1
    f = np.fft.rfft(xc, n=m)
    acov_full = np.fft.irfft(f * np.conj(f), n=m)[:n] / n
    return acov_full[: max_lag + 1]


def tau_int_geyer(
    x: np.ndarray,
    max_lag: Optional[int] = None,
) -> float:
    """Integrated autocorrelation time via Geyer's initial positive sequence.

    Uses the MCMC / turbulence-V&V convention

    .. math::
        \\tau_{\\text{int}} = 1 + 2 \\sum_{k=1}^{\\infty} \\rho_k,

    where :math:`\\rho_k` is the lag-:math:`k` autocorrelation. Geyer's
    *initial positive sequence* (IPS) estimator groups the
    autocorrelations into adjacent pairs starting at lag 0,
    :math:`\\Gamma_m = \\rho_{2m} + \\rho_{2m+1}` (so
    :math:`\\Gamma_0 = 1 + \\rho_1 > 0`), and truncates the sum at the
    first :math:`m` for which :math:`\\Gamma_m` becomes non-positive.
    The estimate is then

    .. math::
        \\hat\\tau_{\\text{int}} = 2 \\sum_{m=0}^{M} \\Gamma_m - 1 .

    This pairing-from-lag-0 is the form Geyer (1992) proves monotone and
    avoids the noise sensitivity of a bare "stop at first negative lag".

    Parameters
    ----------
    x : np.ndarray, shape (N,)
        Time series (assumed stationary on its full extent, and sampled
        at a *uniform* time step — see :func:`window_stats`, which
        resamples before calling this).
    max_lag : int, optional
        Upper cap on lags searched. Defaults to ``N // 4`` and is
        rarely binding because the positive-sequence stop fires first.

    Returns
    -------
    tau : float
        τ_int ≥ 1. White noise → 1. Strong positive autocorrelation →
        arbitrarily large. As a deliberately conservative choice for
        V&V, the estimate is floored at 1.0, i.e. we never report an
        effective sample size larger than the raw sample count even if a
        noisy ACF suggests mild anti-correlation.

    Notes
    -----
    The variance of the sample mean of an autocorrelated series is
    :math:`\\sigma_{\\bar x}^2 = \\tau_{\\text{int}} \\sigma^2 / N`, so the
    effective sample size is :math:`N_{\\text{eff}} = N / \\tau_{\\text{int}}`.
    Verified analytically: iid → 1; AR(1) with parameter φ →
    :math:`(1+\\varphi)/(1-\\varphi)` (see the unit tests).
    """
    x = np.asarray(x, dtype=float).ravel()
    n = x.size
    if n < 4:
        raise ValueError(f"tau_int_geyer needs N ≥ 4, got {n}.")

    if max_lag is None:
        max_lag = max(8, n // 4)
    max_lag = min(max_lag, n - 1)
    # Need complete pairs (2m, 2m+1); make max_lag odd so the last pair
    # is well-defined.
    if max_lag % 2 == 0:
        max_lag -= 1

    acov = _autocov(x, max_lag)
    var = acov[0]
    if var <= 0.0:
        return 1.0  # constant series → iid-equivalent

    rho = acov / var
    # Geyer IPS: pair lags (0,1), (2,3), (4,5), …; stop at first
    # non-positive pair sum Γ_m = ρ_{2m} + ρ_{2m+1}.
    pair_sum = 0.0
    for m in range((max_lag + 1) // 2):
        k0, k1 = 2 * m, 2 * m + 1
        if k1 > max_lag:
            break
        gamma = rho[k0] + rho[k1]
        if gamma <= 0.0:
            break
        pair_sum += gamma
    tau = 2.0 * pair_sum - 1.0
    return float(max(tau, 1.0))  # conservative floor: never sub-iid


# ---------------------------------------------------------------------------
# Autocorrelation-corrected SEM
# ---------------------------------------------------------------------------

def sem_autocorr_corrected(
    x: np.ndarray,
    tau: Optional[float] = None,
) -> tuple[float, float]:
    """Standard error of the mean, corrected for serial correlation.

    Returns
    -------
    sem : float
        Standard error of the mean.
    n_eff : float
        Effective sample size ``N / τ_int``.

    Notes
    -----
    For a stationary time series with sample standard deviation
    :math:`s` and integrated autocorrelation time
    :math:`\\tau_{\\text{int}}` (MCMC convention, iid → 1),

    .. math::
        \\text{SEM} = s \\sqrt{\\tau_{\\text{int}} / N}.

    Setting :math:`\\tau_{\\text{int}} = 1` recovers the textbook
    :math:`s / \\sqrt{N}` (independent samples).
    """
    x = np.asarray(x, dtype=float).ravel()
    n = x.size
    if tau is None:
        tau = tau_int_geyer(x)
    s = float(np.std(x, ddof=1))
    n_eff = n / tau
    sem = s * np.sqrt(tau / n)
    return sem, n_eff


# ---------------------------------------------------------------------------
# KPSS stationarity test (level + trend variants)
# ---------------------------------------------------------------------------

# Critical values from Kwiatkowski et al. (1992), Table 1.
_KPSS_CRIT_LEVEL = {0.10: 0.347, 0.05: 0.463, 0.025: 0.574, 0.01: 0.739}
_KPSS_CRIT_TREND = {0.10: 0.119, 0.05: 0.146, 0.025: 0.176, 0.01: 0.216}


def _bartlett_lrv(e: np.ndarray, lag: int) -> float:
    """Bartlett-kernel long-run variance estimator (Newey-West weights)."""
    n = e.size
    e = e - e.mean()
    s = float(np.dot(e, e) / n)
    for k in range(1, lag + 1):
        w = 1.0 - k / (lag + 1.0)
        s += 2.0 * w * float(np.dot(e[k:], e[:-k])) / n
    return max(s, 1e-300)  # guard against numerical zero


def kpss_test(
    x: np.ndarray,
    regression: str = "c",
    nlags: str | int = "auto",
) -> dict:
    """KPSS test for stationarity.

    Null hypothesis : the series is (level- or trend-) stationary.

    Parameters
    ----------
    x : np.ndarray, shape (N,)
        Time series.
    regression : {"c", "ct"}
        ``"c"`` : test for level (mean) stationarity (regress on constant).
        ``"ct"`` : test for trend stationarity (regress on constant + time).
    nlags : "auto" | int
        Truncation lag for the Bartlett kernel long-run variance.
        ``"auto"`` uses the Schwert rule of thumb
        :math:`\\lfloor 12 \\cdot (N/100)^{1/4} \\rfloor`.

    Returns
    -------
    dict with keys
        ``statistic``  — KPSS LM statistic.
        ``p_value``    — interpolated p-value (≥ 0.10 ⇒ "≥ 0.10").
        ``lag``        — truncation lag used.
        ``critical``   — dict of critical values at 10/5/2.5/1 %.
        ``stationary_5pct`` — bool, statistic < 5 % critical value ⇒ True.

    Notes
    -----
    Implemented from first principles (no statsmodels dependency) to keep
    foamgci a single-purpose, minimum-dependency utility. Matches
    statsmodels.tsa.stattools.kpss to within numerical noise on standard
    test cases.
    """
    x = np.asarray(x, dtype=float).ravel()
    n = x.size
    if n < 10:
        raise ValueError(f"KPSS needs N ≥ 10, got {n}.")

    if regression == "c":
        e = x - x.mean()
        crit = _KPSS_CRIT_LEVEL
    elif regression == "ct":
        t = np.arange(n, dtype=float)
        # Least-squares detrending against [1, t]
        A = np.column_stack([np.ones(n), t])
        coef, *_ = np.linalg.lstsq(A, x, rcond=None)
        e = x - A @ coef
        crit = _KPSS_CRIT_TREND
    else:
        raise ValueError("regression must be 'c' or 'ct'.")

    if nlags == "auto":
        lag = int(np.floor(12.0 * (n / 100.0) ** 0.25))
    else:
        lag = int(nlags)
    lag = max(0, min(lag, n - 1))

    s = _bartlett_lrv(e, lag)
    csum = np.cumsum(e)
    eta = float(np.sum(csum * csum)) / (n * n * s)

    # Linear interpolation of p-value on the (cv, p) grid.
    cvs = np.array(sorted(crit.keys(), reverse=True))           # 0.10, 0.05, 0.025, 0.01
    cv_vals = np.array([crit[c] for c in cvs])                  # ascending in cv_vals
    if eta < cv_vals[0]:
        p = 0.10  # cap; we report it as "≥ 0.10"
    elif eta > cv_vals[-1]:
        p = 0.01
    else:
        p = float(np.interp(eta, cv_vals, cvs))

    return {
        "statistic": eta,
        "p_value": p,
        "lag": lag,
        "critical": dict(crit),
        "stationary_5pct": eta < crit[0.05],
        "regression": regression,
    }


# ---------------------------------------------------------------------------
# Window statistics aggregator
# ---------------------------------------------------------------------------

@dataclass
class WindowStats:
    """Time-averaged statistics on a stationary window of a time series."""

    n: int
    mean: float
    std: float
    tau_int: float
    sem: float
    n_eff: float
    kpss_stat: float
    kpss_p: float
    kpss_stationary_5pct: bool
    t_start: float
    t_end: float


def window_stats(
    time: np.ndarray,
    x: np.ndarray,
    t_start: float,
    t_end: float,
    kpss_regression: str = "c",
    resample: bool = True,
) -> WindowStats:
    """Compute time-averaged stats on the stationary window ``[t_start, t_end]``.

    Parameters
    ----------
    time : np.ndarray, shape (N,)
        Sample times. Need not be uniformly spaced (OpenFOAM with
        adaptive ``adjustTimeStep`` produces non-uniform output).
    x : np.ndarray, shape (N,)
        Series values, paired with ``time``.
    t_start, t_end : float
        Window bounds, inclusive.
    kpss_regression : {"c", "ct"}
        Passed to :func:`kpss_test`.
    resample : bool
        If True (default) and the in-window samples are non-uniformly
        spaced, linearly resample onto a uniform grid at the median
        cadence before computing the autocorrelation-based quantities.
        This is required for τ_int/KPSS to be meaningful: with adaptive
        Δt the time step shrinks near shock passages and shear-layer
        events, so a fixed *lag* does not map to a fixed *time* interval
        unless the series is first put on a uniform grid. The reported
        ``mean`` is always the trapezoidal time-average of the original
        samples, which is unbiased regardless of spacing.

    Returns
    -------
    WindowStats
    """
    time = np.asarray(time, dtype=float).ravel()
    x = np.asarray(x, dtype=float).ravel()
    if time.shape != x.shape:
        raise ValueError("time and x must have the same shape.")
    mask = (time >= t_start) & (time <= t_end)
    if mask.sum() < 10:
        raise ValueError(
            f"Window [{t_start}, {t_end}] contains only "
            f"{mask.sum()} samples (need ≥ 10)."
        )
    tw = time[mask]
    xw = x[mask]
    # Defensive: sort by time and drop duplicate timestamps.
    order = np.argsort(tw)
    tw, xw = tw[order], xw[order]
    uniq = np.concatenate(([True], np.diff(tw) > 0))
    tw, xw = tw[uniq], xw[uniq]
    if tw.size < 10:
        raise ValueError(
            f"Window [{t_start}, {t_end}] has only {tw.size} distinct "
            "timestamps after de-duplication (need ≥ 10)."
        )

    span = tw[-1] - tw[0]
    # Trapezoidal time-average: the physically correct ⟨x⟩ over the
    # window irrespective of sample spacing. This is the QoI fed to GCI.
    # Computed manually (np.trapz was removed in NumPy 2.0) for
    # version-independence.
    if span > 0:
        trap = float(np.sum(0.5 * (xw[:-1] + xw[1:]) * np.diff(tw)))
        mean = trap / span
    else:
        mean = float(xw.mean())

    # Build the series used for the autocorrelation-based quantities.
    dt = np.diff(tw)
    nonuniform = float(dt.std() / dt.mean()) if dt.mean() > 0 else 0.0
    if resample and nonuniform > 1e-6:
        dt_med = float(np.median(dt))
        n_uniform = max(int(np.floor(span / dt_med)) + 1, 10)
        t_uniform = tw[0] + dt_med * np.arange(n_uniform)
        t_uniform = t_uniform[t_uniform <= tw[-1] + 1e-12]
        xr = np.interp(t_uniform, tw, xw)
    else:
        xr = xw

    std = float(np.std(xr, ddof=1))
    tau = tau_int_geyer(xr)
    sem, n_eff = sem_autocorr_corrected(xr, tau=tau)
    kp = kpss_test(xr, regression=kpss_regression)
    return WindowStats(
        n=int(xr.size),
        mean=mean,
        std=std,
        tau_int=tau,
        sem=sem,
        n_eff=n_eff,
        kpss_stat=kp["statistic"],
        kpss_p=kp["p_value"],
        kpss_stationary_5pct=kp["stationary_5pct"],
        t_start=float(t_start),
        t_end=float(t_end),
    )
