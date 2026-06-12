# Changelog

All notable changes to **foamgci** are documented here. Versions follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.0] — 2026-06-12

**Correctness + single-source-of-truth release.**

### Fixed
- `tau_int_geyer`: MCMC convention (τ=1+2Σρ, iid→1, AR(1)→(1+φ)/(1−φ))
  with true Geyer lag-0 pairing; SEM/N_eff and docstrings reconciled.
- `window_stats`: resamples non-uniformly spaced (adaptive-Δt) series
  onto a uniform grid before τ_int/KPSS, and reports a trapezoidal
  time-average.
- `roache_gci`: classifies the convergence regime (monotonic /
  oscillatory / divergent / degenerate) and returns NaN GCI outside the
  asymptotic range instead of a misleading finite value.

### Changed
- The forward-step example runs entirely through `foamgci`
  (`gci/analyze.py`); the duplicate standalone `gci/gci.py` and
  `gci/parse_fieldminmax.py` are removed. The example now reports
  τ_int, SEM, N_eff and KPSS, and frames Rayleigh–Pitot as a lower-bound
  check rather than the convergence target.
- `submit.sh` and `gci/run_all.sh` reduced to the essential steps.

## [0.2.0] — 2026-05-19

**Reproducibility release.** Reorients the package toward an
end-to-end V&V workflow that ingests the user's actual OpenFOAM
output and emits a paper-ready Table 1 directly.

### Added
- `foamgci.report.full_report` — single-call orchestrator that
  reads N `fieldMinMax.dat` files, computes window statistics on
  a stationary window, applies Roache GCI on every consecutive
  triplet, and returns a `ReportTable` with both plain-text and
  LaTeX renderings.
- `foamgci.report.rayleigh_pitot` — analytical Mach-3 normal-shock
  stagnation pressure ratio for cross-checking Richardson
  extrapolation independently of the GCI machinery.
- `foamgci.stats.kpss_test` — implemented from first principles
  (Bartlett-kernel long-run variance, Kwiatkowski et al. 1992
  critical values). No `statsmodels` dependency.
- `foamgci.plot.plot_convergence` — optional two-panel
  convergence figure. Pulls matplotlib only when invoked.

### Changed
- `foamgci.reader.read_fieldminmax` rewritten to handle both
  OpenFOAM `fieldMinMax` dialects (combined multi-field and
  per-field) with `field=` selection.
- `foamgci.gci.roache_gci` switched to the Celik et al. (2008)
  iterative apparent-order solve (handles non-uniform refinement
  ratios) and now reports the asymptotic ratio R as a
  diagnostic in the same call.

### Removed
- Synthetic-data scaffolding in `examples/`. The example now points
  at the user's real cases and ships only the driver script,
  not fake data.

### Verified
- Test suite anchored to Celik et al. (2008) Example 1: apparent
  order 1.53, GCI_21 ≈ 2.2 %, GCI_32 ≈ 5.71 %.
- Synthetic exact-second-order series recovers p̂ = 2.000 and
  R_asymptotic = 1.000.
- KPSS rejects a random walk and accepts iid Gaussian noise; the
  trend variant accepts a deterministic linear trend + noise.

## [0.1.0] — 2026-05-19

Initial release. Package skeleton, `reader.py`, `gci.py`, `stats.py`,
`report.py`, CLI, tests, CI, MIT license, and a worked-synthetic
forward-step example.
