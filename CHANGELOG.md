# Changelog

All notable changes to **foamgci** are documented here. Versions follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.0] — 2026-06-12

**Correctness + single-source-of-truth release.**

### Fixed
- Version is now defined once in `foamgci/_version.py` and read from
  there by the package, the CLI `--version`, the LaTeX header, and
  `pyproject.toml` (`dynamic = ["version"]`). No more drift.
- CI workflow moved to `.github/workflows/tests.yml` (it was at the repo
  root and never ran). Suite now genuinely runs on every push/PR.
- `tests/test_report.py` made deterministic: replaced `hash(label)`
  seeding (Python salts str hashing per process) with fixed seeds, and
  replaced the statistically-fragile "KPSS never rejects" assertion with
  stable bounds. Passes under any `PYTHONHASHSEED`.
- README: corrected the SEM formula to `sigma*sqrt(tau_int/N)` (matching
  the code and the `tau_int = 1 + 2*sum(rho)` convention), removed the
  broken `--coarse/--medium` CLI example (the CLI is `report --case ...`),
  and rewrote the stale repository-layout block.
- `tau_int_geyer`: MCMC convention (tau=1+2*sum(rho), iid->1,
  AR(1)->(1+phi)/(1-phi)) with true Geyer lag-0 pairing; SEM/N_eff and
  docstrings reconciled.
- `window_stats`: resamples non-uniformly spaced (adaptive-dt) series
  onto a uniform grid before tau_int/KPSS, and reports a trapezoidal
  time-average.
- `roache_gci`: classifies the convergence regime (monotonic /
  oscillatory / divergent / degenerate) and returns NaN GCI outside the
  asymptotic range instead of a misleading finite value.

### Changed
- One canonical function-object name everywhere: `fieldMinMax`, writing
  to `postProcessing/fieldMinMax/0/fieldMinMax.dat` (previously the repo
  mixed `fieldMinMax`, `fieldRange`, and `fieldMinMax1`).
- The forward-step example reads four committed KB-sized inputs from
  `gci/data/` instead of from gitignored full-run directories, so the
  figures reproduce from a clean clone. `submit.sh` copies each grid's
  output into `gci/data/` automatically and now also runs `checkMesh`.
- The example `controlDict` is trimmed to the single `fieldMinMax`
  function object needed for GCI; the modal-analysis field writes are
  commented out as clearly-labelled optional.
- The example runs entirely through `foamgci` (`gci/analyze.py`); the old
  duplicate `gci/gci.py` / `gci/parse_fieldminmax.py` and the stale
  `examples/.../data/` folder are gone.
- `gci_over_hierarchy` is now part of the public API.

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
