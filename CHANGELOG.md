# Changelog

All notable changes to **foamgci** are documented here. Versions follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [3.2.1] — 2026-06-14

**Solver-agnostic reader foundation and reproducible example inputs.**

### Added
- Added `read_timeseries()` for generic scalar QoI histories from CSV or whitespace-delimited files.
- Exported `read_timeseries` through the public package API.
- Added reader test coverage for generic `(time, value[, x, y, z])` inputs.
- Committed the four forward-step `fieldMinMax.dat` inputs required by the worked example.

### Changed
- The reader layer now supports both OpenFOAM `fieldMinMax.dat` and solver-independent scalar time histories.
- Updated the forward-step example data status from “rerun required” to reproducible-from-clone.

## [3.2.0] — 2026-06-14

**Pressure-density multi-QoI diagnostics.**

### Added

* Added formal multi-QoI diagnostics to the Mach-3 forward-step example.
* The generated `gci_summary.json` will now contain a `qoi_results` block for:

  * `p_max` — primary maximum-pressure QoI;
  * `rho_max` — secondary maximum-density diagnostic QoI.
* Added separate Rayleigh--Pitot reference-error reporting for the
  primary pressure QoI through `error_table_vs_rayleigh_pitot`.
* Added density-specific statistics to the legacy top-level `cases`
  block:

  * `rho_max_sem`;
  * `rho_max_tau_int`;
  * `rho_max_n_eff`;
  * `rho_kpss_stat`;
  * `rho_kpss_p`;
  * `rho_kpss_stationary`.

### Changed

* Updated the forward-step statistical window to `t in [6, 10]`.
* Updated `analyze.py` and `make_figures.py` 
* Reframed the example from a pressure-only GCI study to a
  pressure-density multi-QoI verification demonstration.
* Kept `p_max` as the primary QoI 
* Reclassified `rho_max` as a secondary diagnostic QoI
  

## [3.1.1] — 2026-06-12

**Stabilisation patch — packaging and documentation only; no library changes.**

### Fixed
- The `v3.1.0` git tag pointed at a commit that predated the 3.1.0
  version bump, so checking out that tag delivered code reporting
  `__version__ = "0.3.1"`. The tag has been removed; `v3.1.1` is the
  first tag of the 3.x line guaranteed to match the code it points to.
  Pin `v3.1.1` (not `v3.1.0`) in any citation.
- READMEs now invoke the example driver as `bash run_all.sh` instead of
  `./run_all.sh`: the executable bit is not stored in the repository,
  so the documented command failed with "Permission denied" on a fresh
  clone. (`submit.sh` is unaffected — it is launched via `sbatch`.)

## [3.1.0] — 2026-06-12

**Correctness of the asymptotic diagnostic + honest reporting.**

### Fixed
- `asymptotic_ratio` was only valid for constant refinement ratio
  (r_21 = r_32): a perfectly asymptotic phi(h) = phi_e + C h^p hierarchy
  with r_21 = 1.5, r_32 = 1.333 returned R ~ 1.61 instead of 1. The
  diagnostic is now `r21^p * eps21 * (r32^p - 1) / (eps32 * (r21^p - 1))`,
  which is exact for non-constant ratios and reduces to the old formula
  when r_21 = r_32. Regression test added (Celik-style non-uniform grids).
- README: removed the false claim that the four `fieldMinMax.dat` inputs
  are committed (they are regenerated with the final workflow and will be
  committed then); fixed the illustrative output table so N scales with
  the CFL-limited sampling cadence (it showed identical N on all grids,
  which the shipped controlDict cannot produce); KPSS p shown as >=0.100.
- Example `controlDict`: added `location true` to `fieldMinMax` (README
  and analyze.py rely on extremum locations); removed stale modal-analysis
  (POD/DMD) comments left over from an unrelated campaign.

### Added
- Celik (2008) oscillatory-convergence uncertainty
  `GCIResult.u_oscillatory_pct` = half the solution span on the triplet,
  reported in the text block instead of bare NaNs.
- Reference cross-check verdict in `ReportTable.as_text()`: explicitly
  states whether |phi_ext - phi_ref| falls inside the GCI_21 band, and
  warns that GCI likely understates total uncertainty when it does not.
- KPSS p-values rendered as `>=0.100` / `<=0.010` at the clamped ends
  (text and LaTeX) so they are not mistaken for exact probabilities.
- `WindowStats.resampled` flag + docstring documenting the interpolation
  bias of uniform-grid resampling (deflated std, inflated tau_int).
- LIMITATIONS.md sections: space-time error confounding under CFL-limited
  time stepping; non-smoothness of extremum QoIs (and inviscid statistical
  non-convergence on shear-layer benchmarks); estimator caveats.
- Optional CI cross-validation of `kpss_test` against statsmodels
  (`tests/test_stats_statsmodels.py`, skipped if statsmodels absent).
- CI: ruff lint step; statsmodels installed on one matrix leg so the
  cross-validation actually runs.

## [3.0.0] — 2026-06-12

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

## [2.0.0] — 2026-05-19

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

## [1.0.0] — 2026-05-19

Initial release. Package skeleton, `reader.py`, `gci.py`, `stats.py`,
`report.py`, CLI, tests, CI, MIT license, and a worked-synthetic
forward-step example.
