# Changelog

All notable changes to **foamgci** are documented here. Versions follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [3.6.0] - 2026-07-21

**Diamond case, plus an entropy QoI that tells the truth.**

The diamond airfoil finally runs on all five grids. It adds a *volume* QoI
(entropy) the other two cases never had, and that's the one that keeps us
honest: the surface pressures match theory to ~0.01%, but the entropy integral
shows the convergence is really first order, with ~11% numerical entropy still
left at 287k cells. We report it as-is.

### Added
- `diamond2D_Ma2` example, five grids (1120 -> 286720). QoIs: `Cd`, both facet
  pressures, the entropy volume integral, and `max(p)` as a mesh-lock diagnostic.
- `compute_entropy.py` -- entropy integral straight from the native `t=10`
  cells. No in-solver function object, so nothing to compile.
- `extract_snapshot.py`, `make_contour_figure.py`, `make_aux_figures.py`.
- Figures now also write `.png`, straight from matplotlib (no `pdftoppm`).

### Changed
- Diamond `controlDict` drops the coded entropy FO (didn't build on ESI v2306)
  and sets `rho rho` on `forceCoeffs`, so the two drag routes agree to ~1e-5%.

## [3.5.0] - 2026-06-23

**Honesty check additions.**

Two new honesty checks and a couple of plots, living in the example drivers
for now. The wedge ran on all four grids and the numbers killed our old
`max(p)` story, so the tooling learned to catch that on its own.

### Added
- `order_consistency` — does the apparent order `p̂` hold steady across
  overlapping triplets? The asymptotic-range signal a single triplet never
  gave us, reported as necessary-not-sufficient.
- `extremum_location_drift` — does an extremum's *location* settle, or just
  march with the mesh? The wedge `max(p)` sits a constant ~13 cells from the
  ramp tip, so it's flagged mesh-locked even though each grid alone looks
  "localized."
- A `fig_diagnostics` plot in both examples: convergence, error-vs-`h`, and a
  location-drift panel that makes the mesh-locking obvious.

### Changed
- `R_asym` reads as `R_fit` everywhere a human sees it, labelled "not an
  asymptotic-range test" (it's ≈1 by construction). API field name unchanged.
- The wedge secondary QoI is relabelled honestly: `max(p)` tracks a near-tip
  mesh feature ~42% above the analytical `p2`, not a plateau pressure.
- Both examples print and store the new verdicts per QoI. On the forward step
  this flags `rho_max` as order-drifting (p̂ 1.49 → 2.53).

## [3.4.0] - 2026-06-21

**Wedge15Ma5 case completed and smoke-tested; first-class surface-region area-average reader.**

This release makes the `wedge15Ma5` example actually run. The 3.3.0 commit
shipped a wedge `analyze.py` that imported `read_surface_field_value` from
`foamgci.reader`, but that function did not exist, the case failed at import.
3.4.0 implements the reader, adds the `fieldMinMax` output alongside the
the area-average output, and verifies the whole pipeline end to end.

### Added
- `foamgci.reader.read_surface_field_value()`, a first-class reader for
  the OpenFOAM surface-region area-average output (`surfaceRegion.dat` in
  OpenFOAM-4.x, `surfaceFieldValue.dat` in v5.0+; integrated surface
  functionals such as `areaAverage(p)`). Auto-detects the value column from the self-describing
  `# Time ...` header; selects by field token, full label, or integer index;
  raises clearly on ambiguity, missing columns, or missing files. Exported,
  with a `SurfaceFieldValueData` container, through the public API.
- Reader test coverage for the area-average reader (header auto-detect, multi-column
  selection, ambiguity and range guards, headerless single column, restrict).
- The `wedge15Ma5` example now reads **two** OpenFOAM outputs per grid and
  reports two QoIs: a primary, reference-anchored ramp wall pressure from
  the `surfaceRegion` area-average (compared to exact oblique-shock `p2/p1`) and a secondary
  `fieldMinMax` `max(p)` carried for cross-case consistency. The secondary QoI
  reports the in-window extremum wander and a localization verdict, surfacing
  the contrast between a well-posed surface integral and a degenerate pointwise
  extremum that target the same physical pressure.

### Changed
- `wedge15Ma5/system/controlDict` writes both a `surfaceRegion`
  (`wallPressure`) and a `fieldMinMax` function object, at a matched sampling
  cadence (`writeInterval 20`).
- `wedge15Ma5/gci/data.py` now carries two input filenames per grid; the
  analysis guards on all eight files and exits with a clear message (and a
  non-zero status, so `run_all.sh` aborts) if any is missing.
- The `wedge15Ma5` README and `gci/data/README.md` document the dual-output
  design and the eight-file input set.

## [3.3.0] - 2026-06-21: superseded by 3.4.0 (do not use)

**Added the wedge15Ma5 case as a work-in-progress.** This commit was never
release-validated: the wedge `analyze.py` referenced a `read_surface_field_value`
reader that had not been committed, so the example raised `ImportError` on run.
Superseded by 3.4.0, which implements the reader and smoke-tests the case.

### Added
- New validation example: `wedge15Ma5` (Mach 5, 15-degree wedge oblique shock).
- Demonstration of integrated/smoother Quantities of Interest (e.g., surface forces or averaged pressures) to provide a formally cleaner target for GCI than pointwise extrema.


## [3.2.2]: 2026-06-20

**Zenodo release**

### Added
- The `v3.2.2` is just for zenodo release
- Added a citation file

## [3.2.1]: 2026-06-14

**Solver-agnostic reader foundation and reproducible example inputs.**

### Added
- Added `read_timeseries()` for generic scalar QoI histories from CSV or whitespace-delimited files.
- Exported `read_timeseries` through the public package API.
- Added reader test coverage for generic `(time, value[, x, y, z])` inputs.
- Committed the four forward-step `fieldMinMax.dat` inputs required by the worked example.

### Changed
- The reader layer now supports both OpenFOAM `fieldMinMax.dat` and solver-independent scalar time histories.
- Updated the forward-step example data status from “rerun required” to reproducible-from-clone.

## [3.2.0]: 2026-06-14

**Pressure-density multi-QoI diagnostics.**

### Added

* Added formal multi-QoI diagnostics to the Mach-3 forward-step example.
* The generated `gci_summary.json` will now contain a `qoi_results` block for:

  * `p_max`, primary maximum-pressure QoI;
  * `rho_max`, secondary maximum-density diagnostic QoI.
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
  

## [3.1.1]: 2026-06-12

**Stabilisation patch, packaging and documentation only; no library changes.**

### Fixed
- The `v3.1.0` git tag pointed at a commit that predated the 3.1.0
  version bump, so checking out that tag delivered code reporting
  `__version__ = "0.3.1"`. The tag has been removed; `v3.1.1` is the
  first tag of the 3.x line guaranteed to match the code it points to.
  Pin `v3.1.1` (not `v3.1.0`) in any citation.
- READMEs now invoke the example driver as `bash run_all.sh` instead of
  `./run_all.sh`: the executable bit is not stored in the repository,
  so the documented command failed with "Permission denied" on a fresh
  clone. (`submit.sh` is unaffected, it is launched via `sbatch`.)

## [3.1.0]: 2026-06-12

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

## [3.0.0]: 2026-06-12

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

## [2.0.0]: 2026-05-19

**Reproducibility release.** Reorients the package toward an
end-to-end V&V workflow that ingests the user's actual OpenFOAM
output and emits a paper-ready Table 1 directly.

### Added
- `foamgci.report.full_report`, single-call orchestrator that
  reads N `fieldMinMax.dat` files, computes window statistics on
  a stationary window, applies Roache GCI on every consecutive
  triplet, and returns a `ReportTable` with both plain-text and
  LaTeX renderings.
- `foamgci.report.rayleigh_pitot`, analytical Mach-3 normal-shock
  stagnation pressure ratio for cross-checking Richardson
  extrapolation independently of the GCI machinery.
- `foamgci.stats.kpss_test`, implemented from first principles
  (Bartlett-kernel long-run variance, Kwiatkowski et al. 1992
  critical values). No `statsmodels` dependency.
- `foamgci.plot.plot_convergence`, optional two-panel
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

## [1.0.0]: 2026-05-19

Initial release. Package skeleton, `reader.py`, `gci.py`, `stats.py`,
`report.py`, CLI, tests, CI, MIT license, and a worked-synthetic
forward-step example.
