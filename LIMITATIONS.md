# Limitations

`foamgci` is intentionally a small, single-purpose tool. The list
below is what it does **not** do, so prospective users can decide
quickly whether it fits their workflow.

## Scope

- **Single scalar QoI per call.** `foamgci report` operates on the
  time-series of one extremum (min or max over space) of one field
  per invocation. Vector and tensor fields are handled by component
  selection at the OpenFOAM side (e.g. write `mag(U)` rather than
  `U`). For coupled QoIs across fields, invoke the report API
  multiple times and combine results yourself.

- **Surface area-averages, but not yet volume integrals or forces.**
  `read_surface_field_value` ingests the `surfaceRegion` /
  `surfaceFieldValue` area-average output (used by the wedge example), so
  integrated surface functionals like `areaAverage(p)` are supported.
  Volume integrals (`volFieldValue`) and force coefficients (`forces`) still
  need their own parsers; a follow-up release will add `volFieldValue` ingest
  if there is interest.

- **No spatial discretisation diagnostics.** The package assumes the
  user has already produced converged time-histories on each mesh
  in a refinement hierarchy. It does not run OpenFOAM, does not
  inspect `checkMesh` logs, and does not compute mesh-quality
  metrics.

- **No iterative-convergence verification.** Roache GCI quantifies
  *discretization* uncertainty only. Iterative-solver convergence
  (residual drop) and round-off uncertainty must be assessed
  separately.

## Statistical assumptions

- **The stationarity window is user-supplied.** `foamgci` does not
  detect when the transient ends; it tests the user-supplied
  `[t0, t1]` window with KPSS. If the test rejects, narrow the
  window manually and re-run.

- **Geyer's τ_int is a conservative estimator.** It can over-
  estimate τ for short series with sharp tails in the
  autocorrelation function. For unambiguous critical applications,
  cross-check with batch-means or block-bootstrap on the same data.

- **Single-window analysis.** No multi-window or sliding-window
  variants are implemented. Users running pre- and post-transient
  comparisons should invoke the report twice with different
  `--window` arguments.

## Refinement-study assumptions

- **Quasi-uniform refinement.** The apparent-order solve handles
  non-uniform refinement ratios formally but degrades in accuracy
  when ratios deviate strongly from constancy. Aim for `r_21 ≈ r_32`
  whenever possible.

- **Three- or four-grid hierarchies.** GCI estimates from fewer
  than three grids are not produced. Studies on a single grid get
  the time-statistics block but no GCI. Hierarchies with more than
  four grids work but only consecutive triplets are reported.

- **Asymptotic-range diagnostic, not enforcement.** The asymptotic
  ratio `R_asym` is reported; the package does **not** refuse to
  emit a GCI when `R_asym` strays from unity. Interpretation is
  the user's responsibility.

## Space-time error confounding (important)

- **CFL-limited time stepping couples Δt to h.** With
  `adjustTimeStep yes` and fixed `maxCo`, halving the mesh spacing
  roughly halves the time step. The "apparent order" from the
  refinement hierarchy therefore reflects a **mixture of spatial and
  temporal discretization error** (rhoCentralFoam time integration is
  first-order). Roache GCI formally quantifies spatial error only.
  Mitigation: run at least one grid pair at a fixed, conservative Δt
  (small enough for the finest grid) and compare the apparent order
  against the CFL-limited result. `foamgci` does not do this for you
  and does not warn about it.

## QoI smoothness

- **A spatial extremum is a non-smooth functional.** The Richardson
  expansion `phi(h) = phi_exact + C h^p + ...` assumed by GCI has no
  formal basis for the max/min of a shock-captured field, whose
  location can jump between cells under refinement. Treat
  `fieldMinMax`-based GCI as a heuristic indicator. An integrated QoI
  (a force, a surface or volume average) is better posed: it has a
  definite continuum value and is insensitive to per-cell grid noise, so the
  Richardson expansion at least applies to a well-defined limit. But its
  *order* of convergence is still capped by the shock-capturing scheme
  (first order at captured discontinuities), so do not expect the formal
  second order even for a smooth functional, measure the apparent order,
  do not assume it. For inviscid shock/shear-layer benchmarks
  (e.g. the forward-facing step), statistics may not converge under
  refinement at all, because resolved Kelvin-Helmholtz structure on
  slip lines grows without a viscous cutoff, a divergent or
  oscillatory regime flag on the coarsest triplet can be physics, not
  a bug.

## Statistical estimator caveats

- **Resampling bias.** Non-uniformly sampled windows are linearly
  interpolated onto a uniform grid before tau_int/KPSS. Interpolation
  low-pass filters the series: it slightly deflates the sample std and
  inflates serial correlation. The net SEM bias is small for mildly
  non-uniform spacing. `WindowStats.resampled` records whether
  resampling occurred.

- **Geyer IPS is proven monotone for reversible Markov chains.** For
  general stationary time series (CFD signals) it is a well-behaved
  heuristic, not a theorem.

- **KPSS p-values are clamped.** Only the published critical values at
  10/5/2.5/1 % are tabulated, so reported p-values are linear
  interpolations clamped to [0.01, 0.10] and rendered as `>=0.100` /
  `<=0.010` at the ends. Do not propagate them as exact probabilities.

## Reproducibility

- Results depend on the OpenFOAM solver, scheme, time step, and
  function-object configuration in `system/controlDict`. Pin
  these explicitly in your case repository alongside the
  `foamgci report` output.
