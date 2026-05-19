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

- **No surface- or volume-integrated quantities.** `fieldMinMax`
  records pointwise extrema, not surface integrals. Drag,
  lift, mass-flow, and other integrated QoIs require separate
  OpenFOAM function objects (`forces`, `volFieldValue`,
  `surfaceFieldValue`) and a separate parser. A follow-up release
  will add `volFieldValue` ingest if there is interest.

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

## Reproducibility

- Results depend on the OpenFOAM solver, scheme, time step, and
  function-object configuration in `system/controlDict`. Pin
  these explicitly in your case repository alongside the
  `foamgci report` output.
