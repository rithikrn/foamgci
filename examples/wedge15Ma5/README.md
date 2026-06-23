# Example: Mach-5 oblique shock over a 15-degree wedge

Grid-convergence verification on a four-grid hierarchy, checked against exact
oblique-shock theory and analysed with the `foamgci` library. This is the
second example. It reuses the library unchanged and reads **two** OpenFOAM
outputs per grid, so it doubles as the portability demonstration: the same GCI
engine handles a new data source (the `surfaceRegion` area-average) *and* the one the
forward step already used (`fieldMinMax`), in a single analysis.

- **Primary QoI** (`surfaceRegion` area-average): the area-averaged ramp-surface
  pressure, an integrated functional, reference-anchored against the exact
  post-shock `p2/p1`. This is a *new* output type vs the forward step and
  carries the verification verdict.
- **Secondary QoI** (`fieldMinMax`): the global `max(p)`, the *same* output
  type the forward step used, kept here for cross-case consistency. Its value
  also tends to `p2`, but the post-shock plateau makes the extremum location
  degenerate, so foamgci's localization check flags it. The contrast between a
  well-posed surface integral and a degenerate pointwise extremum that target
  the same physical pressure is the methodological point.

## What lives where

- `0/`, `constant/`, `system/`: the committed OpenFOAM case (this copy is the
  **fine** grid). These are the OpenFOAM-4.x `wedge15Ma5` tutorial files,
  unchanged except for two case-specific additions in `system/`: two function
  objects in `controlDict` (an area-average on the ramp surface, `surfaceRegion`
  in OpenFOAM-4.x, named `wallPressure`, and a `fieldMinMax`) and a
  `decomposeParDict`.
- `submit.sh`: SLURM runner, identical to the forward-step example. Mesh,
  decompose, run, reconstruct.
- `gci/`: the analysis driver for THIS case. It `import`s `foamgci` (the
  library) and applies it here. You edit `gci/data.py`; you do not touch the
  library.
- `gci/oblique_shock.py`: the analytical reference. An exact theta-beta-M solve
  for the post-shock state, NumPy only.
- `gci/data/`: the **eight** inputs the analysis reads, a
  `surfaceRegion.dat` (the area-average) and a `fieldMinMax.dat` per grid. They
  ship empty; produce them by running the four cases (see `gci/data/README.md`).

## Physics

Inviscid (Euler) Mach-5 flow over a 15-degree compression ramp, solved with
`rhoCentralFoam` (Kurganov fluxes, vanLeer limiters, fixed time step). Target
solver version: OpenFOAM-4.x, which is where these tutorial files come from. The
only version-sensitive piece is the ramp-pressure function object, `surfaceRegion`
in 4.x, `surfaceFieldValue` in v5.0+; the controlDict comment notes the swap.
Normalised perfect gas: a_inf = 1 at T = 1, so M = U = 5, p_inf = 1. The
`obstacle` patch is the ramp surface; its slope, atan(0.08167/0.3048), is 15
degrees exactly, so an attached oblique shock forms at the ramp foot. QoI:
the area-averaged static pressure on the ramp, time-averaged over the steady
window t in [0.15, 0.20].

**Reference value.** For inviscid flow the ramp-surface pressure equals the
post-shock static pressure p2 along the straight ramp, so the area-averaged
wall pressure converges to the analytical oblique-shock value as h goes to 0.
For M1 = 5, theta = 15 degrees, gamma = 1.4 the exact solution
(`gci/oblique_shock.py`) gives shock angle beta = 24.322 degrees, p2/p1 =
4.7808, rho2/rho1 = 2.7535, T2/T1 = 1.7363, M2 = 3.504. The shock leaves
through the outlet at y = 0.138, below the top boundary at 0.152, so it is a
single oblique shock with no reflection inside the domain. These numbers match
NACA Report 1135 (1953) and any compressible-flow text. The reference is exact
and independent of mesh, solver version, and time step, which makes it a
stronger target than another code's result.

## The four grids

The committed case is the **fine** grid. The four grids differ ONLY in the
two block counts in `system/blockMeshDict`; everything else is identical. The
mesh is uniform with square cells, so h = 0.1524 / Ny.

| label                  | cells  | h         | block 1     | block 2      |
|------------------------|-------:|-----------|-------------|--------------|
| coarse                 |  1,200 | 0.0076200 | (20 20 1)   | (40 20 1)    |
| medium                 |  4,800 | 0.0038100 | (40 40 1)   | (80 40 1)    |
| **fine** (committed)   | 19,200 | 0.0019050 | (80 80 1)   | (160 80 1)   |
| extra-fine             | 76,800 | 0.0009525 | (160 160 1) | (320 160 1)  |

i.e. in the `blocks (...)` list, set the two `(nx ny 1)` entries to the row for
that grid.

## Run from scratch (on HPC)

1. Make the four case directories as siblings of `gci/`:
   ```bash
   cd examples/wedge15Ma5
   for d in coarse_grid medium_grid fine_grid extrafine_grid; do
     mkdir -p $d && cp -r 0 constant system submit.sh $d/
   done
   # then edit each $d/system/blockMeshDict block counts per the table above
   ```
2. Submit each. `submit.sh` runs blockMesh -> checkMesh -> decomposePar ->
   rhoCentralFoam -> reconstructPar:
   ```bash
   for d in coarse_grid medium_grid fine_grid extrafine_grid; do
     ( cd $d && sbatch submit.sh ); done
   ```
   (If you run interactively instead of via SLURM, copy the files by hand.)
3. Copy each run's two function-object outputs into `gci/data/`, renamed for
   the grid (example for the coarse run):
   ```bash
   cp postProcessing/wallPressure/0/surfaceRegion.dat  coarse_surfaceFieldValue.dat
   cp postProcessing/fieldMinMax/0/fieldMinMax.dat         coarse_fieldMinMax.dat
   ```
   Do the same for `medium_`, `fine_`, `extrafine_` (eight files total).
4. Analyse:
   ```bash
   pip install -e ../..            # install the foamgci library (once)
   cd gci && bash run_all.sh
   ```
   Writes `gci/gci_summary.json`: per-grid mean, sigma, tau_int, SEM, N_eff,
   KPSS; regime-aware GCI on both triplets; and the error against the
   analytical reference.
5. Commit the eight small `gci/data/*.dat` files. That is what makes the study
   reproducible for the next reader, who then runs step 4 alone.

## Reading the output

- Per grid: `mean` (the time-averaged wall pressure, equal to p2/p1 here since
  p_inf = 1), the autocorrelation-corrected `sem`, `tau_int`, and `n_eff`.
- The case is steady, so the wall-pressure signal is nearly constant in the
  window. Expect tau_int near 1, KPSS stationary, and `sem` close to the naive
  sigma/sqrt(N). The KPSS and tau_int path runs and should pass, but it does
  not change any verdict here. A case where the autocorrelation correction
  changes a conclusion needs a stationary unsteady QoI; this steady case is
  not that one. The forward step and a future vortex-shedding case carry that
  load.
- `triplet_A_CMF` and `triplet_B_MFXF` are the two consecutive triplets; the
  deeper one (B) supplies `phi_star`.
- `reference.p2_p1` (4.7808) is the verification anchor. The finest-grid `mean`
  should approach it, and `phi_star_rel_err_pct` reports the Richardson
  extrapolate's distance from it. `reference_covered_by_gci` is true when that
  distance falls inside the fine-grid GCI band. Unlike the forward step's
  pointwise pressure maximum, this integrated QoI is smooth, so the Richardson
  extrapolate is meaningful rather than a diagnostic.
- If KPSS rejects stationarity on [0.15, 0.20] for any grid, widen the run or
  narrow the window in `gci/data.py` (`T_STAT`) and re-run `analyze.py`. Note
  that `endTime = 0.2` is about 2.2 domain flow-through times (the window starts
  near 1.6), which is on the short side for a startup transient to fully wash
  out. The KPSS check is what tells you whether it has; if it has not, the
  simplest fix is to raise `endTime` and rerun.
- The **secondary** `max(p)` block lives under `secondary_qoi_fieldminmax`.
  Read it as a diagnostic, not a verdict: its `localization.localized` is
  expected to be `false` (the plateau makes the extremum location wander past
  the threshold). That `false` is the *intended* result here, it shows the
  toolkit distinguishing a degenerate pointwise extremum from the well-posed
  surface integral, both of which numerically approach the same `p2`.

This case carries two QoIs from two different OpenFOAM outputs:

| QoI            | Source         | Role                 | Interpretation                                                                                                                                                                                                 |
| -------------- | -------------- | -------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `p_wall_ratio` | `surfaceRegion`| Primary verified QoI | Area-averaged ramp-surface pressure: an integrated functional with a definite continuum value, insensitive to per-cell grid noise. That makes it a better-posed Richardson/GCI target than a pointwise extremum. Its observed order is still set by the shock-capturing scheme (see "Expected order" below), so it is measured, not assumed. |
| `p_max`        | `fieldMinMax`  | Secondary diagnostic | Global `max(p)`. The value tends to p2 as well, but the post-shock plateau makes the extremum location degenerate, so the localization (wander) check flags it. Reported as a diagnostic, and as the `fieldMinMax` output shared with the forward-step example. |

Expected order. `rhoCentralFoam` here uses Kurganov fluxes with vanLeer
reconstruction (second order in smooth regions, first order at the shock) and
Euler time stepping (first order). For a shock-set quantity like the post-shock
pressure, the observed order `p_obs` will realistically sit between 1 and ~1.5,
not 2. The GCI-with-measured-order approach is the right tool precisely because
the order cannot be assumed; the analysis reports the `p_obs` it finds and
classifies the convergence regime.

The analysis driver:

1. Parse the `surfaceRegion` area-average (primary) **and** `fieldMinMax`
   (secondary), two
   OpenFOAM sources, both ingested by the same foamgci readers.
2. Check stationarity using KPSS on each QoI's window.
3. Estimate autocorrelation-corrected SEM and effective sample size.
4. Compute GCI on every consecutive triplet, for each QoI.
5. Extrapolate to h -> 0 and compare against exact oblique-shock theory; for
   the secondary QoI, also report the in-window extremum wander and the
   localization verdict.

### JSON summary structure

`gci_summary.json` contains:

* `qoi`: name, description, and the OpenFOAM source of the **primary** quantity;
* `stationary_window`: the averaging window used for the QoI;
* `reference`: the analytical oblique-shock state (beta, p2/p1, rho2/rho1,
  T2/T1, M2);
* `cases`: per-grid mean, sigma, tau_int, SEM, N_eff, KPSS (primary QoI);
* `triplet_A_CMF`, `triplet_B_MFXF`: regime-aware GCI on both triplets;
* `phi_star`, `phi_star_source`, `phi_star_rel_err_pct`: the extrapolate and
  its error against the reference;
* `reference_covered_by_gci`: whether the reference falls inside the
  fine-grid GCI band;
* `error_table_vs_reference`: per-grid error against p2/p1;
* `secondary_qoi_fieldminmax`: the full GCI block for the `fieldMinMax`
  `max(p)` QoI (its own `cases`, triplets, `phi_star`, error table) plus a
  `localization` sub-object (`median_wander_cells`, `wander_threshold_cells`,
  `localized`, `note`);
* `inputs_per_grid`: a record of the two input files consumed per grid.

The top-level keys describe the primary QoI and are kept stable so the figure
scripts read them directly; the secondary QoI lives entirely under
`secondary_qoi_fieldminmax`.

## Case directory tree

```
wedge15Ma5/
├── 0/                          # initial + boundary conditions (upstream tutorial)
│   ├── T
│   ├── U
│   └── p
├── constant/
│   ├── thermophysicalProperties
│   └── turbulenceProperties
├── system/
│   ├── blockMeshDict           # edit the two block counts per grid (see table)
│   ├── controlDict             # + surfaceRegion AND fieldMinMax (two QoI sources)
│   ├── decomposeParDict
│   ├── fvSchemes
│   └── fvSolution
├── gci/                        # per-case driver; imports foamgci, never edits it
│   ├── analyze.py              # reads BOTH outputs -> GCI(both QoIs) -> vs analytical
│   ├── data.py                 # grid metadata (two files/grid), window, free-stream
│   ├── oblique_shock.py        # exact theta-beta-M reference (NumPy only)
│   ├── run_all.sh
│   └── data/                   # eight inputs (surfaceRegion + fieldMinMax per grid; ship empty)
│       └── README.md
├── README.md
└── submit.sh                   # identical to the forward-step runner
```

## References

- Anderson, J.D. (2003), *Modern Compressible Flow*, 3rd ed., McGraw-Hill: oblique-shock relations.
- Ames Research Staff (1953), NACA Report 1135, *Equations, Tables, and Charts for Compressible Flow*: the tabulated reference.
- Greenshields, C.J. et al. (2010), *Int. J. Numer. Meth. Fluids* **63**(1), 1–21: the rhoCentralFoam central-upwind solver.
- Kurganov, A., Noelle, S., Petrova, G. (2001), *SIAM J. Sci. Comput.* **23**(3), 707–740: the KNP flux scheme.
