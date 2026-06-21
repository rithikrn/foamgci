# Example: Mach-5 oblique shock over a 15-degree wedge

Grid-convergence verification of the area-averaged ramp-surface pressure on a
four-grid hierarchy, checked against exact oblique-shock theory and analysed
with the `foamgci` library. This is the second example. It reuses the library
unchanged and reads a different OpenFOAM output from the forward step, so it
doubles as the portability demonstration: same GCI engine, new data source.

## What lives where

- `0/`, `constant/`, `system/`: the committed OpenFOAM case (this copy is the
  **fine** grid). These are the OpenFOAM-4.x `wedge15Ma5` tutorial files,
  unchanged except for two case-specific additions in `system/` (a
  `surfaceFieldValue` function object in `controlDict` and a
  `decomposeParDict`).
- `submit.sh`: SLURM runner, identical to the forward-step example. Mesh,
  decompose, run, reconstruct.
- `gci/`: the analysis driver for THIS case. It `import`s `foamgci` (the
  library) and applies it here. You edit `gci/data.py`; you do not touch the
  library.
- `gci/oblique_shock.py`: the analytical reference. An exact theta-beta-M solve
  for the post-shock state, NumPy only.
- `gci/data/`: the four `surfaceFieldValue.dat` inputs the analysis reads, one
  per grid. They ship empty; produce them by running the four cases (see
  `gci/data/README.md`).

## Physics

Inviscid (Euler) Mach-5 flow over a 15-degree compression ramp, solved with
`rhoCentralFoam` (Kurganov fluxes, vanLeer limiters, fixed time step).
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
3. Copy each run's `postProcessing/wallPressure/0/surfaceFieldValue.dat` into
   `gci/data/`, renamed `coarse.dat`, `medium.dat`, `fine.dat`,
   `extrafine.dat`.
4. Analyse:
   ```bash
   pip install -e ../..            # install the foamgci library (once)
   cd gci && bash run_all.sh
   ```
   Writes `gci/gci_summary.json`: per-grid mean, sigma, tau_int, SEM, N_eff,
   KPSS; regime-aware GCI on both triplets; and the error against the
   analytical reference.
5. Commit the four small `gci/data/*.dat` files. That is what makes the study
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
  narrow the window in `gci/data.py` (`T_STAT`) and re-run `analyze.py`.

This case carries one formal QoI:

| QoI            | Role                 | Interpretation                                                                                                                                                                                                 |
| -------------- | -------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `p_wall_ratio` | Primary verified QoI | Area-averaged ramp-surface pressure, an integrated functional. Steady and stationary on the window, converges monotonically to the analytical p2/p1. Richardson extrapolation is well founded for a surface integral, unlike for a pointwise extremum. |

The analysis driver:

1. Parse the `surfaceFieldValue` output, a different OpenFOAM source from the
   forward step's `fieldMinMax`.
2. Check stationarity using KPSS.
3. Estimate autocorrelation-corrected SEM and effective sample size.
4. Compute GCI on every consecutive triplet.
5. Extrapolate to h -> 0 and compare against exact oblique-shock theory.

### JSON summary structure

`gci_summary.json` contains:

* `qoi`: name, description, and the OpenFOAM source of the quantity;
* `stationary_window`: the averaging window used for the QoI;
* `reference`: the analytical oblique-shock state (beta, p2/p1, rho2/rho1,
  T2/T1, M2);
* `cases`: per-grid mean, sigma, tau_int, SEM, N_eff, KPSS;
* `triplet_A_CMF`, `triplet_B_MFXF`: regime-aware GCI on both triplets;
* `phi_star`, `phi_star_source`, `phi_star_rel_err_pct`: the extrapolate and
  its error against the reference;
* `reference_covered_by_gci`: whether the reference falls inside the
  fine-grid GCI band;
* `error_table_vs_reference`: per-grid error against p2/p1.

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
│   ├── controlDict             # + surfaceFieldValue function object (QoI source)
│   ├── decomposeParDict
│   ├── fvSchemes
│   └── fvSolution
├── gci/                        # per-case driver; imports foamgci, never edits it
│   ├── analyze.py              # reads surfaceFieldValue -> GCI -> vs analytical
│   ├── data.py                 # grid metadata, window, free-stream conditions
│   ├── oblique_shock.py        # exact theta-beta-M reference (NumPy only)
│   ├── run_all.sh
│   └── data/                   # four surfaceFieldValue.dat inputs (ship empty)
│       └── README.md
├── README.md
└── submit.sh                   # identical to the forward-step runner
```

## References

- Anderson, J.D. (2003), *Modern Compressible Flow*, 3rd ed., McGraw-Hill: oblique-shock relations.
- Ames Research Staff (1953), NACA Report 1135, *Equations, Tables, and Charts for Compressible Flow*: the tabulated reference.
- Greenshields, C.J. et al. (2010), *Int. J. Numer. Meth. Fluids* **63**(1), 1–21: the rhoCentralFoam central-upwind solver.
- Kurganov, A., Noelle, S., Petrova, G. (2001), *SIAM J. Sci. Comput.* **23**(3), 707–740: the KNP flux scheme.
