# Example: Mach-3 Woodward-Colella forward-facing step

Grid-convergence verification of the time-averaged peak pressure on a
four-grid hierarchy, analysed with the `foamgci` library. This directory
is the **template** for adding further cases.

## What lives where

- `0/`, `constant/`, `system/` — the committed OpenFOAM case (this copy
  is the **fine** grid).
- `submit.sh` — SLURM runner: mesh, decompose, run, reconstruct, and
  copy the result into `gci/data/`.
- `gci/` — the analysis driver for THIS case. It `import`s `foamgci`
  (the library) and applies it here. You edit `gci/data.py`; you do not
  touch the library.
- `gci/data/` — the four small `fieldMinMax.dat` inputs the analysis
  reads. Committed so the figures reproduce without rerunning OpenFOAM.

## Physics

Inviscid (Euler) Mach-3 flow over a forward-facing step, solved with
`rhoCentralFoam` (Kurganov/KNP fluxes, vanLeer limiters, CFL <= 0.2).
Normalised perfect gas: a_inf = 1 at T = 1, so M = U = 3, p_inf = 1.
QoI: <max p> over the stationary window t in [6, 10].

**Reference value.** The global pressure maximum sits at the bow-shock
stagnation foot, where the post-shock flow comes to rest behind the
locally normal portion of the bow shock. The inviscid Rayleigh–Pitot
stagnation pressure p02/p1 = 12.061 is therefore a physical **ceiling**:
numerical dissipation smears the captured shock, so the time-mean
<max p> approaches 12.061 **from below** as h → 0 and cannot exceed it.
All four grids confirm this (11.996 → 12.014 → 12.041 → 12.057, all
below 12.061). The verification rests on the finest-grid agreement with
this reference (within 0.029%), **not** on the Richardson extrapolate:
phi_star ≈ 12.080 overshoots the ceiling, the expected signature of
applying Richardson extrapolation to a non-smooth pointwise extremum.

## The four grids

The committed case is the **fine** grid. The four grids differ ONLY in
the three block counts in `system/blockMeshDict`; everything else is
identical.

| label                 | cells   | h        | block 1 | block 2  | block 3   |
|-----------------------|--------:|----------|---------|----------|-----------|
| coarse                |   4,032 | 0.025    | (24 8)  | (24 32)  | (96 32)   |
| medium                |  16,128 | 0.0125   | (48 16) | (48 64)  | (192 64)  |
| **fine** (committed)  |  64,512 | 0.00625  | (96 32) | (96 128) | (384 128) |
| extra-fine            | 258,048 | 0.003125 | (192 64)| (192 256)| (768 256) |

i.e. in the `blocks (...)` list, set the three `(nx ny 1)` entries to the
row for that grid.

## Run from scratch (on HPC)

1. Make the four case directories as siblings of `gci/`:
   ```bash
   cd examples/forwardstep_mach3
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
   (If you run interactively instead of via SLURM, copy the files by hand —
   see `gci/data/README.md`.)
3. Analyse + plot:
   ```bash
   pip install -e ../..            # install the foamgci library (once)
   cd gci && bash run_all.sh
   ```
   Writes `gci/gci_summary.json` (per-grid mean, sigma, tau_int, SEM,
   N_eff, KPSS; regime-aware GCI on both triplets) and the figures.
4. Commit the four KB-sized `gci/data/*.dat` files. That is what makes the
   study reproducible for the next reader — they just run step 3.

## Reading the output

- Per grid: `p_max_mean`, the autocorrelation-corrected `p_max_sem`, and
  `p_max_n_eff` (effective sample count after the tau_int correction —
  well below the raw count).
- `triplet_A_CMF` may be **divergent** (R > 1, pre-asymptotic);
  `triplet_B_MFXF` is the deeper triplet used for `phi_star`.
- `phi_star` is the converged value; `rayleigh_pitot_p02` (12.061) is the
  lower-bound check.
- If KPSS rejects stationarity on `[3, 10]` for any grid, narrow the
  window in `gci/data.py` (`T_STAT`) and re-run `analyze.py`.

### Multi-QoI forward-step diagnostics (latest)

`foamgci` is not restricted to pressure. The current OpenFOAM
`fieldMinMax.dat` reader can extract any scalar field written by the
function object, including `p`, `rho`, and velocity magnitude when
available. The Mach-3 forward-step example currently promotes two formal
QoIs:

| QoI       | Role                     | Interpretation                                                                                                                          |
| --------- | ------------------------ | --------------------------------------------------------------------------------------------------------------------------------------- |
| `p_max`   | Primary verified QoI     | Stationary on all four grids, valid deepest GCI triplet, and extra-fine value agrees with Rayleigh--Pitot within approximately `0.03%`. |
| `rho_max` | Secondary diagnostic QoI | Very clean spatial convergence, but refined-grid density histories fail KPSS stationarity in the current window.                        |

The intended workflow is diagnostic:

1. Parse multiple scalar QoIs from the same OpenFOAM output.
2. Check stationarity using KPSS.
3. Estimate autocorrelation-corrected SEM and effective sample size.
4. Compute GCI on every consecutive triplet.
5. Promote only the QoIs that are stationary, physically meaningful, and
   asymptotically reliable.

### JSON summary structure

`gci_summary.json` contains:

* `stationary_window` — the averaging window used for all QoIs;
* `cases` — backward-compatible top-level pressure-centered output, now
  also carrying density statistics;
* `qoi_results` — full per-QoI statistics and triplet GCI diagnostics;
* `included_qois` — QoIs promoted to the formal table;
* `excluded_qois` — parsed or considered QoIs excluded from formal GCI
  interpretation, with a reason;
* `error_table_vs_rayleigh_pitot` — pressure-only analytical reference
  comparison.

## Adding a new case

Copy this whole directory, swap `system/blockMeshDict` + `0/`, set the
grid metadata and reference in `gci/data.py`, populate `gci/data/`, and
write the new case's README. The library is never edited.
