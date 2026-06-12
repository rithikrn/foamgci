# Example: Mach-3 Woodward–Colella forward-facing step

Grid-convergence verification of the time-averaged peak pressure on a
four-grid hierarchy, analysed with the `foamgci` package. This is the
template for adding further cases.

## Physics

Inviscid (Euler) Mach-3 flow over a forward-facing step, solved with
`rhoCentralFoam` (KNP fluxes, vanLeer limiters, CFL ≤ 0.2). Normalised
perfect gas: a∞ = 1 at T = 1, so M = U = 3, p∞ = 1. QoI: ⟨max p⟩ over
the stationary window t ∈ [3, 10].

**Reference.** The global pressure maximum sits at the bow-shock
stagnation foot, so the inviscid normal-shock value p₀₂/p₁ = 12.061 is a
physical **lower-bound check** — *not* the convergence target. ⟨max p⟩ is
the time-mean of a fluctuating spatial maximum, which exceeds the steady
stagnation value and rises slightly under refinement as finer grids
resolve more transient over-pressure. The convergence target is the
Richardson value φ★ from the deepest monotonic triplet (≈ 12.10).

## The committed case is the FINE grid

`0/`, `constant/`, and `system/` here are the **fine** grid. To run the
study, make four copies and change only the three block counts in
`system/blockMeshDict`:

| label      | n_cells | h        | block 1 | block 2  | block 3   |
|------------|--------:|----------|---------|----------|-----------|
| coarse     |   4 032 | 0.025    | (24 8)  | (24 32)  | (96 32)   |
| medium     |  16 128 | 0.0125   | (48 16) | (48 64)  | (192 64)  |
| **fine** (committed) | 64 512 | 0.00625 | (96 32) | (96 128) | (384 128) |
| extra-fine | 258 048 | 0.003125 | (192 64)| (192 256)| (768 256) |

i.e. in the `blocks (...)` list, set the three `(nx ny 1)` entries to the
row for that grid; everything else stays identical.

## Run from scratch

1. Make the four case directories as siblings of `gci/`:
```bash
   for d in coarse_grid medium_grid fine_grid extrafine_grid; do
     mkdir -p $d && cp -r 0 constant system submit.sh $d/
   done
   # then edit each $d/system/blockMeshDict block counts per the table above
```
2. Submit each (writes `postProcessing/fieldRange/0/fieldMinMax.dat`):
```bash
   for d in coarse_grid medium_grid fine_grid extrafine_grid; do
     ( cd $d && sbatch submit.sh ); done
```
3. Analyse + plot (uses the `foamgci` package — single source of truth):
```bash
   pip install -e ../..
   cd gci && ./run_all.sh
```
   Writes `gci/gci_summary.json` (per-grid mean, σ, τ_int, SEM, N_eff,
   KPSS; regime-aware GCI on both triplets) and the figures.

## Reading the output

- Per grid: `p_max_mean`, the autocorrelation-corrected `p_max_sem`, and
  `p_max_n_eff` (effective sample count after τ_int correction — well
  below the raw count).
- `triplet_A_CMF` is expected **divergent** (R > 1, pre-asymptotic);
  `triplet_B_MFXF` monotonic, p̂ ≈ 1.5, GCI ≈ 0.04 %.
- `phi_star` (≈ 12.10) is the converged value; `rayleigh_pitot_p02`
  (12.061) is the lower-bound check.

## Notes

- The `fieldRange` function object samples min/max every 100 steps. The
  raw `fieldMinMax.dat` files are small (KB) — commit them for
  reproducibility if you wish.
- If KPSS rejects stationarity on `[3, 10]` for any grid, tighten the
  window in `gci/data.py` (`T_STAT`) and re-run.

## Adding a new case

Copy this directory, swap `system/blockMeshDict` + `0/`, set the QoI
field and reference in `gci/analyze.py`, and write the case's README.
