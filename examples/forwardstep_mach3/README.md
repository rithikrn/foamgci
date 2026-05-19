# Worked example — Mach-3 Woodward–Colella forward-facing step

This directory contains the verification driver for the four-grid
hierarchy whose results populate Table 1 of the AIAA SciTech 2027
extended abstract.

## Cases

Four identical OpenFOAM cases differ only in the `system/blockMeshDict`
spacing. All other settings come from the shared `controlDict` listed
in the abstract appendix.

| label       | cell count | $h$       | refinement ratio vs prev |
| ----------- | ----------:| --------- | ------------------------ |
| coarse      |      4,032 | 0.025     | —                        |
| medium      |     16,128 | 0.0125    | 2                        |
| fine        |     64,512 | 0.00625   | 2                        |
| extra-fine  |    258,048 | 0.003125  | 2                        |

Solver: `rhoCentralFoam` (KNP fluxes, vanLeer / vanLeerV limiters,
CFL ≤ 0.2). `t ∈ [0, 10]` flow-through times; the first ≈ 2 are
transient and the time-averaging window is `t ∈ [3, 10]`.

## What this example does **not** ship

To keep the repository small, the **`fieldMinMax.dat` outputs of the
four runs are not included** — they are tens of MB each. Run your
own cases to produce them. The pipeline is:

1. Add the `fieldMinMax` function object to each case's
   `system/controlDict`:

       functions {
           fieldMinMax {
               type            fieldMinMax;
               libs            ("fieldFunctionObjects");
               writeControl    timeStep;
               writeInterval   100;
               fields          (p U rho);
               location        true;
               mode            magnitude;
           }
       }

   This will produce `postProcessing/fieldMinMax/0/fieldMinMax.dat`
   in each case directory.

2. Run all four cases (see `submit.sh` in your case repository).

3. Drop the four files into `data/` of this directory (or pass full
   paths via `--coarse`, `--medium`, `--fine`, `--extra-fine`).

## Running the verification

```bash
# from the foamgci repo root
python examples/forwardstep_mach3/verify_abstract.py \
    --coarse     /path/to/coarse/postProcessing/fieldMinMax/0/fieldMinMax.dat \
    --medium     /path/to/medium/postProcessing/fieldMinMax/0/fieldMinMax.dat \
    --fine       /path/to/fine/postProcessing/fieldMinMax/0/fieldMinMax.dat \
    --extra-fine /path/to/extra-fine/postProcessing/fieldMinMax/0/fieldMinMax.dat \
    --window 3 10
```

Outputs go to `./out/`:

- `report.txt` — human-readable Table 1 + analytical cross-check
- `table1.tex` — drop-in LaTeX block for the abstract
- `fig_convergence.pdf` — two-panel convergence figure

## What to compare against the abstract

The script prints the percent error of two quantities relative to
the analytical Rayleigh-Pitot reference $p_{02}/p_1 = 12.0610$ at
$M=3$, $\gamma=1.4$:

1. The **finest-grid mean** $\langle p_{\max}\rangle_{\text{XF}}$
2. The **Richardson-extrapolated value** $\phi_{\text{exact}}$ from
   the medium → fine → extra-fine triplet

Both should be below ~ 0.1 % for a properly resolved case. The
abstract's specific claim is that the extra-fine result agrees with
Rayleigh-Pitot to within **0.03 %** — re-run this script and update
the abstract figure if your number differs from 0.03 %.

## What to put in Table 1

Copy `out/table1.tex` verbatim into your manuscript and replace the
existing Table 1 environment. The script writes the table in the
exact format used in the abstract (per-grid block + GCI block).

## Notes on the time-averaging window

If KPSS rejects stationarity on the `[3, 10]` window (`KPSS_p < 0.05`)
for any grid, the window contains residual transient. Re-run with
`--window 5 10` or another tightened range and re-verify.
