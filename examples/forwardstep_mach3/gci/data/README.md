# Grid-convergence inputs (`fieldMinMax.dat`)

`analyze.py` expects four files in this folder, one per grid:

| file | grid |
|---|---|
| `coarse.dat` | Coarse |
| `medium.dat` | Medium |
| `fine.dat` | Fine |
| `extrafine.dat` | Extra-fine |

Each file is the OpenFOAM `fieldMinMax` function-object output from one grid-refinement run.

## Current status

The `.dat` files are intentionally not committed yet because the simulation campaign will be rerun. To reproduce the workflow now, run the OpenFOAM cases first and copy the generated files into this folder.

## How to populate this folder

After each grid finishes, copy:

```bash
cp <coarse_case>/postProcessing/fieldMinMax/0/fieldMinMax.dat     coarse.dat
cp <medium_case>/postProcessing/fieldMinMax/0/fieldMinMax.dat     medium.dat
cp <fine_case>/postProcessing/fieldMinMax/0/fieldMinMax.dat       fine.dat
cp <extrafine_case>/postProcessing/fieldMinMax/0/fieldMinMax.dat  extrafine.dat
