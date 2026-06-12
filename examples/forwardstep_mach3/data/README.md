# Grid-convergence inputs (`fieldMinMax.dat`)

`analyze.py` reads four files from this folder, one per grid:

| file            | grid       |
|-----------------|------------|
| `coarse.dat`    | Coarse     |
| `medium.dat`    | Medium     |
| `fine.dat`      | Fine       |
| `extrafine.dat` | Extra-fine |

Each is the `fieldMinMax` function-object output of one OpenFOAM run.
They are **KB-sized and committed to the repo** so anyone can reproduce
the figures without rerunning OpenFOAM.

## How to populate them from your HPC run

After a grid finishes, its min/max history is at
`<grid>/postProcessing/fieldMinMax/0/fieldMinMax.dat`. Copy each into
this folder under the name above:

```bash
cp ../../coarse_grid/postProcessing/fieldMinMax/0/fieldMinMax.dat     coarse.dat
cp ../../medium_grid/postProcessing/fieldMinMax/0/fieldMinMax.dat     medium.dat
cp ../../fine_grid/postProcessing/fieldMinMax/0/fieldMinMax.dat       fine.dat
cp ../../extrafine_grid/postProcessing/fieldMinMax/0/fieldMinMax.dat  extrafine.dat
```

(`submit.sh` does this copy automatically for you at the end of each
run.) Then commit the four `.dat` files — that is what makes the study
reproducible for the next reader.
