# Grid-convergence inputs (`surfaceFieldValue.dat`)

`analyze.py` expects four files in this folder, one per grid:

| file | grid | block1 / block2 cells | total | h = 0.1524/Ny |
|---|---|---|---|---|
| `coarse.dat`    | Coarse     | (20 20 1) / (40 20 1)   | 1 200  | 0.0076200 |
| `medium.dat`    | Medium     | (40 40 1) / (80 40 1)   | 4 800  | 0.0038100 |
| `fine.dat`      | Fine       | (80 80 1) / (160 80 1)  | 19 200 | 0.0019050 |
| `extrafine.dat` | Extra-fine | (160 160 1) / (320 160 1) | 76 800 | 0.0009525 |

Each file is the OpenFOAM `surfaceFieldValue` function-object output
(`areaAverage(p)` on patch `obstacle`) from one grid-refinement run — a
**different** output type from the forward-step example, which used
`fieldMinMax`. The same foamgci GCI/statistics engine consumes both.

## How to produce each file (the only thing you change between runs)

From a copy of this case directory, for each grid:

1. Edit the two `(nx ny 1)` tuples in `system/blockMeshDict` per the table
   above. Change nothing else — same `0/`, `constant/`, `system/fvSchemes`,
   `system/fvSolution`, `system/controlDict`.
2. Run the case (locally or via `submit.sh`):
   ```bash
   blockMesh
   rhoCentralFoam            # or: decomposePar && mpirun -np 16 rhoCentralFoam -parallel && reconstructPar
   ```
3. Copy the result here, renamed for the grid:
   ```bash
   cp postProcessing/wallPressure/0/surfaceFieldValue.dat  <grid>.dat
   ```

## Current status

These `.dat` files are **not** committed: this case ships ready-to-run but
without precomputed CFD output. Populate `data/` from your own four runs,
then:

```bash
pip install -e ".[dev]"          # from the repository root
cd examples/wedge15Ma5/gci
bash run_all.sh                  # analyze.py -> gci_summary.json
```

`analyze.py` will report, per grid, the time-averaged wall pressure, the
Roache GCI on each refinement triplet, the Richardson-extrapolated
`p_wall/p_inf`, and its error against the exact analytical post-shock
pressure `p2/p1` from `oblique_shock.py`.
