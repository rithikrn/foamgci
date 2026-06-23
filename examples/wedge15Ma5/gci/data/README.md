# Grid-convergence inputs (two OpenFOAM outputs per grid)

`analyze.py` reads **two files per grid**, one from each function object in
`system/controlDict`. Eight files total, one pair per grid:

| grid       | block1 / block2 cells     | total  | h = 0.1524/Ny | surfaceFieldValue (primary)      | fieldMinMax (secondary)      |
|------------|---------------------------|-------:|---------------|----------------------------------|------------------------------|
| Coarse     | (20 20 1) / (40 20 1)     |  1 200 | 0.0076200     | `coarse_surfaceFieldValue.dat`   | `coarse_fieldMinMax.dat`     |
| Medium     | (40 40 1) / (80 40 1)     |  4 800 | 0.0038100     | `medium_surfaceFieldValue.dat`   | `medium_fieldMinMax.dat`     |
| Fine       | (80 80 1) / (160 80 1)    | 19 200 | 0.0019050     | `fine_surfaceFieldValue.dat`     | `fine_fieldMinMax.dat`       |
| Extra-fine | (160 160 1) / (320 160 1) | 76 800 | 0.0009525     | `extrafine_surfaceFieldValue.dat`| `extrafine_fieldMinMax.dat`  |

**Why two outputs.** This case deliberately writes both, to (1) stay
consistent with the forward-step example, which used `fieldMinMax`, and (2)
show that the *same* foamgci engine ingests a new output type and multiple
files in one analysis:

- `surfaceFieldValue.dat`, **primary**, reference-anchored QoI: the
  area-averaged ramp-surface pressure (`areaAverage(p)` on patch `obstacle`),
  an integrated functional for which Richardson/GCI is well founded. It is
  compared to the exact oblique-shock `p2/p1`.
- `fieldMinMax.dat`, **secondary**, diagnostic QoI: the global `max(p)`. Its
  value also approaches `p2`, but the post-shock field is a near-uniform
  plateau, so the extremum's *location* is degenerate and foamgci's
  localization check flags it. That contrast, a well-posed surface integral
  vs a degenerate pointwise extremum targeting the *same* pressure, is the
  point. Filenames must contain a field column (`p`) the reader can select.

The same foamgci GCI/statistics engine consumes both.

## How to produce the files (the only thing you change between runs)

From a copy of this case directory, for each grid:

1. Edit the two `(nx ny 1)` tuples in `system/blockMeshDict` per the table
   above. Change nothing else, same `0/`, `constant/`, `system/fvSchemes`,
   `system/fvSolution`, `system/controlDict`.
2. Run the case (locally or via `submit.sh`):
   ```bash
   blockMesh
   rhoCentralFoam            # or: decomposePar && mpirun -np 16 rhoCentralFoam -parallel && reconstructPar
   ```
3. Copy BOTH outputs here, renamed for the grid (example for the coarse run):
   ```bash
   cp postProcessing/wallPressure/0/surfaceRegion.dat  coarse_surfaceFieldValue.dat
   cp postProcessing/fieldMinMax/0/fieldMinMax.dat     coarse_fieldMinMax.dat
   ```
   OpenFOAM names each function object's output folder after the object name
   (`wallPressure/`, `fieldMinMax/`). The file inside is named after the object
   *type*: in OpenFOAM-4.x the area-average file is `surfaceRegion.dat`; in
   OpenFOAM v5.0+ it is `surfaceFieldValue.dat` (and you would set
   `type surfaceFieldValue` in controlDict). The local name on the right
   (`coarse_surfaceFieldValue.dat`) is just our label for the analysis; the
   reader picks the column from the file header, not the filename.

## Current status

These `.dat` files are **not** committed: the case ships ready-to-run but
without precomputed CFD output. Populate `data/` from your own four runs
(eight files), then:

```bash
pip install -e ".[dev]"          # from the repository root
cd examples/wedge15Ma5/gci
bash run_all.sh                  # analyze.py -> gci_summary.json + figures/
```

`analyze.py` will report, per grid and per QoI, the time-averaged value, the
autocorrelation-corrected SEM, the Roache GCI on each refinement triplet, the
Richardson extrapolate, and its error against the exact analytical post-shock
pressure `p2/p1` from `oblique_shock.py`. The secondary QoI additionally
reports the in-window extremum wander (cells) and the localization verdict.

If `analyze.py` exits with `ERROR: missing input file(s)` it lists exactly
which of the eight files is absent; it never runs on a partial input set.
