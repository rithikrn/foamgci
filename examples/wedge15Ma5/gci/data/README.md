# gci/data

`analyze.py` reads **two files per grid** (eight total), one per function
object in `controlDict`:

| grid | surfaceFieldValue (primary) | fieldMinMax (secondary) |
|---|---|---|
| Coarse | `coarse_surfaceFieldValue.dat` | `coarse_fieldMinMax.dat` |
| Medium | `medium_surfaceFieldValue.dat` | `medium_fieldMinMax.dat` |
| Fine | `fine_surfaceFieldValue.dat` | `fine_fieldMinMax.dat` |
| Extra-fine | `extrafine_surfaceFieldValue.dat` | `extrafine_fieldMinMax.dat` |

Copy from each run, renamed for the grid:

```
cp postProcessing/wallPressure/0/surfaceRegion.dat  coarse_surfaceFieldValue.dat
cp postProcessing/fieldMinMax/0/fieldMinMax.dat     coarse_fieldMinMax.dat
```

Naming note: in OpenFOAM-4.x the area-average file is `surfaceRegion.dat`; in
v5.0+ it's `surfaceFieldValue.dat`. The local name is just our label, the reader
picks the column from the header, not the filename.

These eight are committed, so the study runs from a clone:

```bash
pip install -e ".[dev]"
cd examples/wedge15Ma5/gci && bash run_all.sh
```

Missing a file? `analyze.py` lists exactly which and won't run on a partial set.
