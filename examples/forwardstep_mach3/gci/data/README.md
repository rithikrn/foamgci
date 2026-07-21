# gci/data

`analyze.py` reads one `fieldMinMax.dat` per grid:

| file | grid |
|---|---|
| `coarse.dat` | Coarse |
| `medium.dat` | Medium |
| `fine.dat` | Fine |
| `extrafine.dat` | Extra-fine |

Each is the `fieldMinMax` output from one run
(`postProcessing/fieldMinMax/0/fieldMinMax.dat`).

These four are committed, so the study runs from a clone:

```bash
pip install -e ".[dev]"
cd examples/forwardstep_mach3/gci && bash run_all.sh
```
