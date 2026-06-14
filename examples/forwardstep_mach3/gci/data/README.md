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

The four `.dat` files are committed, so the forward-step GCI workflow can be reproduced from a fresh clone without rerunning OpenFOAM.

From the repository root:

```bash
pip install -e ".[dev]"
cd examples/forwardstep_mach3/gci
bash run_all.sh
