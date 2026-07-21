# wedge15Ma5

Mach-5 flow over a 15-degree wedge, inviscid, rhoCentralFoam (OpenFOAM-4.x
tutorial). The **portability** case: same engine, two OpenFOAM outputs per grid.
Ships with its eight `.dat` committed, so it runs from a clone.

## QoIs

| QoI | source | what it's for |
|---|---|---|
| `p_wall_ratio` | `surfaceRegion` area-average | the verified QoI, vs oblique-shock theory |
| `max(p)` | `fieldMinMax` | diagnostic -- same output the step used, comes out mesh-locked |

## The finding

The ramp pressure equals post-shock `p2`, so it converges to the exact
oblique-shock value (`gci/oblique_shock.py`: beta=24.32 deg, p2/p1=4.7808). A
clean, well-posed target. `max(p)` chases the same `p2` but sits on a flat
plateau, so its location is degenerate and the localization check flags it.
Well-posed integral vs degenerate extremum, same pressure, that's the point.
Don't expect second order: it's first order at the shock, so measure it.

## The four grids

Uniform square cells (h = 0.1524/Ny). Only the two block counts change.

| label | cells | h | block 1 | block 2 |
|---|--:|---|---|---|
| coarse | 1,200 | 0.0076200 | (20 20 1) | (40 20 1) |
| medium | 4,800 | 0.0038100 | (40 40 1) | (80 40 1) |
| **fine** (committed) | 19,200 | 0.0019050 | (80 80 1) | (160 80 1) |
| extra-fine | 76,800 | 0.0009525 | (160 160 1) | (320 160 1) |

## Run it

Data's committed, so just:

```bash
pip install -e ../..
cd gci && bash run_all.sh        # -> gci_summary.json + figures
```

To regenerate: make the four `<prefix>_grid/` dirs as siblings of `gci/`, edit
each `blockMeshDict` per the table, run each, then copy both outputs per grid
into `gci/data/` (see `gci/data/README.md` for the names).
