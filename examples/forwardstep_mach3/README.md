# forwardstep_mach3

Mach-3 Woodward-Colella forward-facing step, inviscid, rhoCentralFoam. This is
the **template** case, copy it to start a new one. Ships with its four `.dat`
committed, so it runs from a clone.

## QoIs

| QoI | source | what it's for |
|---|---|---|
| `p_max` | `fieldMinMax` | the verified QoI, vs Rayleigh-Pitot |
| `rho_max` | `fieldMinMax` | diagnostic only -- looks great, fails the gates |

## The finding

`p_max` sits at the bow-shock stagnation foot, so Rayleigh-Pitot (12.061) is a
ceiling it approaches from below (11.996 -> 12.051 across the grids, finest
within 0.085%). `rho_max` is the cautionary twin: its GCI looks lovely
(p_obs~2.1, GCI~0.035%), but KPSS fails on every grid and the max wanders ~79
cells between the triple point and the stagnation foot, so it's comparing
different physical maxima. The toolkit flags it; we keep it as a diagnostic.

## The four grids

Only the three block counts in `system/blockMeshDict` change between grids.

| label | cells | h | block 1 | block 2 | block 3 |
|---|--:|---|---|---|---|
| coarse | 4,032 | 0.025 | (24 8) | (24 32) | (96 32) |
| medium | 16,128 | 0.0125 | (48 16) | (48 64) | (192 64) |
| **fine** (committed) | 64,512 | 0.00625 | (96 32) | (96 128) | (384 128) |
| extra-fine | 258,048 | 0.003125 | (192 64) | (192 256) | (768 256) |

## Run it

Data's committed, so just:

```bash
pip install -e ../..
cd gci && bash run_all.sh        # -> gci_summary.json + figures
```

To regenerate from scratch: make `coarse_grid/` ... `extrafine_grid/` as
siblings of `gci/`, edit each `blockMeshDict` per the table, `sbatch submit.sh`
in each, then copy each grid's `fieldMinMax.dat` into `gci/data/`.
