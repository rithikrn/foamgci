# diamond2D_Ma2 — supersonic diamond airfoil (M=2), multi-QoI verification

A second steady-supersonic benchmark alongside `wedge15Ma5`, built to exercise
the QoI types the wedge and forward-step cases do not: a **surface force**
(drag), **two coupled** quantities that must agree, and a **volume integral**.
The OpenFOAM setup reuses the `wedge15Ma5` rhoCentralFoam configuration
verbatim (same `thermophysicalProperties`, `fvSchemes`, `fvSolution`); only the
geometry and the function-object QoIs are case-specific.

## Flow

Symmetric diamond (double-wedge) airfoil, half-angle 10°, chord 1, in a Mach-2
inviscid stream. The front facets compress the flow through an attached oblique
shock; the rear facets expand it through a Prandtl-Meyer fan. Shock-expansion
theory gives the closed-form references the QoIs are checked against:

| quantity            | value   |
|---------------------|---------|
| wave angle β        | 39.31°  |
| p_front / p∞        | 1.7066  |
| p_rear / p∞         | 0.5508  |
| C_d (shock-expansion)| 0.0728 |
| C_d (linear/Ackeret)| 0.0704  |

`gci/shock_expansion.py` computes these from scratch (numpy bisection, run it
directly to print them).

## Half-domain

The airfoil is symmetric about y=0, so only the **upper half** is meshed
(`symmetry` plane along y=0, slip airfoil surface). Every airfoil coefficient
is therefore the upper-surface contribution and is **doubled** in
post-processing (`CD_HALF_TO_FULL = 2` in `data.py`). This mirrors how
`wedge15Ma5` uses a symmetry plane.

## QoI hierarchy

| QoI       | function object        | type             | reference        |
|-----------|------------------------|------------------|------------------|
| `Cd`      | `forceCoeffs`          | surface force    | 0.0728           |
| `p_front` | `surfaceFieldValue`    | surface average  | 1.7066           |
| `p_rear`  | `surfaceFieldValue`    | surface average  | 0.5508           |
| `Cd_press`| derived (coupling)     | coupled check    | = Cd             |
| `S_vol`   | `volFieldValue`        | volume integral  | converged value  |
| `max(p)`  | `fieldMinMax`          | pointwise (diag) | mesh-locked      |

The **coupling** is the cross-check `Cd_pressure = (p_front − p_rear)·(t/c) /
(½γM²)` against the independent `forceCoeffs` drag. The two routes agree only
if the force normalisation and the pressure integration are both right, so a
disagreement flags a setup error rather than passing silently.

The **smoothness contrast** is deliberate: `Cd`, the facet pressures, and the
entropy integral are smooth functionals that converge with a clean apparent
order; `max(p)` is a pointwise extremum that sits ~42% above the post-shock
pressure and whose location tracks the leading-edge tip at a fixed cell offset
(mesh-locked). The analysis accepts the former and flags the latter.

Five grids (1120 → 286 720 cells, r=2) let the Eça-Hoekstra least-squares fit
run alongside the triplet Roache GCI.

## Run it

1. Edit the four `(nx ny 1)` tuples in `system/blockMeshDict` per the table in
   its header to select a grid level.
2. Submit one job per grid (`submit.sh coarse`, `... medium`, …); each runs
   `blockMesh → rhoCentralFoam` and copies the QoI `.dat` into `gci/data/` with
   the matching prefix.
3. `bash gci/run_all.sh` → `gci_summary.json` and `figures/`.

This example ships **without** precomputed `.dat` (unlike `wedge15Ma5`): the
five grids have not been run here. `analyze.py` and the references are verified
for structure and physics; the mesh and the runs need your cluster.

## Verify on your cluster

- **`blockMesh` + `checkMesh`.** The mesh is hand-written (4-column H-grid,
  airfoil surface along the lower boundary). Confirm no negative volumes and
  that the `airfoilFront`/`airfoilRear`/`symmetry`/`top` faces are assigned as
  intended before trusting any QoI.
- **`forceCoeffs` normalisation.** For compressible rhoCentralFoam the pressure
  is static; confirm the reported `Cd` has the expected magnitude. The coupling
  check (`Cd_force` vs `Cd_pressure`) will catch a wrong `rho`/`rhoInf` setting.
- **Coded entropy FO.** `entropyField` (in `controlDict`) is a `coded` object
  and is version-sensitive. If it does not compile, comment out `entropyField`
  and `sVol`; the other four QoIs are unaffected. The entropy integral can also
  be formed in post from the written `p`/`rho` fields.
- **Top boundary.** `top` is a symmetry plane forming a channel (as in the
  wedge), placed so the front-facet shock leaves downstream of the airfoil. A
  quick pressure contour should confirm no reflected shock strikes the airfoil.
