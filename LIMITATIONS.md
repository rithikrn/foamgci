# Limitations

Small tool, on purpose. Here's what it doesn't do, so you can decide fast.

- **It reads your data; it doesn't run OpenFOAM.** You bring converged
  time-histories on a refinement hierarchy. No meshing, no `checkMesh`, no
  solver.
- **Space and time refine together** unless you fix `deltaT`. With adaptive
  `deltaT` the apparent order mixes space and time error. Use a fixed step
  small enough for the finest grid (the examples do).
- **First order at shocks is normal.** rhoCentralFoam drops to first order
  across a captured shock, so shock-dominated QoIs converge ~first order even
  though the scheme is nominally second. Measure the order, don't assume it.
  The diamond entropy integral is the clearest example: surface pressures look
  great, entropy shows the truth.
- **A max/min isn't a smooth QoI.** GCI on `fieldMinMax` is a heuristic; the
  extremum can hop between cells. Prefer an integral (force, surface, volume).
  The diamond `max(p)` is kept only as a mesh-lock diagnostic.
- **You pick the averaging window.** KPSS checks it; if it rejects, narrow it.
- **Three grids minimum.** Fewer gets you stats but no GCI. `R_fit` is a
  diagnostic, not an asymptotic-range test.
- **Core readers cover surface averages and `fieldMinMax`.** Forces
  (`coefficient.dat`) and the entropy volume integral live in the diamond
  driver for now, not the library.

Results depend on your solver, scheme, time step, and `controlDict`. Pin them.
