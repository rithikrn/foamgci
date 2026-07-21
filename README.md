
![foamgci_logo](media/images/foamgci_logo.png)

# foamgci
 
Grid-convergence (GCI) checks for OpenFOAM, with the statistics done properly:
Roache GCI on every refinement triplet, an autocorrelation-aware standard error,
a KPSS stationarity test, and an analytical cross-check where one exists.
 
## Install
 
```bash
pip install -e ".[dev]"
pytest -q
```
 
Needs NumPy (Python >= 3.10). matplotlib is optional, for figures.
 
## See it work
 
No OpenFOAM needed for the first two, their inputs are committed:
 
```bash
cd examples/forwardstep_mach3/gci && bash run_all.sh   # writes gci_summary.json + figures
```
 
Then open `gci_summary.json` and the case README.
 
## The three cases
 
- **`forwardstep_mach3`** -- Mach-3 step. The template. A pointwise `max(p)` QoI.
- **`wedge15Ma5`** -- Mach-5 wedge. Adds a surface average, checked against exact
  oblique-shock theory.
- **`diamond2D_Ma2`** -- Mach-2 diamond. Five grids, and a volume QoI (entropy)
  that shows the convergence is really first order. 
## How it's split
 
Two things, kept apart:
 
- `foamgci/` is the **library** (case-agnostic, the thing you `import`). Don't edit it.
- `examples/<case>/gci/` is the **driver** for one case. Edit `data.py` there.
To add a case, copy an example folder and point its `data.py` at your files.
 
## Run on your own case
 
1. Add the right function object to `system/controlDict` (see the examples).
2. Run your grids, collect one output file per grid.
3. Copy an example's `gci/`, edit `data.py`, run `analyze.py`.
The Python API and CLI live behind `from foamgci import full_report` and
`foamgci report --help`. The examples show both in context.

## Repository layout

```
foamgci/
├── foamgci/                       # THE LIBRARY (pip-installable, case-agnostic)
│   ├── _version.py                #   single source of the version string
│   ├── __init__.py                #   public API
│   ├── __main__.py                #   CLI: foamgci report ...
│   ├── reader.py                  #   OpenFOAM fieldMinMax + surfaceRegion area-average + generic readers
│   ├── stats.py                   #   Geyer tau_int, KPSS, window stats
│   ├── gci.py                     #   Roache GCI on triplets
│   ├── report.py                  #   end-to-end driver + Rayleigh-Pitot
│   └── plot.py                    #   optional matplotlib figure
├── tests/                         # pytest suite, anchored to Celik (2008)
├── examples/
│   ├── forwardstep_mach3/         
|   ├── wedge15Ma5/                
│   └── diamond2D_Ma2/             
├── README.md  LIMITATIONS.md  CONTRIBUTING.md  CHANGELOG.md  LICENSE
├── pyproject.toml
└── .github/workflows/tests.yml    # CI: Linux + macOS, Py 3.10-3.12
```

## Citing

If you use `foamgci` in published work, please cite it. A plain-text
acknowledgement:

> Grid-convergence study performed with `foamgci` (Roache GCI + Geyer
> tau_int + KPSS stationarity), https://github.com/rithikrn/foamgci.

Or as BibTeX:

```bibtex
@software{nambiar_2026_20946450,
  author       = {Nambiar, Rithik R},
  title        = {foamgci},
  month        = jun,
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v3.5.0},
  doi          = {10.5281/zenodo.20946450},
  url          = {https://doi.org/10.5281/zenodo.20946450},
}
```

Pin the version (release tag or commit hash) you actually used.

## References

**Verification & uncertainty methodology**
- Roache, P.J. (1994), *J. Fluids Eng.* **116**(3), 405–413.
- Celik, I.B. et al. (2008), *J. Fluids Eng.* **130**(7), 078001.
- Eça, L., Hoekstra, M. (2014), *J. Comput. Phys.* **262**, 104–130.
- Oberkampf, W.L., Roy, C.J. (2010), *Verification and Validation in
  Scientific Computing*, Cambridge Univ. Press.

**Autocorrelation-corrected sampling statistics**
- Geyer, C.J. (1992), *Statistical Science* **7**(4), 473–483.
- Flyvbjerg, H., Petersen, H.G. (1989), *J. Chem. Phys.* **91**(1), 461–466.
- Kwiatkowski, D. et al. (1992), *J. Econometrics* **54**(1–3), 159–178.

**Time-averaging uncertainty in CFD (recent context)**
- Rezaeiravesh, S., Vinuesa, R., Schlatter, P. (2022), *J. Comput. Sci.* **62**, 101688.
- Xavier, D., Rezaeiravesh, S., Schlatter, P. (2024), *Phys. Fluids* **36**(10), 105122.
- Related software: Rezaeiravesh, S. et al. (2021), *UQit*, JOSS **6**(60), 2871.

**Solver & benchmark**
- Woodward, P., Colella, P. (1984), *J. Comput. Phys.* **54**(1), 115–173.
- Kurganov, A., Noelle, S., Petrova, G. (2001), *SIAM J. Sci. Comput.* **23**(3), 707–740.
- Greenshields, C.J. et al. (2010), *Int. J. Numer. Meth. Fluids* **63**(1), 1–21.

## License

MIT, see `LICENSE`.
