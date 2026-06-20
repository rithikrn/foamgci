
![foamgci_logo](media/images/foamgci_logo.png)

# foamgci - Autocorrelation-aware grid-convergence verification for OpenFOAM CFD studies.

`foamgci` reads OpenFOAM output and produces a complete
grid-convergence-index (GCI) report: Roache GCI on every refinement
triplet, an autocorrelation-corrected standard error of the mean, a KPSS
stationarity check, an analytical cross-check, and a paper-ready LaTeX
table.

---

## How this repository is organised (read this first)

There are two things in this repo, and keeping them straight removes
most of the confusion:

| | what it is | analogy |
|---|---|---|
| **`foamgci/`** (repo root) | The **library**. Generic, case-agnostic, pip-installable. All the reusable math lives here: read a scalar QoI time series, compute tau_int / SEM / KPSS, run Roache GCI, classify convergence, and render text & LaTeX reports. The built-in readers support OpenFOAM `fieldMinMax.dat` and generic scalar time-series files, so the same verification core can be used with solver-independent QoI histories. | this is `numpy` |
| **`examples/<case>/gci/`** | A **per-case driver** that *imports* the library and applies it to one case. Holds that case's grid metadata, the analysis script, and the figure scripts. | this is your `analysis.py` that does `import foamgci` |

Rule of thumb: anything reusable across cases belongs in `foamgci/`;
anything specific to one case (paths, mesh spacing, reference value,
figure styling) belongs in that example's `gci/` folder. The package
directory shares the repo name on purpose — `foamgci/` is exactly what
you `import foamgci`, which is the standard Python convention.

To add your own case later, copy `examples/forwardstep_mach3/`, swap the
mesh and boundary conditions, and edit only that example's `gci/data.py`.
The library stays untouched.

---

## What the library computes

- **Roache GCI** on every consecutive refinement triplet, using the
  Celik et al. (2008) iterative apparent-order solve, with an explicit
  convergence-regime classification (monotonic / oscillatory /
  divergent / degenerate) and the asymptotic-range diagnostic
  $R_{\mathrm{asym}}$. GCI is returned as `NaN` outside the asymptotic
  range rather than as a misleading finite number.
- **Geyer's integrated autocorrelation time $\tau_{\mathrm{int}}$**
  (initial-positive-sequence estimator) and the autocorrelation-corrected
  standard error of the mean
  $\mathrm{SEM} = \sigma\sqrt{\tau_{\mathrm{int}}/N}$,
  using the convention $\tau_{\mathrm{int}} = 1 + 2\sum_{k\ge1}\rho_k$
  (so iid data give $\tau_{\mathrm{int}}=1$ and recover $\sigma/\sqrt N$).
- **KPSS test** for stationarity of the time-averaging window (level and
  trend variants), implemented from first principles — no `statsmodels`
  dependency.
- **Extremum-localization check** for pointwise QoIs: the in-window
  spread of the extremum *location* (5th–95th percentile, in cell
  widths). A localized QoI stays within a few cells; a maximum that
  migrates between flow features is flagged as not pointwise-localized
  and demoted to a diagnostic — independent of the KPSS value check.
- **Analytical Rayleigh-Pitot reference** for cross-checking the
  Richardson-extrapolated maximum pressure independently of the GCI
  machinery.
- **LaTeX `tabular` output**, drop-in for a paper Table 1.

The motivating finding: unsteady shock-dominated CFD needs more than a single GCI number. The naive $\sigma/\sqrt N$ standard error can understate temporal sampling uncertainty when samples are serially correlated, and different extrema can behave differently under refinement. In the forward-step example, maximum pressure is the primary reference-anchored QoI, while maximum density is retained as a diagnostic QoI because its stationarity and localization behavior reveal additional shock/contact-line dynamics.

Two caveats apply when interpreting GCI on unsteady, CFL-limited runs,
spelled out in `LIMITATIONS.md`: (1) refining the mesh also refines the
time step, so the apparent order mixes spatial and temporal error
unless a fixed-Δt control run is performed; (2) a spatial extremum is a
non-smooth functional, so the Richardson expansion is heuristic for
`fieldMinMax` QoIs — integrated QoIs (forces, surface averages) are the
formally cleaner target and are on the roadmap.

## Installation

```bash
pip install git+https://github.com/rithikrn/foamgci.git
```

Or, for development:

```bash
git clone https://github.com/rithikrn/foamgci.git
cd foamgci
pip install -e ".[dev]"
pytest -v
```

Dependencies: NumPy (Python >= 3.10). Optional matplotlib for plotting.

## Quick start (CLI)

```bash
foamgci report \
    --case coarse:case_C/postProcessing/fieldMinMax/0/fieldMinMax.dat:0.025:4032 \
    --case medium:case_M/postProcessing/fieldMinMax/0/fieldMinMax.dat:0.0125:16128 \
    --case fine:case_F/postProcessing/fieldMinMax/0/fieldMinMax.dat:0.00625:64512 \
    --case extra-fine:case_XF/postProcessing/fieldMinMax/0/fieldMinMax.dat:0.003125:258048 \
    --field p --quantity max --window 3 10 \
    --reference rayleigh-pitot --mach 3 --gamma 1.4 \
    --text out/report.txt --latex out/table1.tex
```

Each `--case` is `label:path:h[:n_cells]`. List cases coarse-to-fine
(`h` strictly decreasing). `--text` / `--latex` are optional; without
them the report just prints to the terminal. (The CLI does not draw
figures — use the Python API or an example's figure scripts for that.)

## Quick start (Python API)

```python
from foamgci import GridCase, full_report, rayleigh_pitot

rep = full_report(
    cases=[
        GridCase("coarse",     "case_C/.../fieldMinMax.dat", h=0.025,    n_cells=4032),
        GridCase("medium",     "case_M/.../fieldMinMax.dat", h=0.0125,   n_cells=16128),
        GridCase("fine",       "case_F/.../fieldMinMax.dat", h=0.00625,  n_cells=64512),
        GridCase("extra-fine", "case_XF/.../fieldMinMax.dat",h=0.003125, n_cells=258048),
    ],
    field="p", quantity="max", window=(3.0, 10.0),
    reference_value=rayleigh_pitot(3.0, 1.4),
    reference_label="Rayleigh-Pitot M=3",
)
print(rep.as_text())
print(rep.as_latex())
```

## The `fieldMinMax` function object

Add this block to each case's `system/controlDict`. The function-object
name (`fieldMinMax`) sets the output folder, so the file lands at
`postProcessing/fieldMinMax/0/fieldMinMax.dat`:

```cpp
functions
{
    fieldMinMax
    {
        type            fieldMinMax;
        libs            ("libfieldFunctionObjects.so");
        mode            magnitude;
        location        true;        // record (x y z) of each extremum
        writeControl    timeStep;
        writeInterval   100;
        fields          (p rho);
    }
}
```

This writes one row per sampled timestep with the spatial min/max of each
field and the location of each extremum. `foamgci.reader.read_fieldminmax`
handles both OpenFOAM dialects (combined multi-field and per-field).

## Repository layout

```
foamgci/
├── foamgci/                       # THE LIBRARY (pip-installable, case-agnostic)
│   ├── _version.py                #   single source of the version string
│   ├── __init__.py                #   public API
│   ├── __main__.py                #   CLI: foamgci report ...
│   ├── reader.py                  #   OpenFOAM fieldMinMax + generic scalar time-series readers
│   ├── stats.py                   #   Geyer tau_int, KPSS, window stats
│   ├── gci.py                     #   Roache GCI on triplets
│   ├── report.py                  #   end-to-end driver + Rayleigh-Pitot
│   └── plot.py                    #   optional matplotlib figure
├── tests/                         # pytest suite, anchored to Celik (2008)
├── examples/
│   └── forwardstep_mach3/         # ONE worked case (the template)
│       ├── 0/ constant/ system/   #   the committed OpenFOAM case (fine grid)
│       ├── submit.sh              #   SLURM runner
│       ├── README.md              #   how to run THIS case end-to-end
│       └── gci/                   #   THIS case's analysis driver
│           ├── data.py            #     grid metadata (edit this per case)
│           ├── analyze.py         #     reads gci/data/*.dat -> gci_summary.json
│           ├── make_*.py          #     paper figures
│           ├── run_all.sh         #     analyze + figures
│           └── data/              #     expected coarse/medium/fine/extrafine QoI inputs
├── README.md  LIMITATIONS.md  CONTRIBUTING.md  CHANGELOG.md  LICENSE
├── pyproject.toml
└── .github/workflows/tests.yml    # CI: Linux + macOS, Py 3.10-3.12
```

## Worked example

`examples/forwardstep_mach3/` is the Mach-3 Woodward-Colella
forward-facing step, used as the template for every future case. It
ships the full OpenFOAM case (the fine grid) and the analysis driver.
The four `fieldMinMax.dat` inputs under `gci/data/` are committed, so the analysis can be reproduced from a fresh clone without rerunning OpenFOAM. Running the OpenFOAM cases is only needed if you want to regenerate the input data.

See **`examples/forwardstep_mach3/README.md`** for step-by-step run
instructions. The short version:

```bash
pip install -e .                       # install the library
cd examples/forwardstep_mach3/gci
bash run_all.sh                           # -> gci_summary.json + figures/
```

| label       | cells   | $h$       | $r$ |
| ----------- | ------: | --------- | --- |
| coarse      |   4,032 | 0.025     | --  |
| medium      |  16,128 | 0.0125    | 2   |
| fine        |  64,512 | 0.00625   | 2   |
| extra-fine  | 258,048 | 0.003125  | 2   |

## Output format

`foamgci report` prints a per-grid statistics table and the GCI block.
The numbers below are a synthetic illustration (not solver output); your
run's exact values depend on the solver, scheme, and time step you used.
Note that `N` roughly doubles per refinement level: with
`writeControl timeStep` the sampling cadence follows the CFL-limited
time step, which halves with `h`.

```
========================================================================
foamgci V&V report — field 'p', quantity 'max', window [3, 10]
========================================================================

Per-grid time-averaged statistics:
  label        N_cells       h     N    mean     std  tau_int     SEM  N_eff   KPSS_p
  coarse          4032   0.025   110  11.986  0.0205     3.37  0.0036     33  >=0.100
  medium         16128  0.0125   210  12.045  0.0205     2.21  0.0021     95  >=0.100
  fine           64512 0.00625   420  12.074  0.0192     2.13  0.0014    197  >=0.100
  extra-fine    258048 0.003125  840  12.083  0.0203     3.70  0.0013    227  >=0.100

Roache GCI (triplet medium, fine, extra-fine):
      regime                   = monotonic
      apparent order p-hat     = 1.72
      Richardson phi_exact     = 12.0871
      GCI_fine_21              = 0.0401 %
      asymptotic ratio (~1)    = 1.000
========================================================================
```

## Citing

If you use `foamgci` in published work, please cite it. A plain-text
acknowledgement:

> Grid-convergence study performed with `foamgci` (Roache GCI + Geyer
> tau_int + KPSS stationarity), https://github.com/rithikrn/foamgci.

Or as BibTeX:

```bibtex
@software{rithik_r_nambiar_2026_20778444,
  author       = {Rithik R Nambiar},
  title        = {rithikrn/foamgci},
  month        = jun,
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v3.2.2},
  doi          = {10.5281/zenodo.20778444},
  url          = {https://doi.org/10.5281/zenodo.20778444},
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

MIT — see `LICENSE`.
