# foamgci

**Autocorrelation-aware verification for unsteady OpenFOAM CFD.**

`foamgci` reads OpenFOAM `fieldMinMax.dat` output and produces a
complete grid-convergence-index (GCI) report including:

- **Roache GCI** on every consecutive refinement triplet, with the
  Celik et al. (2008) iterative apparent-order solve and the
  asymptotic-range diagnostic $R_{\mathrm{asym}}$.
- **Geyer's integrated autocorrelation time $\tau_{\mathrm{int}}$**
  and the autocorrelation-corrected standard error of the mean
  $\sigma\sqrt{2\tau_{\mathrm{int}}/N}$.
- **KPSS test** for stationarity of the time-averaging window
  (level and trend variants), implemented from first principles —
  no `statsmodels` dependency.
- **Analytical Rayleigh-Pitot reference** for cross-checking the
  Richardson-extrapolated maximum pressure independently of the
  GCI machinery.
- **LaTeX `tabular` output** drop-in for a paper Table 1.

The motivating finding is that the naïve $\sigma/\sqrt{N}$ standard
error understates temporal sampling uncertainty for shock-dominated
unsteady flows by 2–3×. Once $\tau_{\mathrm{int}}$ is accounted for,
the apparent grid-refinement increment between the two finest grids
is often smaller than the sampling noise — i.e. the residual is
temporal, not spatial.

## Installation

```bash
pip install git+https://github.com/USERNAME/foamgci.git
```

Or, for development:

```bash
git clone https://github.com/USERNAME/foamgci.git
cd foamgci
pip install -e ".[dev]"
pytest -v
```

Dependencies: NumPy. Optional matplotlib for the convergence plot.
Python ≥ 3.10.

## Quick start

```bash
foamgci report \
    --case coarse:case_C/postProcessing/fieldMinMax/0/fieldMinMax.dat:0.025:4032 \
    --case medium:case_M/postProcessing/fieldMinMax/0/fieldMinMax.dat:0.0125:16128 \
    --case fine:case_F/postProcessing/fieldMinMax/0/fieldMinMax.dat:0.00625:64512 \
    --case extra-fine:case_XF/postProcessing/fieldMinMax/0/fieldMinMax.dat:0.003125:258048 \
    --field p --quantity max --window 3 10 \
    --reference rayleigh-pitot --mach 3 --gamma 1.4 \
    --latex out/table1.tex
```

Each `--case` is `label:path:h[:n_cells]`. Cases must be listed
coarse-to-fine (`h` strictly decreasing) with a constant refinement
ratio.

The same workflow from Python:

```python
from foamgci import GridCase, full_report, rayleigh_pitot

rep = full_report(
    cases=[
        GridCase("coarse",     "case_C/.../fieldMinMax.dat", h=0.025,    n_cells=4032),
        GridCase("medium",     "case_M/.../fieldMinMax.dat", h=0.0125,   n_cells=16128),
        GridCase("fine",       "case_F/.../fieldMinMax.dat", h=0.00625,  n_cells=64512),
        GridCase("extra-fine", "case_XF/.../fieldMinMax.dat",h=0.003125, n_cells=258048),
    ],
    field="p",
    quantity="max",
    window=(3.0, 10.0),
    reference_value=rayleigh_pitot(3.0, 1.4),
    reference_label="Rayleigh-Pitot M=3",
)
print(rep.as_text())
print(rep.as_latex())
```

## What `fieldMinMax` should look like

Add this function-object block to each case's `system/controlDict`:

```cpp
functions
{
    fieldMinMax
    {
        type            fieldMinMax;
        libs            ("fieldFunctionObjects");
        writeControl    timeStep;
        writeInterval   100;
        fields          (p U rho);
        location        true;
        mode            magnitude;
    }
}
```

This produces `postProcessing/fieldMinMax/0/fieldMinMax.dat` containing
one row per write timestep with the spatial min/max of each field and
the location at which the extremum occurs. `foamgci.reader.read_fieldminmax`
handles both OpenFOAM dialects (combined multi-field and per-field).

## Repository layout

```
foamgci/
├── foamgci/                          # the package
│   ├── __init__.py                   # public API
│   ├── __main__.py                   # CLI: python -m foamgci report ...
│   ├── reader.py                     # fieldMinMax.dat parser
│   ├── stats.py                      # Geyer τ_int, KPSS, window stats
│   ├── gci.py                        # Roache GCI on triplets
│   ├── report.py                     # end-to-end driver + Rayleigh-Pitot
│   └── plot.py                       # optional matplotlib figure
├── tests/                            # pytest suite, anchored to Celik 2008
│   ├── test_reader.py
│   ├── test_stats.py
│   ├── test_gci.py
│   ├── test_report.py
│   └── test_cli.py
├── examples/
│   └── forwardstep_mach3/            # worked example for the abstract
│       ├── README.md
│       ├── verify_abstract.py        # regenerates Table 1 from real data
│       └── data/README.md
├── README.md                         # this file
├── LIMITATIONS.md                    # what foamgci does not do
├── CONTRIBUTING.md
├── CHANGELOG.md
├── LICENSE                           # MIT
├── pyproject.toml
├── .gitignore
└── .github/workflows/tests.yml       # CI on Linux + macOS, Py 3.10–3.12
```

## Worked example: Mach-3 forward-facing step

The `examples/forwardstep_mach3/` directory contains the verification
driver for the four-grid Woodward–Colella case used in the AIAA
SciTech 2027 abstract. Cell counts (from `log.checkMesh`):

| label       | cell count | $h$       | $r$  |
| ----------- | ----------:| --------- | ---- |
| coarse      |      4,032 | 0.025     | —    |
| medium      |     16,128 | 0.0125    | 2    |
| fine        |     64,512 | 0.00625   | 2    |
| extra-fine  |    258,048 | 0.003125  | 2    |

To regenerate the abstract's Table 1 against your own runs:

```bash
python examples/forwardstep_mach3/verify_abstract.py \
    --coarse     case_C/postProcessing/fieldMinMax/0/fieldMinMax.dat \
    --medium     case_M/postProcessing/fieldMinMax/0/fieldMinMax.dat \
    --fine       case_F/postProcessing/fieldMinMax/0/fieldMinMax.dat \
    --extra-fine case_XF/postProcessing/fieldMinMax/0/fieldMinMax.dat
```

This writes `out/report.txt`, `out/table1.tex`, and
`out/fig_convergence.pdf`. Paste `table1.tex` directly into the
abstract.

The example does **not** ship the `fieldMinMax.dat` files themselves;
they are tens of MB per case and the abstract's Table 1 is regenerated
from your runs, not cached.

## Output format

`foamgci report` prints something like:

```
========================================================================
foamgci V&V report — field 'p', quantity 'max', window [3, 10]
========================================================================

Per-grid time-averaged statistics:

  label        N_cells         h      N        mean       std    τ_int       SEM   N_eff   KPSS_p
  ----------------------------------------------------------------------------------------------
  coarse          4032     0.025    175     11.9950    0.0020     1.84    0.000218  47.6    0.150
  medium         16128    0.0125    175     12.0150    0.0015     2.21    0.000189  39.6    0.143
  fine           64512   0.00625    175     12.0500    0.0020     2.85    0.000287  30.7    0.165
  extra-fine    258048  0.003125    175     12.0580    0.0030     3.42    0.000462  25.6    0.158

Roache GCI on consecutive triplets (coarse → medium → fine):

  triplet (coarse, medium, fine):
      apparent order p̂        = 0.823
      Richardson φ_exact       = 12.0857
      GCI_fine_21              = 0.4318 %
      GCI_medium_32            = 0.2071 %
      asymptotic ratio (≈1)    = 0.937

  triplet (medium, fine, extra-fine):
      apparent order p̂        = 2.131
      Richardson φ_exact       = 12.0603
      GCI_fine_21              = 0.0207 %
      GCI_medium_32            = 0.0902 %
      asymptotic ratio (≈1)    = 1.001

Analytical reference (Rayleigh-Pitot p_02/p_1 (M=3)):
      reference value         = 12.061
      Richardson extrapolation = 12.0603  (error = 0.0058 %)
      finest-grid mean         = 12.058   (error = 0.0249 %)
========================================================================
```

The above numbers are illustrative; your run's exact values will
differ slightly from cell-by-cell variability in the random initial
phase of the shear-layer roll-up.

## Citing

If `foamgci` contributed to a publication, a citation in the
verification section is appreciated:

> Mach-3 forward-step grid-convergence study performed with
> `foamgci` v0.2.0 (Roache GCI + Geyer τ_int + KPSS stationarity),
> https://github.com/rithikrn/foamgci.

A JOSS paper is planned. The current preferred citation form is the
GitHub URL + commit hash for the release used.

## References

The implementation is anchored to:

- Roache, P.J. (1994), "Perspective: A method for uniform reporting
  of grid refinement studies", *J. Fluids Eng.*, **116**, 405–413.
- Celik, I.B., Ghia, U., Roache, P.J., Freitas, C.J., Coleman, H.,
  Raad, P.E. (2008), "Procedure for estimation and reporting of
  uncertainty due to discretization in CFD applications",
  *J. Fluids Eng.*, **130**, 078001.
- Geyer, C.J. (1992), "Practical Markov chain Monte Carlo",
  *Statistical Science*, **7**, 473–483.
- Kwiatkowski, D., Phillips, P.C.B., Schmidt, P., Shin, Y. (1992),
  "Testing the null hypothesis of stationarity against the
  alternative of a unit root", *J. Econometrics*, **54**, 159–178.
- Woodward, P., Colella, P. (1984), "The numerical simulation of
  two-dimensional fluid flow with strong shocks", *J. Comput. Phys.*,
  **54**, 115–173.
- Greenshields, C.J., Weller, H.G., Gasparini, L., Reese, J.M.
  (2010), "Implementation of semi-discrete, non-staggered central
  schemes in a colocated, polyhedral, finite volume framework",
  *Int. J. Numer. Meth. Fluids*, **63**, 1–21.

## License

MIT — see `LICENSE`.
