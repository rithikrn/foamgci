# Contributing

Bug reports and PRs welcome. The project is deliberately small —
new functionality should be argued from a concrete CFD V&V use
case, not from generality.

## Development setup

```bash
git clone https://github.com/USERNAME/foamgci.git
cd foamgci
python -m pip install -e ".[dev]"
pytest -v
```

## Running tests locally

The suite is pure-Python NumPy and runs in a few seconds:

```bash
pytest -v                # full suite
pytest tests/test_gci.py # one file
pytest -k celik          # by keyword
```

CI runs on Linux + macOS, Python 3.10–3.12, on every push and PR.

## Style

- PEP 8 with line length 88 (Black-compatible).
- Type-annotated public APIs.
- Docstrings cite the V&V or statistics literature for every
  non-obvious formula. If you are adding a new method, also add a
  reference in the relevant module docstring at the top of the file.
- Tests anchored to a published numerical benchmark whenever
  possible (e.g. Celik et al. 2008 Example 1 for GCI).

## Filing an issue

Include:
- A short `fieldMinMax.dat` excerpt (the first ~20 lines is plenty)
- The command line you ran
- The OpenFOAM version, solver name, and `system/controlDict` excerpt
  (`startTime`, `endTime`, `deltaT`, `writeControl`, `writeInterval`,
  and the `fieldMinMax` function-object block)

## Pull-request checklist

- [ ] Tests added or updated and the suite passes locally.
- [ ] `CHANGELOG.md` entry under an `## [Unreleased]` heading.
- [ ] No new runtime dependencies without discussion (the
      `numpy`-only core is a deliberate design choice).
