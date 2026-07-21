# Contributing

PRs and bug reports welcome. It's a small tool on purpose, so new features
should come from a real CFD V&V need, not from "it'd be more general."

## Setup

```bash
pip install -e ".[dev]"
pytest -q
```

## Style

- PEP 8, line length 88 (Black-friendly).
- Type-annotate public functions.
- If a formula isn't obvious, cite the paper in the module docstring.
- New methods get a test anchored to a published benchmark where possible
  (e.g. Celik 2008 for GCI).

## Filing an issue

Paste the first ~20 lines of your `.dat`, the command you ran, and your
OpenFOAM version + the relevant `controlDict` block.

## PR checklist

- [ ] Tests pass locally.
- [ ] `CHANGELOG.md` entry.
- [ ] No new runtime deps without discussing it (NumPy-only core is deliberate).
