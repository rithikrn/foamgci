#!/bin/bash
# Post-processing: parse the per-grid QoI .dat files in data/ into a GCI
# summary and figures. Run after the five grid cases finish and their .dat
# files are copied into gci/data/ as <prefix>_<qoi>.dat (see data.py for the
# exact names). Requires numpy + matplotlib.
set -e
cd "$(dirname "$0")"
python3 analyze.py        # -> gci_summary.json
python3 make_figures.py   # -> figures/fig_convergence,_coupling,_meshlock.pdf
echo "Done. See gci_summary.json and figures/."
