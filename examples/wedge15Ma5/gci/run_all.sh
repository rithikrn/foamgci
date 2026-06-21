#!/bin/bash
# Post-processing: parse surfaceFieldValue -> GCI vs analytical oblique shock,
# then figures. Run after the four grid cases finish and their .dat files are
# in data/. Requires numpy + matplotlib (and pyvista only if regenerating
# snapshots for the contour figure).
set -e
cd "$(dirname "$0")"

python3 analyze.py            # -> gci_summary.json
python3 make_figures.py       # -> figures/fig_grid_convergence.{pdf,png}
python3 make_aux_figures.py   # -> figures/fig_domain.pdf
[ -d snapshots ] && python3 make_contour_figure.py   # -> fig_pressure_contours.pdf
echo "Done. See gci_summary.json and figures/."
