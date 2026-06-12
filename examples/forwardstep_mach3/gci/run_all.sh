#!/bin/bash
# Post-processing: parse fieldMinMax -> GCI -> figures.
# Run after the four cases finish. Requires numpy + matplotlib
# (and pyvista only if regenerating the t=4 snapshots/contours).
set -e
cd "$(dirname "$0")"

python3 analyze.py            # -> gci_summary.json
python3 make_figures.py       # -> figures/fig_grid_convergence.{pdf,png}
python3 make_aux_figures.py   # -> figures/fig_domain.pdf, fig_peak_location.pdf
[ -d snapshots ] && python3 make_contour_figure.py   # -> fig_density_contours.pdf
echo "Done. See gci_summary.json and figures/."
