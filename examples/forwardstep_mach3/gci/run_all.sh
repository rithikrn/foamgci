#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16          # must equal numberOfSubdomains in decomposeParDict
#SBATCH --time=01:00:00               # generous; this case should finish in 15-45 min
#SBATCH --mem=16G                     # 12k-cell 2D Euler: this is already generous
# #SBATCH --gres=gpu:v100:1

#SBATCH --job-name="gci_study"
#SBATCH --mail-user=rithikrn@iastate.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output="gci.out"
#SBATCH --error="gci.err"
# =====================================================================
# run_all.sh -- Slurm-submittable four-grid GCI pipeline for Nova HPC.
#
# Submit:        sbatch run_all.sh
# Interactive:   bash run_all.sh            (SBATCH lines are comments)
#
# Flags:
#   --no-snapshots   skip PyVista extraction + contour figure (fast)
#   --only-figures   regen figures from existing gci_summary.json
#   --only-tex       rebuild verification.pdf only
#   --skip-setup     don't touch modules or pip; use current env as-is
#
# Layout assumed (data.py is the authoritative spec):
#     forwardstep_gridstudy/
#         coarse_grid/coarsecase.foam        + postProcessing/fieldRange/0/
#         medium_grid/mediumcase.foam        + ...
#         fine_grid/finecase.foam            + ...
#         extrafine_grid/extrafinecase.foam  + ...
#         gci/                               <-- this script lives here
# =====================================================================

set -euo pipefail

# -- Locate ourselves --------------------------------------------------
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    cd "$SLURM_SUBMIT_DIR"
else
    cd "$(dirname "${BASH_SOURCE[0]}")"
fi

echo "==========================================================="
echo "  Job started: $(date)"
echo "  Host:        $(hostname)"
echo "  Working dir: $(pwd)"
[[ -n "${SLURM_JOB_ID:-}" ]] && echo "  Slurm job:   $SLURM_JOB_ID"
echo "==========================================================="

# -- Argument parsing --------------------------------------------------
WITH_SNAPSHOTS=1        # ON by default when batch-submitted
ONLY_FIGURES=0
ONLY_TEX=0
SKIP_SETUP=0
for arg in "$@"; do
    case "$arg" in
        --no-snapshots) WITH_SNAPSHOTS=0 ;;
        --only-figures) ONLY_FIGURES=1 ;;
        --only-tex)     ONLY_TEX=1 ;;
        --skip-setup)   SKIP_SETUP=1 ;;
        -h|--help)
            sed -n '/^# ===/,/^# ===$/p' "$0" | sed 's/^# \?//'
            exit 0 ;;
        *)
            echo "Unknown argument: $arg" >&2
            exit 2 ;;
    esac
done

# ---------------------------------------------------------------------
# 0.  Module + Python environment
# ---------------------------------------------------------------------
load_first() {
    # Try each candidate module name; load the first that succeeds.
    local mod
    for mod in "$@"; do
        if module load "$mod" 2>/dev/null; then
            echo "      loaded: $mod"
            return 0
        fi
    done
    return 1
}

if [[ $SKIP_SETUP -eq 0 ]]; then
    echo
    echo ">>>  Step 0/4: environment setup"

    # `module` is a shell function on Lmod systems; use `type` not `command`.
    if [[ -n "$(type -t module 2>/dev/null)" ]]; then
        module purge 2>/dev/null || true
        # ----- Adjust these candidate lists if Nova's module names differ -----
        load_first python/3.11.5 python/3.11 python/3.10 python/3.9 python \
                   miniconda3 anaconda3 \
            || echo "      WARN: no python module loaded; relying on PATH"
        load_first texlive texlive/2023 texlive/2022 \
            || echo "      WARN: no texlive module loaded"
        echo "      modules currently loaded:"
        module list 2>&1 | sed 's/^/        /'
    else
        echo "      no 'module' command found -- assuming Python in PATH"
    fi

    PY=${PY:-python3}
    if ! command -v "$PY" >/dev/null 2>&1; then
        echo "ERROR: $PY not found in PATH after module load." >&2
        echo "       Try: module load python/3.11 (or similar)" >&2
        exit 1
    fi
    echo "      Python: $($PY --version) at $(command -v $PY)"

    # Ensure required packages are available.
    REQ=(numpy matplotlib scipy)
    if [[ $WITH_SNAPSHOTS -eq 1 && $ONLY_FIGURES -eq 0 && $ONLY_TEX -eq 0 ]]; then
        REQ+=(pyvista)
    fi
    MISSING=()
    for pkg in "${REQ[@]}"; do
        if ! $PY -c "import $pkg" 2>/dev/null; then
            MISSING+=("$pkg")
        fi
    done
    if [[ ${#MISSING[@]} -gt 0 ]]; then
        echo "      installing missing packages: ${MISSING[*]}"
        $PY -m pip install --user --quiet --upgrade pip 2>/dev/null || true
        if ! $PY -m pip install --user --quiet "${MISSING[@]}"; then
            echo "ERROR: pip install failed. Install manually with:" >&2
            echo "    $PY -m pip install --user ${MISSING[*]}" >&2
            exit 1
        fi
    else
        echo "      all required packages present"
    fi
fi

PY=${PY:-python3}

# Headless rendering (compute nodes have no X11).
export PYVISTA_OFF_SCREEN=true
export MPLBACKEND=Agg
unset DISPLAY 2>/dev/null || true

# Cap thread oversubscription on shared nodes.
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS
export MKL_NUM_THREADS=$OMP_NUM_THREADS

# ---------------------------------------------------------------------
# 1.  Parse fieldMinMax + Roache GCI
# ---------------------------------------------------------------------
if [[ $ONLY_FIGURES -eq 0 && $ONLY_TEX -eq 0 ]]; then
    echo
    echo ">>>  Step 1/4: parse fieldMinMax + run GCI"
    "$PY" analyze.py
fi

# ---------------------------------------------------------------------
# 2.  PyVista snapshot extraction at t = 4
# ---------------------------------------------------------------------
if [[ $WITH_SNAPSHOTS -eq 1 && $ONLY_TEX -eq 0 && $ONLY_FIGURES -eq 0 ]]; then
    echo
    echo ">>>  Step 2/4: extract t=4 snapshots (PyVista)"
    "$PY" extract_snapshot.py --time 4.0 --nx 600 --ny 200
fi

# ---------------------------------------------------------------------
# 3.  Figures
# ---------------------------------------------------------------------
if [[ $ONLY_TEX -eq 0 ]]; then
    echo
    echo ">>>  Step 3/4: generate figures"
    "$PY" make_figures.py
    "$PY" make_aux_figures.py
    if compgen -G "snapshots/snap_*.npz" > /dev/null; then
        "$PY" make_contour_figure.py
    else
        echo "      (no snapshots cached; skipping contour figure)"
    fi
fi

# ---------------------------------------------------------------------
# 4.  Compile LaTeX
# ---------------------------------------------------------------------
echo
echo ">>>  Step 4/4: compile verification.tex"
if command -v pdflatex >/dev/null 2>&1; then
    pdflatex -interaction=nonstopmode verification.tex > pdflatex.log 2>&1
    pdflatex -interaction=nonstopmode verification.tex > pdflatex.log 2>&1
    if [[ -f verification.pdf ]]; then
        echo "      wrote verification.pdf  ($(du -h verification.pdf | cut -f1))"
    else
        echo "      pdflatex failed; last 50 lines of pdflatex.log:" >&2
        tail -50 pdflatex.log >&2
        exit 3
    fi
else
    echo "      pdflatex not found; skipping LaTeX compile"
    echo "      (load a texlive module to enable, then re-run with --only-tex)"
fi

# ---------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------
echo
echo "==========================================================="
echo "  Outputs:"
[[ -f gci_summary.json ]] && \
    echo "    gci_summary.json  ($(du -h gci_summary.json | cut -f1))"
for f in figures/*.pdf; do
    [[ -e "$f" ]] && echo "    $f  ($(du -h "$f" | cut -f1))"
done
[[ -d snapshots ]] && \
    echo "    snapshots/        ($(ls snapshots/ 2>/dev/null | wc -l) NPZ files)"
[[ -f verification.pdf ]] && \
    echo "    verification.pdf  ($(du -h verification.pdf | cut -f1))"
echo
echo "  Finished: $(date)"
echo "==========================================================="