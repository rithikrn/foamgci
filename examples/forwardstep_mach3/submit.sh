#!/bin/bash
# Submit with: sbatch submit.sh
# Monitor with: squeue -u $USER  ;  tail -f fs_base-<jobid>.out

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16          # must equal numberOfSubdomains in decomposeParDict
#SBATCH --time=01:00:00               # generous; this case should finish in 15-45 min
#SBATCH --mem=16G                     # 12k-cell 2D Euler: this is already generous
# #SBATCH --gres=gpu:v100:1

#SBATCH --job-name="fs_base"
#SBATCH --mail-user=rithikrn@iastate.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output="fs_base-%j.out"
#SBATCH --error="fs_base-%j.err"

# ----- safety: abort on any error, undefined var, or pipe failure -----
set -euo pipefail

# ----- environment -----
module purge
module load openfoam

# ----- sanity checks before we burn allocation -----
echo "=== Job $SLURM_JOB_ID on $SLURM_JOB_NODELIST ==="
echo "Tasks: $SLURM_NTASKS  |  Case: $(pwd)"
which rhoCentralFoam || { echo "ERROR: rhoCentralFoam not in PATH"; exit 1; }

# ----- defensive cleanup in case of re-submission over old run -----
rm -rf processor* constant/polyMesh 0.0* 0.[1-9]* [1-9]* \
       log.blockMesh log.checkMesh log.decomposePar log.run log.reconstruct || true

# ----- workflow -----
blockMesh        2>&1 | tee log.blockMesh
checkMesh        2>&1 | tee log.checkMesh
decomposePar     2>&1 | tee log.decomposePar
mpirun -np "$SLURM_NTASKS" rhoCentralFoam -parallel 2>&1 | tee log.run
reconstructPar   2>&1 | tee log.reconstruct
rm -rf processor*

echo "=== DONE: $(date) ==="
