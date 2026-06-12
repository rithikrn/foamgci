#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16          # must equal numberOfSubdomains in decomposeParDict
#SBATCH --time=01:00:00               # generous; this case should finish in 15-45 min
#SBATCH --mem=16G                     # 12k-cell 2D Euler: this is already generous

#SBATCH --job-name="foamgci"
#SBATCH --output="foamgci-%j.out"
#SBATCH --error="foamgci-%j.err"

set -e
module purge
module load openfoam

blockMesh            2>&1 | tee log.blockMesh
checkMesh -allGeometry -allTopology 2>&1 | tee log.checkMesh
decomposePar         2>&1 | tee log.decompose
mpirun -np "$SLURM_NTASKS" rhoCentralFoam -parallel  2>&1 | tee log.run
reconstructPar       2>&1 | tee log.reconstruct
rm -rf processor*
