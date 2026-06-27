#!/bin/bash
#SBATCH --job-name=diamond2D
#SBATCH --nodes=1
#SBATCH --ntasks=16            # must match numberOfSubdomains in decomposeParDict
#SBATCH --time=08:00:00
#SBATCH --partition=short
# Mach-2 diamond airfoil, one grid level. Set the four (nx ny 1) tuples in
# system/blockMeshDict per the table in its header, then submit one job per
# grid (coarse..ultrafine). After each run, copy the function-object .dat into
# gci/data/ with the matching <prefix> (see gci/data.py).

set -e
module load openfoam/4.x || true     # adjust to your site
. "$WM_PROJECT_DIR/etc/bashrc" 2>/dev/null || true

blockMesh
checkMesh                 # VERIFY: no negative volumes, faces match patches
decomposePar -force
mpirun -np 16 rhoCentralFoam -parallel
reconstructPar -latestTime

# collect the QoI time-series (paths vary slightly by OpenFOAM version)
PREFIX=${1:-coarse}
mkdir -p gci/data
cp postProcessing/forceCoeffs1/0/forceCoeffs.dat   gci/data/${PREFIX}_forceCoeffs.dat
cp postProcessing/pFront/0/surfaceFieldValue.dat   gci/data/${PREFIX}_pFront.dat
cp postProcessing/pRear/0/surfaceFieldValue.dat    gci/data/${PREFIX}_pRear.dat
cp postProcessing/fieldMinMax/0/fieldMinMax.dat    gci/data/${PREFIX}_fieldMinMax.dat
cp postProcessing/sVol/0/volFieldValue.dat         gci/data/${PREFIX}_sVol.dat 2>/dev/null || \
  echo "  (no sVol.dat - coded entropy FO did not run; other QoIs unaffected)"
echo "grid '${PREFIX}' done; .dat copied to gci/data/"
