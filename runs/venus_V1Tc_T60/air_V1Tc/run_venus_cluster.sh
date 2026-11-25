#!/bin/bash
#SBATCH --job-name=venus_v1tc
#SBATCH --nodes=5
#SBATCH --ntasks=168
#SBATCH --partition=t2standard
#SBATCH --time=12:00:00
#SBATCH --output=log.txt
# submit from inside this folder:  sbatch run_venus_cluster.sh

set -e
if [ ! -f DATA/Par_file ]; then
  echo "run this script from inside the run folder (DATA/ not found)"
  exit 1
fi
mkdir -p OUTPUT_FILES MESH bin
rm -rf OUTPUT_FILES/*
ln -sf /home/irseppi/REPOS/specfem2d/bin/xmeshfem2D bin/xmeshfem2D
ln -sf /home/irseppi/REPOS/specfem2d/bin/xspecfem2D bin/xspecfem2D

# the solver must be built with fortran/define_external_model.f90
if [ ! -x bin/xspecfem2D ]; then
  echo "ERROR: /home/irseppi/REPOS/specfem2d/bin/xspecfem2D not found"
  exit 1
fi
if ! grep -aq "Venus 1-D table" bin/xspecfem2D; then
  echo "ERROR: /home/irseppi/REPOS/specfem2d/bin/xspecfem2D was not built with the Venus model"
  echo "reader. Copy fortran/define_external_model.f90 into"
  echo "src/specfem2D/ of that SPECFEM2D and run 'make xspecfem2D'"
  echo "(step 1 of the README)."
  exit 1
fi

# partition for the number of processes used on this machine
sed "s/^NPROC *=.*/NPROC = 168/" DATA/Par_file > DATA/Par_file.tmp
mv DATA/Par_file.tmp DATA/Par_file

echo "building the mesh: $(date)"
gfortran -O3 -o bin/xcreate_mesh create_mesh_venus.F90
./bin/xcreate_mesh > OUTPUT_FILES/output_create_mesh.txt

cp DATA/Par_file DATA/SOURCE DATA/STATIONS DATA/venus_model.txt \
   OUTPUT_FILES/
echo "partitioning: $(date)"
mpirun -np 168 ./bin/xmeshfem2D
echo "solver: $(date)"
mpirun -np 168 ./bin/xspecfem2D
echo "done: $(date)"
