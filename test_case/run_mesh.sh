# Script to run MESH-SVS
# Usage: ./run_mesh.sh PATH_TO_CODE

here=`pwd -P`

# Create symbolic link for mpi_sa_mesh
ln -sf $1/MESH_SVS/mpi_sa_mesh . 

# Load ECCC environment
source $1/sps/.eccc_setup_intel

# Back to the run directory and run mpi_sa_mesh
cd $here
./mpi_sa_mesh
