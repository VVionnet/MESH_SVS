
#tag_sps=6.3.0-a17
branch_sps=630a18_fora20

# Change dir
cd ../

# Create code structure
mkdir sps sps_build

# Extract sps code from gitlab
cd sps
git clone --no-checkout git@gitlab.science.gc.ca:continental-surface-hydrology/sps-dev.git .
#git checkout $tag_sps
git checkout $branch_sps
git submodule update --init --recursive

# Load compiler
. .eccc_setup_intel

# Go to build dir
cd ../sps_build

# Compile rpn physics
cmake ../sps
make rpnphy -j4

# Compile MESH
cd ../MESH_SVS
make clean
make mpi_intel debug


