
tag_sps=630a20_vvi001_surface_fora22
system=Science # Science | GPSCC
activate_debug=true

# Change dir
cd ../

# Create code structure
mkdir sps sps_build

# Extract sps code from gitlab
cd sps

# Clone based on system
if [ "$system" = "Science" ]; then
    git clone --no-checkout git@gitlab.science.gc.ca:continental-surface-hydrology/sps-dev.git .
else
    git clone --branch 6.3 git@github.com:VVionnet/sps_dev.git .
fi

git checkout $tag_sps
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

if [ "$activate_debug" = true ]; then
   make mpi_intel debug
else
   make mpi_intel
fi

