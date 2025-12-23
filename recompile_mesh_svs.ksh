# Load compiler
cd ../sps
. .eccc_setup_intel_2022.1.2

# Go to sps build dir
cd ../sps_build

# ReCompile rpn physics
make rpnphy -j4

# Compile MESH
cd ../MESH_SVS
make mpi_intel debug
