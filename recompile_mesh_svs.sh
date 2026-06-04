# Name of the system where MESH-SVS is compiled
# Science: internal ECCC network (use ECCC Gitlab to retrieve SVS and SVS2 code)
# GPSCC: ECCC collaboration server (use Github to retrieve SVS and SVS2 code)
# Other: other machine (use Github to retrieve SVS and SVS2 code)
system=Science # Science | GPSCC | Other

# Load compiler
cd ../sps
# Load compiler
if [ "$system" = "Science" ] || [ "$system" = "GPSCC"  ]; then	
   . .eccc_setup_intel
elif [ "$system" = "Other" ]; then
   . .common_setup gnu	
fi

# Go to sps build dir
cd ../sps_build

# ReCompile rpn physics
make rpnphy -j4

# Compile MESH
cd ../MESH_SVS
if [ "$system" = "Science" ] || [ "$system" = "GPSCC"  ]; then	
   make mpi_intel debug
elif [ "$system" = "Other" ]; then
   make mpi_gcc debug
fi
