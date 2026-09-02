# Name of the system where MESH-SVS is compiled
# Science: internal ECCC network (use ECCC Gitlab to retrieve SVS and SVS2 code)
# GPSCC: ECCC collaboration server (use Github to retrieve SVS and SVS2 code)
# Other: other machine (use Github to retrieve SVS and SVS2 code)
system=GPSCC # Science | GPSCC | Other

# Tag of SPS version or name of SPS branch to be extracted from reference SPS repository on Gitlab or Gitbub
# If tag_sps_user is not speficied, the most recent branch is used as a default. 
#tag_sps_user=630a20_vvi001_surface_fora22

# Compile in debug mode
activate_debug=true

# Extract sps code from gitlab
cd ../sps

# Load compiler
if [ "$system" = "Science" ] || [ "$system" = "GPSCC"  ]; then	
   . .eccc_setup_intel
elif [ "$system" = "Other" ]; then
   . .common_setup gnu	
fi

# Go to build dir
cd ../sps_build

# Compile rpn physics
cmake --debug-trycompile ../sps
make rpnphy -j4

# Compile MESH
cd ../MESH_SVS
make clean

if [ "$activate_debug" = true ]; then
   if [ "$system" = "Science" ] || [ "$system" = "GPSCC"  ]; then		
       make mpi_intel debug
   elif [ "$system" = "Other" ]; then
       make mpi_gcc debug
   fi
else
   if [ "$system" = "Science" ] || [ "$system" = "GPSCC"  ]; then		
      make mpi_intel
   elif [ "$system" = "Other" ]; then
       make mpi_gcc
   fi
fi

