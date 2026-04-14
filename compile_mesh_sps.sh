# Name of the system where MESH-SVS is compiled
# Science: internal ECCC network (use ECCC Gitlab to retrieve SVS and SVS2 code)
# GPSCC: ECCC collaboration server (use Github to retrieve SVS and SVS2 code)
# Other: other machine (use Github to retrieve SVS and SVS2 code)
system=Other # Science | GPSCC | Other

# Tag of SPS version or name of SPS branch to be extracted from reference SPS repository on Gitlab or Gitbub
# If tag_sps_user is not speficied, the most recent branch is used as a default. 
#tag_sps_user=630a20_vvi001_surface_fora22
tag_sps_user=6.3

# Compile in debug mode
activate_debug=true

###### No changes required below this line

# Select tag of branch name to be used
if [[ -n $tag_sps_user ]]; then
   tag_sps=$tag_sps_user
else
  if [ "$system" = "Science" ]; then
    tag_sps=630-a22
  elif [ "$system" = "Other" ] || [ "$system" = "GPSCC"  ]; then	
    tag_sps=6.3 
  fi
fi

# Change dir
cd ../

# Create code structure
mkdir sps sps_build

# Extract sps code from gitlab
cd sps

# Clone based on system
if [ "$system" = "Science" ]; then
    git clone --no-checkout git@gitlab.science.gc.ca:continental-surface-hydrology/sps-dev.git .
elif [ "$system" = "Other" ] || [ "$system" = "GPSCC"  ]; then	
    git clone --branch 6.3 git@github.com:VVionnet/sps_dev.git .
else
    echo "$system is an unvalid machine name. Please choose among: 'Scicence', 'GPSCC' and 'Other'"	
    exit
fi

git checkout $tag_sps

# Adjust relative path in .gitsubmodule
if [ "$system" = "Other" ] || [ "$system" = "GPSCC"  ]; then
    git config submodule.cmake_rpn.url https://github.com/ECCC-ASTD-MRD/cmake_rpn
    git config submodule."src/rpn-si/vgrid".url https://github.com/ECCC-ASTD-MRD/vgrid
    git config submodule."src/rpn-si/rpncomm".url https://github.com/ECCC-ASTD-MRD/rpncomm
    git config submodule."src/rpn-si/rmn".url https://github.com/ECCC-ASTD-MRD/librmn
    git config submodule."src/rpn-si/tdpack".url https://github.com/ECCC-ASTD-MRD/tdpack
fi	

git submodule update --init --recursive

# Load compiler
if [ "$system" = "Science" ] || [ "$system" = "GPSCC"  ]; then	
   . .eccc_setup_intel
elif [ "$system" = "Other" ]; then
   . .common_setup gnu	
fi

# Go to build dir
cd ../sps_build

# Compile rpn physics
cmake ../sps
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

