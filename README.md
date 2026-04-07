# MESH model with SVS 1.0 and 2.0 

This repository contains the version of the MESH model compatible with the latest versions of the land surface scheme Soil Vegetation and Snow (SVS) version 1.0 and 2.0.  

# Installation and compilation

Create first a repository where the necessary code will be downloaded and compiled: 

```
mkdir example_install
cd example_install
```
To get the code, use `git clone`:

```
git clone https://github.com/VVionnet/MESH_SVS.git
```

To compile the code, the user needs to edit the script `compile_mesh_sps.sh` in the `MESH_SVS` directory. 

First, the user needs to specify the type of machine used to run MESH-SVS:

- `Science`: internal ECCC network
- `GPSCC`: ECCC collaboration server
- `Other`: other machine

When specifying `GPSCC` or `Other`, the SVS code is obtained from the developement repository of the ECCC Surface Prediction System on Github (https://github.com/VVionnet/sps_dev). This repository is a fork from the main official SPS repo on Github (https://github.com/ECCC-ASTD-MRD/sps). When specifying `Science`, the SVS code is obtained from the SPS developement repository on the ECCC internal Gitlab.  

The user can also specify the specific branch or tag that they want to compile from the SPS repository (key `tag_sps_user`). If this branch is not specified, the more recent branch is used as a default.  

Finally, the user can choose to compile the MESH-SVS code with or without a debug option by setting the key `activate_debug` to True or False. 

Once `compile_mesh_sps.sh` has been edited, the user can run the script. It will create two repositories: 

- `sps`: it contains the routines of the ECCC Surface Prediction System, including the SVS code. 
- `sps_build`: used when compiling SPS.

  The compilation will generate the executable: `mpi_sa_mesh`

# Code modification and new compilation 

The code can be modified at several places: 

- `/MESH_SVS/LSS_Model/SVS/runsvs_mesh.F90` contains the interfaces routines between the MESH code and the SVS code (useful to modify the outputs)
- `/sps/src/rpnphy/src/surface` contains the SVS 1.0 and SVS 2.0 code.
- `/sps/src/rpnphy/src/surface/from_surfex* includes the part of the code in common with the SURFEX platform, including the detailed snowpack scheme Crocus.

Once the code has been modified, it needs to be recompiled using the script `recompile_mesh_svs.sh` located in the `MESH_SVS` directory. The user needs to edit this script to specify the type of machine used to run MESH-SVS (see above for `compile_mesh_svs.sh`). 

# More information

Information about MESH are provided on the [MESH wiki](https://mesh-model.atlassian.net/wiki/spaces/USER/overview?mode=global). Specific information on the use of SVS 1.0 and 2.0 in MESH are detailed [here](https://mesh-model.atlassian.net/wiki/spaces/USER/pages/6390037/Soil-Vegetation-Snow+SVS). In particular, the instructions to configure the model in point-scale mode are given [here](https://mesh-model.atlassian.net/wiki/spaces/USER/pages/6390475/How+to+configure+MESH-SVS+for+point+mode+1D+including+SVS2)

# Test case 
The directory *test_case* contains an example of a MESH-SVS experiment in point-scale mode. 
