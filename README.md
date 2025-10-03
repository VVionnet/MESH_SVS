# MESH model with SVS 1.0 and 2.0 

This repository contains the version of the MESH model that includes the land surface scheme SVS2. 

# Installation

To get the code, use `git clone`:

```
git clone https://github.com/VVionnet/MESH_SVS.git
```
Infomation about the compilation of MESH are provided in file README.txt. Several compilers are supported (ifort, gfortran, ..).

To compile with ifort, execute the command `make ifort` in the main directory. A debug option is availalbe using the command: `make ifort debug`. 

To compile with gfortran, execute the command `make gfortran` in the main directory. A debug option is availalbe using the command: `make gfortran debug`. 

On the contained of ECCC collaboration server **GPSCC** (`inter-c-eccc-ubuntu2204.collab.science.gc.ca`), load the Intel compiler ifort using the following command: 

```
. r.load.dot  /fs/ssm/eccc/mrd/rpn/code-tools/20250826/env/ubuntu-22.04-amd64-64@inteloneapi-2023.2.0

```

and then compile the code with `make ifort`. 

# General information

Information about MESH are provided on the [MESH wiki](https://mesh-model.atlassian.net/wiki/spaces/USER/overview?mode=global). Specific information on the use of SVS 1.0 and 2.0 in MESH are detailed [here](https://mesh-model.atlassian.net/wiki/spaces/USER/pages/6390037/Soil-Vegetation-Snow+SVS). In particular, the instructions to configure the model in point-scale mode are given [here](https://mesh-model.atlassian.net/wiki/spaces/USER/pages/6390475/How+to+configure+MESH-SVS+for+point+mode+1D+including+SVS2)

# Code organization
* */Modules/rpnphy/6.1.0/src/surface* contains the codes of SVS 1.0 and 2.0 (including the code of the snowpack models Crocus and ES). 
* *LSS_Model/SVS/svs1/src* contains the interface routine between MESH and SVS (useful to modify the outputs)

# Test case 
The directory *test_case* contains an example of a MESH-SVS experiment in point-scale mode. 
