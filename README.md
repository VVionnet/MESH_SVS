# MESH model with SVS 2.0 

This repository contains the version of the MESH model that includes the land surface scheme SVS2. 

# Installation

Infomation about the compilation of MESH are provided in file README.txt. Several compilers are supported (ifort, gfortran, ..).

To compile with ifort, execute the command `make ifort` in the main directory. A debug option is availalbe using the command: `make ifort debug`. 

To compile with gfort, execute the command `make gfortran` in the main directory. A debug option is availalbe using the command: `make gfortran debug`. 

# General information

Information about MESH are provided on the [MESH wiki](https://wiki.usask.ca/pages/viewpage.action?pageId=220332269). Specific information on the use of SVS 1.0 and 2.0 in MESH are detailed [here](https://wiki.usask.ca/pages/viewpage.action?pageId=1303674916). In particular, the instructions to configure the model in point mode are given [here](https://wiki.usask.ca/pages/viewpage.action?pageId=1716094475). 

# Code organization
* */Modules/rpnphy/6.1.0/src/surface* contains the codes of SVS 1.0 and 2.0 (including the code of the snowpack models Crocus and ES). 
* *LSS_Model/SVS/svs1/src* contains the interface routine between MESH and SVS (useful to modify the outputs)

# Test case 
The directory *test_case* contains an example of a MESH-SVS experiment in point scale mode. 
