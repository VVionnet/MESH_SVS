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

On the ECCC collaboration server **GPSCC**, load the Intel compiler ifort using the following command: 

```
. r.load.dot /fs/ssm/eccc/mrd/rpn/code-tools/ENV/cdt-1.5.7-inteloneapi-2022.1.2
```
and then compile the code with `make ifort`. 

# General information

Information about MESH are provided on the [MESH wiki](https://wiki.usask.ca/pages/viewpage.action?pageId=220332269). Specific information on the use of SVS 1.0 and 2.0 in MESH are detailed [here](https://wiki.usask.ca/pages/viewpage.action?pageId=1303674916). In particular, the instructions to configure the model in point mode are given [here](https://wiki.usask.ca/display/MESH/How+to+configure+MESH-SVS+for+point+mode+%281D%29%2C+including+SVS2). 

# Code organization
* */Modules/rpnphy/6.1.0/src/surface* contains the codes of SVS 1.0 and 2.0 (including the code of the snowpack models Crocus and ES). 
* *LSS_Model/SVS/svs1/src* contains the interface routine between MESH and SVS (useful to modify the outputs)

# Test case 
The directory *test_case* contains an example of a MESH-SVS experiment in point scale mode. 
