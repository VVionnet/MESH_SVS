# Test case for MESH-SVS in point mode

This folder contains a setup for a test of MESH-SVS in  point mode at the Col de Porte experimental catchment in the French Alps. Detailed instructions to configure the model in point mode are given [here](https://wiki.usask.ca/pages/viewpage.action?pageId=1716094475). 

# Configuration files

**MESH_input_run_options.ini**  MESH configuration file including information about the meteorologicla forcing

**MESH_input_soil_levels.txt** MESH configuration file about the layering of the soil 

**MESH_parameters_cdp_svs2.txt** Configuration file for SVS2 at Col de Porte 

**MESH_parameters_cdp_svs1.txt** Configuration file for SVS1 at Col de Porte used in [Leonadini et al. (2021)](https://journals.ametsoc.org/view/journals/hydr/aop/JHM-D-20-0249.1/JHM-D-20-0249.1.xml) 

# Input files 

**basin_forcing.met** Meteorological driving data used in [Leonadini et al. (2021)](https://journals.ametsoc.org/view/journals/hydr/aop/JHM-D-20-0249.1/JHM-D-20-0249.1.xml) including a 5-y spin up (first year of the CDP dataset is repeated 5 times) 

**obs_insitu_cdp_1994_2014.nc** Evaluation data from the [ESM-Snow MIP dataset](https://doi.pangaea.de/10.1594/PANGAEA.897575)

# Routines 

