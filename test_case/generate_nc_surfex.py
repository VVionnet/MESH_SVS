import os,shutil,pdb,glob
import pandas as pd
import xarray as xr
import numpy as np


# Open Netcdf file from MESH-SVS
mod = xr.open_dataset('output/out_svs2.nc')

# Define fill value
fill_val = -999.

# Select time when profiles are defined
mask_undef = np.isnan(mod.SNOMA_ML.values[:,0]) # Identify date when profiles are not defined
time_def = mod.time[~mask_undef]   # Extract time when profiles are defined
mod = mod.loc[dict(time=time_def)]

# Add layer thickness and SSA
mod['SNOWDZ'] = mod['SNOMA_ML']/mod['SNODEN_ML'] # in m 
mod['SNOWSSA'] =  6./( mod['SNODOPT_ML'] * 917) # in m2 kg-1
mask=mod['SNOWSSA'].values==np.inf
mod['SNOWSSA'].values[mask] = np.nan

# Rename variables to match the SURFEX terminology
mod = mod.rename({'SNODEN_ML':'SNOWRO','SNOAGE_ML':'SNOWAGE','SNOHIST_ML':'SNOWHIST','WSNOW_ML':'SNOWLIQ','SNOSPH_ML':'SNOWSPHER','TSNOW_ML':'SNOWTEMP','SNOMA_ML':'WSN_VEG','SNOTYPE_ML':'SNOWTYPE'})

# Convert from 'SNOWLIQ' from m to kg m-3
mod['SNOWLIQ'] = 1000.*mod['SNOWLIQ']/mod['SNOWDZ'] 

# Sent missing values to 
mask_novalue=mod['WSN_VEG'].values[:,:]==0

for var in  ['SNOWSSA','SNOWRO','SNOWAGE','SNOWTEMP','SNOWLIQ','SNOWSPHER','SNOWTYPE']:
    print(var)
    mod[var].values[mask_novalue] = fill_val


# Extract variables for the generation of the SURFEX-like netcdf file
list_var = ['SNOWDZ','SNOWSSA','SNOWRO','SNOWAGE','SNOWTEMP','SNOMA','SNODP','SNODEN','SNOALB','SNOWLIQ','SNOWSPHER','WSN_VEG','SNOWTYPE']
sub = mod[list_var].copy()

# Add the dimension 'Number_of_points' present in the SURFEX netcdf
sub = sub.expand_dims(dim='Number_of_points',axis=-1)

# Remove the snow layer variable
sub = sub.drop_vars('snow_layer')

# Add the variables metadata
dic_variables = {'SNOWRO':{'long_name':'Snow_density','units':'Kg/m3'},
                 'SNOWAGE':{'long_name':'Snow_age_param','units':'days'},
                 'SNOWHIST':{'long_name':'Snow_historical_param','units':'-'},
                 'SNOWLIQ':{'long_name':'Snow_liquid_water_layer','units':'m'},
                 'SNOWTEMP':{'long_name':'Snow_Temp_layer','units':'K'},
                 'SNOWSPHER':{'long_name':'Snow_layer_sphericity','units':'-'},
                 'SNOWSSA':{'long_name':'Snow_layer_specific_surface_area','units':'m2 kg-1'},
                 'SNOWDZ':{'long_name':'Snow_layer_thickness','units':'m'},
                 'SNODP':{'long_name':'Total_snow_thickness','units':'m'},
                 'SNOMA':{'long_name':'Total_snow_water_equivalent','units':'kg m-2'},
                 'WSN_VEG':{'long_name':'Snow_layer_water_equivalent','units':'kg m-2'},
                 'SNODEN':{'long_name':'Bulk_snowpack_density','units':'kg m-3'},
                 'SNOALB':{'long_name':'Snow_surface_albedo','units':'-'},
                 'SNOWTYPE':{'long_name':'Snow_grain_type','units':'-'},

               }

for var in list_var:
	if var in dic_variables.keys():
		sub[var].attrs['long_name'] = dic_variables[var]['long_name']
		sub[var].attrs['units'] = dic_variables[var]['units']

# functions and encoding info
DEFAULT_ENCODING = {
    'zlib': True,
    'shuffle': True,
    'complevel': 9,
    'fletcher32': False,
    'contiguous': False,
    '_FillValue': fill_val,
}

def generate_encodings(data):
    encoding = {}
    for var in data.data_vars:
        encoding[var] = DEFAULT_ENCODING.copy()
    return encoding

# Write netcdf   
encoding = generate_encodings(sub) 
netcdf_file_out = 'output/pro_svs.nc'
sub.to_netcdf(netcdf_file_out,unlimited_dims='time',encoding=encoding)
