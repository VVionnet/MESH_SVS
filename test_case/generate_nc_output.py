import os,shutil,pdb,glob
import pandas as pd
import xarray as xr
import numpy as np
import sys


if len(sys.argv) != 2:
    raise ValueError('Requires one argument [svs1 or svs2]')

# Land surface scheme : select among svs1 and svs2
lss = sys.argv[-1]

# Go to exp directory
os.chdir('output')

# Get list of files containing SVS outputs
list_fic = glob.glob(lss+'*.csv')

for ific,fic in enumerate(list_fic):

	# Type of file
	fic_type = fic.split('_')[1]

	df= pd.read_csv(fic,skipinitialspace=True)

	# Create date block to generate time (see https://stackoverflow.com/questions/48587595/convert-julian-dates-to-normal-dates-in-a-dataframe)
	df['date_block'] = df['YEAR'].astype(str)+"-"+df['JDAY'].astype(str)+"-"+df['HOUR'].astype(str)
	df['time'] = pd.to_datetime(df['date_block'], format='%Y-%j-%H')
	#Create time index
	df.index = pd.DatetimeIndex(df.time)
	#Remove unnecessary variables
	df = df.drop(columns=['YEAR', 'JDAY','HOUR','MINS','date_block','time'])

        # get rid of last columns 
	df = df.iloc[:, :-1]

        # Get list of columns 
	col = list(df.columns)
	col_short = []
	for nam in col:
		ss = nam.split('_')
		if(ss[-1].isnumeric()):
			col_short.append('_'.join(ss[:-1]))
		else:
			col_short.append(nam)

	var_ref = set(col_short)

	if(len(var_ref) == len(col_short)): # File containing multilayer variables
		da = df.to_xarray()
	else:
		da_int = df.to_xarray()

		for ivar, var in enumerate(var_ref):

			mask = [var == tt for tt in col_short]
			lvar = np.array(col)[mask].tolist()
			da_sel = da_int[lvar]  # Extract relevant data

			nlayer = len(da_sel)

			if(nlayer>1):
				lis_layer = np.arange(0,nlayer)
				yy=np.zeros((len(da_sel.time),nlayer))

				if(var in ['TGROUND','TVEG']):
					name_layer = 'fr_layer'
				else:
					name_layer = fic_type+'_layer'
				da_var = xr.DataArray(yy,coords={'time':da_sel.time,name_layer:lis_layer},dims=('time',name_layer))

				for ilayer,var_sel in enumerate(lvar):
					da_var[:,ilayer] = da_sel[var_sel].values
				da_var.name = var
				da_var.to_dataset()
			else: 
				da_var = da_sel 


			if(ivar ==0):
				da=da_var
			else:
				da=xr.merge([da,da_var])

	if(ific ==0):
		ref=da
	else:
		ref=xr.merge([ref,da])

# List of cumulated variable that need to be reprocessed to handle the fact that they are reset to zero 
# every day at 12 UTC by MESH-SVS (to mimic GEM-Hydro daily integration cycle)
var_cum=['RSNOW_AC']

for var in var_cum:
	if var in ref.data_vars:
		#Get the hourly increase from the cumulated values
		ext = ref[var].diff(dim='time', label='upper')

		# Extract data at 13 UTC 	
		mm=ext.time.dt.hour==13
		ref_var = ref[var][1:]
		
		# Adjust the houlry increase at 13 UTC
		ext.loc[dict(time=ext.time[mm])] = ref_var[mm].values

		# Compute cumulated values
		ext = ext.cumsum()

		# Update value in NetCDf file so that the correct cumulated values is used. 
		ref[var][1:] = ext[:]


# functions and encoding info
DEFAULT_ENCODING = {
    'zlib': True,
    'shuffle': True,
    'complevel': 4,
    'fletcher32': False,
    'contiguous': False,
}

def generate_encodings(data):
    encoding = {}
    for var in data.data_vars:
        encoding[var] = DEFAULT_ENCODING.copy()
    return encoding

# Write netcdf   
encoding = generate_encodings(ref) 
netcdf_file_out = 'out_'+lss+'.nc'
ref.to_netcdf(netcdf_file_out)
