from argparse import ArgumentParser
import os,shutil,pdb,glob
import pandas as pd
import xarray as xr
import numpy as np


# Go to exp directory
os.chdir('output')

# Get list of files containing SVS outputs
list_fic = glob.glob('svs2*.csv')

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
			col_short.append(ss[:-1][0])
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
			lis_layer = np.arange(0,nlayer)
			yy=np.zeros((len(da_sel.time),nlayer))

			da_var = xr.DataArray(yy,coords={'time':da_sel.time,fic_type+'_layer':lis_layer},dims=('time',fic_type+'_layer'))

			for ilayer,var_sel in enumerate(lvar):
				da_var[:,ilayer] = da_sel[var_sel].values
			da_var.name = var
			da_var.to_dataset()

			if(ivar ==0):
				da=da_var
			else:
				da=xr.merge([da,da_var])

	if(ific ==0):
		ref=da
	else:
		ref=xr.merge([ref,da])


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
netcdf_file_out = 'out_svs.nc'
ref.to_netcdf(netcdf_file_out)
