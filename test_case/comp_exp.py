from argparse import ArgumentParser
import os,shutil,pdb,glob
import pandas as pd
import xarray as xr
import numpy as np
import pdb
import matplotlib.pyplot as plt


fic1 = '/home/vvi001/Codes/MESH/mesh_svs2_6-1/test_case/output/out_svs2.nc'
ds1 = xr.open_dataset(fic1)

fic2 = '/home/vvi001/Codes/MESH/Dev_SVS2/MESH_SVS1_Freezing/MESH_SVS/test_case/output/out_svs2.nc'
ds2 = xr.open_dataset(fic2)


list_var = ['SNODEN']
#for var in list(ds1.data_vars):
for var in list_var:

	val1 = ds1[var]
	val2 = ds2[var]

        
	diff = val1-val2
        
	#diff.plot()
	#plt.show()

	#diff_var =  np.isclose(val1.values,val2.values,equal_nan=True)




list_var_2d = ['WSOIL']
nlayer=0


for var in list_var_2d:

	val1 = ds1[var]
	val2 = ds2[var]
	#pdb.set_trace()

	diff = (val2-val1)/val1 *100. 
        
	diff.plot()
        
	plt.show()

	#diff_var =  np.isclose(val1.values,val2.values,equal_nan=True)


	#pdb.set_trace()


