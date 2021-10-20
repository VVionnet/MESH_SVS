from argparse import ArgumentParser
import os,shutil,pdb,glob
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt


# Beginning and end date
date_begp = '1994-9-1'
date_endp = '2014-9-1'

# Min and max for temperature on y-axis
tmin =265
tmax=300


# Read Observations using xarray
ds = xr.open_dataset('obs_insitu_cdp_1994_2014.nc')
print('read obs')

#fig=plt.figure(figsize=(5,10))

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True)

print('1')
#ax1=fig.add_subplot(4,1,1)
#ax2=fig.add_subplot(4,1,2)
#ax3=fig.add_subplot(4,1,3)
#ax4=fig.add_subplot(4,1,4)

ax1.set_ylabel('SD (m)', fontsize=15)

ax1.plot(ds.time.loc[dict(time=slice(date_begp,date_endp))].values,ds.snd_auto.loc[dict(time=slice(date_begp,date_endp))].values,'ko',label='Obs. (Auto)',linewidth=2.0)

ax1.set_ylim([0,2.2])
ax1.grid(True)

ax2.set_ylabel('T 10 cm (K)', fontsize=15)

ax2.plot(ds.time.loc[dict(time=slice(date_begp,date_endp))].values,ds.tsl.loc[dict(time=slice(date_begp,date_endp))].values[:,0]+273.15,'k-',label='Obs.',linewidth='2.0')
ax2.axhline(273.15,color='darkgrey',ls='--')
ax2.set_ylim([tmin,tmax])
ax2.grid(True)

ax3.set_ylabel('T 20 cm (K)', fontsize=15)
ax3.plot(ds.time.loc[dict(time=slice(date_begp,date_endp))].values,ds.tsl.loc[dict(time=slice(date_begp,date_endp))].values[:,1]+273.15,'k-',label='Obs.',linewidth='2.0')
ax3.axhline(273.15,color='darkgrey',ls='--')
ax3.set_ylim([tmin,tmax])
ax3.grid(True)

ax4.set_ylabel('T 50 cm (K)', fontsize=15)
ax4.plot(ds.time.loc[dict(time=slice(date_begp,date_endp))].values,ds.tsl.loc[dict(time=slice(date_begp,date_endp))].values[:,2]+273.15,'k-',label='Obs.',linewidth='2.0')
ax4.axhline(273.15,color='darkgrey',ls='--')
ax4.set_ylim([tmin,tmax])
ax4.grid(True)

#Extract model data
mod = xr.open_dataset('output/out_svs.nc')

ax1.plot(mod.time.loc[dict(time=slice(date_begp,date_endp))].values,mod['SNODP'].loc[dict(time=slice(date_begp,date_endp))].values,label='Obs. (Auto)',linewidth=2.0,ls='solid',color='red',alpha=0.8)

t10cm = (2*mod['TPSOIL'].loc[dict(time=slice(date_begp,date_endp))].values[:,1]+mod['TPSOIL'].loc[dict(time=slice(date_begp,date_endp))].values[:,2])/3.
ax2.plot(mod.time.loc[dict(time=slice(date_begp,date_endp))].values,t10cm,label='Obs. (Auto)',linewidth=2.0,ls='solid',color='red',alpha=0.8)

t20cm = (2*mod['TPSOIL'].loc[dict(time=slice(date_begp,date_endp))].values[:,2]+mod['TPSOIL'].loc[dict(time=slice(date_begp,date_endp))].values[:,3])/3.
ax3.plot(mod.time.loc[dict(time=slice(date_begp,date_endp))].values,t20cm,label='Obs. (Auto)',linewidth=2.0,ls='solid',color='red',alpha=0.8)

t50cm = (mod['TPSOIL'].loc[dict(time=slice(date_begp,date_endp))].values[:,3]+2*mod['TPSOIL'].loc[dict(time=slice(date_begp,date_endp))].values[:,4])/3.
ax4.plot(mod.time.loc[dict(time=slice(date_begp,date_endp))].values,t50cm,label='Obs. (Auto)',linewidth=2.0,ls='solid',color='red',alpha=0.8)


namePlot='Eval_CDP_'+date_begp+'_'+date_endp+'.png'
plt.savefig(namePlot,format="png",bbox_inches = "tight")

fig, (ax1, ax2, ax3,) = plt.subplots(3, 1, sharex=True)

print('1')


ax1.set_ylabel('SD (m)', fontsize=15)
ax1.plot(ds.time.loc[dict(time=slice(date_begp,date_endp))].values,ds.snd_auto.loc[dict(time=slice(date_begp,date_endp))].values,'ko',label='Obs. (Auto)',linewidth=2.0)

ax1.plot(ds.time.loc[dict(time=slice(date_begp,date_endp))].values,ds.snd_man.loc[dict(time=slice(date_begp,date_endp))].values,'ks',label='Obs. (Auto)',linewidth=2.0)

ax1.set_ylim([0,2.3])
ax1.grid(True)

ax2.set_ylabel('SWE (kg/m2)', fontsize=15)
ax2.plot(ds.time.loc[dict(time=slice(date_begp,date_endp))].values,ds.snw_auto.loc[dict(time=slice(date_begp,date_endp))].values,'ko',label='Obs. (Auto)',linewidth=2.0)

ax2.plot(ds.time.loc[dict(time=slice(date_begp,date_endp))].values,ds.snw_man.loc[dict(time=slice(date_begp,date_endp))].values,'-k',label='Obs. (Auto)',linewidth=2.0)


ax2.set_ylim([0,750])
ax2.grid(True)

ax3.set_ylabel('T_surf (K)', fontsize=15)
ax3.plot(ds.time.loc[dict(time=slice(date_begp,date_endp))].values,ds.ts.loc[dict(time=slice(date_begp,date_endp))].values+273.15,'k-',label='Obs.',linewidth='2.0')
ax3.axhline(273.15,color='darkgrey',ls='--')
ax3.set_ylim([245,275])
ax3.grid(True)

mod = xr.open_dataset('output/out_svs.nc')

ax1.plot(mod.time.loc[dict(time=slice(date_begp,date_endp))].values,mod['SNODP'].loc[dict(time=slice(date_begp,date_endp))].values,label='Obs. (Auto)',linewidth=2.0,ls='solid',color='red',alpha=0.8)

ax2.plot(mod.time.loc[dict(time=slice(date_begp,date_endp))].values,mod['SNOMA'].loc[dict(time=slice(date_begp,date_endp))].values,label='Obs. (Auto)',linewidth=2.0,ls='solid',color='red',alpha=0.8)

ax3.plot(mod.time.loc[dict(time=slice(date_begp,date_endp))].values,mod['TSNO_SURF'].loc[dict(time=slice(date_begp,date_endp))].values,label='Obs. (Auto)',linewidth=2.0,ls='solid',color='red',alpha=0.8)


namePlot='Eval_CDP_snow_'+date_begp+'_'+date_endp+'.png'
plt.savefig(namePlot,format="png",bbox_inches = "tight")
plt.show()


pdb.set_trace()
