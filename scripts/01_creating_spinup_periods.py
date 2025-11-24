#!/usr/bin/env python3

# Loading libraries
import os
from glob import glob
import numpy as np
import xarray as xr
import useful_functions as uf
import pandas as pd
import dask
from distributed import Client
from multiprocessing import Process, freeze_support

#Start cluster
if __name__ == '__main__':
    freeze_support()

    #Start a cluster
    client = Client(threads_per_worker = 1, memory_limit = 0)

    #Base folder where GFDL outputs are stored 
    base_dir = '/g/data/vf71/fishmip_inputs/ISIMIP3a'

    #Define resolutions
    resolutions = ['1deg', '025deg']

    #Define variables for which a spinup period will be created
    dynamic_vars = ['er', 'intercept', 'slope', 'tob', 'tos']

    #Defining stable spin and spinup periods
    stable_spin = pd.date_range('1741-01', end = '1840-12', freq = 'MS')
    spinup_period = pd.date_range('1841-01', end = '1960-12', freq = 'MS')

    # Choose whether smoothing of inputs will be performed by LOESS (smoothed) or
    # deseasoning data (deseasoned). Select None for no smoothing.
    smoothing = 'deseasoned'
    #Loop through experiments and resolutions
    for res in resolutions:
        #Define GFDL folder where sea ice masks are stored
        gfdl_folder = os.path.join(base_dir, 'global_gridded_zarr', res)
        #Create files for spinup period
        [si_files] = glob(os.path.join(gfdl_folder, f'gfdl-*ctrlclim_simask_*'))
        f_out = (si_files.replace('ctrlclim', 'spinup').replace('1961', '1841').
            replace('2010', '1960'))
        uf.gridded_spinup(si_files, '1961-01', '1980-12', spinup_period, file_out = f_out)
        #Create files for the stable spinup period from sea ice concentration
        #Load area of grid cells
        area = xr.open_zarr(glob(os.path.join(gfdl_folder, '*areacello*'))[0])['cellareao']
        #Load sea ice concentration
        sic = (xr.open_zarr(glob(os.path.join(
            gfdl_folder, f'gfdl-*ctrlclim_siconc_*'))[0])['siconc']).where(np.isfinite(area))
        #Getting SIC spinup
        da = uf.gridded_spinup(sic, '1961-01', '1961-12', stable_spin, mean_spinup = True)
        da_mask = xr.where(da >= 15, True, False)
        #Split into northern and southern hemispheres
        #Sea ice kept from 42N towards the north pole as the Sea of Okhotsk (45N) is the
        #lowest latitude area where sea ice forms each winter according to NASA's 
        #Earth Observatory
        da_mask_north = (xr.where(da_mask.lat > 42, da_mask, False).
            isel(lat = slice(None, None, -1)).cumsum('lat'))
        # da_mask_north.name = 'simask'
        #Sea ice kept from 52S towards the south pole as 55S is the lowest latitude area
        # where sea ice forms each winter according to NASA's Earth Observatory
        da_mask_south = xr.where(da_mask.lat <= -52, da_mask, False).cumsum('lat')
        # da_mask_south.name = 'simask'
        da_mask = (da_mask_north+da_mask_south)
        da_mask = xr.where(da_mask > 0, np.nan, 1).where(np.isfinite(area))
        #Rechunk data
        da_mask = da_mask.chunk({'time': '500MB'})
        #Update data array variable name
        da_mask.name = 'simask'
        #Create file name before saving
        fout_stable = (si_files.replace('ctrlclim', 'stable-spin').
            replace('1961', '1741').replace('2010', '1840'))
        da_mask.to_zarr(fout_stable, consolidated = True, mode = 'w')

        #Change GFDL folder if needed
        if smoothing is not None:
            gfdl_folder = os.path.join(base_dir, f'global_gridded-{smoothing}_zarr', res)
        
        #The spinup period goes from 1841 and 1960. It is created by repeating inputs
        #from "ctrlclim" experiment between 1961 and 1980
        #The stable spinup period goes from 1741 to 1840. It is created by repeating
        #the mean for the year 1841 in the spinup period
        for dv in dynamic_vars:
            [f_in] = glob(os.path.join(gfdl_folder, f'gfdl-*ctrlclim_{dv}_*monthly*'))
            f_out = (f_in.replace('ctrlclim', 'spinup').replace('1961', '1841').
                replace('2010', '1960'))
            uf.gridded_spinup(f_in, '1961-01', '1980-12', spinup_period, 
                              file_out = f_out)
            #The stable spinup period 
            fout_stable = (f_in.replace('ctrlclim', 'stable-spin').
                replace('1961', '1741').replace('2010', '1840'))
            uf.gridded_spinup(f_out, '1841-01', '1841-12', stable_spin,
                              mean_spinup = True, file_out = fout_stable)
