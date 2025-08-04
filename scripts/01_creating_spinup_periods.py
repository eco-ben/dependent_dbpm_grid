#!/usr/bin/env python3

# Loading libraries
import os
from glob import glob
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

    #Define experiments and resolution
    exp_name = ['obsclim', 'ctrlclim']
    resolutions = ['1deg']#, '025deg']

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
        #Define GFDL folder
        if smoothing is None:
            gfdl_folder = os.path.join(base_dir, 'global_gridded_zarr', res)
        else:
            gfdl_folder = os.path.join(base_dir, f'global_gridded-{smoothing}_zarr', res)
        
        #The spinup period goes from 1841 and 1960. It is created by repeating inputs
        #from "ctrlclim" experiment between 1961 and 1980
        #The stable spinup period goes from 1741 to 1840. It is created by repeating
        #the mean for the year 1841 in the spinup period
        for dv in dynamic_vars:
            [f_in] = glob(os.path.join(gfdl_folder, f'gfdl-*ctrlclim_{dv}_*monthly*'))
            f_out = (f_in.replace('ctrlclim', 'spinup').
                replace('1961', '1841').replace('2010', '1960'))
            uf.gridded_spinup(f_in, '1961-01', '1980-12', spinup_period,
                              file_out = f_out)
            #The stable spinup period 
            fout_stable = (f_in.replace('ctrlclim', 'stable-spin').
                replace('1961', '1741').replace('2010', '1840'))
            uf.gridded_spinup(f_out, '1841-01', '1841-12', stable_spin,
                              mean_spinup = True, file_out = fout_stable)
