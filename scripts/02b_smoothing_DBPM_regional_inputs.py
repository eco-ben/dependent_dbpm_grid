#!/usr/bin/env python3

# Loading libraries
import os
import pandas as pd
import useful_functions as uf
from glob import glob
import shutil

# Defining base folder
base_dir = '/g/data/vf71/fishmip_inputs/ISIMIP3a/fao_inputs/'
# Getting list of FAO regions
# fao_code = [f for f in os.listdir(base_dir) if 'fao' in f]
fao_code = ['fao-21']
# Define resolutions
# resolutions = ['1deg', '025deg']
resolutions = ['1deg']
# Choose whether smoothing of inputs will be performed by LOESS (smoothed) or
# deseasoning data (deseasoned)
smoothing = 'deseasoned'

#Defining stable spin and spinup periods
stable_spin = pd.date_range('1741-01', end = '1840-12', freq = 'MS')
spinup_period = pd.date_range('1841-01', end = '1960-12', freq = 'MS')

for fao, res in zip(fao_code, resolutions):
    #Gridded folder
    gridded_folder = os.path.join(base_dir, fao, 'gridded', res)
    # Define output folder
    folder_out = os.path.join(base_dir, fao, f'gridded-{smoothing}', res)
    # Ensure output folder exists
    os.makedirs(folder_out, exist_ok = True)
    
    # Get list of files to be copied across without any smoothing
    fixed_vars = glob(os.path.join(gridded_folder, 'gfdl-*fixed*'))
    for fv in fixed_vars:
        f_out = os.path.join(folder_out, os.path.basename(fv))
        try:
            shutil.copytree(fv, f_out)
        except:
            pass
    
    # Get a list of variables to be smoothed
    file_list = glob(os.path.join(gridded_folder, 'gfdl-*clim_*monthly*'))
    # Apply smoothing to all files
    for f in file_list:
        #Create file path to save outputs
        f_out = os.path.basename(f).replace('_monthly_', '_monthly-smoothed_')
        path_out = os.path.join(folder_out, f_out)
        if smoothing == 'deseasoned':
            # Seasonal decomposition 
            uf.seasonal_decomposition(f, path_out, period = 12, 
                                      component = 'seasonal')
        elif smoothing == 'smoothed':
            # LOESS smoothing
            uf.smoothing_loess(f, path_out, frac = 12, it = 0)
            
        #Create spinup and stable-spin periods
        if '_ctrlclim_' in f_out:
            # Spinup period
            fn_spin = (f_out.replace('ctrlclim', 'spinup').
                       replace('1961', '1841').replace('2010', '1960'))
            fn_spin = os.path.join(folder_out, fn_spin)
            uf.gridded_spinup(path_out, '1961-01', '1980-12', spinup_period,
                              file_out = fn_spin)
            
            # Stable spinup period
            fn_stable_spin = (fn_spin.replace('spinup', 'stable-spin').
                              replace('1841', '1741').replace('1960', '1840'))
            uf.gridded_spinup(fn_spin, '1841-01', '1841-12', stable_spin,
                              mean_spinup = True, file_out = fn_stable_spin)
    