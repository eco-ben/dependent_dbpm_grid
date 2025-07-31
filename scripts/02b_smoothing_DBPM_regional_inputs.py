#!/usr/bin/env python3

# Loading libraries
import os
import useful_functions as uf
from glob import glob
import shutil

# Defining base folder
base_dir = '/g/data/vf71/fishmip_inputs/ISIMIP3a/fao_inputs/'
# Getting list of FAO regions
# fao_code = [f for f in os.listdir(base_dir) if 'fao' in f]
fao_code = [34]
# Define resolutions
# resolutions = ['1deg', '025deg']
resolutions = ['1deg']
# Choose whether smoothing of inputs will be performed by LOESS (smoothed) or
# deseasoning data (deseasoned)
smoothing = 'deseasoned'

for fao, res in zip(fao_code, resolutions):
    file_list = glob(os.path.join(base_dir, f'fao-{fao}', 'gridded', res,
                                  'gfdl-*'))
    folder_out = os.path.join(base_dir, f'fao-{fao}', f'gridded-{smoothing}', res)
    os.makedirs(folder_out, exist_ok = True)

    #Extracting data for FAO area
    for f in file_list:
        if 'monthly' in f and 'stable-spin' not in f:
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
            
        else:
            f_out = f.replace('gridded/', f'gridded-{smoothing}/')
            try:
                shutil.copytree(f, f_out)
            except:
                pass
    