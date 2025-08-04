#!/usr/bin/env python3

# Loading libraries
import os
from glob import glob
import useful_functions as uf
import dask
from distributed import Client
from multiprocessing import Process, freeze_support
import shutil

########################### THIS STEP IS OPTIONAL ###########################
# Only run this script if inputs need to be smoothed out before running 
# DBPM (gridded or non-spatial)

#Start cluster
if __name__ == '__main__':
    freeze_support()

    #Start a cluster
    client = Client(threads_per_worker = 1, memory_limit = 0)

    #Base folder where GFDL outputs are stored 
    base_dir = '/g/data/vf71/fishmip_inputs/ISIMIP3a/global_gridded_zarr'

    #Define experiments and resolution
    exp_name = ['obsclim', 'ctrlclim']
    resolutions = ['025deg']#, '025deg']

    #Define variables of interest
    dynamic_vars = ['er', 'intercept', 'slope']#, 'tob', 'tos']

    # Choose whether smoothing of inputs will be performed by LOESS (smoothed) or
    # deseasoning data (deseasoned)
    smoothing = 'deseasoned'
    #Create folder for smoothed data if needed
    smooth_folder = base_dir.replace('_gridded_', f'_gridded-{smoothing}_')
    # Ensure output folder exists
    os.makedirs(smooth_folder, exist_ok = True)

    for res in resolutions:
        #Input folder
        f_in = os.path.join(base_dir, res)
        #Output folder
        f_out = os.path.join(smooth_folder, res)
        os.makedirs(f_out, exist_ok = True)
        #Identify area file
        [area_file] = glob(os.path.join(f_in, f'*_areacello_*fixed*'))
        area_out = os.path.join(f_out, os.path.basename(area_file))
        try:
            shutil.copytree(area_file, area_out)
        except:
            pass
        #Identify depth files
        for exp in exp_name:
            [depth_file] = glob(os.path.join(f_in, f'*{exp}_deptho_*fixed*'))
            depth_out = os.path.join(f_out, os.path.basename(depth_file))
            try:
                shutil.copytree(depth_file, depth_out)
            except:
                pass

            #Find files to be smoothed
            for dv in dynamic_vars:
                [file_in] = glob(os.path.join(f_in, f'*{exp}_{dv}_*'))
                file_out = os.path.join(
                    f_out, os.path.basename(file_in).replace('_monthly_',
                                                             '_monthly-smoothed_'))
                if smoothing == 'deseasoned':
                    
                    # Seasonal decomposition 
                    uf.seasonal_decomposition(file_in, file_out, period = 12, 
                                              component = 'seasonal')
                elif smoothing == 'smoothed':
                    # LOESS smoothing
                    uf.smoothing_loess(file_in, file_out, frac = 12, it = 0)
