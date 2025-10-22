#!/usr/bin/env python3

# Loading libraries
import os
import useful_functions as uf
import xarray as xr
import numpy as np
from glob import glob
from distributed import Client
from multiprocessing import Process, freeze_support

#Start cluster
if __name__ == '__main__':
    freeze_support()

    #Start a cluster
    client = Client(threads_per_worker = 1, memory_limit = 0)

    # Defining folder where masks with all FAO regions are stored
    mask_folder = '/g/data/vf71/shared_resources/fao_fishing_areas'
    
    # Define location of DBPM inputs
    base_dir = '/g/data/vf71/fishmip_inputs/ISIMIP3a/'
    
    # Define variables for which data will be extracted
    vars_int = ['simask', 'er', 'intercept', 'slope', 'tob', 'tos', 'deptho', 'areacello']
    
    # Define resolutions
    resolutions = ['1deg', '025deg']
    
    # Choose whether smoothing of inputs will be performed by LOESS (smoothed) or
    # deseasoning data (deseasoned). Select None for no smoothing.
    smoothing = 'deseasoned'
    
    for res in resolutions:
        if res == '1deg':
            mask_all = xr.open_dataarray(
                os.path.join(mask_folder,
                             'gfdl-mom6-cobalt2_FAO_MajorAreas_60arcmin_global_fixed.nc'))
        elif res == '025deg':
            mask_all = xr.open_dataarray(
                os.path.join(mask_folder, 
                             'gfdl-mom6-cobalt2_FAO_MajorAreas_15arcmin_global_fixed.nc'))
        
        # Getting region codes included in mask
        fao_code = np.unique(mask_all.values[np.isfinite(mask_all.values)]).astype(int)
        
        #Define GFDL folder
        if smoothing is None:
            file_list = glob(os.path.join(base_dir, 'global_gridded_zarr', res, '*'))
            out_name = 'gridded'
        else:
            file_list = glob(
                os.path.join(base_dir, f'global_gridded-{smoothing}_zarr', res, '*'))
            si_list = glob(os.path.join(base_dir, 'global_gridded_zarr', res, '*_simask_*'))
            file_list = file_list+si_list
            out_name = f'gridded-{smoothing}'
            

        #List all files to be extracted
        for fao in fao_code:
            gfdl_out = os.path.join(base_dir, 'fao_inputs', f'fao-{fao}', 
                                    out_name, res)
            os.makedirs(gfdl_out, exist_ok = True)
            mask = xr.where(mask_all == fao, 1, np.nan)
            for dv in vars_int:
                #Ignoring files not needed
                file_dv = [f for f in file_list if f'_{dv}_' in f]
                
                #Extracting data for FAO area
                for f in file_dv:
                    #Create file path to save outputs
                    f_out = os.path.basename(f).replace('global', f'fao-{fao}')
                    f_out = os.path.join(gfdl_out, f_out)
                    #Apply function
                    if fao in [18, 61, 71, 81, 88]:
                        cross_dateline = True
                    else:
                        cross_dateline = False
                    uf.extract_gfdl(f, mask, f_out, cross_dateline = cross_dateline)
    