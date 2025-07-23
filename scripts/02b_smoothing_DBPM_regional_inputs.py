#!/usr/bin/env python3

# Loading libraries
import os
import useful_functions as uf
from glob import glob
import shutil

# Gaussian smoothing
# import xarray as xr
# from scipy.ndimage import gaussian_filter1d
# def apply_gaus_fil(da, sigma):
#     return xr.apply_ufunc(
#         gaussian_filter1d,
#         da,
#         input_core_dims = [['time']],
#         output_core_dims = [['time']],
#         kwargs = {'sigma': sigma},
#         dask = 'parallelized', # For larger datasets, enables Dask
#         output_dtypes = [da.dtype]
#     )

# Defining base folder
base_dir = '/g/data/vf71/fishmip_inputs/ISIMIP3a/fao_inputs/'
# Getting list of FAO regions
# fao_code = [f for f in os.listdir(base_dir) if 'fao' in f]
# fao_code = [48, 51, 57, 58, 61, 67, 71, 77, 81, 87, 88]
fao_code = [31]
# Define resolutions
# resolutions = ['1deg', '025deg']
resolutions = ['1deg']

for fao, res in zip(fao_code, resolutions):
    file_list = glob(os.path.join(base_dir, f'fao-{fao}', 'gridded', res,
                                  'gfdl-*'))
    folder_out = os.path.join(base_dir, f'fao-{fao}', 'gridded-smoothed', res)
    os.makedirs(folder_out, exist_ok = True)

    #Extracting data for FAO area
    for f in file_list:
        if 'monthly' in f:
            #Create file path to save outputs
            f_out = os.path.basename(f).replace('_monthly_', 
                                                '_monthly-smoothed_')
            path_out = os.path.join(folder_out, f_out)
            uf.smoothing_loess(f, path_out, frac = 12, it = 0)
            # # Loading data
            # da = xr.open_zarr(f)
            # #Getting name of variable contained in dataset
            # [var] = list(da.keys())
            # da = da[var]
            # da_gf = apply_gaus_fil(da, 1)
            # da_gf.to_zarr(path_out, consolidated = True, mode = 'w')
        else:            
            f_out = f.replace('gridded/', 'gridded-smoothed/')
            try:
                shutil.copytree(f, f_out)
            except:
                pass
    