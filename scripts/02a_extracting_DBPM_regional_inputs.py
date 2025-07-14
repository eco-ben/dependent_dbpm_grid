#!/usr/bin/env python3

# Loading libraries
import os
import useful_functions as uf
import xarray as xr
import numpy as np
from glob import glob

# Loading mask with all FAO regions
mask_folder = '/g/data/vf71/shared_resources/fao_fishing_areas'
mask_1deg = xr.open_dataarray(os.path.join(mask_folder, 
                         'gfdl-mom6-cobalt2_FAO_MajorAreas_60arcmin_global_fixed.nc'))
mask_025deg = xr.open_dataarray(os.path.join(mask_folder, 
                         'gfdl-mom6-cobalt2_FAO_MajorAreas_15arcmin_global_fixed.nc'))

# Getting region codes included in mask
fao_code = np.unique(mask_025deg.values[np.isfinite(mask_025deg.values)]).astype(int)

# Choose a resolution
resolutions = ['1deg', '025deg']

for fao in fao_code:
    reg = f'fao-{fao}'
    for res in resolutions:
        if res == '1deg':
            mask = xr.where(mask_1deg == fao, 1, np.nan)
        elif res == '025deg':
            mask = xr.where(mask_025deg == fao, 1, np.nan)
        file_list = glob(f'/g/data/vf71/fishmip_inputs/ISIMIP3a/global_gridded_zarr/{res}/*')
        #Ignoring files not needed
        file_list = [f for f in file_list if '-vint_' not in f]
        gfdl_out = f'/g/data/vf71/fishmip_inputs/ISIMIP3a/fao_inputs/{reg}/gridded/{res}'
        os.makedirs(gfdl_out, exist_ok = True)
    
        #Extracting data for FAO area
        for f in file_list:
            #Create file path to save outputs
            f_out = os.path.basename(f).replace('global', reg)
            f_out = os.path.join(gfdl_out, f_out)
            #Apply function
            if fao in [18, 61, 71, 81, 88]:
                cross_dateline = True
            else:
                cross_dateline = False
            uf.extract_gfdl(f, mask, f_out, cross_dateline = cross_dateline)
    