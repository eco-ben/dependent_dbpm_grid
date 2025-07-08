#!/usr/bin/env python3

# Loading libraries
import os
import useful_functions as uf
import xarray as xr
import numpy as np
import os
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
            if 'monthly' in f:
                loess_smooth = True
                weighted_search = '_monthly-smoothed_'
                weighted_fn = '-smoothed'
            else:
                loess_smooth = False
                weighted_search = '_monthly_'
                weighted_fn = ''
            uf.extract_gfdl(f, mask, f_out, cross_dateline = cross_dateline,
                            loess_smooth = loess_smooth)
    
        #Load area of grid cell to be used as weights
        weights = xr.open_zarr(glob(os.path.join(gfdl_out, '*area*'))[0]).cellareao
        #Areas outside LME need to be given a value of 0
        weights = weights.fillna(0)
        
        region_int = reg.replace('-', ' ').upper()
        
        obs_fn = (glob(os.path.join(gfdl_out, f'*gfdl*obsclim*{weighted_search}*'))+
                  glob(os.path.join(gfdl_out, f'*obsclim_deptho*')))
        weighted_obs_df = uf.weighted_mean_timestep(obs_fn, weights, region_int)
        
        ctrl_fn = (glob(os.path.join(gfdl_out, f'gfdl*ctrlclim*{weighted_search}*'))+ 
                   glob(os.path.join(gfdl_out, f'*ctrlclim_deptho*')))
        weighted_ctrl_df = uf.weighted_mean_timestep(ctrl_fn, weights, region_int)
        
        spinup_fn = (glob(os.path.join(gfdl_out, f'gfdl*spinup*{weighted_search}*'))+
                     glob(os.path.join(gfdl_out, f'*ctrlclim_deptho*')))
        weighted_spinup_df = uf.weighted_mean_timestep(spinup_fn, weights, region_int)
        
        stable_fn = (glob(os.path.join(gfdl_out, f'*gfdl*stable-spin*{weighted_search}*'))+
                     glob(os.path.join(gfdl_out, f'*ctrlclim_deptho*')))
        weighted_stable_spin_df = uf.weighted_mean_timestep(stable_fn, weights, region_int)
        
        #Defining output folder
        gfdl_out = f'/g/data/vf71/fishmip_inputs/ISIMIP3a/fao_inputs/{reg}/monthly_weighted/{res}'
        os.makedirs(gfdl_out, exist_ok = True)

        #Saving data
        weighted_obs_df.to_parquet(
            os.path.join(gfdl_out, 
                         f'obsclim_dbpm_clim-inputs{weighted_fn}_{reg}_1961-2010.parquet'), 
            index = False)
        
        weighted_ctrl_df.to_parquet(
            os.path.join(gfdl_out,
                         f'ctrlclim_dbpm_clim-inputs{weighted_fn}_{reg}_1961-2010.parquet'), 
            index = False)
        
        weighted_spinup_df.to_parquet(
            os.path.join(gfdl_out,
                         f'spinup_dbpm_clim-inputs{weighted_fn}_{reg}_1841-1960.parquet'), 
            index = False)
        
        weighted_stable_spin_df.to_parquet(
            os.path.join(gfdl_out, 
                         f'stable-spin_dbpm_clim-inputs{weighted_fn}_{reg}_1741-1840.parquet'), 
            index = False)
    