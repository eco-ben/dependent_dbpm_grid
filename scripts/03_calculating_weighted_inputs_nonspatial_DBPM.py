#!/usr/bin/env python3

# Loading libraries
import os
import useful_functions as uf
import xarray as xr
import numpy as np
import os
from glob import glob
import shutil

# Defining base folder
base_dir = '/g/data/vf71/fishmip_inputs/ISIMIP3a/fao_inputs/'
# Getting list of FAO regions
fao_code = [f for f in os.listdir(base_dir) if 'fao' in f]
# Define resolutions
# resolutions = ['1deg', '025deg']
resolutions = ['025deg']
# Process smoothed inputs?
smoothed = True

for fao, res in zip(fao_code, resolutions):
    if not smoothed:
        gridded_folder = os.path.join(base_dir, fao, 'gridded', res, '*')
        weighted_fn = ''
    else:
        gridded_folder =os.path.join(base_dir, fao, 'gridded-smoothed',
                                      res, '*')
        weighted_fn = '-smoothed'

    #Extracting data for FAO area
    #Load area of grid cell to be used as weights
    weights = xr.open_zarr(glob(os.path.join(gridded_folder, '*area*'))[0]).cellareao
    #Areas outside LME need to be given a value of 0
    weights = weights.fillna(0)
    
    region_int = fao.replace('-', ' ').upper()
    
    obs_fn = (glob(os.path.join(gridded_folder, f'*gfdl*obsclim*'))+
              glob(os.path.join(gridded_folder, f'*obsclim_deptho*')))
    weighted_obs_df = uf.weighted_mean_timestep(obs_fn, weights, region_int)
    
    ctrl_fn = (glob(os.path.join(gridded_folder, f'gfdl*ctrlclim*'))+ 
               glob(os.path.join(gridded_folder, f'*ctrlclim_deptho*')))
    weighted_ctrl_df = uf.weighted_mean_timestep(ctrl_fn, weights, region_int)
    
    spinup_fn = (glob(os.path.join(gridded_folder, f'gfdl*spinup*'))+
                 glob(os.path.join(gridded_folder, f'*ctrlclim_deptho*')))
    weighted_spinup_df = uf.weighted_mean_timestep(spinup_fn, weights, region_int)
    
    stable_fn = (glob(os.path.join(gridded_folder, f'*gfdl*stable-spin*'))+
                 glob(os.path.join(gridded_folder, f'*ctrlclim_deptho*')))
    weighted_stable_spin_df = uf.weighted_mean_timestep(stable_fn, weights, region_int)

    #Defining output folder
    gfdl_out = os.path.join(base_dir, f'{fao}/monthly_weighted/{res}')
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
    