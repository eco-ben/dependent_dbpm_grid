#!/usr/bin/env python3

# Loading libraries
import os
import useful_functions as uf
import xarray as xr
from glob import glob

# Defining base folder
base_dir = '/g/data/vf71/fishmip_inputs/ISIMIP3a/fao_lme_inputs/'
# Getting list of FAO regions
fao_lme_code = [f for f in os.listdir(base_dir) if 'fao_lme' in f]
# Using higher resolution inputs only - Non-spatial calibration is done only once
res = '025deg'
# Choose whether smoothing of inputs will be performed by LOESS (smoothed) or
# deseasoning data (deseasoned)
smoothing = None

# Set variables to find correct input files based on 'smoothing' variable - Also
# used to name processed inputs
if smoothing != None:
    smoothing = f'-{smoothing}'
    weighted_fn = '-smoothed'
else:
    smoothing = ''
    weighted_fn = ''

# Apply data processing workflow to all regions - High resolution only
for aoi in fao_lme_code:
    gridded_folder = os.path.join(base_dir, aoi, f'gridded{smoothing}', res)
    
    #Extracting data for FAO-LME area
    #Load area of grid cell to be used as weights
    weights = xr.open_zarr(
        glob(os.path.join(gridded_folder, '*area*'))[0]).cellareao
    #Areas outside LME need to be given a value of 0
    weights = weights.fillna(0)
    
    region_int = aoi.replace('-', ' ').upper()
    
    obs_fn = glob(os.path.join(gridded_folder, f'gfdl*obsclim*'))
    weighted_obs_df = uf.weighted_mean_timestep(obs_fn, weights, region_int)
    
    ctrl_fn = glob(os.path.join(gridded_folder, f'gfdl*ctrlclim*'))
    weighted_ctrl_df = uf.weighted_mean_timestep(ctrl_fn, weights, region_int)
    
    spinup_fn = (glob(os.path.join(gridded_folder, f'gfdl*spinup*'))+
                 glob(os.path.join(gridded_folder, f'*ctrlclim_deptho*')))
    weighted_spinup_df = uf.weighted_mean_timestep(spinup_fn, weights, region_int)
    
    stable_fn = (glob(os.path.join(gridded_folder, f'gfdl*stable-spin*'))+
                 glob(os.path.join(gridded_folder, f'*ctrlclim_deptho*')))
    weighted_stable_spin_df = uf.weighted_mean_timestep(stable_fn, weights, region_int)

    #Defining output folder
    gfdl_out = os.path.join(base_dir, f'{aoi}/monthly_weighted{smoothing}')
    os.makedirs(gfdl_out, exist_ok = True)

    #Saving data
    weighted_obs_df.to_parquet(
        os.path.join(gfdl_out, 
                     f'obsclim_dbpm_clim-inputs{weighted_fn}_{aoi}_1961-2010.parquet'), 
        index = False)
    
    weighted_ctrl_df.to_parquet(
        os.path.join(gfdl_out,
                     f'ctrlclim_dbpm_clim-inputs{weighted_fn}_{fao}_1961-2010.parquet'), 
        index = False)
    
    weighted_spinup_df.to_parquet(
        os.path.join(gfdl_out,
                     f'spinup_dbpm_clim-inputs{weighted_fn}_{aoi}_1841-1960.parquet'), 
        index = False)
    
    weighted_stable_spin_df.to_parquet(
        os.path.join(gfdl_out, 
                     f'stable-spin_dbpm_clim-inputs{weighted_fn}_{aoi}_1741-1840.parquet'), 
        index = False)
