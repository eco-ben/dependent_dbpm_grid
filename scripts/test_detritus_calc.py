#!/usr/bin/env python3

#Loading libraries
import os
from glob import glob
import pandas as pd
import xarray as xr
import useful_functions as uf
import dask
from distributed import Client
from multiprocessing import Process, freeze_support
import json
import numpy as np

if __name__ == '__main__':
    freeze_support()

    #Start a cluster
    client = Client(threads_per_worker = 1, memory_limit = 0)
    
    ## Name of region and model resolution ----
    reg = 'fao-34'
    res = '1deg'

    ## Define if run should include fishing ----
    fishing = True
    # Choose whether smoothing of inputs will be performed by LOESS (smoothed) or
    # deseasoning data (deseasoned)
    smoothing = 'deseasoned'
    ## If starting DBPM run from a specific time step ----
    # Character: Year and month from when DBPM initialisation values should be loaded
    # If starting model for the first time, it should be set to None
    init_time = None

    # Set variables to find correct input files based on 'smoothing' variable - Also
    # used to name processed inputs
    # Create variables based on 'smoothing' selection
    if smoothing != None:
        smoothing = f'-{smoothing}'
    else:
        smoothing = ''

    ## Defining input and output folders ----
    base_folder = '/g/data/vf71/fishmip_inputs/ISIMIP3a/fao_inputs'
    base_out_folder = '/g/data/vf71/fishmip_outputs/ISIMIP3a/fao_outputs'
    
    #Location of gridded inputs
    gridded_folder = os.path.join(base_folder, reg, f'gridded{smoothing}', res)

    #Folder where outputs will be stored 
    out_folder = f'fao-34_problem_cell_fishing{smoothing}_full_dbpm_euler'
    # out_folder = f'fao-34_problem_cell_no-fishing{smoothing}_full_dbpm_euler'
    # out_folder = f'fao-34_good_cell_fishing{smoothing}_full_dbpm_euler'
    # out_folder = f'fao-34_good_cell_no-fishing{smoothing}_full_dbpm_rk4'
    #If output folder does not exist, it will create it
    os.makedirs(out_folder, exist_ok = True) 
    
    ## Loading fixed DBPM parameters ----
    ds_fixed = uf.loading_dbpm_fixed_inputs(gridded_folder)
    #Adding additional fixed DBPM parameters to dataset
    #Size bins in log10
    log10_size_bins_mat = xr.open_zarr(os.path.join(
        base_folder, 'log10_size_bins_matrix.zarr/'))['size_bins']
    ds_fixed['log10_size_bins'] = log10_size_bins_mat
    ds_fixed['size_bin_vals'] = 10**log10_size_bins_mat
    #Removing datarrays added to fixed inputs
    del log10_size_bins_mat

    ## Loading predator, detritivores and detritus initialisation data ----
    if init_time is None:
        ds_init = uf.loading_dbpm_biomass_inputs(gridded_folder)
    else:
        ds_init = uf.loading_dbpm_biomass_inputs(out_folder, init_time)
    
    ## Loading dynamic data ----
    ds_dynamic = uf.loading_dbpm_dynamic_inputs(gridded_folder, init_time, 
                                                fishing = fishing)
  
    if init_time is not None and fishing:
        #Timestep from when to restart DBPM 
        subset_time = (pd.Timestamp(init_time)+
                       pd.DateOffset(months = 1)).strftime('%Y-%m')
        #Timestep from when to add init effort data
        effort_time = (pd.Timestamp(init_time)+
                       pd.DateOffset(months = 2)).strftime('%Y-%m')
        #Load effort for time step DBPM starts
        e_start = xr.open_dataarray(glob(os.path.join(
            out_folder, f'effort_*_{subset_time}.nc'))[0])
        #Subset effort data from the timestep after DBPM restart 
        effort = ds_dynamic['effort'].sel(time = slice(effort_time, None))
        #Combine both data arrays
        effort = xr.concat([e_start, effort], dim = 'time')
        effort = effort.chunk({'lat': -1, 'lon': -1, 'time': -1})
    
        #Creating a single dataset for all dynamic inputs
        ds_dynamic['effort'] = effort.load()
    
    #Gridded parameters
    gridded_params = json.load(open(
        os.path.join(base_folder, reg, 'init_fish_vals', res, 
                     f'best_fish_vals{smoothing}', 
                     f'dbpm_gridded_size_params_{reg}_python.json')))

    ds_init = ds_init.sel(lon = 11.5, lat = -4.5)
    ds_fixed = ds_fixed.sel(lon = 11.5, lat = -4.5)
    ds_dynamic = ds_dynamic.sel(time = slice(None, '1840'), lon = 11.5, lat = -4.5)

    ## Running spatial DBPM ----
    for t in range(0, len(ds_dynamic.time)):
        ds_dyn = ds_dynamic.isel(time = t)
        # Redistribute total fishing effort across grid cells (if fishing is on)
        if fishing:
            try:
                eff_short = uf.effort_calculation(ds_init['predators'], 
                                                  ds_init['detritivores'], 
                                                  ds_dynamic['effort'].isel(time = t+1), 
                                                  ds_fixed['depth'], 
                                                  ds_fixed['log10_size_bins'], 
                                                  gridded_params)
                #Getting year and month 
                dt_eff = pd.to_datetime(eff_short.time.values[0]).strftime('%Y-%m')
                # Creating file name
                fn = f'effort_{res}_{reg}_{dt_eff}.nc'
                eff_short.to_netcdf(os.path.join(out_folder, fn))
                ds_dynamic['effort'] = xr.where(ds_dynamic.time == ds_dynamic.time[t+1], 
                                                eff_short.values, ds_dynamic['effort'])
                #Remove variables not needed
                del dt_eff, fn
            except:
                dt = pd.to_datetime(ds_dyn.time.values).strftime('%Y-%m')
                eff_short = xr.open_dataarray(glob(os.path.join(out_folder, 
                                                                f'effort*{dt}*'))[0])
                ds_dynamic['effort'] = xr.where(ds_dynamic.time == ds_dynamic.time[t], 
                                                eff_short.values, ds_dynamic['effort'])
        
        ds_init = uf.gridded_sizemodel(gridded_params, ds_fixed, ds_init, 
                                       ds_dyn, region = reg, model_res = res, 
                                       out_folder = out_folder)
        # ds_init = uf.gridded_sizemodel_rk4(gridded_params, ds_fixed, ds_init, 
        #                                    ds_dyn, region = reg, model_res = res, 
        #                                    out_folder = out_folder, rk4_substeps = 4)
        # ds_init = uf.gridded_sizemodel_solve_ivp(gridded_params, ds_fixed, ds_init, 
        #                                          ds_dyn, region = reg, model_res = res, 
        #                                          out_folder = out_folder)

