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
    client = Client(threads_per_worker = 1)

    ## Name of region and model resolution ----
    region = 'fao-48'
    model_res = '025deg'
    
    ## Defining input and output folders ----
    base_folder = '/g/data/vf71/la6889/dbpm_inputs/weddell/'
    gridded_folder = os.path.join(base_folder, 'gridded_params', model_res)
    out_folder = os.path.join(base_folder, 'run_fishing', model_res)
    #If output folder does not exist, it will create it
    os.makedirs(out_folder, exist_ok = True) 

    ## If starting DBPM run from a specific time step ----
    # Character: Year and month from when DBPM initialisation values should be loaded
    # If starting model for the first time, it should be set to None
    init_time = '2010-10'
    
    ## Loading fixed DBPM parameters ----
    ds_fixed = uf.loading_dbpm_fixed_inputs(gridded_folder)
    #Adding additional fixed DBPM parameters to dataset
    #Depth
    depth = xr.open_zarr(glob(os.path.join(base_folder, 'gridded',
                                           model_res, '*obsclim_deptho_*'))[0])['deptho']
    ds_fixed['depth'] = depth
    #Create a land mask from depth data
    ds_fixed['mask'] = np.isfinite(depth)
    #Size bins in log10
    log10_size_bins_mat = xr.open_zarr('outputs/log10_size_bins_matrix.zarr/')['size_bins']
    ds_fixed['log10_size_bins'] = log10_size_bins_mat
    ds_fixed['size_bin_vals'] = 10**log10_size_bins_mat
    #Removing datarrays added to fixed inputs
    del depth, log10_size_bins_mat
    
    ## Loading predator, detritivores and detritus initialisation data ----
    if init_time is None:
        ds_init = uf.loading_dbpm_biomass_inputs(gridded_folder)
    else:
        ds_init = uf.loading_dbpm_biomass_inputs(out_folder, init_time)

    ## Loading dynamic data ----
    gridded_esm = os.path.join(base_folder, 'gridded', model_res)
    ds_dynamic = uf.loading_dbpm_dynamic_inputs(gridded_esm, gridded_folder,
                                                init_time = init_time)
  
    if init_time is not None:
        init_yr = pd.Timestamp(init_time).year
        #Timestep from when to restart DBPM 
        subset_time = (pd.Timestamp(init_time)+pd.DateOffset(months = 1)).strftime('%Y-%m')
        #Timestep from when to add init effort data
        effort_time = (pd.Timestamp(init_time)+pd.DateOffset(months = 2)).strftime('%Y-%m')
        
    #Spinup data is loaded if init_time is None or if the init_time year is less than 1960
    if init_time is None or init_yr < 1959:
        # Loading effort
        effort = xr.open_zarr(glob(os.path.join(gridded_folder, 'effort_spinup*'))[0])['effort']
    #Spinup data plus obsclim are loaded if init_time is 1960
    elif init_yr >= 1959 and init_yr < 1961:
        # Loading effort
        effort = xr.open_mfdataset(glob(os.path.join(gridded_folder, 'effort_*')),
                                   engine = 'zarr')['effort']
    #Obsclim data loaded from 1961 onwards:
    else:
        # Loading effort
        effort = xr.open_zarr(glob(os.path.join(gridded_folder, 'effort_[0-9]*'))[0])['effort']
    
    if init_time is not None:
        #Load effort for time step DBPM starts
        e_start = xr.open_dataarray(glob(os.path.join(out_folder, 
                                                      f'effort_*_{subset_time}.nc'))[0])
        #Subset effort data from the timestep after DBPM restart 
        effort = effort.sel(time = slice(effort_time, None))
        #Combine both data arrays
        effort = xr.concat([e_start, effort], dim = 'time')
        effort = effort.chunk({'lat': len(effort.lat), 
                               'lon': len(effort.lon),
                               'time': -1})
    
    #Creating a single dataset for all dynamic inputs
    ds_dynamic['effort'] = effort
    ds_dynamic = ds_dynamic.chunk({'time': 1})

    #Gridded parameters
    gridded_params = json.load(open(
        glob(os.path.join(gridded_folder, 'dbpm_gridded_*_python.json'))[0]))

    ## Running spatial DBPM ----
    for t in range(0, len(ds_dynamic.time)):
        ds_dyn = ds_dynamic.isel(time = t)
        # Redistribute total effort across grid cells 
        try:
            eff_short = uf.effort_calculation(ds_init['predators'], ds_init['detritivores'], 
                                              ds_dynamic['effort'].isel(time = t+1), 
                                              ds_fixed['depth'], 
                                              ds_fixed['log10_size_bins'], 
                                              gridded_params)
            # Saving predation mortality
            #Getting year and month 
            dt_eff = pd.to_datetime(eff_short.time.values[0]).strftime('%Y-%m')
            # Creating file name
            fn = f'effort_{model_res}_{region}_{dt_eff}.nc'
            eff_short.to_netcdf(os.path.join(out_folder, fn))
            ds_dynamic['effort'] = xr.where(ds_dynamic.time == ds_dynamic.time[t+1], 
                                            eff_short.values, ds_dynamic['effort'])
            #Remove variables not needed
            del dt_eff, fn
        except:
            dt = pd.to_datetime(ds_dyn.time.values).strftime('%Y-%m')
            eff_short = xr.open_dataarray(glob(os.path.join(out_folder, f'effort*{dt}*'))[0])
            ds_dynamic['effort'] = xr.where(ds_dynamic.time == ds_dynamic.time[t], 
                                            eff_short.values, ds_dynamic['effort'])
        
        ds_init = uf.gridded_sizemodel(gridded_params, ds_fixed, ds_init, 
                                       ds_dyn, region = region, model_res = model_res, 
                                       out_folder = out_folder)

