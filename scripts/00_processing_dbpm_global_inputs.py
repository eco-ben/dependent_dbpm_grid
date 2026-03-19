#!/usr/bin/env python3

# Loading libraries
import os
import numpy as np
from glob import glob
import useful_functions as uf
import xarray as xr
import dask
from distributed import Client
from multiprocessing import Process, freeze_support

# This script pre-processes outputs from GFDL-MOM6-COBALT2 (referred to as GFDL 
# from hereon) prior to their use as forcings in DBPM. GFDL outputs are  
# available at two horizontal resolutions (1 deg and 0.25 deg).
# Here, GFDL outputs are stored as zarr files to process data faster and 
# phytoplankton intercept and slope, as well as export ratio are calculated.

#Start cluster
if __name__ == '__main__':
    freeze_support()

    #Start a cluster
    client = Client(threads_per_worker = 1, memory_limit = 0)

    #Base folder where GFDL outputs are stored 
    base_dir = '/g/data/vf71/fishmip_inputs/ISIMIP3a/global_inputs'

    #Define location of area of grid cell
    grid_dir = '/g/data/vf71/shared_resources/grid_cell_vars_ESMs/isimip3a'
    
    #Define experiments and resolution
    exp_name = ['obsclim', 'ctrlclim']
    resolutions = ['1deg', '025deg']

    #Define variables of interest
    dbpm_var = ['deptho', 'phyc', 'phypico', 'siconc', 'deptho', 'expc-bot', 'tob', 'tos']
    
    #Loop through experiments and resolutions
    for res in resolutions:
        #Define output folder
        gfdl_out = f'/g/data/vf71/fishmip_inputs/ISIMIP3a/global_gridded_zarr/{res}'
        os.makedirs(gfdl_out, exist_ok = True)
        
        #Store degrees to arcmin
        if res == '1deg':
            arc_res = '60arcmin'
        elif res == '025deg':
            arc_res = '15arcmin'

        # Process area of grid cell files
        [area_file] = glob(os.path.join(grid_dir, f'gfdl*_areacello_{arc_res}*.nc'))
        f_out = os.path.basename(area_file).replace('.nc', '.zarr')
        f_out = os.path.join(gfdl_out, f_out)
        uf.netcdf_to_zarr(area_file, f_out)
        #Load file area to use as land mask for sea ice processing step
        area = xr.open_zarr(f_out)['cellareao']

        # Process thickness of grid cell files
        [thick_file] = glob(os.path.join(grid_dir, f'gfdl*_thkcello_{arc_res}*.nc'))
        f_out = os.path.basename(thick_file).replace('.nc', '.zarr')
        f_out = os.path.join(gfdl_out, f_out)
        uf.netcdf_to_zarr(thick_file, f_out)

        #Process all other variables
        for exp in exp_name:
            #Define folder containing netCDF files
            gfdl_folder = os.path.join(base_dir, exp, res)
            base_fn = f'gfdl-mom6-cobalt2_{exp}_var_{arc_res}_global_monthly_1961_2010.zarr'
            for var in dbpm_var:
                [gfdl_file] = glob(os.path.join(gfdl_folder, f'*clim_{var}_*'))
                f_out = os.path.basename(gfdl_file).replace('.nc', '.zarr')
                f_out = os.path.join(gfdl_out, f_out)
                #Apply function
                uf.netcdf_to_zarr(gfdl_file, f_out)

            # Vertically integrate phytoplankton inputs up to threshold depth
            phyc, phypico = uf.integrating_phyto(gfdl_out, exp, thresh_depth = 200)
            #Save outputs
            phyc.to_zarr(
                os.path.join(gfdl_out, base_fn.replace('_var_', '_phyc-vint200m_')), 
                consolidated = True, mode = 'w')
            phypico.to_zarr(
                os.path.join(gfdl_out, base_fn.replace('_var_', '_phypico-vint200m_')), 
                consolidated = True, mode = 'w')

            #Calculate phytoplankton size distribution and export ratio
            sphy, lphy, er = uf.getExportRatio(gfdl_out, exp)
            #Save outputs
            sphy.to_zarr(
                os.path.join(gfdl_out, base_fn.replace('_var_', '_sphy_')), 
                consolidated = True, mode = 'w')
            lphy.to_zarr(
                os.path.join(gfdl_out, base_fn.replace('_var_', '_lphy_')),
                consolidated = True, mode = 'w')
            er.to_zarr(
                os.path.join(gfdl_out, base_fn.replace('_var_', '_er_')),
                consolidated = True, mode = 'w')
            
            #Calculate intercept and slope
            intercept, slope = uf.GetPPIntSlope(gfdl_out, exp)
            #Save outputs
            intercept.to_zarr(
                os.path.join(gfdl_out, base_fn.replace('_var_', '_intercept_')), 
                consolidated = True, mode = 'w')
            slope.to_zarr(
                os.path.join(gfdl_out, base_fn.replace('_var_', '_slope_')), 
                consolidated = True, mode = 'w')

            #Create sea ice masks 
            #Any grid cells with sea ice concentration of 15% or above will not be
            #available for fishing
            #Identify sea ice files
            si_files = glob(os.path.join(gfdl_out, '*siconc*'))

            #Create masks for all sea ice files
            for f in si_files:
                #Load files
                da = xr.open_zarr(f)['siconc'].where(np.isfinite(area))
                #Identify grid cells with at least 15% SIC
                da_mask = xr.where(da >= 15, True, False)
                #Split into northern and southern hemispheres before calculating cumulative sum
                
                #Sea ice kept from 42N towards the north pole as the Sea of Okhotsk (45N) is the
                #lowest latitude area where sea ice forms each winter according to NASA's 
                #Earth Observatory
                da_mask_north = (xr.where(da_mask.lat > 42, da_mask, False).
                    isel(lat = slice(None, None, -1)).cumsum('lat')))
                #Sea ice kept from 52S towards the south pole as 55S is the lowest latitude area
                # where sea ice forms each winter according to NASA's Earth Observatory
                da_mask_south = xr.where(da_mask.lat <= -52, da_mask, False).cumsum('lat')
                #Create a single global mask
                da_mask = (da_mask_north+da_mask_south)
                #Remove any grid cells that are on land or have sea ice within their boundaries
                da_mask = xr.where(da_mask > 0, np.nan, 1).where(np.isfinite(area))
                #Rechunk data
                da_mask = da_mask.chunk({'time': '500MB'})
                #Update data array variable name
                da_mask.name = 'simask'
                #Create file path to save outputs
                f_out = f.replace('siconc', 'simask')
                #Save results
                da_mask.to_zarr(f_out, consolidated = True, mode = 'w')
