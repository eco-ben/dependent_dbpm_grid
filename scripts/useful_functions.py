#!/usr/bin/env python3

# Library of useful functions
# Authors: Denisse Fierro Arcos
# Date last update: 2025-03-25

# Loading libraries
import xarray as xr
import re
import numpy as np
import os
from glob import glob
import pandas as pd
import json
from statsmodels.nonparametric.smoothers_lowess import lowess
from statsmodels.tsa.seasonal import STL

#Transforming netCDF files to zarr
def netcdf_to_zarr(file_path, path_out):
    '''
    Inputs:
    - file_path (character) File path where GFDL output is located
    - path_out (character) File path where outputs should be stored as zarr files

    Outputs:
    - None. This function saves results as zarr files in the path provided.
    '''

    # Loading and rechunking data
    da = xr.open_dataarray(file_path)

    # Rechunking before saving 
    if 'time' in da.dims:
        da = da.chunk({'time': '500MB', 'lat': -1, 'lon': -1})
        # Change date format to DateTimeIndex if needed
        try:
            da['time'] = pd.DatetimeIndex(da.indexes['time'].to_datetimeindex(),
                                          dtype = 'datetime64[ns]')
        except:
            pass
    else:
        da = da.chunk({'lat': -1, 'lon': -1})

    #Save results
    da.to_zarr(path_out, consolidated = True, mode = 'w')


# Apply LOESS smoothing to gridded data
def loess_xarray(da, time_dim, **kwargs):
    '''
    Inputs:
    - da (Data Array) Three dimensional xarray data array (time, lat, lon) to be smoothed.
    - time_dim (Data Array) Time dimension from "da" input.
    **Optional**: 
    - frac (numeric) Proportion of timesteps to be used in data smoothing. If not provided,
    it defaults to 0.5

    Outputs:
    smooth (1D numpy array) Containing smoothed "da" values.
    '''
    # If "frac" is not provided, default to 0.5
    if kwargs.get('frac') is not None:
        frac = kwargs['frac']
    else:
        frac = 0.5

    # If "it" is not provided, default to 0
    if kwargs.get('it') is not None:
        it = kwargs['it']
    else:
        it = 0

    #Turn time_sim and da to numpy arrays
    data_vals = np.asarray(da)
    time_vals = np.asarray(time_dim)

    #Applying LOESS smoother - Keep only smoothed values (i.e., ignore time column)
    smooth = lowess(data_vals, time_vals, is_sorted = True, frac = frac, it = it, 
                    missing = 'none')[:, 1]
    return smooth


# Smoothing DBPM inputs using LOESS
def smoothing_loess(file_path, path_out, **kwargs):
    '''
    - file_path (character) File path where GFDL zarr file is located
    - path_out (character) File path where outputs should be stored as zarr files
    **Optional**: 
    - mth_smooth (numeric) Number of timesteps to be used in data smoothing. If not
    provided, it defaults to 4 timesteps.

    Outputs:
    - None. This function saves results as zarr files in the path provided.    
    '''
    # Loading data
    da = xr.open_zarr(file_path)
    #Getting name of variable contained in dataset
    [var] = list(da.keys())
    da = da[var].chunk({'time': -1, 'lat': '500MB', 'lon': '500MB'})
    
    # If "frac" is not provided, default to 0.5
    if kwargs.get('mth_smooth') is not None:
        frac = (kwargs['mth_smooth']/len(da.time))
    else:
        frac = (4/len(da.time))
    
    # If "it" is not provided, default to 0
    if kwargs.get('it') is not None:
        it = kwargs['it']
    else:
        it = 0
    
    # Define function being used first
    da_smooth = xr.apply_ufunc(
        loess_xarray, 
        # Pass arguments needed by function being applied
        da, da.time, kwargs = {'frac': frac, 'it': it},
        # Dimensions the function needs for each argument (i.e.,
        # the function will be applied along the time dimension)
        input_core_dims = [['time'], ['time']], 
        # Dimensions of each output from the function
        output_core_dims = [['time']],
        # Allows function needs to applied along time dimension
        # (i.e., one lat and lon at a time)
        vectorize = True,
        # Parallise function using Dask
        dask = 'parallelized',
        # Allow rechunking of data if needed
        dask_gufunc_kwargs = {'allow_rechunk': True},
        # Define the type of object the function will return
        output_dtypes = [float]).load()
    # Add metadata and record smoothing step
    da_smooth = da_smooth.assign_attrs(da.attrs)
    da_smooth.attrs['fishmip_postprocess'] = f'LOESS smooth applied to data, span = {frac}, it = {it}'
    # Add variable name
    da_smooth.name = var
    # Save results
    da_smooth.to_zarr(path_out, consolidated = True, mode = 'w')


# Apply seasonal trend decomposition using LOESS (STL) to gridded data
def stl_xarray(da, **kwargs):
    '''
    Inputs:
    - da (Data Array) - 1 dimensional vector from a data array
    **Optional**
    - period (integer) - Default is 12 (annual). The periodicity of the time series
    (e.g., 12 for monthly, 365 for daily).
    - component (character) - Default is seasonal. Returns specified STL component

    Outputs:
    res - Specified STL component
    '''
    #Get data into correct format to use an input in STL function
    data_vals = np.asarray(da)

    #If more than half the values are NAs, then return NAs
    if np.isnan(data_vals).sum() > len(data_vals)*0.5:
        return np.full_like(data_vals, np.nan)
    
    if kwargs.get('period') is not None:
        period = kwargs['period']
    else:
        period = 12

    if kwargs.get('component') is not None:
        component = kwargs['component']
    else:
        component = 'seasonal'
    
    #Seasonal trend decomposition
    res = STL(data_vals, period = period).fit()
    if component == 'trend':
        return res.trend
    elif component == 'seasonal':
        return res.seasonal
    elif component == 'resid':
        return res.resid
    else:
        raise ValueError("component must be 'trend', 'seasonal', or 'resid'")


# Applying STL decomposition to an xarray DataArray along the 'time' dimension
def seasonal_decomposition(file_path, path_out, period, component):
    '''
    - file_path (character) File path where GFDL zarr file is located
    - period (integer) The periodicity of the time series (e.g., 12 for monthly, 
    365 for daily).
    - component (character) Default is seasonal. Returns specified STL component

    Outputs:
    - None. A new DataArray containing the original data minus the specified 
    decomposed component is saved in the file path provided in path_out.  
    '''
    # Loading data
    da = xr.open_zarr(file_path)
    #Getting name of variable contained in dataset
    [var] = list(da.keys())
    da = da[var].chunk({'time': -1, 'lat': '500MB', 'lon': '500MB'})

    # Use apply_ufunc to apply the wrapper function across spatial dimensions
    da_comp = xr.apply_ufunc(stl_xarray, da, 
                             kwargs = {'period': period, 'component': component},
                             input_core_dims = [['time']], 
                             output_core_dims = [['time']],
                             dask = 'parallelized', vectorize = True, 
                             output_dtypes = [da.dtype],
                             # Allow rechunking of data if needed
                             dask_gufunc_kwargs = {'allow_rechunk': True}).load()

    #Remove calculated component from original data
    da_decomp = da-da_comp
    # Add metadata and record smoothing step
    da_decomp = da_decomp.assign_attrs(da.attrs)
    da_decomp.attrs['fishmip_postprocess'] = f'Original data minus {component} component'
    # Add variable name
    da_decomp.name = var
    # Save results
    da_decomp.to_zarr(path_out, consolidated = True, mode = 'w')

    
## Extracting GFDL outputs for region of interest using boolean mask
def extract_gfdl(file_path, mask, path_out, cross_dateline = False):
    '''
    Inputs:
    - file_path (character) File path where GFDL zarr file is located
    - mask (boolean data array) Grid cells within region of interest should be identified
    as 1.
    - path_out (character) File path where outputs should be stored as zarr files
    - cross_dateline (boolean) Default is False. If set to True, it will convert longitudes 
    from +/-180 to 0-360 degrees before extracting data for region of interest
    
    Outputs:
    - None. This function saves results as zarr files in the path provided.
    '''

    #Loading and rechunking data
    da = xr.open_zarr(file_path)

    #Fix time format if needed
    try:
        new_time = da.indexes['time'].to_datetimeindex()
        da['time'] = new_time
    except:
        pass

    #Getting name of variable contained in dataset
    [var] = list(da.keys())
    da = da[var]
    
    #Apply mask and remove rows where all grid cells are empty to 
    #reduce data array size
    if cross_dateline:
        da = da.where(mask == 1)
        da['lon'] = da.lon%360
        da = da.sortby('lon')
        da = (da.dropna(dim = 'lon', how = 'all').
              dropna(dim = 'lat', how = 'all')).load()
    else:
        da = da.where(mask == 1, drop = True).load()
    
    #Rechunking data
    if 'time' in da.dims:
        da = da.chunk({'time': -1, 'lat': -1, 'lon': -1})
    else:
        da = da.chunk({'lat':-1, 'lon': -1})
    
    #Save subsetted data
    da.to_zarr(path_out, consolidated = True, mode = 'w')


## Calculating area weighted means
def weighted_mean_timestep(file_paths, weights, region):
    '''
    Inputs:
    - file_paths (list) File paths pointing to zarr files from which weighted means 
    will be calculated and stored in a single data frame
    - weights (data array) Contains the weights to be used when calculating weighted mean.
    It should NOT include NaN, zeroes (0) should be used instead.
    - region (character) Name of the region to be recorded in data frame

    Outputs:
    - df (pandas data frame) containing area weighted mean per time step
    '''
    
    #Loading all zarr files into a single dataset
    da = xr.open_mfdataset(file_paths, engine = 'zarr')
    #Fix time format if needed
    try:
        new_time = da.indexes['time'].to_datetimeindex()
        da['time'] = new_time
    except:
        pass

    #Apply weights
    da_weighted = da.weighted(weights)
    #Calculate weighted mean
    da_weighted_mean = da_weighted.mean(('lat', 'lon'))
    
    #Transform to data frame
    df = da_weighted_mean.to_dataframe().reset_index()

    #Getting names of variables in data frame
    col_names = [i for i in df.columns if i != 'time']
    #Getting name of experiment from file path
    [exp] = re.findall('cobalt2_(.*?)_', file_paths[0])

    #If no depth file is included, add variable and leave it empty
    if 'deptho' not in col_names:
        df['depth_m'] = np.nan
    else:
        col_names.remove('deptho')
        df = df.rename(columns = {'deptho': 'depth_m'})
    
    #Add metadata to data frame
    df['tot_area_m2'] = weights.values.sum()
    df['year'] = df.apply(lambda x: x.time.year, axis = 1)
    df['month'] = df.apply(lambda x: x.time.strftime('%B'), axis = 1)
    df['region'] = region
    df['scenario'] = exp

    #Rearrange columns 
    names = ['region', 'scenario', 'time', 'year', 'month', 
             'depth_m', 'tot_area_m2'] + col_names
    df = df[names]

    return df


#Calculating export ratio
def get_threshold_depth(folder_gridded_data, gfdl_exp, max_depth = 200):
    '''
    Inputs:
    - folder_gridded_data (character) File path pointing to folder containing
    zarr files with GFDL data for the region of interest
    - gfdl_exp (character) Select GFDL experiment 'ctrl_clim' or 'obs_clim'
    - max_depth (numeric) Default is 200 (m) Maximum depth in meters to be considered when
    processing DBPM phytoplankton inputs
    
    Outputs:
    - thresh_depth (data array) Contains the maximum depth to integrate phytoplankton
    inputs
    '''
    
    #load depth
    depth = xr.open_zarr(glob(
        os.path.join(folder_gridded_data, f'*{gfdl_exp}_deptho_*'))[0])['deptho']
    
    # Identify areas where depth is 'max_depth' or less
    depth_n = depth.where(depth <= max_depth)

    #Load mixed layer depth
    mld = xr.open_zarr(glob(
        os.path.join(folder_gridded_data, f'*{gfdl_exp}_mlotst-0125_*'))[0])['mlotst-0125']

    # Check if MLD depth in shallow areas (<= 100 m depth) is less 100 m
    # Keep the largest depth (water column depth vs MLD)
    thresh_depth = np.maximum(depth_n, mld)

    # For all other areas use MLD
    thresh_depth = xr.where(np.isnan(thresh_depth), mld, thresh_depth)

    # Store description for 'long name' attribute
    # Updating metadata
    thresh_depth.name = 'threshold_depth'
    thresh_depth.attrs = {'standard_name': f'threshold_depth_mld_{max_depth}m',
                          'long_name':
                          f'Threshold depth based on MLD and at least {max_depth}m in shallow waters', 
                          'units': 'm'}
    return thresh_depth


# Integrating (mean) phytoplankton inputs to threshold depth
def integrating_phyto(folder_gridded_data, gfdl_exp, thresh_depth = 200):
    '''
    Inputs:
    - folder_gridded_data (character) File path pointing to folder containing
    zarr files with vertically resolved phytopklankton outputs (`phyc` and `phypico`) 
    from GFDL
    - gfdl_exp (character) Select GFDL experiment 'ctrl_clim' or 'obs_clim'
    - thresh_depth (numeric) Default is 200 (m) Maximum depth in meters to be considered 
    when processing DBPM phytoplankton inputs

    Outputs:
    - phyc (data array) Contains integrated phytoplankton values up to threshold depth.
    - phypico (data array) Contains integrated picophytoplankton values up to threshold
    depth.
    '''

    #load depth
    depth = (xr.open_zarr(glob(os.path.join(folder_gridded_data, 
                                            f'*_thkcello_*'))[0])['thkcello'].
        drop_vars('time').squeeze().sel(lev = slice(None, thresh_depth)).fillna(0))
    
    #Load phypico
    phypico = (xr.open_zarr(glob(
        os.path.join(folder_gridded_data, f'*{gfdl_exp}_phypico_*'))[0])['phypico'].
        sel(lev = slice(None, thresh_depth)))
    phypico_weighted = phypico.weighted(depth).mean('lev')
    phypico_weighted = phypico_weighted.assign_attrs(
        {'standard_name': 'mean_mole_concentration_of_picophytoplankton_expressed_as_carbon_in_sea_water',
         'long_name': 'Depth Weighted Mean of Picophytoplankton Carbon Concentration', 
         'units': 'mol m-3'})
    
    #Load phyc
    phyc = (xr.open_zarr(glob(
        os.path.join(folder_gridded_data, f'*{gfdl_exp}_phyc_*'))[0])['phyc'].
        sel(lev = slice(None, thresh_depth)))
    phyc_weighted = phyc.weighted(depth).mean('lev')
    phyc_weighted = phyc_weighted.assign_attrs(
        {'standard_name': 'mean_mole_concentration_of_phytoplankton_expressed_as_carbon_in_sea_water',
         'long_name': 'Depth Weighted Mean of Phytoplankton Carbon Concentration', 
         'units': 'mol m-3'})
    
    return phyc_weighted, phypico_weighted


#Calculating export ratio
def getExportRatio(folder_gridded_data, gfdl_exp):
    '''
    Inputs:
    - folder_gridded_data (character) File path pointing to folder containing
    zarr files with GFDL data for the region of interest
    - gfdl_exp (character) Select GFDL experiment 'ctrl_clim' or 'obs_clim'

    Outputs:
    - sphy (data array) Contains small phytoplankton. This is data frame simply
    renames 'phypico-vint' data
    - lphy (data array) Contains large phytoplankton: difference between 'phyc_vint'
    and 'phypico_vint'
    - er (data array) Contains export ratio
    '''

    #Load sea surface temperature
    tos = xr.open_zarr(glob(
        os.path.join(folder_gridded_data, f'*{gfdl_exp}_tos_*'))[0])['tos']
    
    #load depth
    depth = xr.open_zarr(glob(
        os.path.join(folder_gridded_data, f'*{gfdl_exp}_deptho_*'))[0])['deptho']
    
    #Load phypico-vint
    sphy = xr.open_zarr(glob(
        os.path.join(folder_gridded_data, 
                     f'*{gfdl_exp}_phypico-vint200m_*'))[0])['phypico-vint']
    #Rename phypico-vint to sphy
    sphy.name = 'sphy'
    sphy = sphy.assign_attrs({'short_name': 'sphy',
                              'long_name': 'Small phytoplankton carbon content', 
                              'units': 'mol m-3'})

    #Load phyc-vint
    ptotal = xr.open_zarr(glob(
        os.path.join(folder_gridded_data, f'*{gfdl_exp}_phyc-vint200m_*'))[0])['phyc-vint']

    #Calculate large phytoplankton
    lphy = ptotal-sphy
    #Give correct name to large phytoplankton dataset
    lphy.name = 'lphy'
    lphy = lphy.assign_attrs(ptotal.attrs)
    lphy = lphy.assign_attrs({'short_name': 'lphy',
                              'long_name': 'Large phytoplankton carbon content', 
                              'units': 'mol m-3'})

    #Calculate phytoplankton size ratios
    plarge = (lphy/ptotal)
    psmall = (sphy/ptotal)

    #Calculate export ratio
    er = (np.exp(-0.032*tos)*((0.14*psmall)+(0.74*(plarge)))+
          (0.0228*(plarge)*(depth*0.004)))/(1+(depth*0.004))
    #If values are negative, assign a value of 0
    er = xr.where(er < 0, 0, er)
    #If values are above 1, assign a value of 1
    er = xr.where(er > 1, 1, er)
    er.name = 'export_ratio'
    er = er.assign_attrs({'short_name': 'export_ratio_expressed_as_carbon_in_sea_water',
                          'long_name': 'Export ratio of organic matter',
                          'units': 'mol m-3'})
    
    return sphy, lphy, er


#Calculate slope and intercept
def GetPPIntSlope(gfdl_folder, gfdl_exp, mmin = 10**(-14.25), mmid = 10**(-10.184), 
                  mmax = 10**(-5.25), output = 'both'):
    '''
    Inputs:
    - gfdl_folder (character) File path pointing to folder containing
    zarr files with GFDL data for the region of interest
    - gfdl_exp (character) Select GFDL experiment 'ctrl_clim' or 'obs_clim'
    - mmin (numeric)  Default is 10**(-14.25). Minimum body mass of primary
    producers representing phycophytoplankton
    - mmid (numeric)  Default is 10**(-10.184). Body mass of primary producers
    representing limit between phycophytoplankton (small phytoplankton) and 
    microphytoplankton (large phytoplankton)
    - mmax (numeric)  Default is 10**(-5.25). Maximum body mas of primary producers
    representing microphytoplankton
    - output (character) Default is 'both'. Select what outputs should be returned. 
    Choose from 'both', 'slope', or 'intercept'

    Outputs:
    - (Data array) - Depends on value of 'output' parameter. 
    '''
    
    #load large phytoplankton
    lphy = xr.open_zarr(glob(os.path.join(gfdl_folder, 
                                          f'*{gfdl_exp}_lphy_*'))[0])['lphy']
    
    #Load small phytoplankton
    sphy = xr.open_zarr(glob(os.path.join(gfdl_folder, 
                                          f'*{gfdl_exp}_sphy_*'))[0])['sphy']
    
    #Convert sphy and lphy from mol C / m^3 to g C / m^3
    sphy = sphy*12.0107
    lphy = lphy*12.0107

    #From Appendix of Barnes 2010 JPR, used in Woodworth-Jefcoats et al 2013
    #GCB. The scaling of biomass with body mass can be described as B=aM^b
    #the exponent b (also equivalent to the slope b in a log B vs log M 
    #relationship) can be assumed:
    #0.25 (results in N ~ M^-3/4) or 0 (results in N ~ M^-1)
    #most studies seem to suggest N~M^-1, so can assume that and test 
    #sensitivity of our results to this assumption. 
      
    #Calculate a and b in log B (log10 abundance) vs. log M (log10 gww)
    #in log10 gww
    midsmall = np.log10((mmin+mmid)/2)
    midlarge = np.log10((mmid+mmax)/2)

    #convert to log10 (gww/size class median size) for log10 abundance
    small = (np.log10((sphy*10)/10**midsmall))
    #convert to log10 (gww/size class median size) for log10 abundance
    large = (np.log10((lphy*10)/10**midlarge))

    #Calculating lope
    slope = ((small-large)/(midsmall-midlarge))
    slope.name = 'slope'
    slope = slope.assign_attrs({
        'short_name': 'slope',
        'long_name': 'Slope of primary producer spectrum',
        'comment': ('Calculations described in full in Woodworth-'+ 
                    'Jefcoats et al 2013 (DOI: 10.1111/gcb.12076)')})

    #a is really log10(a), same a when small, midsmall are used
    intercept = (large-(slope*midlarge))
    intercept.name = 'intercept'
    intercept = intercept.assign_attrs({
        'short_name': 'intercept',
        'long_name': 'Intercept of primary producer spectrum',
        'comment': ('Calculations described in full in Woodworth-'+ 
                    'Jefcoats et al 2013 (DOI: 10.1111/gcb.12076)')})

    # a could be used directly to replace 10^pp in sizemodel()
    if output == 'slope':
        return slope
    if output == 'intercept':
        return intercept
    if output == 'both':
        return intercept, slope


# Calculate spinup from gridded data
def gridded_spinup(file_path_or_data_array, start_spin, end_spin, spinup_period,
                   mean_spinup = False, **kwargs):
    '''
    Inputs:
    - file_path_or_data_array (character or Data Array) Provide either file path 
    where GFDL zarr file is located or Data Array
    - start_spin (character or numeric) Year or date spinup calculation starts
    - end_spin (character or numeric) Year or date spinup calculation ends
     - spinup_period (pandas Datetime array) New time labels for spinup period.
    Must be a multiple of spinup range (i.e., difference between start and end
    spin)
    - mean_spinup (boolean) Default is False. If set to True, then the spinup
    period will be based on the mean value over the spin period.
    **Optional**: 
    - file_out (character) File path to save results

    Outputs:
    spinup_da (data array) Spinup data containing information within spinup
    range
    '''

    #Loading data
    if isinstance(file_path_or_data_array, str):
        da = xr.open_zarr(file_path_or_data_array)
        #Getting name of variable contained in dataset
        [var] = list(da.keys())
        da = da[var]
    elif isinstance(file_path_or_data_array, xr.DataArray):
        da = file_path_or_data_array
    else:
        print('Please provide either a file path or data array to run function.')

    #Select period to be used for spinup
    da = da.sel(time = slice(str(start_spin), str(end_spin)))
    #If spinup should be created based on mean values over spinup period
    if mean_spinup:
        da = da.mean('time')
        spinup_da = [da] * len(spinup_period)
    else:
        spinup_da = [da] * int(len(spinup_period)/len(da.time))

    #Create spinup data array
    spinup_da = xr.concat(spinup_da, dim = 'time')
    spinup_da['time'] = spinup_period
    
    #Updating chunks
    spinup_da = spinup_da.chunk({'time': '500MB', 'lat': -1, 'lon': -1})
    spinup_da = spinup_da.drop_encoding()

    #Save result if path is provided
    if kwargs.get('file_out', None):
        f_out = kwargs.get('file_out')
        #Ensure folder in file path exists
        os.makedirs(os.path.dirname(f_out), exist_ok = True)
        spinup_da.to_zarr(f_out, consolidated = True, mode = 'w')
    
    return spinup_da


#Format gridded_params in Python friendly way
def gridded_param_python(gridded_params):
    '''
    Inputs:
    - gridded_params (dictionary) Containing the output from the `sizeparam` function in R.
    Stored as a `json` file.

    Outputs:
    griddded_python (dictionary) Containing gridded_params in a Python friendly format. New
    entries calculated as needed, and entries not used in the DBPM gridded model were removed.
    '''

    gridded_python = {'timesteps_years': gridded_params['timesteps_years'][0],
                      'numb_time_steps': gridded_params['numb_time_steps'][0],
                      'effort': gridded_params['effort'],
                      'fish_mort_pred': gridded_params['fish_mort_pred'][0],
                      'fish_mort_detritivore': gridded_params['fish_mort_detritivore'][0],
                      'hr_volume_search': gridded_params['hr_volume_search'][0],
                      'detritus_coupling': gridded_params['detritus_coupling'][0],
                      'log10_pred_prey_ratio': gridded_params['log10_pred_prey_ratio'][0],
                      'log_prey_pref': gridded_params['log_prey_pref'][0],
                      'hr_vol_filter_benthos': gridded_params['hr_vol_filter_benthos'][0],
                      'metabolic_req_pred': gridded_params['metabolic_req_pred'][0],
                      'metabolic_req_detritivore': gridded_params['metabolic_req_detritivore'][0],
                      'defecate_prop': gridded_params['defecate_prop'][0],
                      'growth_prop': 1-(gridded_params['defecate_prop'][0]),
                      'def_low': gridded_params['def_low'][0],
                      'high_prop': 1-(gridded_params['def_low'][0]),
                      'growth_pred': gridded_params['growth_pred'][0],
                      'growth_detritivore': gridded_params['growth_detritivore'][0],
                      'growth_detritus': gridded_params['growth_detritus'][0],
                      'energy_pred': gridded_params['energy_pred'][0],
                      'energy_detritivore': gridded_params['energy_detritivore'][0],
                      'handling': gridded_params['handling'][0],
                      'dynamic_reproduction': gridded_params['dynamic_reproduction'][0],
                      'c1': gridded_params['c1'][0],
                      'activation_energy': gridded_params['activation_energy'][0],
                      'boltzmann': gridded_params['boltzmann'][0],
                      'natural_mort': gridded_params['natural_mort'][0],
                      'size_senescence': gridded_params['size_senescence'][0],
                      'exp_senescence_mort': gridded_params['exp_senescence_mort'][0],
                      'const_senescence_mort': gridded_params['const_senescence_mort'][0],
                      'log_size_increase': gridded_params['log_size_increase'][0],
                      'log10_size_bins': gridded_params['log10_size_bins'],
                      'numb_size_bins':  gridded_params['numb_size_bins'][0],
                      'ind_min_pred_size': (gridded_params['ind_min_pred_size'][0])-1,
                      'ind_min_detritivore_size': (gridded_params['ind_min_detritivore_size'][0])-1,
                      'idx_new': list(range(gridded_params['ind_min_detritivore_size'][0],
                                            gridded_params['numb_size_bins'][0])),
                      'ind_min_fish_pred': int(gridded_params['ind_min_fish_pred'][0]-1),
                      'ind_min_fish_det': int(gridded_params['ind_min_fish_det'][0]-1),
                      'idx': (np.array(gridded_params['idx'])-1).tolist(),
                      'init_pred': gridded_params['init_pred'],
                      'init_detritivores': gridded_params['init_detritivores'],
                      'init_detritus': gridded_params['init_detritus'][0]
                     }
    
    return gridded_python


### DPBM functions ----
# Build a lookup table for diet preference. Looks at all combinations of predator 
# and prey body size: diet preference (in the predator spectrum only)
def phi_f(q, log10_pred_prey_ratio, log_prey_pref):
    phi = np.where(q > 0, 
                   np.exp(-(q-log10_pred_prey_ratio)*(q-log10_pred_prey_ratio)/
                          (2*log_prey_pref**2))/(log_prey_pref*np.sqrt(2.0*np.pi)),
                   0) 
    return phi


# Function to build lookup tables for (constant) growth
# Considers components which remain constant
def gphi_f(pred_prey_matrix, log10_pred_prey_ratio, log_prey_pref):
    gphi = 10**(-pred_prey_matrix)*phi_f(pred_prey_matrix, log10_pred_prey_ratio, 
                                         log_prey_pref)
    return gphi


# Function to build lookup tables for (constant) mortality
def mphi_f(rev_pred_prey_matrix, log10_pred_prey_ratio, log_prey_pref, 
           metabolic_req_pred):
    mphi = (10**(metabolic_req_pred*rev_pred_prey_matrix)*
            phi_f(rev_pred_prey_matrix, log10_pred_prey_ratio, log_prey_pref))
    return mphi


# Function to build lookup table for components of 10^(alpha*x)
def expax_f(log10_size_bins, metabolic_req_pred):
    expax = 10**(log10_size_bins*metabolic_req_pred)
    return expax


# Create predator-prey combination matrix
# The square matrix holds the log(predatorsize/preysize) for all combinations of
# predator-prey sizes
def pred_prey_matrix(log10_size_bins):
    '''
    Inputs:
    - log10_size_bins (numpy array) Containing the log of all sizes of predator and
    prey items
    
    Outputs:
    ppm (numpy array) Contains all combinations of predator-prey sizes
    '''
    
    ppm = np.empty([len(log10_size_bins), len(log10_size_bins)])
    for i, j in enumerate(log10_size_bins):
        ppm[:, i] = log10_size_bins[i] - log10_size_bins
    return ppm


# Initialising matrices using size classes and time as dimensions
def init_da(log10_size_bins, time):
    '''
    Inputs:
    - log10_size_bins (numpy array) Containing the log of all sizes of predator and
    prey items
    - time (numpy array) Containing dates to be included in final data array
    
    Outputs:
    da (data array) Containing zeroes with size classes and time as dimensions
    '''
    data_start = np.zeros((len(log10_size_bins), len(time)))
    da = xr.DataArray(data = data_start, dims = ['size_class', 'time'], 
                      coords = {'size_class': log10_size_bins, 'time': time})
    return da


# Gravity model ----
# Redistribute total effort across grid cells according to proportion of biomass in that 
# grid cell using gravity model, Walters & Bonfil, 1999, Gelchu & Pauly 2007 ideal free
# distribution - Blanchard et al 2008
def gravitymodel(effort, prop_b, depth, n_iter, sea_ice_mask):
    '''
    Inputs:
    - effort (Data array) Fishing effort for a single time step
    - prop_b (Data array) Proportion of total fishable biomass for each grid cell at a 
    single time step
    - depth (Data array) Bathymetry of the area of interest
    - n_iter (integer) Number of iterations needed to redistribute fishing effort
    - sea_ice_mask (Data array) Mask describing sea ice covered areas inaccessible to 
    fishers
    
    Outputs:
    eff (data array) Containing redistributed fishing effort
    '''

    eff = effort
    d = (1-(depth/depth.max()))

    #Initialise loop
    i = 1
    while(i <= n_iter):
        suit = prop_b*d*sea_ice_mask
        rel_suit = suit/(suit.sum())
        neweffort = eff+(rel_suit*eff)
        mult = (eff.sum())/(neweffort.sum())
        eff = neweffort*mult
        i += i

    #Ensuring NA values are transformed to zeroes
    eff = xr.where(np.isfinite(eff), eff, 0)

    return eff


# Calculate fishing mortality and effort ------
def effort_calculation(predators, detritivores, effort, depth, size_bin_vals, 
                       gridded_params, sea_ice_mask):
    '''
    Inputs:
    - predators (2D data array) Pelagic predator biomass
    - detritivores (2D data array) Bethic detritivore biomass
    - effort (2D data array) Fishing effort
    - depth (2D data array) Bathymetry of the area of interest
    - size_bins_vals (1D data array) Size classes in grams
    - gridded_params (dictionary) DBPM parameters
    - sea_ice_mask (Data array) Mask describing sea ice covered areas inaccessible 
    to fishers

    Outputs:
    - new_effort (2D data array) Fishing effort calculated for next time step
    '''
    if effort.sum() != 0:
        pred_bio = ((predators*size_bin_vals*gridded_params['log_size_increase']).
                    isel(size_class = slice(gridded_params['ind_min_fish_pred'],
                                            -1))).sum('size_class')
    
        det_bio = ((detritivores*size_bin_vals*gridded_params['log_size_increase']).
                   isel(size_class = slice(gridded_params['ind_min_fish_det'],
                                           -1))).sum('size_class')
    
        sum_bio = pred_bio+det_bio
    
        prop_b = sum_bio/sum_bio.sum()
    
        #Calculate new effort
        new_effort = gravitymodel(effort, prop_b, depth, 1, sea_ice_mask)
        new_effort = xr.where(np.isfinite(depth), new_effort, np.nan)
    else:
        new_effort = effort

    # Adjusting time stamp
    new_effort['time'] = [effort.time.values]
    try:
        new_effort = new_effort.transpose('time', 'lat', 'lon')
    except:
        pass
    #Adding name
    new_effort.name = 'effort'

    #Return effort
    return new_effort


# Loading fixed gridded DBPM inputs ------
def loading_dbpm_fixed_inputs(base_folder):
    '''
    Inputs:
    - base_folder (character) Full path to folder where fixed DBPM gridded inputs
    are stored

    Outputs:
    - dbpm_inputs (xarray Dataset) Contains fixed gridded inputs needed to run DBPM
    '''
    #Loading data from base folder
    #Depth
    depth = xr.open_zarr(glob(os.path.join(
        base_folder, '*obsclim_deptho_*'))[0])['deptho']
    
    #Preference benthic
    pref_benthos = xr.open_zarr(glob(
        os.path.join(base_folder, 
                     'pref-benthos_*'))[0])['pref_benthos']
    #Preference pelagic
    pref_pelagic = xr.open_zarr(glob(
        os.path.join(base_folder, 
                     'pref-pelagic_*'))[0])['pref_pelagic']
    #Constant growth
    constant_growth = xr.open_zarr(glob(
        os.path.join(base_folder, 
                     'const-growth_*'))[0])['constant_growth']
    #Constant mortality
    constant_mortality = xr.open_zarr(glob(
        os.path.join(base_folder, 
                     'const-mort_*'))[0])['constant_mortality']
    #Consumption pelagic
    consume_pelagic = xr.open_zarr(glob(
        os.path.join(base_folder, 
                     'consume-pelagic_*'))[0])['consume_pelagic']
    #Consumption benthos
    consume_benthos = xr.open_zarr(glob(
        os.path.join(base_folder,
                     'consume-benthos_*'))[0])['consume_benthos']
    #Predator mortality from fishing
    fish_mort_pred = xr.open_zarr(glob(
        os.path.join(base_folder,
                     'fish-mort-pred_*'))[0])['fish_mort_pred']
    #Detritivores mortality from fishing
    fish_mort_det = xr.open_zarr(glob(
        os.path.join(base_folder,
                     'fish-mort-det_*'))[0])['fish_mort_det']
    #Predator mortality from others
    other_mort_pred = xr.open_zarr(glob(
        os.path.join(base_folder,
                     'other-mort-pred_*'))[0])['other_mort_pred']
    #Detritivore mortality from others
    other_mort_det = xr.open_zarr(glob(
        os.path.join(base_folder,
                     'other-mort-det_*'))[0])['other_mort_det']
    #Predator mortality from senescence
    senes_mort_pred = xr.open_zarr(glob(
        os.path.join(base_folder,
                     'senes-mort-pred_*'))[0])['senes_mort_pred']
    #Detritivore mortality from senescence
    senes_mort_det = xr.open_zarr(glob(
        os.path.join(base_folder,
                     'senes-mort-det_*'))[0])['senes_mort_det']

    dbpm_inputs = xr.Dataset(data_vars = {'depth': depth,
                                          'pref_benthos': pref_benthos,
                                          'pref_pelagic': pref_pelagic,
                                          'constant_growth': constant_growth,
                                          'constant_mortality': constant_mortality,
                                          'consume_pelagic': consume_pelagic,
                                          'consume_benthos': consume_benthos,
                                          'fish_mort_pred': fish_mort_pred,
                                          'fish_mort_det': fish_mort_det,
                                          'other_mort_pred': other_mort_pred,
                                          'other_mort_det': other_mort_det,
                                          'senes_mort_pred': senes_mort_pred,
                                          'senes_mort_det': senes_mort_det,
                                          'mask': np.isfinite(depth)})

    return dbpm_inputs


# Loading initialising values for gridded DBPM biomass ------
def loading_dbpm_biomass_inputs(data_folder, init_time = None):
    '''
    Inputs:
    - data_folder (character) Full path to folder where initialising DBPM gridded inputs
    are stored
    - init_time (character) Default is None. Year and month from when to restart gridded 
    DBPM runs. If set to None, DBPM is run from the beginning

    Outputs:
    - ds_init (xarray Dataset) Contains initialising gridded inputs needed to run DBPM
    '''
    if init_time is None:
        predators = xr.open_zarr(
            glob(os.path.join(data_folder, 'predators_*'))[0])['predators'] 
        detritivores = xr.open_zarr(
            glob(os.path.join(data_folder, 'detritivores_*'))[0])['detritivores']
        detritus = xr.open_zarr(
            glob(os.path.join(data_folder, 'detritus_*'))[0])['detritus']
    #If 'init_time' is provided, model runs from 'init_time'
    else:
        predators = xr.open_dataarray(
            glob(os.path.join(data_folder, f'predators_*_{init_time}.nc'))[0])
        detritivores = xr.open_dataarray(
            glob(os.path.join(data_folder, f'detritivores_*_{init_time}.nc'))[0])
        detritus = xr.open_dataarray(
            glob(os.path.join(data_folder, f'detritus_*_{init_time}.nc'))[0])
    
    #Create dataset for predator, detritivores and detritus initialisation data
    ds_init = xr.Dataset(data_vars = {'predators': predators, 
                                      'detritivores': detritivores, 
                                      'detritus': detritus})

    return ds_init


# Loading initialising values for gridded DBPM biomass ------
def loading_dbpm_dynamic_inputs(gridded_folder, init_time = None, fishing = True):
    '''
    Inputs:
    - gridded_esm (character) Full path to folder where processed DBPM inputs are stored
    - init_time (character) Default is None. Year and month from when to restart gridded 
    DBPM runs. If set to None, DBPM is run from the beginning
    - fishing (boolean) Default is True. It indicates whether or not fishing effort 
    should be included in the dynamic inputs dataset
    
    Outputs:
    - ds_dynamic (xarray Dataset) Contains dynamic gridded inputs needed to run DBPM
    '''

    #Load data
    [ui0_search] = glob(os.path.join(gridded_folder, 'ui0_spinup_obsclim*'))
    ui0 = xr.open_zarr(ui0_search, chunks = {'time': -1, 'lon': -1, 
                                             'lat': -1})['ui0']

    slope_search = [f for f in glob(os.path.join(gridded_folder, '*_slope_*')) 
                    if 'ctrlclim' not in f]
    slope = xr.open_mfdataset(
        slope_search, engine = 'zarr', 
        parallel = True)['slope'].chunk({'time': -1, 'lon': -1, 'lat': -1})
    
    [pel_temp_search] = glob(os.path.join(gridded_folder,
                                          'pel-temp-eff_spinup_obsclim*'))
    pel_tempeffect = xr.open_zarr(
        pel_temp_search, chunks = {'time': -1, 'lon': -1, 
                                   'lat': -1})['pel_temp_eff']
    
    [ben_temp_search] = glob(os.path.join(gridded_folder, 
                                          'ben-temp-eff_spinup_obsclim*'))
    ben_tempeffect = xr.open_zarr(
        ben_temp_search, chunks = {'time': -1, 'lon': -1, 
                                   'lat': -1})['ben_temp_eff']

    sinking_rate_search = [f for f in glob(os.path.join(gridded_folder, '*_er_*')) 
                           if 'ctrlclim' not in f]
    sinking_rate = xr.open_mfdataset(
        sinking_rate_search, engine = 'zarr', 
        parallel = True)['export_ratio'].chunk({'time': -1, 'lon': -1, 'lat': -1})

    if fishing:
        [effort_search] = glob(os.path.join(gridded_folder, 'effort_spinup_obsclim*'))
        effort = xr.open_zarr(
            effort_search, chunks = {'time': -1, 'lon': -1, 'lat': -1})['effort']
        simask_search = [f for f in glob(os.path.join(gridded_folder, '*_simask_*')) 
                    if 'ctrlclim' not in f]
        simask = xr.open_mfdataset(
            simask_search, engine = 'zarr', 
            parallel = True)['simask'].chunk({'time': -1, 'lon': -1, 'lat': -1})
    else:
        effort = xr.where(np.isfinite(ui0), 0, np.nan)
        effort.name = 'effort'
        simask = xr.where(np.isfinite(ui0), 0, np.nan)
        simask.name = 'simask'

    #Subset data
    if init_time is not None:
        #Timestep from when to restart DBPM 
        subset_time = (pd.Timestamp(init_time)+
                       pd.DateOffset(months = 1)).strftime('%Y-%m')
        
        #Subset data from timestep above until the end of the available data
        ui0 = ui0.sel(time = slice(subset_time, None))
        slope = slope.sel(time = slice(subset_time, None))
        pel_tempeffect = pel_tempeffect.sel(time = slice(subset_time, None))
        ben_tempeffect = ben_tempeffect.sel(time = slice(subset_time, None))
        sinking_rate = sinking_rate.sel(time = slice(subset_time, None))
        effort = effort.sel(time = slice(subset_time, None))
        simask = simask.sel(time = slice(subset_time, None))

    ds_dynamic = xr.Dataset(data_vars = {'ui0': ui0, 'slope': slope,
                                         'pel_tempeffect': pel_tempeffect,
                                         'ben_tempeffect': ben_tempeffect, 
                                         'sinking_rate': sinking_rate,
                                         'effort': effort,
                                         'sea_ice_mask': simask})

    return ds_dynamic


# Calculating feeding and satiation rates ------
def feeding_satiation_rates(gridded_params, dbpm_fixed_inputs, dbpm_init_inputs, 
                            dbpm_dynamic_inputs, group):
    '''
    Inputs:
    - gridded_params (dictionary) Gridded DBPM parameters obtained in step 04.
    - dbpm_fixed_inputs (xarray Dataset) Contains fixed gridded inputs needed to run 
    DBPM
    - dbpm_init_inputs (xarray Dataset) Contains gridded inputs with initialisation 
    values for predators, detritivores and detritus
    - dbpm_dynamic_inputs (xarray Dataset) Contains dynamic gridded inputs need to run
    DBPM
    - group (character): Choose the group for which feeding and satiation rates will 
    be calculated for. Choices: predators, detritivores, detritus.

    Outputs:
    - feed_sat_rates (xarray Dataset). Contains feeding and satiation rates for chosen
    group
    '''
    
    if group == 'predators':
        #To be applied to feeding rates for pelagics and benthic groups
        feed_mult_pel = (gridded_params['hr_volume_search']*
                         (10**(dbpm_fixed_inputs['log10_size_bins']*
                               gridded_params['metabolic_req_pred']))*
                         dbpm_fixed_inputs['pref_pelagic'])
        # Predators (f_pel)
        pred_growth = (((dbpm_init_inputs['predators']*
                         gridded_params['log_size_increase']).
                        dot(dbpm_fixed_inputs['constant_growth']).
                       rename({'sc': 'size_class'}))*feed_mult_pel)
        feed_rate_pel = (dbpm_dynamic_inputs['pel_tempeffect']*
                         (pred_growth/(1+gridded_params['handling']*pred_growth)))
        # Satiation level of predator for pelagic prey (sat.pel)
        sat_pel = xr.where(feed_rate_pel > 0, feed_rate_pel/pred_growth, 0)

        #Create dataset with feeding and satiation outputs for predators
        feed_sat_rates = xr.Dataset(data_vars = {'f_pel': feed_rate_pel,
                                                 'sat_pel': sat_pel})
    
    elif group == 'detritivores':
        #To be applied to feeding rates for benthic group
        feed_mult_ben = (gridded_params['hr_volume_search']*
                         (10**(dbpm_fixed_inputs['log10_size_bins']*
                               gridded_params['metabolic_req_pred']))*
                        dbpm_fixed_inputs['pref_benthos'])
        # Detritivores (f_ben)
        detrit_growth = (((dbpm_init_inputs['detritivores']*
                           gridded_params['log_size_increase']).
                          dot(dbpm_fixed_inputs['constant_growth']).
                          rename({'sc': 'size_class'}))*feed_mult_ben)
        feed_rate_bent = (dbpm_dynamic_inputs['pel_tempeffect']*
                          (detrit_growth/(1+gridded_params['handling']*detrit_growth)))
        # Satiation level of predator for benthic prey
        calc_growth = ((gridded_params['hr_volume_search']*
                        (10**(dbpm_fixed_inputs['log10_size_bins']*
                              gridded_params['metabolic_req_detritivore']))*
                        dbpm_fixed_inputs['pref_benthos'])*
                       ((dbpm_init_inputs['detritivores']*
                         gridded_params['log_size_increase']).
                        dot(dbpm_fixed_inputs['constant_growth'])).
                       rename({'sc': 'size_class'}))
        #Satiation level of predator for benthic prey (sat.ben)
        sat_ben = xr.where(feed_rate_bent > 0, feed_rate_bent/calc_growth, 0)

        #Create dataset with feeding and satiation outputs for detritivores
        feed_sat_rates = xr.Dataset(data_vars = {'f_ben': feed_rate_bent,
                                                 'sat_ben': sat_ben})

    elif group == 'detritus':
        # Feeding rates
        detritus_multiplier = ((1/dbpm_fixed_inputs['size_bin_vals'])*
                               gridded_params['hr_vol_filter_benthos']*
                               10**(dbpm_fixed_inputs['log10_size_bins']*
                                    gridded_params['metabolic_req_detritivore'])*
                               dbpm_init_inputs['detritus'])
        feed_rate_det = (dbpm_dynamic_inputs['ben_tempeffect']*detritus_multiplier/
                         (1+gridded_params['handling']*detritus_multiplier))
        
        #Create dataset with feeding and satiation outputs
        feed_sat_rates = xr.Dataset(data_vars = {'f_det': feed_rate_det})

    else:
        print('the "group" parameter must be "predators", "detritivores", or "detritus"')

    #Reorganise output dimensions
    try:
        feed_sat_rates = feed_sat_rates.transpose('time', 'size_class', 'lat', 'lon')
    except:
        feed_sat_rates = feed_sat_rates.transpose('time', 'size_class')
    #Apply spatial mask
    feed_sat_rates = feed_sat_rates.where(dbpm_fixed_inputs['mask'])
    #Ensure correct time is applied
    feed_sat_rates['time'] = dbpm_init_inputs.time
    
    return feed_sat_rates


def mortality_calc(gridded_params, dbpm_fixed_inputs, dbpm_init_inputs, 
                   dbpm_dynamic_inputs, feed_sat_rates, group):
    '''
    Inputs:
    - gridded_params (dictionary) Gridded DBPM parameters obtained in step 04.
    - dbpm_fixed_inputs (xarray Dataset) Contains fixed gridded inputs needed to run 
    DBPM
    - dbpm_init_inputs (xarray Dataset) Contains gridded inputs with initialisation 
    values for predators, detritivores and detritus
    - dbpm_dynamic_inputs (xarray Dataset) Contains dynamic gridded inputs need to run
    DBPM
    - feed_sat_rates (xarray Dataset) Contains feeding and satiation rates for the 
    group for which mortality rates will be calculated. This is the output of the
    "feeding_satiation_rates" function
    - group (character): Choose the group for which mortality rates will be calculated.
    Choices: predators, detritivores.

    Outputs:
    - mortality_terms (xarray Dataset). Contains mortality terms for all groups modelled
    '''

    if group == 'predators': 
        # Fishing mortality (FVec.u, FVec.v) from Benoit & Rochet 2004. Here fish_mort_pred 
        # and fish_mort_pred= fixed catchability term for predators and detritivores to be 
        # estimated along with ind_min_det and ind_min_fish_pred
        fishing_mort_pred = (dbpm_fixed_inputs['fish_mort_pred']*
                             dbpm_dynamic_inputs['effort']).drop_vars('time')
        #Time dimension matches feed rates time stamp
        fishing_mort_pred = fishing_mort_pred.expand_dims({'time': dbpm_init_inputs.time})
        # Predator death integrals
        # Predation mortality
        # Predators (PM.u)
        pred_mort_pred = (dbpm_fixed_inputs['consume_pelagic']*
                          ((dbpm_init_inputs['predators']*feed_sat_rates['sat_pel']*
                             gridded_params['log_size_increase']).
                           dot(dbpm_fixed_inputs['constant_mortality'])).
                          rename({'sc': 'size_class'}))
        # Total mortality
        # Predators (Z.u)
        tot_mort_pred = (pred_mort_pred+dbpm_dynamic_inputs['pel_tempeffect']*
                         dbpm_fixed_inputs['other_mort_pred']+
                         dbpm_fixed_inputs['senes_mort_pred']+fishing_mort_pred)
        
        # Create dataset with mortality for predators
        mortality_terms = xr.Dataset(data_vars = {'Fvec_u': fishing_mort_pred,
                                                  'PM_u': pred_mort_pred,
                                                  'Z_u': tot_mort_pred})
    elif group == 'detritivores':
        # Fishing mortality (FVec.u, FVec.v) from Benoit & Rochet 2004. Here 
        # fish_mort_pred and fish_mort_pred= fixed catchability term for predators and 
        # detritivores to be estimated along with ind_min_det and ind_min_fish_pred
        fishing_mort_det = (dbpm_fixed_inputs['fish_mort_det']*
                            dbpm_dynamic_inputs['effort']).drop_vars('time')
        #Time dimension matches feed rates time stamp
        fishing_mort_det = fishing_mort_det.expand_dims({'time': dbpm_init_inputs.time})
        # Predation mortality
        # Detritivores (PM.v)
        pred_mort_det = xr.where(feed_sat_rates['sat_ben'] > 0, 
                                 (dbpm_fixed_inputs['consume_benthos']*
                                  ((dbpm_init_inputs['predators']*
                                    feed_sat_rates['sat_ben']*
                                    gridded_params['log_size_increase']).
                                   dot(dbpm_fixed_inputs['constant_mortality'])).
                                  rename({'sc': 'size_class'})), 0)
        # Total mortality
        # Detritivores (Z.v)
        tot_mort_det = (pred_mort_det+dbpm_dynamic_inputs['ben_tempeffect']*
                        dbpm_fixed_inputs['other_mort_det']+
                        dbpm_fixed_inputs['senes_mort_det']+fishing_mort_det)

        # Create dataset with mortality for detritivores
        mortality_terms = xr.Dataset(data_vars = {'Fvec_v': fishing_mort_det,
                                                  'PM_v': pred_mort_det,
                                                  'Z_v': tot_mort_det})
    else:
        print('the "group" parameter must be either "predators" or "detritivores"')
        
    #Reorganise dimensions
    try:
        mortality_terms = mortality_terms.transpose('time', 'size_class', 'lat', 'lon')
    except:
        mortality_terms = mortality_terms.transpose('time', 'size_class')
    #Apply spatial mask
    mortality_terms = mortality_terms.where(dbpm_fixed_inputs['mask'])
    #Ensure correct time is applied
    mortality_terms['time'] = dbpm_init_inputs.time
    
    return mortality_terms


def growth_rates_calc(gridded_params, feed_sat_rates, dbpm_fixed_inputs):
    '''
    Inputs:
    - gridded_params (dictionary) Gridded DBPM parameters obtained in step 04.
    - feed_sat_rates (xarray Dataset) Contains feeding rates calculated for predators,
    detritivores and detritus calculated with the "feeding_satiation_rates" function
    - dbpm_fixed_inputs (xarray Dataset) Contains fixed gridded inputs needed to run 
    DBPM
    
    Outputs:
    - growth_rates (xarray Dataset). Contains growth and estimates for detritivores 
    and predators
    '''
    # Predator growth integral (GG_u)
    growth_int_pred = (gridded_params['growth_prop']*gridded_params['growth_pred']*
                       feed_sat_rates['f_pel']+gridded_params['high_prop']*
                       gridded_params['growth_detritivore']*feed_sat_rates['f_ben'])
    # Benthos growth integral yr-1 (GG_v)
    growth_int_det = (gridded_params['high_prop']*gridded_params['growth_detritus']*
                      feed_sat_rates['f_det'])
    # Create dataset with growth rates for detritivores and predators
    growth_reprod = xr.Dataset(data_vars = {'GG_u': growth_int_pred,
                                            'GG_v': growth_int_det})

    #Reorganise dimensions
    try:
        growth_reprod = growth_reprod.transpose('time', 'size_class', 'lat', 'lon')
    except:
        growth_reprod = growth_reprod.transpose('time', 'size_class')
    #Apply spatial mask
    growth_reprod = growth_reprod.where(dbpm_fixed_inputs['mask'])
    #Ensure correct time is applied
    growth_reprod['time'] = feed_sat_rates.time
    
    return growth_reprod


def reproduction_rates(gridded_params, feed_sat_rates, dbpm_fixed_inputs, group):
    '''
    Inputs:
    - gridded_params (dictionary) Gridded DBPM parameters obtained in step 04.
    - feed_sat_rates (xarray Dataset) Contains feeding rates calculated for predators,
    detritivores and detritus calculated with the "feeding_satiation_rates" function
    - dbpm_fixed_inputs (xarray Dataset) Contains fixed gridded inputs needed to run 
    DBPM
    - group (character): Choose the group for which mortality rates will be 
    calculated. Choices: predators, detritivores.
    
    Outputs:
    - growth_rates (xarray Dataset). Contains growth and estimates for detritivores 
    and predators
    '''

    if group == 'predators':
        #Keep relevant size classes
        feed_sat_rates = (feed_sat_rates.
                          isel(size_class = 
                               slice(gridded_params['ind_min_pred_size']+1, None)))
        #Predators
        reprod_pred = (gridded_params['growth_prop']*
                       (1-(gridded_params['growth_pred']+
                           gridded_params['energy_pred']))*
                       feed_sat_rates['f_pel']+gridded_params['growth_prop']*
                       (1-(gridded_params['growth_detritivore']+
                           gridded_params['energy_detritivore']))*
                       feed_sat_rates['f_ben'])
        #Adding to output dataset
        reprod_rates = xr.Dataset(data_vars = {'R_u': reprod_pred})
        
    elif group == 'detritivores':
        #Keep relevant size classes
        feed_sat_rates = (feed_sat_rates.
                          isel(size_class = 
                               slice(gridded_params['ind_min_detritivore_size']+1, 
                                     None)))
        #Detritivores
        reprod_det = (gridded_params['high_prop']*
                      (1-(gridded_params['growth_detritus']+
                          gridded_params['energy_detritivore']))*
                      feed_sat_rates['f_det'])
        #Adding to output dataset
        reprod_rates = xr.Dataset(data_vars = {'R_v': reprod_det})

    else:
        print('the "group" parameter must be either "predators" or "detritivores"')

    #Reorganise dimensions
    try:
        reprod_rates = reprod_rates.transpose('time', 'size_class', 'lat', 'lon')
    except:
        reprod_rates = reprod_rates.transpose('time', 'size_class')
    #Apply spatial mask
    reprod_rates = reprod_rates.where(dbpm_fixed_inputs['mask'])
    #Ensure correct time is applied
    reprod_rates['time'] = feed_sat_rates.time
    
    return reprod_rates


def detritus_output(gridded_params, dbpm_fixed_inputs, dbpm_init_inputs, 
                    feed_detritus):
    '''
    Inputs:
    - gridded_params (dictionary) Gridded DBPM parameters obtained in step 04.
    - dbpm_fixed_inputs (xarray Dataset) Contains fixed gridded inputs needed to run 
    DBPM
    - dbpm_init_inputs (xarray Dataset) Contains gridded inputs with initialisation 
    values for predators, detritivores and detritus
    - feed_detritus (xarray Data Array) Contains detritus rates calculated with the 
    "feeding_satiation_rates" function
    
    Outputs:
    - output_w (xarray Data Array). Contains detritus output (g.m-2.yr-1), which are
    defined as losses from detritivore scavenging/filtering only.
    '''
    # Detritus output (g.m-2.yr-1)
    # losses from detritivore scavenging/filtering only:
    output_w = (dbpm_fixed_inputs['size_bin_vals']*
                gridded_params['log_size_increase']*
                dbpm_init_inputs['detritivores']*
                feed_detritus).sum('size_class')

    # Adding variable name
    output_w.name = 'output_w'
    #Reorganise dimensions
    try:
        output_w = output_w.transpose('time', 'lat', 'lon')
    except:
        pass
    #Apply spatial mask
    output_w = output_w.where(dbpm_fixed_inputs['mask'])
    
    return output_w
    

def defecation_predators(gridded_params, dbpm_fixed_inputs, dbpm_init_inputs, 
                         feed_sat_rates):
    '''
    Inputs:
    - gridded_params (dictionary) Gridded DBPM parameters obtained in step 04.
    - dbpm_fixed_inputs (xarray Dataset) Contains fixed gridded inputs needed to run 
    DBPM
    - dbpm_init_inputs (xarray Dataset) Contains gridded inputs with initialisation 
    values for predators, detritivores and detritus
    - feed_sat_rates (xarray Dataset) Contains feeding rates for predators and 
    detritivores calculated with the "feeding_satiation_rates" function
    
    Outputs:
    - defbypred (xarray Data Array). Contains total biomass defecated by predators
    '''
    # Select relevant predator size classes
    pred = (dbpm_init_inputs['predators'].
            isel(size_class = slice(gridded_params['ind_min_pred_size'], None)))

    # Total biomass density defecated by pred (g.m-2.yr-1)
    defbypred = ((gridded_params['defecate_prop']*feed_sat_rates['f_pel']*
                  dbpm_fixed_inputs['size_bin_vals']*
                  pred+gridded_params['def_low']*feed_sat_rates['f_ben']*
                  dbpm_fixed_inputs['size_bin_vals']*pred)*
                 gridded_params['log_size_increase']).sum('size_class')  

    # Adding variable name
    defbypred.name = 'defbypred'
    #Reorganise dimensions
    try:
        defbypred = defbypred.transpose('time', 'lat', 'lon')
    except:
        pass
    #Apply spatial mask
    defbypred = defbypred.where(dbpm_fixed_inputs['mask'])
    
    return defbypred


def detritus_pool(gridded_params, dbpm_fixed_inputs, dbpm_init_inputs, 
                  dbpm_dynamic_inputs, defbypred, output_w):
    '''
    Inputs:
    - gridded_params (dictionary) Gridded DBPM parameters obtained in step 04.
    - dbpm_fixed_inputs (xarray Dataset) Contains fixed gridded inputs needed to 
    run DBPM
    - dbpm_init_inputs (xarray Dataset) Contains gridded inputs with initialisation 
    values for predators, detritivores and detritus
    - dbpm_dynamic_inputs (xarray Dataset) Contains dynamic gridded inputs need to
    run DBPM
     - defbypred (xarray Data Array). Contains total biomass defecated by predators 
     calculated with the "defecation_predators" function
     - output_w (xarray Data Array). Contains detritus output (g.m-2.yr-1) 
     calculated with the "detritus_output" function
    
    Outputs:
    - det_pool (xarray Dataset). Contains detritus biomass density pool fluxes
    (g.m-2.yr-1)
    '''

    # Considering pelagic faeces as input as well as dead bodies 
    # from both pelagic and benthic communities and phytodetritus
    # (dying sinking phytoplankton)
    if gridded_params['detritus_coupling']:
        # Pelagic spectrum inputs (sinking dead bodies and faeces): 
        # Export ratio used for "sinking rate" + benthic spectrum 
        # inputs (dead stuff already on/in seafloor)
        input_w = (dbpm_dynamic_inputs['sinking_rate']*
                   (defbypred+
                    ((dbpm_fixed_inputs['other_mort_pred']*
                      dbpm_dynamic_inputs['pel_tempeffect']*
                      dbpm_init_inputs['predators']*
                      dbpm_fixed_inputs['size_bin_vals']*
                      gridded_params['log_size_increase'])+
                     (dbpm_fixed_inputs['senes_mort_pred']*
                      dbpm_dynamic_inputs['pel_tempeffect']*
                      dbpm_init_inputs['predators']*
                      dbpm_fixed_inputs['size_bin_vals']*
                      gridded_params['log_size_increase'])+
                     (dbpm_fixed_inputs['other_mort_det']*
                      dbpm_dynamic_inputs['ben_tempeffect']*
                      dbpm_fixed_inputs['size_bin_vals']*
                      gridded_params['log_size_increase']*
                      dbpm_init_inputs['detritivores'])+
                    (dbpm_fixed_inputs['senes_mort_det']*
                      dbpm_dynamic_inputs['ben_tempeffect']*
                      dbpm_fixed_inputs['size_bin_vals']*
                      gridded_params['log_size_increase']*
                      dbpm_init_inputs['detritivores']))
                    .sum('size_class')))
    else:
        input_w = (((dbpm_fixed_inputs['other_mort_det']*
                     dbpm_dynamic_inputs['ben_tempeffect']*
                     dbpm_fixed_inputs['size_bin_vals']*
                     gridded_params['log_size_increase']*
                     dbpm_init_inputs['detritivores'])).
                   sum('size_class'))
   
   # Get burial rate from Dunne et al. 2007 equation 3
    burial = input_w*(0.013+0.53*input_w**2/(7+input_w)**2)
    output_w = output_w+burial
    
    # Losses from detritivory + burial rate (not including 
    # remineralisation because that goes to p.p. after sediment, 
    # we are using realised p.p. as inputs to the model) 
    dW = input_w-output_w
   
    # Create dataset with changes in detritus biomass
    det_pool = xr.Dataset(data_vars = {'input_w': input_w, 
                                       'dW': dW})
    #Reorganise dimensions
    try:
        det_pool = det_pool.transpose('time', 'lat', 'lon')
    except:
        pass
    #Apply spatial mask
    det_pool = det_pool.where(dbpm_fixed_inputs['mask'])
    
    return det_pool


def biomass_density(gridded_params, dbpm_fixed_inputs, growth_rate, 
                     total_mortality, biomass, group):
    '''
    Inputs:
    - gridded_params (dictionary) Gridded DBPM parameters obtained in step 04.
    - dbpm_fixed_inputs (xarray Dataset) Contains fixed gridded inputs needed to 
    run DBPM
    - growth_rates (xarray Dataset). Contains growth and estimates for detritivores 
    and predators estimated by the "growth_rate" function
    - total_mortality (xarray Data Array). Contains total mortality for group of
    interest estimated by the "mortality_calc" function
    - biomass (xarray Data Array). Contains total biomass for group of interest
    - group (character): Choose the group for which mortality rates will be 
    calculated. Choices: predators, detritivores.
    
    Outputs:
    - density (xarray Dataset). Contains biomass density (nos.m-2) for the group
    of interest
    '''
    
    if group == 'predators':
        index_sc = np.array(gridded_params['idx'])
        ref_sc =  gridded_params['log10_ind_min_pred_size']
    elif group == 'detritivores':
        index_sc = np.array(gridded_params['idx_new'])
        ref_sc =  gridded_params['log10_ind_min_detritivore_size']
    else:
        print('the "group" parameter must be either "predators" or ' + 
              '"detritivores"')
    
    # Pelagic Predator Density (nos.m-2) ----
    # Solve for time + delta.t using implicit time Euler upwind finite 
    # difference (help from Ken Andersen and Richard Law)

    # Shifting size classes by one before calculating Ai
    ggp_shift = growth_rate.shift({'size_class': 1})
    ggp_shift, range_sc = xr.align(ggp_shift, 
                                   dbpm_fixed_inputs['log10_size_bins'].
                                   isel(size_class = index_sc))
      
    # Matrix setup for implicit differencing 
    Ai = (((1/np.log(10))*-ggp_shift*gridded_params['timesteps_years']/
           gridded_params['log_size_increase']).
          reindex(size_class = dbpm_fixed_inputs['log10_size_bins'], 
                  fill_value = 0))     
    
    Bi = ((1+(1/np.log(10))*growth_rate.isel(size_class = index_sc)*
           gridded_params['timesteps_years']/
           gridded_params['log_size_increase']+
           total_mortality.isel(size_class = index_sc)*
             gridded_params['timesteps_years']).
            reindex(size_class = dbpm_fixed_inputs['log10_size_bins'],
                    fill_value = 0))
    
    Si = (biomass.isel(size_class = index_sc).
          reindex(size_class = dbpm_fixed_inputs['log10_size_bins'],
                  fill_value = 0))

    # Boundary condition at upstream end
    Ai = xr.where(Ai.size_class == ref_sc, 0, Ai)   
    Bi = xr.where(Bi.size_class == ref_sc, 1, Bi)
    if ref_sc not in growth_rate.size_class[index_sc]:
        Si = xr.where(Si.size_class == ref_sc, biomass, Si)

    density = xr.Dataset(data_vars = {'Ai': Ai, 'Bi': Bi, 'Si': Si})
    #Reorganise dimensions
    try:
        density = density.transpose('time', 'size_class', 'lat', 'lon')
    except:
        density = density.transpose('time', 'size_class')
    #Apply spatial mask
    density = density.where(dbpm_fixed_inputs['mask'])

    return density


def tot_biomass_calc(gridded_params, dbpm_fixed_inputs, group, biomass_current,
                     biomass_next, biomass_density, reprod_rate = None, 
                     growth_rate = None, total_mortality = None):

    if group == 'predators':
        ref_sc =  gridded_params['log10_ind_min_pred_size']
        ref = gridded_params['ind_min_pred_size']
    elif group == 'detritivores':
        ref_sc =  gridded_params['log10_ind_min_detritivore_size']
        ref = gridded_params['ind_min_detritivore_size']
    else:
        print('the "group" parameter must be either "predators" or ' +
              '"detritivores"')

    # Apply transfer efficiency of 10% *plankton density at same size
    # reproduction from energy allocation
    if gridded_params['dynamic_reproduction']:
        if growth_rate is None:
            raise Exception('When "dynamic_reproduction" is set to True, ' + 
                            '"growth_rate" must be provided as input')
        elif total_mortality is None:
            raise Exception('When "dynamic_reproduction" is set to True, ' +
                            '"total_mortality" must be provided as input')
        elif reprod_rate is None:
            raise Exception('When "dynamic_reproduction" is set to True, ' + 
                            '"reprod_rate" must be provided as input')
        else:
            growth_rate = growth_rate.sel(size_class = ref_sc)
            total_mortality = total_mortality.sel(size_class = ref_sc)
        bio_ref = biomass_current.sel(size_class = ref_sc)
        bio_ref_next = (bio_ref+
                        ((reprod_rate*dbpm_fixed_inputs['size_bin_vals']*
                          biomass_current*gridded_params['log_size_increase']).
                         sum('size_class')*gridded_params['timesteps_years'])/
                        (gridded_params['log_size_increase']*
                         dbpm_fixed_inputs['size_bin_vals'].sel(size_class = ref_sc))-
                        (gridded_params['timesteps_years']/
                         gridded_params['log_size_increase'])*(1/np.log(10))*
                        growth_rate*bio_ref-gridded_params['timesteps_years']*
                        total_mortality*bio_ref)
        # Change time stamp to next month
        if isinstance(biomass_next.time.values, np.ndarray):
            bio_ref_next['time'] = biomass_next.time.values
        else:
            bio_ref_next['time'] = [biomass_next.time.values]
        #Replace boundary condition in next timestep
        biomass_next = xr.where(biomass_next.size_class == ref_sc, bio_ref_next,
                                biomass_next)
        del bio_ref, bio_ref_next

    #main loop calculation
    biomass_density = biomass_density.drop_vars('time').squeeze()
    for sc in range(ref+1, gridded_params['numb_size_bins']):
        bio_next_shift = biomass_next.isel(size_class = sc-1).drop_vars('size_class')
        da = ((biomass_density['Si'].isel(size_class = sc)-
               biomass_density['Ai'].isel(size_class = sc)*bio_next_shift)/
              biomass_density['Bi'].isel(size_class = sc))
        biomass_next = xr.where(biomass_next.size_class == da.size_class, da, 
                                biomass_next)
    
    #Adding name to data array
    biomass_next.name = group
    #Reorganise dimensions
    try:
        biomass_next = biomass_next.transpose('time', 'size_class', 'lat', 'lon')
    except:
        biomass_next = biomass_next.transpose('time', 'size_class')
    #Apply spatial mask
    biomass_next = biomass_next.where(dbpm_fixed_inputs['mask'])
        
    return biomass_next


def detritus_rk4_step(gridded_params, dbpm_fixed_inputs, dbpm_init_inputs, 
                      dbpm_dynamic_inputs, detritus_current, h):
    """
    Single RK4 step for detritus calculation
    
    Parameters:
    - h: timestep size (gridded_params['timesteps_years'])
    - detritus_current: current detritus values
    """
    
    def calculate_dW_dt(detritus_state, t_fraction = 0):
        """Calculate dW/dt for given detritus state"""
        # Create temporary inputs with updated detritus
        temp_inputs = dbpm_init_inputs.copy()
        temp_inputs['detritus'] = detritus_state
        
        # Recalculate feeding rates with updated detritus
        feed_detritus = feeding_satiation_rates(gridded_params, dbpm_fixed_inputs, 
                                               temp_inputs, dbpm_dynamic_inputs, 
                                               group = 'detritus')
        
        # Calculate outputs and defecation (these depend on current biomass)
        output_w = detritus_output(gridded_params, dbpm_fixed_inputs, temp_inputs,
                                  feed_detritus.f_det)
        
        feed_sat_pred = feeding_satiation_rates(gridded_params, dbpm_fixed_inputs, 
                                               temp_inputs, dbpm_dynamic_inputs, 
                                               group = 'predators')
        feed_sat_det = feeding_satiation_rates(gridded_params, dbpm_fixed_inputs, 
                                              temp_inputs, dbpm_dynamic_inputs, 
                                              group = 'detritivores')
        
        defbypred = defecation_predators(gridded_params, dbpm_fixed_inputs,
                                        temp_inputs, 
                                        xr.Dataset({'f_pel': feed_sat_pred.f_pel,
                                                   'f_ben': feed_sat_det.f_ben}))
        
        # Calculate detritus pool changes
        det_pool = detritus_pool(gridded_params, dbpm_fixed_inputs, temp_inputs,
                                dbpm_dynamic_inputs, defbypred, output_w)
        
        return det_pool['dW']
    
    # RK4 implementation
    k1 = calculate_dW_dt(detritus_current)
    k2 = calculate_dW_dt(detritus_current + h * k1 / 2)
    k3 = calculate_dW_dt(detritus_current + h * k2 / 2)
    k4 = calculate_dW_dt(detritus_current + h * k3)
    
    # Combine the slopes
    detritus_next = detritus_current + h * (k1 + 2*k2 + 2*k3 + k4) / 6
    
    return detritus_next


# Run model per grid cell or averaged over an area ------
def gridded_sizemodel_rk4(gridded_params, dbpm_fixed_inputs, dbpm_init_inputs, 
                          dbpm_dynamic_inputs, region, model_res, out_folder, 
                          weekly = False, force_finite = True, rk4_substeps = 4):
    """
    Modified gridded_sizemodel with RK4 for detritus calculation
    """
    
    # Storing date from dynamic dataset
    dbpm_time = dbpm_dynamic_inputs.time.values
    if not weekly:
        pred_ts_next = pd.Timestamp(dbpm_time).strftime('%Y-%m')
    else:
        pred_ts_next = pd.Timestamp(dbpm_time).strftime('%Y-%m-%d')

    # Storing date from initialising dataset
    [pred_time] = dbpm_init_inputs.time.values
    if not weekly:
        pred_ts = pd.Timestamp(pred_time).strftime('%Y-%m')
    else:
        pred_ts = pd.Timestamp(pred_time).strftime('%Y-%m-%d')

    # Feeding and satiation rates ----
    # Predators
    feed_sat_pred = feeding_satiation_rates(gridded_params, dbpm_fixed_inputs, 
                                            dbpm_init_inputs, dbpm_dynamic_inputs, 
                                            group = 'predators')
    # Detritivores
    feed_sat_det = feeding_satiation_rates(gridded_params, dbpm_fixed_inputs, 
                                           dbpm_init_inputs, dbpm_dynamic_inputs, 
                                           group = 'detritivores')
    # Detritus
    feed_detritus = feeding_satiation_rates(gridded_params, dbpm_fixed_inputs, 
                                            dbpm_init_inputs, dbpm_dynamic_inputs, 
                                            group = 'detritus')
    
    
    # Mortality terms ----
    # Predators
    mortality_pred = mortality_calc(gridded_params, dbpm_fixed_inputs, 
                                    dbpm_init_inputs, dbpm_dynamic_inputs, 
                                    feed_sat_pred, group = 'predators')
    # Saving mortality rates for predators
    fn = f'pred_mort_pred_{model_res}_{region}_{pred_ts}.nc'
    mortality_pred['PM_u'].to_netcdf(os.path.join(out_folder, fn))

    # Detritivores
    mortality_det = mortality_calc(gridded_params, dbpm_fixed_inputs, 
                                   dbpm_init_inputs, dbpm_dynamic_inputs, 
                                   feed_sat_det, group = 'detritivores')
    # Saving mortality rates for detritivores
    fn = f'pred_mort_det_{model_res}_{region}_{pred_ts}.nc'
    mortality_det['PM_v'].to_netcdf(os.path.join(out_folder, fn))

    
    # Growth rates ----
    # Predators and detritivores
    growth_rates_pred_det = growth_rates_calc(gridded_params, 
                                              xr.Dataset({'f_pel': feed_sat_pred.f_pel,
                                                          'f_ben': feed_sat_det.f_ben,
                                                          'f_det': feed_detritus.f_det}), 
                                              dbpm_fixed_inputs)
    # Saving growth rates for predators
    fn = f'growth_int_pred_{model_res}_{region}_{pred_ts}.nc'
    growth_rates_pred_det['GG_u'].to_netcdf(os.path.join(out_folder, fn))
    # Saving growth rates for detritivores
    fn = f'growth_int_det_{model_res}_{region}_{pred_ts}.nc'
    growth_rates_pred_det['GG_v'].to_netcdf(os.path.join(out_folder, fn))


    # Reproduction rates ----
    if gridded_params['dynamic_reproduction']:
        reprod_pred = reproduction_rates(gridded_params, 
                                         xr.Dataset({'f_pel': feed_sat_pred.f_pel,
                                                     'f_ben': feed_sat_det.f_ben}),
                                         dbpm_fixed_inputs, group = 'predators')
        reprod_det = reproduction_rates(gridded_params, feed_detritus,
                                        dbpm_fixed_inputs, group = 'detritivores')

    # Replace the simple Euler detritus calculation with RK4
    h = gridded_params['timesteps_years'] / rk4_substeps  # Subdivide timestep
    detritus = dbpm_init_inputs['detritus'].copy()
    
    # Perform multiple RK4 substeps for stability
    for substep in range(rk4_substeps):
        detritus = detritus_rk4_step(gridded_params, dbpm_fixed_inputs, 
                                    dbpm_init_inputs, dbpm_dynamic_inputs, 
                                    detritus, h)
        detritus = xr.where(detritus < 0, 0, detritus).load()
        
    # Updating timestamp (results for next time step)
    if force_finite:
        detritus = detritus.where(np.isfinite(detritus))
    detritus['time'] = [dbpm_time]
    detritus.name = 'detritus'
    
    # Save detritus biomass
    fn = f'detritus_{model_res}_{region}_{pred_ts_next}.nc'
    detritus.to_netcdf(os.path.join(out_folder, fn))

    # Pelagic predator density ----
    pred_density = biomass_density(gridded_params, dbpm_fixed_inputs, 
                                   growth_rates_pred_det['GG_u'],
                                   mortality_pred['Z_u'], 
                                   dbpm_init_inputs['predators'], 
                                   group = 'predators')

    
    # Recruitment at smallest consumer mass ----
    # continuation of plankton hold constant  
    pred_next = (dbpm_dynamic_inputs['ui0']*
                 (10**(dbpm_dynamic_inputs['slope']*
                       dbpm_fixed_inputs['log10_size_bins'])))
    pred_next = xr.where((pred_next.size_class <= 
                          gridded_params['log10_ind_min_pred_size']),
                         pred_next, 0)

    predators = tot_biomass_calc(gridded_params, dbpm_fixed_inputs, 
                                 'predators', dbpm_init_inputs['predators'], pred_next, 
                                 pred_density, reprod_rate = reprod_pred['R_u'],
                                 growth_rate = growth_rates_pred_det['GG_u'], 
                                 total_mortality = mortality_pred['Z_u']).load()
    
    #If values are negative, assign a value of 0
    predators = xr.where(predators < 0, 0, predators)
    if force_finite:
        predators = predators.where(np.isfinite(predators))
    
    # Saving predators biomass
    # Creating file name
    fn = f'predators_{model_res}_{region}_{pred_ts_next}.nc'
    predators.to_netcdf(os.path.join(out_folder, fn))

    
    # Benthic detritivore density ----
    det_density = biomass_density(gridded_params, dbpm_fixed_inputs, 
                                  growth_rates_pred_det['GG_v'], 
                                  mortality_det['Z_v'], 
                                  dbpm_init_inputs['detritivores'], 
                                  group = 'detritivores')

    
    # Recruitment at smallest detritivore mass ----
    # hold constant continution of plankton with sinking rate multiplier 
    detriti_next = xr.where((dbpm_init_inputs['detritivores'].size_class <= 
                             gridded_params['log10_ind_min_detritivore_size']),
                            dbpm_init_inputs['detritivores'], 0)
    # Updating timestamp (results for next time step)
    detriti_next['time'] = [dbpm_time]

    detritivores = tot_biomass_calc(gridded_params, dbpm_fixed_inputs, 
                                    'detritivores', dbpm_init_inputs['detritivores'],
                                    detriti_next, det_density, 
                                    reprod_rate = reprod_det['R_v'],
                                    growth_rate = growth_rates_pred_det['GG_v'], 
                                    total_mortality = mortality_det['Z_v']).load()
    #If values are negative, assign a value of 0
    detritivores = xr.where(detritivores < 0, 0, detritivores)
    if force_finite:
        detritivores = detritivores.where(np.isfinite(detritivores))
    
    # Saving detritivores biomass
    fn = f'detritivores_{model_res}_{region}_{pred_ts_next}.nc'
    detritivores.to_netcdf(os.path.join(out_folder, fn))

    #Creating dataset with new values for predators and detritivores biomass, and detritus
    ds_next = xr.Dataset(data_vars = {'predators': predators, 
                                      'detritivores': detritivores, 
                                      'detritus': detritus})

    return ds_next


# Run model per grid cell or averaged over an area ------
def gridded_sizemodel(gridded_params, dbpm_fixed_inputs, dbpm_init_inputs, 
                      dbpm_dynamic_inputs, region, model_res, out_folder, 
                      weekly = False, force_finite = True):
    '''
    Inputs:
    - gridded_params (dictionary) Gridded DBPM parameters obtained in step 04.
    - dbpm_fixed_inputs (xarray Dataset) Contains fixed gridded inputs needed to run 
    DBPM
    - dbpm_init_inputs (xarray Dataset) Contains gridded inputs with initialisation 
    values for predators, detritivores and detritus
    - dbpm_dynamic_inputs (xarray Dataset) Contains dynamic gridded inputs need to run
    DBPM
    - region (character) Name of region being processed (as included in folder and file
    names)
    - model_res (character) Spatial resolution of DBPM inputs (as included in folder
    and file names)
    - out_folder (character) Full path to folder where DBPM outputs will be stored
    - weekly (boolean) Default is False. If set to True, the model will run at weekly
    timesteps

    Outputs:
    - None. This function does not return any objects. Instead outputs are saved in
    the path provided in the "out_folder" argument
    '''
    
    # Storing date from dynamic dataset
    dbpm_time = dbpm_dynamic_inputs.time.values
    if not weekly:
        pred_ts_next = pd.Timestamp(dbpm_time).strftime('%Y-%m')
    else:
        pred_ts_next = pd.Timestamp(dbpm_time).strftime('%Y-%m-%d')

    # Storing date from initialising dataset
    [pred_time] = dbpm_init_inputs.time.values
    if not weekly:
        pred_ts = pd.Timestamp(pred_time).strftime('%Y-%m')
    else:
        pred_ts = pd.Timestamp(pred_time).strftime('%Y-%m-%d')

    # Feeding and satiation rates ----
    # Predators
    feed_sat_pred = feeding_satiation_rates(gridded_params, dbpm_fixed_inputs, 
                                            dbpm_init_inputs, dbpm_dynamic_inputs, 
                                            group = 'predators')
    # Detritivores
    feed_sat_det = feeding_satiation_rates(gridded_params, dbpm_fixed_inputs, 
                                           dbpm_init_inputs, dbpm_dynamic_inputs, 
                                           group = 'detritivores')
    # Detritus
    feed_detritus = feeding_satiation_rates(gridded_params, dbpm_fixed_inputs, 
                                            dbpm_init_inputs, dbpm_dynamic_inputs, 
                                            group = 'detritus')
    
    
    # Mortality terms ----
    # Predators
    mortality_pred = mortality_calc(gridded_params, dbpm_fixed_inputs, 
                                    dbpm_init_inputs, dbpm_dynamic_inputs, 
                                    feed_sat_pred, group = 'predators')
    # Saving mortality rates for predators
    fn = f'pred_mort_pred_{model_res}_{region}_{pred_ts}.nc'
    mortality_pred['PM_u'].to_netcdf(os.path.join(out_folder, fn))

    # Detritivores
    mortality_det = mortality_calc(gridded_params, dbpm_fixed_inputs, 
                                   dbpm_init_inputs, dbpm_dynamic_inputs, 
                                   feed_sat_det, group = 'detritivores')
    # Saving mortality rates for detritivores
    fn = f'pred_mort_det_{model_res}_{region}_{pred_ts}.nc'
    mortality_det['PM_v'].to_netcdf(os.path.join(out_folder, fn))

    
    # Growth rates ----
    # Predators and detritivores
    growth_rates_pred_det = growth_rates_calc(gridded_params, 
                                              xr.Dataset({'f_pel': feed_sat_pred.f_pel,
                                                          'f_ben': feed_sat_det.f_ben,
                                                          'f_det': feed_detritus.f_det}), 
                                              dbpm_fixed_inputs)
    # Saving growth rates for predators
    fn = f'growth_int_pred_{model_res}_{region}_{pred_ts}.nc'
    growth_rates_pred_det['GG_u'].to_netcdf(os.path.join(out_folder, fn))
    # Saving growth rates for detritivores
    fn = f'growth_int_det_{model_res}_{region}_{pred_ts}.nc'
    growth_rates_pred_det['GG_v'].to_netcdf(os.path.join(out_folder, fn))


    # Reproduction rates ----
    if gridded_params['dynamic_reproduction']:
        reprod_pred = reproduction_rates(gridded_params, 
                                         xr.Dataset({'f_pel': feed_sat_pred.f_pel,
                                                     'f_ben': feed_sat_det.f_ben}),
                                         dbpm_fixed_inputs, group = 'predators')
        reprod_det = reproduction_rates(gridded_params, feed_detritus,
                                        dbpm_fixed_inputs, group = 'detritivores')
    
    # Detritus output (g.m-2.yr-1) ----
    # losses from detritivore scavenging/filtering only:
    output_w = detritus_output(gridded_params, dbpm_fixed_inputs, dbpm_init_inputs,
                               feed_detritus.f_det)


    # Defecation by predators ----
    defbypred = defecation_predators(gridded_params, dbpm_fixed_inputs,
                                     dbpm_init_inputs, 
                                     xr.Dataset({'f_pel': feed_sat_pred.f_pel,
                                                 'f_ben': feed_sat_det.f_ben}))

    # Detritus Biomass Density Pool ----
    # Fluxes in and out (g.m-2.yr-1) of detritus pool. Solve for detritus biomass
    # density in next timestep 
    # Increment values of detritus, predators & detritivores for next 
    # timestep
    det_pool = detritus_pool(gridded_params, dbpm_fixed_inputs, dbpm_init_inputs,
                             dbpm_dynamic_inputs, defbypred, output_w)

    # Biomass density of detritus g.m-2
    detritus = (dbpm_init_inputs['detritus']+det_pool['dW']*
                gridded_params['timesteps_years']).load()
    # detritus = (dbpm_init_inputs['detritus']*det_pool['dW']).load()

    #If values are negative, assign a value of 0
    detritus = xr.where(detritus < 0, 0, detritus)
    if force_finite:
        detritus = detritus.where(np.isfinite(detritus))
        
    # Updating timestamp (results for next time step)
    detritus['time'] = [dbpm_time]
    #Adding name
    detritus.name = 'detritus'
    # Saving detritus biomass
    fn = f'detritus_{model_res}_{region}_{pred_ts_next}.nc'
    detritus.to_netcdf(os.path.join(out_folder, fn))

    
    # Pelagic predator density ----
    pred_density = biomass_density(gridded_params, dbpm_fixed_inputs, 
                                   growth_rates_pred_det['GG_u'],
                                   mortality_pred['Z_u'], 
                                   dbpm_init_inputs['predators'], 
                                   group = 'predators')

    
    # Recruitment at smallest consumer mass ----
    # continuation of plankton hold constant  
    pred_next = (dbpm_dynamic_inputs['ui0']*
                 (10**(dbpm_dynamic_inputs['slope']*
                       dbpm_fixed_inputs['log10_size_bins'])))
    pred_next = xr.where((pred_next.size_class <= 
                          gridded_params['log10_ind_min_pred_size']),
                         pred_next, 0)

    predators = tot_biomass_calc(gridded_params, dbpm_fixed_inputs, 
                                 'predators', dbpm_init_inputs['predators'], pred_next, 
                                 pred_density, reprod_rate = reprod_pred['R_u'],
                                 growth_rate = growth_rates_pred_det['GG_u'], 
                                 total_mortality = mortality_pred['Z_u']).load()
    
    #If values are negative, assign a value of 0
    predators = xr.where(predators < 0, 0, predators)
    if force_finite:
        predators = predators.where(np.isfinite(predators))
    
    # Saving predators biomass
    # Creating file name
    fn = f'predators_{model_res}_{region}_{pred_ts_next}.nc'
    predators.to_netcdf(os.path.join(out_folder, fn))

    
    # Benthic detritivore density ----
    det_density = biomass_density(gridded_params, dbpm_fixed_inputs, 
                                  growth_rates_pred_det['GG_v'], 
                                  mortality_det['Z_v'], 
                                  dbpm_init_inputs['detritivores'], 
                                  group = 'detritivores')

    
    # Recruitment at smallest detritivore mass ----
    # hold constant continution of plankton with sinking rate multiplier 
    detriti_next = xr.where((dbpm_init_inputs['detritivores'].size_class <= 
                             gridded_params['log10_ind_min_detritivore_size']),
                            dbpm_init_inputs['detritivores'], 0)
    # Updating timestamp (results for next time step)
    detriti_next['time'] = [dbpm_time]

    detritivores = tot_biomass_calc(gridded_params, dbpm_fixed_inputs, 
                                    'detritivores', dbpm_init_inputs['detritivores'],
                                    detriti_next, det_density, 
                                    reprod_rate = reprod_det['R_v'],
                                    growth_rate = growth_rates_pred_det['GG_v'], 
                                    total_mortality = mortality_det['Z_v']).load()
    #If values are negative, assign a value of 0
    detritivores = xr.where(detritivores < 0, 0, detritivores)
    if force_finite:
        detritivores = detritivores.where(np.isfinite(detritivores))
    
    # Saving detritivores biomass
    fn = f'detritivores_{model_res}_{region}_{pred_ts_next}.nc'
    detritivores.to_netcdf(os.path.join(out_folder, fn))

    #Creating dataset with new values for predators and detritivores biomass, and detritus
    ds_next = xr.Dataset(data_vars = {'predators': predators, 
                                      'detritivores': detritivores, 
                                      'detritus': detritus})

    return ds_next


# Calculate decadal catches ------
def decadal_catches(catches_da, mask_da, **kwargs):
    '''
    Inputs:
    - catches_da (data array) 3D data array containing total catches over time. 
    - mask_da (data array) 2D data array used to mask out any areas outside the
    region of interest. Areas outside boundaries must be identified with NaN.
    **Optional**: 
    - name (character) Name to be applied to decadal data array
    - attrs (dictionary) Containing attributes to be added to data array

    Output:
    - dec_da (data array): Mean catches per decade
    '''

    catches_da['decade'] = (np.floor(catches_da.time.dt.year/10)*10).astype(int)
    dec_da = catches_da.groupby('decade').mean()
    dec_da = dec_da.where(np.isfinite(mask_da))
    if 'name' in kwargs.keys():
        dec_da.name = kwargs.get('name')
    if 'attrs' in kwargs.keys():
        dec_da = dec_da.assign_attrs(kwargs.get('attrs'))
    return dec_da


# Merging DBPM gridded outputs (yearly or decadal) ----
def merge_files(var, folder, merge_by = 'decade', **kwargs):
    #Getting list of files in folder
    file_list = sorted(glob(os.path.join(folder, f'{var}*')))

    #Getting the base file name to save merged data
    base_file = os.path.basename(file_list[-1])

    #Get years if provided
    if 'years' in kwargs:
        yrs = kwargs.get('years')
    else:
        #Get year from file names
        yrs = np.unique([int(re.findall('_(\\d{4})-\\d{2}', 
                                        os.path.basename(f))[0]) for f in file_list])
    
    #Merge by year or decade
    if merge_by == 'decade':
        dec = np.unique([np.floor(y/10).astype('int') for y in yrs])
        for d in dec:
            #Get files included in current decade
            sub_list = [f for f in file_list if str(d) in f]
            #Getting minimum and maximum years in the decade
            ymax = max([y for y in yrs if str(d) in str(y)])
            ymin = min([y for y in yrs if str(d) in str(y)])
            #Create file name to save merged data
            f_out = re.sub('\\d{4}-\\d{2}|\\d{4}', f'{str(ymin)}-{str(ymax)}', base_file)
            #Open dataset
            ds = xr.open_mfdataset(sub_list)
            #Save dataset
            ds.to_netcdf(os.path.join(folder, f_out))
    if merge_by == 'year':
        for y in yrs:
            #Get files included in current decade
            sub_list = [f for f in file_list if str(y) in f]
            #Create file name to save merged data
            f_out = re.sub('\\d{4}-\\d{2}|\\d{4}', f'{str(y)}', base_file)
            #Open dataset
            ds = xr.open_mfdataset(sub_list)
            #Save dataset
            ds.to_netcdf(os.path.join(folder, f_out))
    else:
        print('"merge_by" should be either "decade" or "year"')

