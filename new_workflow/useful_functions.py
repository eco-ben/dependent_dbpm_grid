#!/usr/bin/env python3

# Library of useful functions
# Authors: Denisse Fierro Arcos
# Date last update: 2024-10-10

# Loading libraries
import xarray as xr
import re
import numpy as np
import os
from glob import glob
import pandas as pd
import json


#Transforming netCDF files to zarr
def netcdf_to_zarr(file_path, path_out):
    '''
    Inputs:
    - file_path (character) File path where GFDL output is located
    - path_out (character) File path where outputs should be stored as zarr files

    Outputs:
    - None. This function saves results as zarr files in the path provided.
    '''

    #Loading and rechunking data
    da = xr.open_dataarray(file_path).chunk({'lat': '50MB', 'lon': '50MB'})

    #Change date format
    try:
        da['time'] = pd.DatetimeIndex(da.indexes['time'].to_datetimeindex(),
                                      dtype = 'datetime64[ns]')
    except:
        pass

    #Save results
    da.to_zarr(path_out, consolidated = True, mode = 'w')


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
    
    #Apply mask and remove rows where all grid cells are empty to reduce data array size
    if cross_dateline:
        da = da.where(mask == 1)
        da['lon'] = da.lon%360
        da = da.sortby('lon')
        da = da.dropna(dim = 'lon', how = 'all').dropna(dim = 'lat', how = 'all')
    else:
        da = da.where(mask == 1, drop = True)
    
    #Rechunking data
    if 'time' in da.dims:
        da = da.chunk({'time': '50MB', 'lat': len(da.lat), 'lon': len(da.lon)})
    else:
        da = da.chunk({'lat': len(da.lat), 'lon': len(da.lon)})

    #Save results
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
def getExportRatio(folder_gridded_data, gfdl_exp):
    '''
    Inputs:
    - folder_gridded_data (character) File path pointing to folder containing
    zarr files with GFDL data for the region of interest
    - gfdl_exp (character) Select 'ctrl_clim' or 'obs_clim'

    Outputs:
    - sphy (data array) Contains small phytoplankton. This is data frame simply
    renames 'phypico-vint' data
    - lphy (data array) Contains large phytoplankton: difference between 'phyc_vint'
    and 'phypico_vint'
    - er (data array) Contains export ratio
    '''

    #Get list of files in experiment
    file_list = glob(os.path.join(folder_gridded_data, f'*_{gfdl_exp}_*'))
    
    #Load sea surface temperature
    tos = xr.open_zarr([f for f in file_list if '_tos_' in f][0])['tos']
    
    #load depth
    depth = xr.open_zarr([f for f in file_list if '_deptho_' in f][0])['deptho']
    
    #Load phypico-vint
    sphy = xr.open_zarr([f for f in file_list if '_phypico-vint_' in f][0])['phypico-vint']
    #Rename phypico-vint to sphy
    sphy.name = 'sphy'

    #Load phyc-vint
    ptotal = xr.open_zarr([f for f in file_list if '_phyc-vint_' in f][0])['phyc-vint']

    #Calculate large phytoplankton
    lphy = ptotal-sphy
    #Give correct name to large phytoplankton dataset
    lphy.name = 'lphy'

    #Calculate phytoplankton size ratios
    plarge = lphy/ptotal
    psmall = sphy/ptotal

    #Calculate export ration
    er = (np.exp(-0.032*tos)*((0.14*psmall)+(0.74*(plarge)))+
          (0.0228*(plarge)*(depth*0.004)))/(1+(depth*0.004))
    #If values are negative, assign a value of 0
    er = xr.where(er < 0, 0, er)
    #If values are above 1, assign a value of 1
    er = xr.where(er > 1, 1, er)
    er.name = 'export_ratio'
    
    return sphy, lphy, er


#Calculate slope and intercept
def GetPPIntSlope(file_paths, mmin = 10**(-14.25), mmid = 10**(-10.184), 
                  mmax = 10**(-5.25), output = 'both'):
    '''
    Inputs:
    - file_paths (list) File paths pointing to zarr files from which slope and
    intercept will be calculated
    - mmin (numeric)  Default is 10**(-14.25). ????
    - mmid (numeric)  Default is 10**(-10.184). ????
    - mmax (numeric)  Default is 10**(-5.25). ????
    - output (character) Default is 'both'. Select what outputs should be returned. 
    Choose from 'both', 'slope', or 'intercept'

    Outputs:
    - (Data array) - Depends on value of 'output' parameter. 
    '''
        
    #Load depth
    depth = xr.open_zarr([f for f in file_paths if '_deptho_' in f][0])['deptho']
    
    #load large phytoplankton
    lphy = xr.open_zarr([f for f in file_paths if '_lphy_' in f][0])['lphy']
    
    #Load small phytoplankton
    sphy = xr.open_zarr([f for f in file_paths if '_sphy_' in f][0])['sphy']
    
    #Convert sphy and lphy from mol C / m^3 to g C / m^3
    sphy = sphy*12.0107
    lphy = lphy*12.0107

    #From Appendix of Barnes 2010 JPR, used in Woodworth-Jefcoats et al 2013
    #GCB. The scaling of biomass with body mass can be described as B=aM^b
    #the exponent b (also equivalent to the slope b in a log B vs log M 
    #relationship) can be assumed:
    #0.25 (results in N ~ M^-3/4) or 0 (results in N ~ M^-1)
    # most studies seem to suggest N~M^-1, so can assume that and test 
    #sensitivity of our results to this assumption. 
      
    #Calculate a and b in log B (log10 abundance) vs. log M (log10 gww)
    #in log10 gww
    midsmall = np.log10((mmin+mmid)/2) 
    midlarge = np.log10((mmid+mmax)/2)

    #convert to log10 (gww/size class median size) for log10 abundance
    small = np.log10((sphy*10)/(10**midsmall))
    #convert to log10 (gww/size class median size) for log10 abundance
    large = np.log10((lphy*10)/(10**midlarge))

    #Calculating lope
    b = (small-large)/(midsmall-midlarge)
    b.name = 'slope'

    #a is really log10(a), same a when small, midsmall are used
    a = large-(b*midlarge)
    a.name = 'intercept'

    # a could be used directly to replace 10^pp in sizemodel()
    if output == 'slope':
        return b
    if output == 'intercept':
        return a
    if output == 'both':
        return a, b


# Calculate spinup from gridded data
def gridded_spinup(file_path, start_spin, end_spin, spinup_period, 
                   mean_spinup = False, **kwargs):
    '''
    Inputs:
    - file_path (character) File path pointing to zarr file from which spinup
    will be calculated
    - start_spin (character or numeric) Year or date spinup starts
    - end_spin (character or numeric) Year or date spinup ends
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
    da = xr.open_zarr(file_path)
    #Getting name of variable contained in dataset
    [var] = list(da.keys())
    #Select period to be used for spinup
    da = da[var].sel(time = slice(str(start_spin), str(end_spin)))
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
    spinup_da = spinup_da.chunk({'time': '50MB', 'lat': 100, 'lon': 240})
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
# grid cell using graivity model, Walters & Bonfil, 1999, Gelchu & Pauly 2007 ideal free
# distribution - Blanchard et al 2008
def gravitymodel(effort, prop_b, depth, n_iter):
    '''
    Inputs:
    - effort (Data array) Fishing effort for a single time step
    - prop_b (Data array) Proportion of total fishable biomass for each grid cell at a 
    single time step
    - depth (Data array) Bathymetry of the area of interest
    - n_iter (integer) Number of iterations needed to redistribute fishing effort
    
    Outputs:
    eff (data array) Containing redistributed fishing effort
    '''

    eff = effort
    d = (1-(depth/depth.max()))

    #Initialise loop
    i = 1
    while(i <= n_iter):
        suit = prop_b*d
        rel_suit = suit/(suit.sum())
        neweffort = eff+(rel_suit*eff)
        mult = (eff.sum())/(neweffort.sum())
        eff = neweffort*mult
        i += i

    return eff


# Calculate fishing mortality and effort ------
def effort_calculation(predators, detritivores, effort, depth, size_bin_vals, 
                       gridded_params):
    '''
    Inputs:
    - predators (2D data array) Pelagic predator biomass
    - detritivores (2D data array) Bethic detritivore biomass
    - effort (2D data array) Fishing effort
    - depth (2D data array) Bathymetry of the area of interest
    - size_bins_vals (1D data array) Size classes in grams
    - gridded_params (dictionary) DBPM parameters

    Outputs:
    - new_effort (2D data array) Fishing effort calculated for next time step
    '''
    pred_bio = ((predators*size_bin_vals*gridded_params['log_size_increase']).
                isel(size_class = slice(gridded_params['ind_min_fish_pred'],
                                        -1))).sum('size_class')

    det_bio = ((detritivores*size_bin_vals*gridded_params['log_size_increase']).
               isel(size_class = slice(gridded_params['ind_min_fish_det'],
                                       -1))).sum('size_class')

    sum_bio = pred_bio+det_bio

    prop_b = sum_bio/sum_bio.sum()

    #Calculate new effort
    new_effort = gravitymodel(effort, prop_b, depth, 1)

    #Return effort
    return new_effort


# Loading fixed gridded DBPM inputs ------
def loading_dbpm_inputs(base_folder):
    '''
    Inputs:
    - base_folder (character) Full path to folder where fixed DBPM gridded inputs
    are stored

    Outputs:
    - dbpm_inputs (xarray Dataset) Contains fixed gridded inputs needed to run DBPM
    '''
    #Loading data from base folder
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

    dbpm_inputs = xr.Dataset(data_vars = {'pref_benthos': pref_benthos,
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
                                          'senes_mort_det': senes_mort_det})

    return dbpm_inputs


# Run model per grid cell or averaged over an area ------
def gridded_sizemodel(base_folder, dbpm_fixed_inputs, dbpm_init_inputs, 
                      dbpm_dynamic_inputs, region, model_res, out_folder, 
                      ERSEM_det_input = False):
    '''
    Inputs:
    - base_folder (character) Full path to folder where fixed DBPM gridded inputs
    are stored
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
    - ERSEM_det_input (boolean) Default is False. ?????

    Outputs:
    - None. This function does not return any objects. Instead outputs are saved in
    the path provided in the "out_folder" argument
    '''
    
    #Gridded parameters
    gridded_params = json.load(open(
        glob(os.path.join(base_folder, 'dbpm_gridded_*_python.json'))[0]))

    #Create a land mask from depth data
    mask = np.isfinite(dbpm_fixed_inputs['depth'])
    
    #Preparing most used parameters from gridded parameters
    log10_size_bins = np.array(gridded_params['log10_size_bins'])
    size_bin_vals = 10**dbpm_fixed_inputs['log10_size_bins_mat']
    log_size_increase = gridded_params['log_size_increase']
    ind_min_pred_size = gridded_params['ind_min_pred_size']
    ind_min_detritivore_size = gridded_params['ind_min_detritivore_size']
    
    #Calculating multipliers that do not require iteration
    #Size bin related
    size_bin_multi = size_bin_vals*log_size_increase
    ts_size_inc = gridded_params['timesteps_years']/log_size_increase
    
    #To be applied to feeding rates for pelagics and benthic groups
    feed_mult_pel = (gridded_params['hr_volume_search']*
                     (10**(dbpm_fixed_inputs['log10_size_bins_mat']*
                           gridded_params['metabolic_req_pred']))*
                     dbpm_fixed_inputs['pref_pelagic'])
    
    feed_mult_ben = (gridded_params['hr_volume_search']*
                     (10**(dbpm_fixed_inputs['log10_size_bins_mat']*
                           gridded_params['metabolic_req_pred']))*
                    dbpm_fixed_inputs['pref_benthos'])

    # Choosing first time step for effort to initialise model
    sinking_time = dbpm_dynamic_inputs['sinking_rate'].time.values
    eff_short = dbpm_dynamic_inputs['effort'].sel(time = min(sinking_time))

    #Separating predators, detritivores and 
    predators = dbpm_init_inputs['predators']
    detritivores = dbpm_init_inputs['detritivores']
    detritus = dbpm_init_inputs['detritus']
    
    #Looping through time
    for i in range(0, len(sinking_time)):
        #Saving time stamps for predators
        dbpm_time = sinking_time[i]
        pred_ts_next = pd.Timestamp(dbpm_time).strftime('%Y-%m')
        pred_time = pd.Timestamp(dbpm_time)-pd.DateOffset(months = 1)
        pred_ts = pd.Timestamp(pred_time).strftime('%Y-%m')
        
        #Selecting data for current time step - Removing timestamp to avoid issues
        #due to mismatching timestamps
        ben_temp_short = dbpm_dynamic_inputs['ben_tempeffect'].isel(time = i)
        pel_temp_short = dbpm_dynamic_inputs['pel_tempeffect'].isel(time = i)
        
        # Fishing mortality (FVec.u, FVec.v)
        # from Benoit & Rochet 2004
        # here fish_mort_pred and fish_mort_pred= fixed catchability term for
        # predators and detritivores to be estimated along with ind_min_det 
        # and ind_min_fish_pred
        fishing_mort_pred = (dbpm_fixed_inputs['fish_mort_pred']*
                             eff_short).drop_vars('time')
        fishing_mort_det = (dbpm_fixed_inputs['fish_mort_det']*
                            eff_short).drop_vars('time')
        
        # Feeding rates
        # Predators (f_pel)
        pred_growth = (((predators*log_size_increase).
                        dot(dbpm_fixed_inputs['constant_growth']).
                       rename({'sc': 'size_class'}))*feed_mult_pel)
        feed_rate_pel = (pel_temp_short*
                         (pred_growth/(1+gridded_params['handling']*pred_growth)))
        
        # Detritivores (f_ben)
        detrit_growth = (((detritivores*log_size_increase).
                          dot(dbpm_fixed_inputs['constant_growth']).
                          rename({'sc': 'size_class'}))*feed_mult_ben)
        feed_rate_bent = (pel_temp_short*
                          (detrit_growth/(1+gridded_params['handling']*detrit_growth)))
        #del detrit_growth

        # Feeding rates
        detritus_multiplier = ((1/size_bin_vals)*gridded_params['hr_vol_filter_benthos']*
                               10**(dbpm_fixed_inputs['log10_size_bins_mat']*
                                    gridded_params['metabolic_req_detritivore'])*
                               detritus)
        feed_rate_det = (ben_temp_short*detritus_multiplier/
                         (1+gridded_params['handling']*detritus_multiplier))
        #del detritus_multiplier

        # Calculate Growth and Mortality
        # Predator growth integral (GG_u)
        growth_int_pred = (gridded_params['growth_prop']*gridded_params['growth_pred']*
                           feed_rate_pel+gridded_params['high_prop']*
                           gridded_params['growth_detritivore']*feed_rate_bent)
        
        # Reproduction (R_u)
        if gridded_params['dynamic_reproduction']:
            reprod_pred = (gridded_params['growth_prop']*(1-(gridded_params['growth_pred']+
                                                             gridded_params['energy_pred']))*
                           feed_rate_pel.isel(size_class = slice(ind_min_pred_size+1, None))+
                           gridded_params['growth_prop']*
                           (1-(gridded_params['growth_detritivore']+
                               gridded_params['energy_detritivore']))*
                           feed_rate_bent.isel(size_class = slice(ind_min_pred_size+1, None)))
        
        # Predator death integrals
        # Satiation level of predator for pelagic prey (sat.pel)
        sat_pel = xr.where(feed_rate_pel > 0, feed_rate_pel/pred_growth, 0)
        #del pred_growth

        # Predation mortality
        # Predators (PM.u)
        pred_mort_pred = (dbpm_fixed_inputs['consume_pelagic']*
                          ((predators*sat_pel*log_size_increase).
                           dot(dbpm_fixed_inputs['constant_mortality'])).rename({'sc': 'size_class'}))
        #del sat_pel
        
        # Total mortality
        # Predators (Z.u)
        tot_mort_pred = (pred_mort_pred+pel_temp_short*
                         dbpm_fixed_inputs['other_mort_pred']+
                         dbpm_fixed_inputs['senes_mort_pred']+fishing_mort_pred)
    
        # Saving predation mortality
        # Adding time stamp
        if 'time' in pred_mort_pred.dims:
            pred_mort_pred['time'] = [pred_time]
        else:
            pred_mort_pred = pred_mort_pred.expand_dims({'time': [pred_time]})
        #Apply spatial mask
        pred_mort_pred = pred_mort_pred.where(mask)
        #Adding name
        pred_mort_pred.name = 'pred_mort_pred'
        # Creating file name
        fn = f'pred_mort_pred_{model_res}_{region}_{pred_ts}.nc'
        pred_mort_pred = pred_mort_pred.transpose('time', 'size_class', 'lat', 'lon')
        pred_mort_pred.to_netcdf(os.path.join(out_folder, fn))

        #Remove variables no longer needed
        #del pred_mort_pred, fishing_mort_pred, fn

        # Benthos growth integral yr-1 (GG_v)
        growth_int_det = (gridded_params['high_prop']*gridded_params['growth_detritus']*
                          feed_rate_det)
        
        # Reproduction (R_v)
        if gridded_params['dynamic_reproduction']:
            reprod_det = (gridded_params['high_prop']*
                          (1-(gridded_params['growth_detritus']+
                              gridded_params['energy_detritivore']))*
                          feed_rate_det.isel(size_class = 
                                             slice(ind_min_detritivore_size+1, None)))

        # Benthos death integral (sat_ben)
        # Satiation level of predator for benthic prey
        calc_growth = ((gridded_params['hr_volume_search']*
                        (10**(dbpm_fixed_inputs['log10_size_bins_mat']*
                              gridded_params['metabolic_req_detritivore']))*
                        dbpm_fixed_inputs['pref_benthos'])*
                       ((detritivores*log_size_increase).
                        dot(dbpm_fixed_inputs['constant_growth'])).
                       rename({'sc': 'size_class'}))
        sat_ben = xr.where(feed_rate_bent > 0, feed_rate_bent/calc_growth, 0)
        #del calc_growth

        # Predation mortality
        # Detritivores (PM.v)
        pred_mort_det = xr.where(sat_ben > 0, 
                                 (dbpm_fixed_inputs['consume_benthos']*
                                  ((predators*sat_ben*log_size_increase).
                                   dot(dbpm_fixed_inputs['constant_mortality'])).
                                  rename({'sc': 'size_class'})), 0)
        #del sat_ben

        # Total mortality
        # Detritivores (Z.v)
        tot_mort_det = (pred_mort_det+ben_temp_short*
                        dbpm_fixed_inputs['other_mort_det']+
                        dbpm_fixed_inputs['senes_mort_det']+fishing_mort_det)
        
        # Saving predation mortality
        # Adding time stamp
        if 'time' in pred_mort_det.dims:
            pred_mort_det['time'] = [pred_time]
        else:
            pred_mort_det = pred_mort_det.expand_dims({'time': [pred_time]})
        #Apply spatial mask
        pred_mort_det = pred_mort_det.where(mask)
        #Adding name
        pred_mort_det.name = 'pred_mort_det'
        # Creating file name
        fn = f'pred_mort_det_{model_res}_{region}_{pred_ts}.nc'
        pred_mort_det = pred_mort_det.transpose('time', 'size_class', 'lat', 'lon')
        pred_mort_det.to_netcdf(os.path.join(out_folder, fn))
        #Remove variables not in use
        #del pred_mort_det, fishing_mort_det, fn

        # Redistribute total effort across grid cells 
        try:
            effort_next = dbpm_dynamic_inputs['effort'].isel(time = i+1)
            t_next = dbpm_dynamic_inputs['effort'].time[i+1].values
            eff_short = effort_calculation(predators, detritivores, effort_next,
                                           dbpm_fixed_inputs['depth'], size_bin_vals, 
                                           gridded_params)
            # Saving predation mortality
            # Adding time stamp
            if 'time' in eff_short.dims:
                eff_short['time'] = [t_next]
            else:
                eff_short = eff_short.expand_dims({'time': [t_next]})
            #Apply spatial mask
            eff_short = eff_short.where(mask)
            #Adding name
            eff_short.name = 'effort'
            #Getting year and month 
            dt_eff = pd.to_datetime(t_next).strftime('%Y-%m')
            # Creating file name
            fn = f'effort_{model_res}_{region}_{dt_eff}.nc'
            eff_short = eff_short.transpose('time', 'lat', 'lon')
            eff_short.to_netcdf(os.path.join(out_folder, fn))
            #Remove variables not needed
            #del dt_eff, fn, effort_next
        except:
            pass
        
        # Detritus output (g.m-2.yr-1)
        # losses from detritivore scavenging/filtering only:
        output_w = (size_bin_multi*detritivores*feed_rate_det).sum('size_class')
        #del feed_rate_det

        # Total biomass density defecated by pred (g.m-2.yr-1)
        defbypred_size = (((gridded_params['defecate_prop']*feed_rate_pel*
                            size_bin_vals*predators+gridded_params['def_low']*
                            feed_rate_bent*size_bin_vals*predators).
                           isel(size_class = slice(ind_min_pred_size, None)))*
                         log_size_increase).sum('size_class')
        #del feed_rate_pel, feed_rate_bent
        
        # Increment values of detritus, predators & detritivores for next 
        # timestep
        # Detritus Biomass Density Pool - fluxes in and out (g.m-2.yr-1) of
        # detritus pool. Solve for detritus biomass density in next timestep 
        if not ERSEM_det_input:
            # Considering pelagic faeces as input as well as dead bodies 
            # from both pelagic and benthic communities and phytodetritus
            # (dying sinking phytoplankton)
            if gridded_params['detritus_coupling']:
                # Pelagic spectrum inputs (sinking dead bodies and faeces): 
                # Export ratio used for "sinking rate" + benthic spectrum 
                # inputs (dead stuff already on/in seafloor)
                input_w = (dbpm_dynamic_inputs['sinking_rate'].isel(time = i)* 
                           (defbypred_size+
                            (pel_temp_short*dbpm_fixed_inputs['other_mort_pred']*
                             predators*size_bin_multi).sum('size_class')+ 
                            (pel_temp_short*dbpm_fixed_inputs['senes_mort_pred']*
                             predators*size_bin_multi).sum('size_class'))+
                           ((ben_temp_short*dbpm_fixed_inputs['other_mort_det']*
                             detritivores*size_bin_multi).sum('size_class')+ 
                            (ben_temp_short*dbpm_fixed_inputs['senes_mort_det']*
                             detritivores*size_bin_multi).sum('size_class')))
            else:
                input_w = ((ben_temp_short*dbpm_fixed_inputs['other_mort_det']*
                             detritivores*size_bin_multi).sum('size_class')+ 
                            (ben_temp_short*dbpm_fixed_inputs['senes_mort_det']*
                             detritivores*size_bin_multi).sum('size_class'))
            
            # Get burial rate from Dunne et al. 2007 equation 3
            burial = input_w*(0.013+0.53*input_w**2/(7+input_w)**2)
            # Losses from detritivory + burial rate (not including 
            # remineralisation because that goes to p.p. after sediment, 
            # we are using realised p.p. as inputs to the model) 
            dW = input_w-(output_w+burial)
            # Biomass density of detritus g.m-2
            detritus = (detritus+dW*gridded_params['timesteps_years'])
        else:
            detritus = xr.zeros_like(detritus)
            
        # Saving detritus
        # Adding time stamp
        if 'time' in detritus.dims:
            detritus['time'] = [dbpm_time]
        else:
            detritus = detritus.expand_dims({'time': [dbpm_time]})
        #Apply spatial mask
        detritus = detritus.where(mask)
        #Adding name
        detritus.name = 'detritus'
        # Creating file name
        fn = f'detritus_{model_res}_{region}_{pred_ts_next}.nc'
        detritus = detritus.transpose('time', 'lat', 'lon')
        detritus.to_netcdf(os.path.join(out_folder, fn))
        #Remove variables not in use
        #del output_w, defbypred_size, input_w, burial, dW, fn

        # Pelagic Predator Density (nos.m-2) solved for time using implicit 
        # time Euler upwind finite difference (help from Ken Andersen and 
        # Richard Law)
        ggp_shift = growth_int_pred.shift({'size_class': min(gridded_params['idx'])})
        ggp_shift, range_sc = xr.align(ggp_shift, 
                                       dbpm_fixed_inputs['log10_size_bins_mat'].
                                       isel(size_class = gridded_params['idx']))
        
        Ai_u = (((1/np.log(10))*-ggp_shift*ts_size_inc).
                reindex(size_class = dbpm_fixed_inputs['log10_size_bins_mat'],
                        fill_value = 0))
        Bi_u = ((1+(1/np.log(10))*
                 growth_int_pred.isel(size_class = gridded_params['idx'])*
                 ts_size_inc+
                 tot_mort_pred.isel(size_class = gridded_params['idx'])*
                 gridded_params['timesteps_years']).
                reindex(size_class = dbpm_fixed_inputs['log10_size_bins_mat'],
                        fill_value = 0))
        Si_u = ((predators.isel(size_class = gridded_params['idx'])).
                reindex(size_class = dbpm_fixed_inputs['log10_size_bins_mat'],
                        fill_value = 0))
        #del ggp_shift, range_sc

        # Boundary condition at upstream end
        Ai_u = xr.where(Ai_u.size_class == gridded_params['log10_ind_min_pred_size'],
                        0, Ai_u)
        Bi_u = xr.where(Bi_u.size_class == gridded_params['log10_ind_min_pred_size'],
                        1, Bi_u)

        #Calculate predator biomass for next time step
        pred_next = (dbpm_dynamic_inputs['ui0'].isel(time = i)*
                     (10**(dbpm_dynamic_inputs['slope'].isel(time = i)*
                           dbpm_fixed_inputs['log10_size_bins_mat'])))
        pred_next = xr.where((pred_next.size_class <= 
                              gridded_params['log10_ind_min_pred_size']),
                             pred_next, 0)

        # Apply transfer efficiency of 10% *plankton density at same size
        # reproduction from energy allocation
        if gridded_params['dynamic_reproduction']:
            pred_repro = (predators.isel(size_class = ind_min_pred_size)+
                         ((reprod_pred*size_bin_multi*predators).sum('size_class')*
                          gridded_params['timesteps_years'])/
                          size_bin_multi.isel(size_class = ind_min_pred_size)-
                          ts_size_inc*(1/np.log(10))*
                          (growth_int_pred.isel(size_class = ind_min_pred_size))*
                          predators.isel(size_class = ind_min_pred_size)-
                          gridded_params['timesteps_years']*
                          tot_mort_pred.isel(size_class = ind_min_pred_size)*
                          predators.isel(size_class = ind_min_pred_size))
        
            pred_next = xr.where((pred_next.size_class == 
                                  gridded_params['log10_ind_min_pred_size']), 
                                 pred_repro, pred_next)
            #del pred_repro

        # Saving growth predators
        # Adding time stamp
        if 'time' in growth_int_pred.dims:
            growth_int_pred['time'] = [pred_time]
        else:
            growth_int_pred = growth_int_pred.expand_dims({'time': [pred_time]})
        #Apply spatial mask
        growth_int_pred = growth_int_pred.where(mask)
        #Adding name
        growth_int_pred.name = 'growth_int_pred'
        # Creating file name
        fn = f'growth_int_pred_{model_res}_{region}_{pred_ts}.nc'
        growth_int_pred = growth_int_pred.transpose('time', 'size_class', 'lat', 'lon')
        growth_int_pred.to_netcdf(os.path.join(out_folder, fn))
        #Remove variables not in use
        #del growth_int_pred, reprod_pred, tot_mort_pred, fn

        # Calculating predator biomass for next time step
        for sc in range(ind_min_pred_size+1, len(pred_next.size_class)):
            pred_next_shift = pred_next.isel(size_class = sc-1).drop_vars('size_class')
            da = ((Si_u.isel(size_class = sc)-Ai_u.isel(size_class = sc)*pred_next_shift)/
                  Bi_u.isel(size_class = sc))
            pred_next = xr.where(pred_next.size_class == log10_size_bins[sc], da, pred_next)
        
        predators = pred_next
        
        #Removing variables not needed
        #del pred_next_shift, da

        # Saving predators biomass
        # Adding time stamp
        if 'time' in predators.dims:
            predators['time'] = [dbpm_time]
        else:
            predators = predators.expand_dims({'time': [dbpm_time]})
        #Apply spatial mask
        predators = predators.where(mask)
        #Adding name
        predators.name = 'predators'
        # Creating file name
        fn = f'predators_{model_res}_{region}_{pred_ts_next}.nc'
        predators = predators.transpose('time', 'size_class', 'lat', 'lon')
        predators.to_netcdf(os.path.join(out_folder, fn))
        #Remove variables not in use
        #del Ai_u, Bi_u, Si_u, fn

        # Benthic Detritivore Density (nos.m-2)
        Ai_v = (((1/np.log(10))*
                 -growth_int_det.isel(size_class = np.array(gridded_params['idx_new'])-1)*
                 ts_size_inc).reindex(size_class = dbpm_fixed_inputs['log10_size_bins_mat'],
                                      fill_value = 0)).shift({'size_class': 1}, 
                                                             fill_value = 0)
        Bi_v = ((1+(1/np.log(10))*
                 growth_int_det.isel(size_class = gridded_params['idx_new'])*
                 ts_size_inc+
                 tot_mort_det.isel(size_class = gridded_params['idx_new'])*
                 gridded_params['timesteps_years']).
                reindex(size_class = dbpm_fixed_inputs['log10_size_bins_mat'],
                        fill_value = 0))
        
        Si_v = (detritivores.isel(size_class = gridded_params['idx_new']).
                reindex(size_class = dbpm_fixed_inputs['log10_size_bins_mat'],
                        fill_value = 0))
        
        # Boundary condition at upstream end
        Ai_v = xr.where((Ai_v.size_class == 
                         gridded_params['log10_ind_min_detritivore_size']), 
                        0, Ai_v)
        
        Bi_v = xr.where((Bi_v.size_class == 
                         gridded_params['log10_ind_min_detritivore_size']),
                        1, Bi_v)

        # Recruitment at smallest detritivore mass
        # hold constant continuation of plankton with sinking rate multiplier
        detriti_next = xr.where((detritivores.size_class <= 
                                 gridded_params['log10_ind_min_detritivore_size']),
                                detritivores, 0)

        # Apply a very low of transfer efficiency 1%* total biomass of 
        # detritus divided by minimum size
        if gridded_params['dynamic_reproduction']:
            det_repro = (detritivores.isel(size_class = ind_min_detritivore_size)+
                         ((reprod_det*size_bin_multi*detritivores).sum('size_class')*
                          gridded_params['timesteps_years'])/
                         (size_bin_multi.isel(size_class = ind_min_detritivore_size))-
                         ts_size_inc*(1/np.log(10))*
                         (growth_int_det.isel(size_class = ind_min_detritivore_size))*
                         (detritivores.isel(size_class = ind_min_detritivore_size))-
                         gridded_params['timesteps_years']*
                         (tot_mort_det.isel(size_class = ind_min_detritivore_size))*
                         (detritivores.isel(size_class = ind_min_detritivore_size)))
        
            detriti_next = xr.where((detriti_next.size_class ==
                                     gridded_params['log10_ind_min_detritivore_size']),
                                    det_repro, detriti_next)
            #del det_repro

        # Saving growth detritivores
        # Adding time stamp
        if 'time' in growth_int_det.dims:
            growth_int_det['time'] = [pred_time]
        else:
            growth_int_det = growth_int_det.expand_dims({'time': [pred_time]})
        #Apply spatial mask
        growth_int_det = growth_int_det.where(mask)
        #Adding name
        growth_int_det.name = 'growth_int_det'
        # Creating file name
        fn = f'growth_int_det_{model_res}_{region}_{pred_ts}.nc'
        growth_int_det = growth_int_det.transpose('time', 'size_class', 'lat', 'lon')
        growth_int_det.to_netcdf(os.path.join(out_folder, fn))
        #Remove variables not in use
        #del growth_int_det, reprod_det, tot_mort_det, fn

        # Calculate detritivore biomass for next time step
        for sc in gridded_params['idx_new']:
            detriti_next_shift = detriti_next.isel(size_class = sc-1).drop_vars('size_class')
            da = ((Si_v.isel(size_class = sc)-Ai_v.isel(size_class = sc)*detriti_next_shift)/
                  Bi_v.isel(size_class = sc))
            detriti_next = xr.where(detriti_next.size_class == log10_size_bins[sc], da, 
                                    detriti_next)
        detritivores = detriti_next
        #Removing variables not in use
        #del da, detriti_next_shift, Ai_v, Bi_v, Si_v

        # Saving predation mortality
        # Adding time stamp
        if 'time' in detritivores.dims:
            detritivores['time'] = [dbpm_time]
        else:
            detritivores = detritivores.expand_dims({'time': [dbpm_time]})
        #Apply spatial mask
        detritivores = detritivores.where(mask)
        #Adding name
        detritivores.name = 'detritivores'
        # Creating file name
        fn = f'detritivores_{model_res}_{region}_{pred_ts_next}.nc'
        detritivores = detritivores.transpose('time', 'size_class', 'lat', 'lon')
        detritivores.to_netcdf(os.path.join(out_folder, fn))
