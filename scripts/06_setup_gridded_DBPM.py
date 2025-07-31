#!/usr/bin/env python3

# Loading libraries
import os
from glob import glob
import xarray as xr
import json
import pandas as pd
import numpy as np
import useful_functions as uf

# Defining base variables
base_folder = '/g/data/vf71/fishmip_inputs/ISIMIP3a/fao_inputs'
# regions = [f for f in os.listdir(base_folder) if 'fao' in f]
regions = ['fao-31']
# resolutions = ['1deg', '025deg']
resolutions = ['1deg']
# Choose whether smoothing of inputs will be performed by LOESS (smoothed) or
# deseasoning data (deseasoned)
smoothing = 'deseasoned'

# Set variables to find correct input files based on 'smoothing' variable - Also
# used to name processed inputs
# Create variables based on 'smoothing' selection
if smoothing != None:
    smoothing = f'-{smoothing}'
    fn_search = '-smoothed'
else:
    smoothing = ''
    fn_search = ''

# Looping through all regions and resolutions
for reg in regions:
    for res in resolutions:
        # Base folder for FAO region
        reg_folder = os.path.join(base_folder, reg)

        # Model resolution in arcmin to add to file name
        if res == '1deg':
            res_arc = '60arcmin'
        elif res == '025deg':
            res_arc = '15arcmin'

        # Defining output folder
        out_folder = os.path.join(reg_folder, f'gridded{smoothing}', res)

        # Transforming DBPM parameters to Python-friendly format ----
        # These parameters are the outputs of the `sizeparam`function
        # for R. See script 04_calculating_dbpm_fishing_params for 
        # more details.
        fish_param_file = os.path.join(
            reg_folder, 'init_fish_vals', res, f'best_fish_vals{smoothing}',
                f'dbpm_gridded_size_params_{reg}.json')
    
        gridded_params = json.load(open(fish_param_file))
        
        #Transforming them to Python friendly format
        gridded_python = uf.gridded_param_python(gridded_params)

        #Adding useful entries
        # Size classes
        log10_size_bins = np.array(gridded_python['log10_size_bins'])
        # Minimum (log10) size class for predators
        log10_ind_min_pred_size = log10_size_bins[gridded_python['ind_min_pred_size']]
        ind_min_detritivore_size = gridded_python['ind_min_detritivore_size']
        log10_ind_min_detritivore_size = log10_size_bins[ind_min_detritivore_size]
        log10_ind_min_fish_pred = log10_size_bins[gridded_python['ind_min_fish_pred']]
        log10_ind_min_fish_det = log10_size_bins[gridded_python['ind_min_fish_det']]
        #Adding new variables to gridded parameters
        gridded_python['log10_ind_min_pred_size'] = log10_ind_min_pred_size
        gridded_python['log10_ind_min_detritivore_size'] = log10_ind_min_detritivore_size

        #Saving Python-friendly gridded parameters
        fout = os.path.basename(fish_param_file).replace(reg, f'{reg}_python')
        #Save to disk
        with open(os.path.join(
            os.path.dirname(fish_param_file), fout), 'w') as outfile:
            json.dump(gridded_python, outfile)

        # Loading gridded input data ----
        #Loading depth
        depth = xr.open_zarr(glob(os.path.join(
            out_folder, '*obsclim_deptho_*'))[0])['deptho']
        #Loading intercept data for stable spinup period
        int_phy_zoo_stable_spin = xr.open_zarr(glob(os.path.join(
            out_folder, f'*stable-spin_intercept_*_monthly_*'))[0])['intercept']
        #Loading intercept data for spinup period
        int_phy_zoo_spinup = xr.open_zarr(glob(os.path.join(
            out_folder, f'*spinup_intercept_*_monthly{fn_search}_*'))[0])['intercept']
        #Loading intercept data for model period 
        int_phy_zoo  = xr.open_zarr(glob(os.path.join(
            out_folder, f'*obsclim_intercept_*_monthly{fn_search}_*'))[0])['intercept']
        #Loading sea surface temperature data for stable spinup period
        sea_surf_temp_stable_spin = xr.open_zarr(glob(os.path.join(
            out_folder, f'*stable-spin_tos_*_monthly_*'))[0])['tos']
        #Loading sea surface temperature data for spinup period
        sea_surf_temp_spinup = xr.open_zarr(glob(os.path.join(
            out_folder, f'*spinup_tos_*_monthly{fn_search}_*'))[0])['tos']
        #Loading sea surface temperature data for model period
        sea_surf_temp  = xr.open_zarr(glob(os.path.join(
            out_folder, f'*obsclim_tos_*_monthly{fn_search}_*'))[0])['tos']
        #Loading bottom ocean temperature data for spinup period
        sea_floor_temp_stable_spin = xr.open_zarr(glob(os.path.join(
            out_folder, f'*stable-spin_tob_*_monthly_*'))[0])['tob']
        #Loading bottom ocean temperature data for spinup period
        sea_floor_temp_spinup = xr.open_zarr(glob(os.path.join(
            out_folder, f'*spinup_tob_*_monthly{fn_search}_*'))[0])['tob']
        #Loading bottom ocean temperature data for model period
        sea_floor_temp  = xr.open_zarr(glob(os.path.join(
            out_folder, f'*obsclim_tob_*_monthly{fn_search}_*'))[0])['tob']

        # Timesteps for stable spinup period
        spinup_period = pd.date_range('1741-01-01', end = '1840-12-31', 
                                      freq = 'MS')
        
        # Log10 size of individuals found in the model
        log10_file = os.path.join(base_folder, 'log10_size_bins_matrix.zarr')
        if not os.path.exists(log10_file):
            log10_size_bins_mat = xr.DataArray(
                data = log10_size_bins,
                dims = ['size_class'], 
                coords = {'size_class': log10_size_bins})
            
            # size of individuals found in the model
            size_bin_vals = 10**log10_size_bins_mat
            log10_size_bins_mat.name = 'size_bins'
            log10_size_bins_mat.to_zarr(
                log10_file, consolidated = True, mode = 'w')
        else:
            log10_size_bins_mat = xr.open_dataarray(log10_file)

        # Fishing effort ----       
        effort_spinup = (xr.DataArray(
            gridded_python['effort'][:len(int_phy_zoo_spinup.time)], 
            dims = 'time', 
            coords = {'time': int_phy_zoo_spinup.time}).
                         expand_dims({'lat': depth.lat, 
                                      'lon': depth.lon}).
                         transpose('time', 'lat', 
                                   'lon')).where(np.isfinite(depth))

        #Select first year of spinup data to create stable spinup dataset
        da = (effort_spinup.sel(time = slice('1841-01', '1841-12')).
              mean('time'))
        effort_stable_spin = [da] * len(spinup_period)
        effort_stable_spin = xr.concat(effort_stable_spin, dim = 'time')
        effort_stable_spin['time'] = spinup_period
        
        #Data fro modelling period
        effort = (xr.DataArray(
            gridded_python['effort'][len(int_phy_zoo_spinup.time):],
            dims = 'time', coords = {'time': int_phy_zoo.time}).
                  expand_dims({'lat': depth.lat, 'lon': depth.lon}).
                  transpose('time', 'lat',
                            'lon')).where(np.isfinite(depth))
        
        # Create a single effort file
        effort_out = xr.concat([effort_stable_spin, effort_spinup, effort],
                               dim = 'time').chunk({'lon': -1, 'lat': -1,
                                                    'time': -1}).load()
        effort_out.name = 'effort'
        effort_out.to_zarr(os.path.join(
            out_folder, 
            f'effort_spinup_obsclim_{res_arc}_{reg}_monthly{fn_search}_1741_2010.zarr/'),
                           consolidated = True, mode = 'w')

        # Habitat preferences ----
        pref_benthos = (0.8*np.exp((-1/250*depth)))
        pref_benthos.name = 'pref_benthos'
        pref_benthos.to_zarr(os.path.join(
            out_folder, f'pref-benthos_{res_arc}_{reg}_fixed.zarr/'),
                             consolidated = True, mode = 'w')

        pref_pelagic = (1-pref_benthos)
        pref_pelagic.name = 'pref_pelagic'
        pref_pelagic.to_zarr(os.path.join(
            out_folder, 
            f'pref-pelagic_{res_arc}_{reg}_fixed.zarr/'), 
                             consolidated = True, mode = 'w')

        # Metabolic requirements ----
        met_req_log10_size_bins = uf.expax_f(
            log10_size_bins_mat,
            gridded_python['metabolic_req_pred'])

        consume_pelagic = (pref_pelagic*
                           gridded_python['hr_volume_search']*
                           met_req_log10_size_bins)
        consume_pelagic.name = 'consume_pelagic'
        consume_pelagic.to_zarr(os.path.join(
            out_folder, 
            f'consume-pelagic_{res_arc}_{reg}_fixed.zarr/'),
                                consolidated = True, mode = 'w')
        consume_benthos = (pref_benthos*
                           gridded_python['hr_volume_search']*
                           met_req_log10_size_bins)
        consume_benthos.name = 'consume_benthos'
        consume_benthos.to_zarr(os.path.join(
            out_folder, 
            f'consume-benthos_{res_arc}_{reg}_fixed.zarr/'),
                                consolidated = True, mode = 'w')

        # Constant growth and mortality
        pred_prey_mat = uf.pred_prey_matrix(log10_size_bins)

        constant_growth = xr.DataArray(
            uf.gphi_f(pred_prey_mat, 
                      gridded_python['log10_pred_prey_ratio'],
                      gridded_python['log_prey_pref']),
            dims = ['size_class', 'sc'])
        constant_growth.name = 'constant_growth'
        constant_growth.to_zarr(os.path.join(
            out_folder, f'const-growth_{res_arc}_{reg}_fixed.zarr/'),
                                consolidated = True, mode = 'w')

        constant_mortality = xr.DataArray(
            uf.mphi_f(-pred_prey_mat, 
                      gridded_python['log10_pred_prey_ratio'],
                      gridded_python['log_prey_pref'],
                      gridded_python['metabolic_req_pred']),
            dims = ['size_class', 'sc'])
        constant_mortality.name = 'constant_mortality'
        constant_mortality.to_zarr(os.path.join(
            out_folder, f'const-mort_{res_arc}_{reg}_fixed.zarr/'),
                                consolidated = True, mode = 'w')

        # Gridded intercept plankton spectrum
        ui0 = xr.concat([int_phy_zoo_stable_spin, int_phy_zoo_spinup,
                         int_phy_zoo], 
                        dim = 'time').chunk({'lon': -1, 'lat': -1,
                                             'time': -1}).load()
        ui0_out = 10**ui0
        ui0_out.name = 'ui0'
        ui0_out.to_zarr(os.path.join(
            out_folder, 
            f'ui0_spinup_obsclim_{res_arc}_{reg}_monthly{fn_search}_1741_2010.zarr/'),
                    consolidated = True, mode = 'w')

        # Initial biomass ----
        # Timestep to initialise biomass values
        time_init = (pd.Timestamp(ui0_out.time.min().values)-
                     pd.DateOffset(months = 1))
        # Predators
        predators = (xr.DataArray(
            data = gridded_python['init_pred'], 
            dims = ['size_class'], 
            coords = {'size_class': log10_size_bins}).
                     expand_dims({'time': [time_init], 
                                  'lat': depth.lat, 
                                  'lon': depth.lon}))
        #Applying spatial mask
        predators = predators.where(np.isfinite(depth))
        predators.name = 'predators'
        predators.to_zarr(os.path.join(
            out_folder, 
            f'predators_{res_arc}_{reg}_init_{time_init.year}.zarr/'),
                          consolidated = True, mode = 'w')

        # Detritivores
        init_detritivores = (xr.DataArray(
            data = gridded_python['init_detritivores'], 
            dims = ['size_class'], 
            coords = {'size_class': log10_size_bins}).
                             expand_dims({'lat': depth.lat, 
                                          'lon': depth.lon}))
        # Assigning initial values
        detritivores = xr.where(
            (init_detritivores.size_class >= log10_ind_min_detritivore_size),
            init_detritivores, 0)
        # Adding time dimension
        detritivores = detritivores.expand_dims({'time': [time_init]})
        #Applying spatial mask
        detritivores = detritivores.where(np.isfinite(depth))
        detritivores.name = 'detritivores'
        detritivores.to_zarr(os.path.join(
            out_folder, 
            f'detritivores_{res_arc}_{reg}_init_{time_init.year}.zarr/'),
                             consolidated = True, mode = 'w')

        # Initialising detritus
        detritus = xr.zeros_like(detritivores.
                                 isel(size_class = 0).
                                 drop_vars('size_class'))
        # Assigning initial values 
        detritus = xr.where(detritus == 0,
                            gridded_python['init_detritus'], 
                            detritus)
        #Applying spatial mask
        detritus = detritus.where(np.isfinite(depth))
        detritus.name = 'detritus'
        detritus.to_zarr(os.path.join(
            out_folder, f'detritus_{res_arc}_{reg}_init_{time_init.year}.zarr/'),
                         consolidated = True, mode = 'w')

        # Other intrinsic natural mortality ----
        other_mort = (gridded_python['natural_mort']*
                      10**(-0.25*log10_size_bins_mat))
        other_mort.name = 'other_mort_pred'
        other_mort.to_zarr(os.path.join(
            out_folder, 
            f'other-mort-pred_{res_arc}_{reg}_fixed.zarr/'),
                           consolidated = True, mode = 'w')
        
        other_mort.name = 'other_mort_det'
        other_mort.to_zarr(os.path.join(
            out_folder, 
            f'other-mort-det_{res_arc}_{reg}_fixed.zarr/'),
                           consolidated = True, mode = 'w')
        
        # Senescence mortality ----
        senes_mort = (gridded_params['const_senescence_mort']*
                      10**(gridded_params['exp_senescence_mort']*
                           (log10_size_bins_mat-
                            gridded_params['size_senescence'])))
        senes_mort.name = 'senes_mort_pred'
        senes_mort.to_zarr(os.path.join(
            out_folder, 
            f'senes-mort-pred_{res_arc}_{reg}_fixed.zarr/'),
                           consolidated = True, mode = 'w')
        
        senes_mort.name = 'senes_mort_det'
        senes_mort.to_zarr(os.path.join(
            out_folder, 
            f'senes-mort-det_{res_arc}_{reg}_fixed.zarr/'),
                           consolidated = True, mode = 'w')

        # Fishing mortality ----
        # Using fish_mort_pred as a base for predator mortality
        fish_mort_pred = xr.DataArray(
            np.repeat(gridded_python['fish_mort_pred'], 
                      np.prod(consume_pelagic.shape)).
            reshape(consume_pelagic.shape),
            dims = ['lat', 'lon', 'size_class'],
            coords = {'lat': consume_pelagic.lat, 
                      'lon': consume_pelagic.lon,
                     'size_class': consume_pelagic.size_class})
        # Changing mortality to zero outside the predator size fished
        fishing_mort_pred = xr.where(
            (fish_mort_pred.size_class >= log10_ind_min_fish_pred) &
            (fish_mort_pred.size_class < fish_mort_pred.size_class.max()),
            fish_mort_pred, 0).where(np.isfinite(depth))
        fishing_mort_pred.name = 'fish_mort_pred'
        fishing_mort_pred.to_zarr(os.path.join(
            out_folder,
            f'fish-mort-pred_{res_arc}_{reg}_fixed.zarr/'),
                                  consolidated = True, mode = 'w')

        # Using fish_mort_det as a base for detritivore mortality
        fish_mort_det = xr.DataArray(
            np.repeat(gridded_python['fish_mort_detritivore'], 
                      np.prod(consume_benthos.shape)).
            reshape(consume_benthos.shape),
            dims = ['lat', 'lon', 'size_class'],
            coords = {'lat': consume_benthos.lat, 
                      'lon': consume_benthos.lon,
                      'size_class': consume_benthos.size_class})
        
        # Changing mortality to zero outside the size fished
        fishing_mort_det = xr.where(
            (fish_mort_det.size_class >= log10_ind_min_fish_det) &
            (fish_mort_det.size_class < fish_mort_det.size_class.max()),
            fish_mort_det, 0).where(np.isfinite(depth))
        fishing_mort_det.name = 'fish_mort_det'
        fishing_mort_det.to_zarr(os.path.join(
            out_folder,
            f'fish-mort-det_{res_arc}_{reg}_fixed.zarr/'),
                                 consolidated = True, mode = 'w')

        # Temperature effects ----
        # Pelagics
        sst = xr.concat([sea_surf_temp_stable_spin,
                          sea_surf_temp_spinup, sea_surf_temp],
                        dim = 'time')
        pel_tempeffect = np.exp(
            gridded_python['c1']-gridded_python['activation_energy']/
            (gridded_python['boltzmann']*(sst+273)))
        pel_tempeffect.name = 'pel_temp_eff' 
        pel_tempeffect_out = pel_tempeffect.chunk({'lon': -1, 'lat': -1,
                                                   'time': -1}).load()
        pel_tempeffect_out.to_zarr(os.path.join(
            out_folder, 
            f'pel-temp-eff_spinup_obsclim_{res_arc}_{reg}_monthly{fn_search}_1741_2010.zarr/'),
            consolidated = True, mode = 'w')

        #Benthics
        sbt = xr.concat([sea_floor_temp_stable_spin, sea_floor_temp_spinup,
                         sea_floor_temp], dim = 'time')
        ben_tempeffect = np.exp(
            gridded_python['c1']-gridded_python['activation_energy']/
            (gridded_python['boltzmann']*(sbt+273)))
        ben_tempeffect.name = 'ben_temp_eff'
        ben_tempeffect_out = ben_tempeffect.chunk({'lon': -1, 'lat': -1,
                                                   'time': -1}).load()
        ben_tempeffect_out.to_zarr(os.path.join(
            out_folder, 
            f'ben-temp-eff_spinup_obsclim_{res_arc}_{reg}_monthly{fn_search}_1741_2010.zarr/'),
            consolidated = True, mode = 'w')





