
# Loading libraries -------------------------------------------------------
source("scripts/useful_functions.R")
library(dplyr)
library(arrow)
library(jsonlite)
library(purrr)
library(tibble)

# Loading DBPM climate and fishing inputs ---------------------------------
base_folder <- "/g/data/vf71/fishmip_inputs/ISIMIP3a/fao_inputs"
fao <- list.dirs(base_folder, recursive = F, full.names = F) |> 
  str_subset(pattern = "fao-")

# Resolution
resolutions <- c("025deg", "1deg")

# Setting "smoothed" to TRUE will use 'smoothed' instead of original GFDL 
# outputs to force DBPM
smoothed <- TRUE
if(smoothed){
  fn_search <- "-smoothed"
}else{
  fn_search <- ""
}

for(f in fao){
  for(res in resolutions){
    dbpm_inputs <- file.path(base_folder, f, "monthly_weighted", res,
                             paste0("dbpm_clim-fish-inputs", fn_search, "_", f, 
                                    "_1841-2010.parquet")) |> 
      read_parquet()
    
    # Searching best fishing parameters values for area of interest -----------
    #Path to folder where results will be stored
    results_folder <- file.path(base_folder, f, "fishing_params", res,
                                paste0("best_fish_vals", fn_search))
    #Number of iterations
    no_iter <- 100
    params_calibration <- LHSsearch(num_iter = no_iter,
                                    forcing_file = dbpm_inputs, 
                                    gridded_forcing = NULL, 
                                    best_val_folder = results_folder, 
                                    best_param = F, new_detritus_calc = F) |> 
      rowid_to_column("id")
    
    ## Creating plots with fishing parameters calculated above --------------
    # Calculate errors and correlations with tuned fishing parameters and save 
    # plot
    params_corr <- params_calibration |> 
      split(params_calibration$id) |>
      map_df(\(x) getError(x, dbpm_inputs, corr = T, new_detritus_calc = F)) 
    
    #Adding correlation to fishing parameter data frame
    params_calibration <- params_calibration |> 
      #Removing column to avoid duplication
      select(!region) |> 
      bind_cols(params_corr) |> 
      select(!id) |> 
      #Remove any rows where simulation returned NA values
      filter(catchNA == 0) |> 
      arrange(desc(cor), rmse) |> 
      relocate(region, .before = fmort_u)
    
    #Saving results
    params_calibration |> 
      write_parquet(file.path(
        results_folder,
        paste0("best-fishing-parameters_", f, 
               "_searchvol_estimated_numb-iter_", no_iter,".parquet")))
    
    # Calibration plots to be done after all fishing parameters are calculated
    #Filter best fishing parameters
    good <- params_calibration |> 
      filter(cor >= 0.5) |> 
      filter(rmse == min(rmse))
    #If nothing is returned, then use parameters for lowest rmse
    if(nrow(good) == 0){
      good <- params_calibration |> 
        filter(rmse == min(rmse))
    }
    #Create plot with best performing parameters
    good |> 
      corr_calib_plots(dbpm_inputs, results_folder, new_detritus_calc = F)
    
    # Calibration plots with parameters that had highest correlation values
    params_calibration |> 
      slice(1) |>
      corr_calib_plots(dbpm_inputs, file.path(results_folder, "high_corr"),
                       new_detritus_calc = F)
    
  }
}


## Optimising underperforming regions --------------------------------------
# This section may need to be ran multiple times until all regions have good
# parameters

# Load all fishing parameter files for each resolution
for(res in resolutions){
  # Getting a list of files containing fishing parameters calculated for all 
  # regions
  fish_param <- fao |> 
    map_chr(\(x) list.files(file.path(base_folder, x, "fishing_params", res,
                                      paste0("best_fish_vals", fn_search)),
                            pattern = "best-fishing-parameters", recursive = T,
                            full.names = T)) |> 
    map(~read_parquet(.)) |> 
    bind_rows()
  
  # Find the regions for which fishing parameters perform well
  good_params <- fish_param |> 
    group_by(region) |> 
    # Lowest RMSE and correlation must be 0.5 or higher
    filter(cor >= 0.5) |> 
    filter(rmse == min(rmse)) |> 
    mutate(region = str_replace(str_to_lower(region), " ", "-")) |> 
    pull(region)
  
  # Compare to full region list to find underperforming regions
  bad_params <- fao[!fao %in% good_params]
  
  # Since correlation is below 0.5 and the plots comparing estimates and obs do 
  # not look like a great fit, we will calculate fishing parameters again 
  no_iter <- 1000
  
  for(f in bad_params){
    dbpm_inputs <- file.path(base_folder, f, "monthly_weighted", res,
                             paste0("dbpm_clim-fish-inputs", fn_search, "_", f, 
                                    "_1841-2010.parquet")) |> 
      read_parquet()
    
    # Searching best fishing parameters values for area of interest -----------
    #Path to folder where results will be stored
    results_folder <- file.path(base_folder, f, "fishing_params", res,
                                paste0("best_fish_vals", fn_search))
    
    params_calibration_optim <- LHSsearch(num_iter = no_iter, seed = 1234,
                                          forcing_file = dbpm_inputs, 
                                          gridded_forcing = NULL, 
                                          best_val_folder = results_folder, 
                                          best_param = F, 
                                          new_detritus_calc = F) |> 
      rowid_to_column("id")
    
    params_corr <- params_calibration_optim |> 
      split(params_calibration_optim$id) |>
      map_df(\(x) getError(x, dbpm_inputs, corr = T, new_detritus_calc = F))
    
    #Adding correlation to fishing parameter data frame
    params_calibration_optim <- params_calibration_optim |> 
      select(!region) |> 
      bind_cols(params_corr) |> 
      select(!id) |> 
      filter(catchNA == 0) |> 
      arrange(desc(cor), rmse) |> 
      relocate(region, .before = fmort_u)
    
    #Saving results
    params_calibration_optim |> 
      write_parquet(file.path(
        results_folder,
        paste0("best-fishing-parameters_", f, 
               "_searchvol_estimated_numb-iter_", no_iter, ".parquet")))
    
    #Identifying best performing parameters
    good <- params_calibration_optim |> 
      filter(cor >= 0.5) |> 
      filter(rmse == min(rmse))
    #If nothing is returned, select lowest RMSE
    if(nrow(good) == 0){
      good <- params_calibration |> 
        filter(rmse == min(rmse))
    }
    
    #Create calibration plot
    corr_calib_plots(good, dbpm_inputs, results_folder, new_detritus_calc = F)
    
    #Create plot with parameters that resulted in highest correlation
    params_calibration_optim |>
      filter(cor == max(corr)) |> 
      corr_calib_plots(dbpm_inputs, file.path(results_folder, "high_corr"),
                       new_detritus_calc = F)
    
  }
}


# Getting DBPM parameters -------------------------------------------------
out_folder <- "/g/data/vf71/fishmip_outputs/ISIMIP3a/fao_outputs"

for(res in resolutions){
  # Getting a list of files containing fishing parameters calculated for all 
  # regions
  fishing_params <- fao |> 
    map_chr(\(x) list.files(file.path(base_folder, x, "fishing_params", res,
                                      paste0("best_fish_vals", fn_search)),
                            pattern = "best-fishing-parameters", recursive = T,
                            full.names = T)) |> 
    # Load all fishing parameter files for coarser resolution
    map(~read_parquet(.)) |> 
    bind_rows() |> 
    group_by(region) |> 
    # Find parameters with lowest RMSE and correlation of 0.5 or higher
    filter(cor >= 0.5) |> 
    filter(rmse == min(rmse))
  
  #Filtering best fishing parameters
  #Find regions that did not meet the correlation requirement (FAO 21 and 58)
  # fp <- fishing_params |> 
  #   group_by(region) |> 
  #   # Find parameters with lowest RMSE
  #   filter(rmse == min(rmse) & str_detect(region, "21|58"))

  for(f in fao){
    results_folder <- file.path(base_folder, f, "fishing_params", res,
                                paste0("best_fish_vals", fn_search))
    
    dbpm_inputs <- file.path(base_folder, f, "monthly_weighted", res,
                             paste0("dbpm_clim-fish-inputs", fn_search, "_", f, 
                                    "_1841-2010.parquet")) |> 
      read_parquet()
    
    fish_param <- fishing_params |> 
      filter(region == str_replace(str_to_upper(f), "-", " "))
    
    params <- sizeparam(dbpm_inputs, fish_param, xmin_consumer_u = -3, 
                        xmin_consumer_v = -3)
    
    # Saving non-spatial parameters
    params |> 
      #Ensuring up to 10 decimal places are saved in file
      write_json(file.path(results_folder, 
                           paste0("dbpm_size_params_", f, ".json")), 
                 digits = 10)
    
    # Run non-spatial DBPM.  This step is necessary to get the initial 
    # conditions to be used in the gridded DBPM
    init_results <- run_model(fish_param, dbpm_inputs, withinput = F, 
                              new_detritus_calc = F)
    
    # Prepare fishing parameters for gridded DBPM 
    pred_initial <- rowMeans(init_results$predators)
    detritivore_initial <- rowMeans(init_results$detritivores)
    detritus_initial <- mean(init_results$detritus)
    
    gridded_params <- sizeparam(dbpm_inputs, fish_param, xmin_consumer_u = -3, 
                                xmin_consumer_v = -3, use_init = T, 
                                pred_initial = pred_initial, 
                                detritivore_initial = detritivore_initial, 
                                detritus_initial = detritus_initial,
                                gridded = T)
    
    #Save for use in gridded DBPM (step 05)
    gridded_params |> 
      write_json(file.path(results_folder, 
                           paste0("dbpm_gridded_size_params_", f, ".json")),
                 digits = 10)
    
    # Defining folder to save non-spatial results
    dbpm_out_folder <- file.path(out_folder, f, 
                                 paste0("fishing_runs", fn_search),
                                 "nonspatial")
    if(!dir.exists(dbpm_out_folder)){
      dir.create(dbpm_out_folder, recursive = T)
    }
    
    # Running non-spatial DBPM and saving results - This step is needed only 
    # once
    fout <- file.path(dbpm_out_folder, 
                      paste0("dbpm_nonspatial_", f, "_1841-2010.parquet"))
    if(!file.exists(fout)){
      run_model(fish_param, dbpm_inputs, new_detritus_calc = F) |> 
        write_parquet(fout)
    }
  }
}


    