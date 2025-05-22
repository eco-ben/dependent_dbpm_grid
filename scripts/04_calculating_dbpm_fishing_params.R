
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

for(f in fao){
  for(res in resolutions){
    dbpm_inputs <- file.path(base_folder, f, "monthly_weighted", res,
                             paste0("dbpm_clim-fish-inputs_", f, 
                                    "_1841-2010.parquet")) |> 
      read_parquet()
    
    # Searching best fishing parameters values for area of interest -----------
    #Path to folder where results will be stored
    results_folder <- file.path(base_folder, f, "fishing_params", res,
                                "best_fish_vals")
    #Number of iterations
    no_iter <- 100
    params_calibration <- LHSsearch(num_iter = no_iter, 
                                    forcing_file = dbpm_inputs, 
                                    gridded_forcing = NULL, 
                                    best_val_folder = results_folder, 
                                    best_param = F) |> 
      rowid_to_column("id")
    
    ## Creating plots with fishing parameters calculated above --------------
    # Calculate errors and correlations with tuned fishing parameters and save 
    # plot
    params_corr <- params_calibration |> 
      split(params_calibration$id) |>
      map_df(\(x) getError(x, dbpm_inputs, corr = T)) 
    
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
      filter(rmse == min(rmse) & cor >= 0.5)
    #If nothing is returned, then use parameters for lowest rmse
    if(nrow(good) == 0){
      good <- params_calibration |> 
        filter(rmse == min(rmse))
    }
    #Create plot with best performing parameters
    good |> 
      arrange(rmse) |> 
      slice(1) |> 
      corr_calib_plots(dbpm_inputs, results_folder)
    
    # Calibration plots with parameters that had highest correlation values
    params_calibration |> 
      slice(1) |>
      corr_calib_plots(dbpm_inputs, file.path(results_folder, "high_corr"))
    
  }
}


## Optimising underperforming regions --------------------------------------
# This section may need to be ran multiple times until all regions have good
# parameters

# Getting a list of files containing fishing parameters calculated for all 
# regions
fish_param_files <- list.files(base_folder, pattern = "parquet$", 
                               recursive = T, full.names = T) |> 
  str_subset("fishing_params") 

# Load all fishing parameter files for each resolution
for(res in resolutions){
  fish_param <- fish_param_files |> 
    str_subset(res) |> 
    map(~read_parquet(.)) |> 
    bind_rows()
  
  # Find the regions for which fishing parameters perform well
  good_params <- fish_param |> 
    group_by(region) |> 
    # Lowest RMSE and correlation must be 0.5 or larger
    filter(rmse == min(rmse) & cor >= 0.5) |> 
    mutate(region = str_replace(str_to_lower(region), " ", "-")) |> 
    pull(region)
  
  # Compare to full region list to find underperforming regions
  bad_params <- fao[!fao %in% good_params]
  
  # Since correlation is below 0.5 and the plots comparing estimates and obs do 
  # not look like a great fit, we will calculate fishing parameters again 
  no_iter <- 1000
  
  for(f in bad_params){
    dbpm_inputs <- file.path(base_folder, f, "monthly_weighted", res,
                             paste0("dbpm_clim-fish-inputs_", f, 
                                    "_1841-2010.parquet")) |> 
      read_parquet()
    
    # Searching best fishing parameters values for area of interest -----------
    #Path to folder where results will be stored
    results_folder <- file.path(base_folder, f, "fishing_params", res,
                                "best_fish_vals")
    
    params_calibration_optim <- LHSsearch(num_iter = no_iter, seed = 1234,
                                          forcing_file = dbpm_inputs, 
                                          gridded_forcing = NULL, 
                                          best_val_folder = results_folder, 
                                          best_param = F) |> 
      rowid_to_column("id")
    
    params_corr <- params_calibration_optim |> 
      split(params_calibration_optim$id) |>
      map_df(\(x) getError(x, dbpm_inputs, corr = T))
    
    #Adding correlation to fishing parameter data frame
    params_calibration_optim <- params_calibration_optim |> 
      select(!region) |> 
      bind_cols(params_corr) |> 
      filter(cor >= 0) |> 
      arrange(desc(cor), rmse) |> 
      relocate(region, .before = fmort_u)
    
    #Saving results
    params_calibration_optim |> 
      write_parquet(file.path(
        results_folder,
        paste0("best-fishing-parameters_", f, 
               "_searchvol_estimated_numb-iter_", no_iter, ".parquet")))
    
    params_calibration_optim |>
      slice(1) |> 
      corr_calib_plots(dbpm_inputs, file.path(results_folder, "high_corr"))
    
    params_calibration_optim |> 
      arrange(rmse) |>
      slice(1) |> 
      corr_calib_plots(dbpm_inputs, results_folder)
  }
}


# Getting DBPM parameters -------------------------------------------------
# Getting a list of files containing fishing parameters calculated for all 
# regions
fish_param_files <- list.files(base_folder, 
                               pattern = "best-fishing-parameters", 
                               recursive = T, full.names = T) |> 
  str_subset("/best_fish_vals/") 


out_folder <- "/g/data/vf71/fishmip_outputs/ISIMIP3a/fao_outputs"

# Load all fishing parameter files for coarser resolution
fishing_params <- fish_param_files |> 
  str_subset("1deg") |> 
  map(~read_parquet(.)) |> 
  bind_rows() 

#Filtering best fishing parameters
#Find regions that did not meet the correlation requirement (FAO 21 and 58)
fp <- fishing_params |> 
  group_by(region) |> 
  # Find parameters with lowest RMSE
  filter(rmse == min(rmse) & str_detect(region, "21|58")) |> 
  select(!id)

#Find parameters for all other regions
fishing_params <- fishing_params |> 
  group_by(region) |> 
  # Find parameters with lowest RMSE and correlation of 0.5 or higher
  filter(rmse == min(rmse) & cor >= 0.5) |> 
  # Merge with previous dataframe
  select(!id) |> 
  bind_rows(fp) |> 
  arrange(region)

#Removing variable not needed
rm(fp)

for(f in fao){
  results_folder <- file.path(base_folder, f, "fishing_params", "1deg",
                              "best_fish_vals")
  # if(!dir.exists(results_folder)){
  #   dir.create(results_folder)
  # }
  
  dbpm_inputs <- file.path(base_folder, f, "monthly_weighted", "1deg",
                           paste0("dbpm_clim-fish-inputs_", f, 
                                  "_1841-2010.parquet")) |> 
    read_parquet()
  
  fish_param <- fishing_params |> 
    filter(region == str_replace(str_to_upper(f), "-", " "))
  
  params <- sizeparam(dbpm_inputs, fish_param, xmin_consumer_u = -3, 
                      xmin_consumer_v = -3, tstepspryr = 12)
  
  # Saving non-spatial parameters
  params |> 
    #Ensuring up to 10 decimal places are saved in file
    write_json(file.path(results_folder, 
                         paste0("dbpm_size_params_", f, ".json")), 
               digits = 10)
  
  # Run non-spatial DBPM.  This step is necessary to get the initial 
  # conditions to be used in the gridded DBPM
  init_results <- run_model(fish_param, dbpm_inputs, withinput = F)
  
  # Prepare fishing parameters for gridded DBPM 
  pred_initial <- rowMeans(init_results$predators)
  detritivore_initial <- rowMeans(init_results$detritivores)
  detritus_initial <- mean(init_results$detritus)
  
  gridded_params <- sizeparam(dbpm_inputs, fish_param, xmin_consumer_u = -3, 
                              xmin_consumer_v = -3, tstepspryr = 12, 
                              use_init = T, pred_initial = pred_initial, 
                              detritivore_initial = detritivore_initial, 
                              detritus_initial = detritus_initial,
                              gridded = T)
  
  #Save for use in gridded DBPM (step 05)
  gridded_params |> 
    write_json(file.path(results_folder, paste0("dbpm_gridded_size_params_", 
                                             f, ".json")),
               digits = 10)
  
  
  # Running non-spatial DBPM and saving results
  non_spatial_run <- run_model(fish_param, dbpm_inputs)
  
  # Defining folder to save non-spatial results
  dbpm_out_folder <- file.path(out_folder, f, "fishing_runs", "nonspatial")
  if(!dir.exists(dbpm_out_folder)){
    dir.create(dbpm_out_folder, recursive = T)
  }
  
  non_spatial_run |> 
    write_parquet(
      file.path(dbpm_out_folder, 
                paste0("dbpm_nonspatial_", f, "_1841-2010.parquet")))
}

    