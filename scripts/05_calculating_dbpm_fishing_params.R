
# Loading libraries -------------------------------------------------------
source("scripts/useful_functions.R")
library(dplyr)
library(arrow)
library(jsonlite)
library(purrr)
library(GGally)
library(tibble)


# Loading DBPM climate and fishing inputs ---------------------------------
base_folder <- "/g/data/vf71/fishmip_inputs/ISIMIP3a/fao_inputs"
fao <- list.dirs(base_folder, recursive = F, full.names = F) |> 
  str_subset(pattern = "fao-")

# Define search volume
search_volume <- 0.5

# Calibration runs should include fishing?
fishing <- F

# Fishing parameters already created?
params_ready <- F

# Number of iterations (minimum recommend - 500 iterations)
no_iter <- 500

# "smoothed" can be either NULL to use original inputs, 'smoothed' to use LOESS
# smoothed inputs or 'deseasoned' to use deseasoned inputs
# outputs to force DBPM
smoothed <- NULL
if(!is.null(smoothed)){
  fn_search <- "-smoothed"
  smoothed <- paste0("-", smoothed)
}else{
  fn_search <- ""
  smoothed <- ""
}


# Initial calibration loop ------------------------------------------------
for(f in fao){
  # Load inputs (weighted means) - We will only use 1 deg inputs to calculate 
  # fishing parameters
  dbpm_inputs <- file.path(base_folder, f, paste0("monthly_weighted", smoothed),
                           "025deg", paste0("dbpm_clim-fish-inputs", fn_search,
                                          "_", f, "_1841-2010.parquet")) |> 
    read_parquet()
  
  ## Searching best fishing parameters values for area of interest ----------
  #Path to folder where results will be stored
  results_folder <- file.path(base_folder, f, "fishing_params",
                                paste0("best_fish_params", smoothed))
 
  if(!params_ready & fishing){
    params_calibration <- LHSsearch(num_iter = no_iter,
                                    search_volume = search_volume, 
                                    forcing_file = dbpm_inputs,
                                    gridded_forcing = NULL,
                                    best_val_folder = results_folder,
                                    best_param = F) |>
      rowid_to_column("id")
  }else{
    params_calibration <- read_parquet(
      file.path(results_folder, 
                paste0("best-fishing-parameters_", f, "_searchvol_", 
                       search_volume, "_numb-iter_", no_iter, ".parquet"))) |>
      rowid_to_column("id")
  }

  ## Creating plots with fishing parameters calculated above --------------
  # Calculate errors and correlations with tuned fishing parameters and save 
  # plot
  params_corr <- tryCatch(
    {params_calibration |> 
      split(params_calibration$id) |>
      map_df(\(x) getError(x, dbpm_inputs, corr = T))},
    error = function(e){
      NULL},
    warning = function(w){
      NULL},
    finally = {message("Correlation calculation for region ", f, " completed.")}
  )
    
  #Adding correlation to fishing parameter data frame
  if(!is.null(params_corr)){
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
    if(nrow(params_calibration) > 0){
      params_calibration |> 
        write_parquet(
          file.path(results_folder, 
                    paste0("best-fishing-parameters_", f, "_searchvol_", 
                           search_volume, "_numb-iter_", no_iter, ".parquet")))
      
      # Prepare data to create correlation plot
      p1_data <- params_calibration |> 
        select(fmort_u:fminx_v, rmse:cor) |> 
        filter(cor >= 0.2) |>
        mutate(qc = ifelse(cor >= 0.5, "good", "bad"))
        
      # Apply different design to plot based on data available
      if(sum(table(p1_data$qc) <= 1)){
        p1 <- p1_data |> 
          ggpairs(columns = 1:6, aes(alpha = 0.5))+
          labs(title = str_replace(str_to_upper(f), "-", " "),
               subtitle = "Showing parameter values when correlation >= 0.2")
      }else{
        p1 <- p1_data |> 
          ggpairs(columns = 1:6, aes(color = qc, alpha = 0.5))+
          labs(title = str_replace(str_to_upper(f), "-", " "),
               subtitle = "Showing parameter values when correlation >= 0.2")
      }
      
      # Save correlation plot
      p1 |> 
        ggsave(filename = 
                 file.path(results_folder, 
                           paste0("corr_plot_best-fishing-parameters_", f, 
                                  "_searchvol_", search_volume, "_numb-iter_",
                                  no_iter, ".png")))
      
      # Calibration plots to be done after all fishing parameters are calculated
      #Filter best fishing parameters
      good <- params_calibration |> 
        filter(cor >= 0.5) |> 
        filter(rmse == min(rmse))
      #If nothing is returned, then use parameters for lowest rmse
      if(nrow(good) == 0){
        good <- params_calibration |> 
          filter(cor > 0) |> 
          filter(rmse == min(rmse))
      }
      #Create plot with best performing parameters
      good |> 
        corr_calib_plots(dbpm_inputs, results_folder)
      
      # Calibration plots with parameters that had highest correlation values
      params_calibration |> 
        filter(cor == max(cor)) |> 
        corr_calib_plots(dbpm_inputs, file.path(results_folder, "high_corr"))
    }else{
        print(
          paste0("No fishing parameters could be successfully calculated for ", 
                 f))
      }
  }
}


# Optimising underperforming regions --------------------------------------
# This section may need to be ran multiple times until all regions have good
# parameters

# Getting a list of files containing fishing parameters calculated for all 
# regions
fish_param <- fao |> 
  map_chr(\(x) file.path(base_folder, x, "fishing_params",
                         paste0("best_fish_params", smoothed))) |> 
  list.files(pattern = "^best-fishing-parameters", recursive = T, 
             full.names = T) |> 
  str_subset(paste0("searchvol_", search_volume)) |> 
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
no_iter <- 2000

for(f in bad_params){
  dbpm_inputs <- file.path(base_folder, f, paste0("monthly_weighted", smoothed),
                           "025deg", paste0("dbpm_clim-fish-inputs", fn_search,
                                            "_", f, "_1841-2010.parquet")) |> 
    read_parquet()
  
  ## Searching best fishing parameters values for area of interest ----------
  #Path to folder where results will be stored
  results_folder <- file.path(base_folder, f, "fishing_params",
                              paste0("best_fish_params", smoothed))
  
  params_calibration_optim <- LHSsearch(num_iter = no_iter, 
                                        search_volume = search_volume,
                                        seed = 42,
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
    select(!id) |> 
    filter(catchNA == 0) |> 
    arrange(desc(cor), rmse) |> 
    relocate(region, .before = fmort_u)
  
  #Saving results
  params_calibration_optim |> 
    write_parquet(file.path(
      results_folder,
      paste0("best-fishing-parameters_", f, "_searchvol_", search_volume, 
             "_numb-iter_", no_iter, ".parquet")))
  
  # Prepare data to create correlation plot
  p1_data <- params_calibration_optim |> 
    select(fmort_u:fminx_v, rmse:cor) |> 
    mutate(qc = ifelse(cor >= 0.5, "good", "bad"))
  
  # Apply different design to plot based on data available
  if(sum(table(p1_data$qc) <= 1)){
    p1 <- p1_data |> 
      ggpairs(columns = 1:6, aes(alpha = 0.5))+
      labs(title = str_replace(str_to_upper(f), "-", " "))
  }else{
    p1 <- p1_data |> 
      ggpairs(columns = 1:6, aes(color = qc, alpha = 0.5))+
      labs(title = str_replace(str_to_upper(f), "-", " "))
  }
  
  # Save correlation plot
  p1 |> 
    ggsave(filename = 
             file.path(results_folder, 
                       paste0("corr_plot_best-fishing-parameters_", f, 
                              "_searchvol_", search_volume, "_numb-iter_",
                              no_iter, ".png")))
  
  #Identifying best performing parameters
  good <- params_calibration_optim |> 
    filter(cor >= 0.5) |> 
    filter(rmse == min(rmse))
  #If nothing is returned, select lowest RMSE
  if(nrow(good) == 0){
    good <- params_calibration_optim |> 
      filter(cor > 0) |> 
      filter(rmse == min(rmse))
  }
  
  #Create calibration plot
  corr_calib_plots(good, dbpm_inputs, results_folder)

  #Create plot with parameters that resulted in highest correlation
  params_calibration_optim |>
    filter(cor == max(cor)) |> 
    corr_calib_plots(dbpm_inputs, file.path(results_folder, "high_corr"))
}


# Running DBPM calibration (non-spatial runs) -----------------------------
out_folder <- "/g/data/vf71/fishmip_outputs/ISIMIP3a/fao_outputs"

# Start calibration runs
for(f in fao){
  if(fishing){
    results_folder <- file.path(base_folder, f, "init_fish_vals", 
                                paste0("best_fish_vals", smoothed))
  }else{
    results_folder <- file.path(base_folder, f, "init_fish_vals_no-fishing", 
                                paste0("best_fish_vals", smoothed))
  }
  
  # If the folder does not exist, create a new one
  if(!dir.exists(results_folder)){
    dir.create(results_folder, recursive = T)
  }
  
  dbpm_inputs <- file.path(base_folder, f, paste0("monthly_weighted", smoothed),
                           "025deg", paste0("dbpm_clim-fish-inputs", fn_search,
                                            "_", f, "_1841-2010.parquet")) |> 
    read_parquet()
  
  if(fishing){
    fishing_params <- read_parquet(
      list.files(file.path(base_folder, f, "fishing_params",
                           paste0("best_fish_params", smoothed)), 
                 "^best-fishing-parameters_", full.names = T))
      
    # Find parameters with lowest RMSE and correlation of 0.5 or higher          
    fish_param <- fishing_params |> 
      filter(cor >= 0.5) |> 
      filter(rmse == min(rmse))
    # If above condition is not met, then find the lowest RMSE but consider 
    # results from runs with positive correlations
    if(nrow(fish_param) == 0){
      fish_param <- fishing_params |> 
        filter(cor > 0) |> 
        filter(rmse == min(rmse))
    }
  }else if(!fishing){
    if(is.numeric(search_volume)){
      fish_param <- data.frame("region" = str_replace(str_to_upper(f),
                                                      "-", " "),
                               "fmort_u" = 0, "fmort_v" = 0, "fminx_u" = 0, 
                               "fminx_v" = 0, "search_vol" = search_volume)
    }else if(search_volume == "estimated"){
      fish_param <- data.frame("region" = str_replace(str_to_upper(f), 
                                                      "-", " "),
                               "fmort_u" = 0, "fmort_v" = 0, "fminx_u" = 0, 
                               "fminx_v" = 0, 
                               "search_vol" = lhs::randomLHS(1, 1))
    }
  }
  
  # Get parameters ready for calibration run
  params <- sizeparam(dbpm_inputs, fish_param, xmin_consumer_u = -3, 
                      xmin_consumer_v = -3)
  
  # Saving non-spatial parameters
  params |> 
    #Ensuring up to 10 decimal places are saved in file
    write_json(file.path(results_folder, 
                         paste0("dbpm_size_params_", f, ".json")), 
               digits = 10)
  
  
  # Defining folder to save non-spatial results
  if(fishing){
    dbpm_out_folder <- file.path(out_folder, f, 
                                 paste0("fishing_runs", smoothed),
                                 "nonspatial")
  }else{
    dbpm_out_folder <- file.path(out_folder, f, 
                                 paste0("no-fishing_runs", smoothed),
                                 "nonspatial")
  }
  
  if(!dir.exists(dbpm_out_folder)){
    dir.create(dbpm_out_folder, recursive = T)
  }
  
  # Run non-spatial DBPM.  This step is necessary to get the initial 
  # conditions to be used in the gridded DBPM
  init_results <- run_model(fish_param, dbpm_inputs, withinput = F, 
                            xmin_consumer_u = -3, xmin_consumer_v = -3)
  
  # Saving initial results for non-spatial run
  init_results |> 
    #Ensuring up to 10 decimal places are saved in file
    write_json(file.path(dbpm_out_folder, 
                         paste0("init_dbpm_nonspatial_", f, ".json")), 
               digits = 10)
  
  
  ## Size spectrum plots per group (predators and detritivores) ---------
  # Transform density matrix to data frame to create plots
  density_df <- dbpm_density_mat_to_df(init_results, dbpm_inputs$time)
  
  # Create size spectrum plots
  size_sp_plot <- plotsizespectrum(density_df, init_results$params, f, 
                                   fishing_params = fish_param, mean_decade = T,
                                   nrow = 3)
  
  # Saving size spectrum plot for non-spatial runs
  ggsave(file.path(dbpm_out_folder, paste0("size_spectrum_", f, ".png")), 
         size_sp_plot, bg = "white")
  

  ## Combined growth rate plots -------------------------------------------
  # Getting relevant size bins to be used in plots
  size_bins <- params$log10_size_bins[
    params$ind_min_pred_size:params$numb_size_bins]
  
  # Adding growth rates for predators and detritivores and excluding first 
  # column as it is used for initialisation only
  tot_growth <- as.data.frame((init_results$growth_det+
                                 init_results$growth_int_pred), 
                              row.names = params$log10_size_bins)[,-1]
  # Adding timesteps as column names
  colnames(tot_growth) <- dbpm_inputs$time
  
  # Reorganise data before it can be used in plots
  tot_growth <- tot_growth |>
    rownames_to_column("size_class") |> 
    pivot_longer(!size_class, names_to = "time", values_to = "growth_rate") |> 
    mutate(time = as_date(time), decade = ((year(time)-1) %/% 10)*10,
           size_class = as.numeric(size_class), .before = growth_rate) |> 
    filter(size_class %in% size_bins) |> 
    filter(decade >= 1960) |> 
    group_by(decade, size_class) |> 
    summarise(mean_growth = mean(growth_rate, na.rm = T)) 
  
  # Create growth rate plot
  growth_plot_all <- tot_growth |> 
    ggplot(aes(size_class, mean_growth))+
    geom_line()+
    scale_y_continuous(trans = "log10", name = "Relative growth rate per year")+
    labs(x = expression("" *log[10] ~ "body mass (g)"),
         title = paste0(
           unique(dbpm_inputs$region),
           ": Mean growth rate per decade\n(detritivores+predators)"))+
    facet_wrap(~decade)+
    theme_bw()+
    theme(panel.grid.minor = element_blank(), 
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  
  # Saving size spectrum plot for non-spatial runs
  ggsave(file.path(dbpm_out_folder, 
                   paste0("growth_rates_combined_", f, ".png")), 
         growth_plot_all, bg = "white")
  
  
  ## Prepare fishing parameters for gridded DBPM ---------------------------
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


  ## Running non-spatial DBPM and saving results --------------------------
  fout <- file.path(dbpm_out_folder, 
                    paste0("dbpm_nonspatial_", f, "_1841-2010.parquet"))

  run_model(fish_param, dbpm_inputs, xmin_consumer_u = -3, 
            xmin_consumer_v = -3) |> 
    write_parquet(fout)
}

    