# Running no-fishing simulations for entire globe using equilibrium results
# as initial biomass values

# Choose local R library
.libPaths("/g/data/vf71/la6889/R_personal_lib/")

# Loading libraries -------------------------------------------------------
source("scripts/useful_functions.R")
library(dplyr)
library(arrow)
library(jsonlite)
library(purrr)
# library(GGally)
# library(tibble)


# Loading DBPM climate and fishing inputs ---------------------------------
base_folder <- "/g/data/vf71/fishmip_inputs/ISIMIP3a/fao_lme_inputs"
fao_lme <- list.dirs(base_folder, recursive = F, full.names = F) |> 
  str_subset(pattern = "fao_lme-")

# Define maxium search volume that was successful in all LME-FAO regions
search_volume <- 12.8

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

# Running DBPM calibration (non-spatial runs) -----------------------------
out_folder <- "/g/data/vf71/fishmip_outputs/ISIMIP3a/fao_lme_outputs"

# Start calibration runs
for(f in fao){
  if(f == "fao_lme-64"){
    results_folder <- file.path(out_folder, f, "dynamic_equilibrium_run")
  }else{
    results_folder <- file.path(out_folder, f, "equilibrium_run")
  }
  
  # Print message to console to keep track of region being processed
  print(paste0("Running DBPM for region: ", f))
  
  # Load results from equilibrium run
  init_results <- read_json(list.files(
    results_folder, pattern = paste0("^init_dbpm.*_", search_volume, ".json"),
    full.names = T), simplifyVector = T)
  
  # Load DBPM inputs - Excluding inputs from stable spinup
  dbpm_inputs <- read_parquet(
    file.path(base_folder, f, paste0("monthly_weighted", smoothed),
              paste0("dbpm_clim-fish-inputs", fn_search, "_", f, 
                     "_1641-2010.parquet"))) |> 
    filter(scenario != "stable-spin")
  
  fish_param <- data.frame("region" = str_replace(str_to_upper(f), "-", " "),
                           "fmort_u" = 0, "fmort_v" = 0, "fminx_u" = 0, 
                           "fminx_v" = 0, "search_vol" = search_volume)
  
  ts_stablespin <- init_results$params$numb_time_steps
  
  # Selecting values for last time step in stable spin period
  # Transforming predators and detritivores from gm-2yr-1 to gm-3yr-1
  pred_initial <- init_results$predators[,ts_stablespin]/
    init_results$params$depth
  detritivore_initial <- init_results$detritivores[,ts_stablespin]/20
  detritus_initial <- init_results$detritus[ts_stablespin]
  
  # Set DBPM params for spinup and modelling period
  no_fishing_run <- run_model(fish_param, dbpm_inputs, withinput = F, 
                              xmin_consumer_u = -3, xmin_consumer_v = -3, 
                              pred_initial = pred_initial, 
                              detritivore_initial = detritivore_initial,
                              detritus_initial = detritus_initial, 
                              set_plankton = T, use_init = T)
  
  
  density_df <- dbpm_output_mat_to_df(no_fishing_run, dbpm_inputs$time,
                                      "density")
  
  # Create size spectrum plots
  size_sp_plot <- plotsizespectrum(density_df, no_fishing_run$params, f, 
                                   fishing_params = fish_param, mean_decade = T,
                                   nrow = 3)
  
  size_sp_plot2 <- plotsizespectrum(density_df, no_fishing_run$params, f, 
                                   fishing_params = fish_param, mean_decade = F)
  
  
  ## Creating growth rate plots -------------------------------------------
  dates_model <- c(min(dbpm_inputs$time)%m-% months(1), dbpm_inputs$time)
  growth_df <- dbpm_output_mat_to_df(no_fishing_run, dates_model, "growth") |>
    # Excluding growth value for first time step as it is used for model
    # initialisation only
    filter(time >= min(as_date(dbpm_inputs$time)))
  
  # Create growth rate plot
  growth_plot <- growth_df |> 
    filter(decade == max(decade)) |> 
    plot_growth_rate(no_fishing_run$params, f, fishing_params = fish_param)
  
  growth_df |> 
    filter(time != max(time)) |> 
    filter(time == max(time)) |> 
    plot_growth_rate(no_fishing_run$params, f, fishing_params = fish_param)
  
  no_fishing_run_bio <- run_model(fish_param, dbpm_inputs, withinput = T, 
                                  xmin_consumer_u = -3, xmin_consumer_v = -3, 
                                  include_plankton = T, 
                                  pred_initial = pred_initial,
                                  detritivore_initial = detritivore_initial,
                                  detritus_initial = detritus_initial, 
                                  set_plankton = T, use_init = T)
  
  biomass_data <- no_fishing_run_bio |>
    select(year, ends_with("biomass"), total_detritus) |>
    group_by(year) |>
    summarise(across(starts_with("total"), ~ mean(.x, na.rm = T))) |>
    pivot_longer(!year, names_to = "group", values_to = "values",
                 names_prefix = "total_") |>
    separate_wider_delim(group, delim = "_", names = c("group", "type"),
                         too_few = "align_start") |>
    filter(group != "plankton") |>
    replace_na(list(type = "detritus"))
  
  bio_plot <- biomass_data |>
    ggplot(aes(x = year, y = values, color = group))+
    geom_line()+
    geom_point()+
    facet_grid(type~., scales = "free")+
    labs(title =
           paste0(str_replace_all(str_to_upper(f), "_|-", " "),
                  ": Predator and detritivore biomass, and detritus density"))+
    theme_bw()+
    theme(axis.title.x = element_blank())
}

  