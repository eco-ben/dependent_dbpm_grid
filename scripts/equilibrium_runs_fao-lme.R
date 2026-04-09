# Equilibrium runs

# Choose local R library
.libPaths("/g/data/vf71/la6889/R_personal_lib/")

# Loading libraries -------------------------------------------------------
source("scripts/useful_functions.R")
library(dplyr)
library(arrow)
library(jsonlite)
library(purrr)
library(furrr)
library(GGally)
library(tibble)


# Loading DBPM climate and fishing inputs ---------------------------------
base_folder <- "/g/data/vf71/fishmip_inputs/ISIMIP3a/fao_lme_inputs"
out_folder <- "/g/data/vf71/fishmip_outputs/ISIMIP3a/fao_lme_outputs"

fao_lme <- list.dirs(base_folder, recursive = F, full.names = F) |> 
  str_subset(pattern = "fao_lme-")

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

# Testing no-fishing and search volume from 0.64 to 640
fish_param <- data.frame("region" = NA, "fmort_u" = 0, "fmort_v" = 0,
                         "fminx_u" = 0, "fminx_v" = 0, 
                         "search_vol" = 64*seq(0.1, 10, by = 0.1))

# Parallelise using workers available using parallelly
plan(multisession, workers = availableCores())

# Initial calibration loop ------------------------------------------------
for(f in fao_lme){
  fish_param <- fish_param |> 
    mutate(region = f)
  
  # Load inputs (weighted means) - We will only use 1 deg inputs to calculate 
  # fishing parameters
  dbpm_inputs <- read_parquet(list.files(
    file.path(base_folder, f, paste0("monthly_weighted", smoothed)), 
    paste0("dbpm_clim-fish-inputs", fn_search, "_", f), full.names = T)) |> 
    filter(str_detect(scenario, "stable")) |> 
    replicate(2, expr = _, simplify = F) |> 
    bind_rows() |> 
    mutate(time = seq(as_date("1801-01-01"), as_date("2000-12-31"), 
                      by = "month"), 
           year = year(time), month = month(time, label = T, abbr = F))
  
  results_folder <- file.path(out_folder, f, "equilibrium_run")
  
  # If the folder does not exist, create a new one
  if(!dir.exists(results_folder)){
    dir.create(results_folder, recursive = T)
  }
  
  print(paste0("Running DBPM for region: ", f))
  
  # Parallelising equilibrium runs per region of interest
  future_walk(1:nrow(fish_param), ~{
    i <- .x
    # Get parameters ready for calibration run
    params <- sizeparam(dbpm_inputs, fish_param[i,], xmin_consumer_u = -3,
                        xmin_consumer_v = -3)

    # Saving non-spatial parameters
    params |>
      #Ensuring up to 10 decimal places are saved in file
      write_json(file.path(results_folder,
                           paste0("dbpm_size_params_", f, "_searchvol_",
                                  fish_param[i,]$search_vol, ".json")),
                 digits = 10)

    # Run non-spatial DBPM.  This step is necessary to get the initial
    # conditions to be used in the gridded DBPM
    init_results <- run_model(fish_param[i,], dbpm_inputs, withinput = F,
                              xmin_consumer_u = -3, xmin_consumer_v = -3,
                              include_plankton = T)

    # Saving initial results for non-spatial run
    init_results |>
      #Ensuring up to 10 decimal places are saved in file
      write_json(file.path(results_folder,
                           paste0("init_dbpm_nonspatial_", f, "_searchvol_",
                                  fish_param[i,]$search_vol, ".json")),
                 digits = 10)

    ## Size spectrum plots per group (predators and detritivores) ---------
    # Transform density matrix to data frame to create plots
    density_df <- dbpm_output_mat_to_df(init_results, dbpm_inputs$time, 
                                        "density")

    # Create size spectrum plots
    size_sp_plot <- plotsizespectrum(density_df, init_results$params, f,
                                     fishing_params = fish_param[i,],
                                     mean_decade = T, nrow = 3)

    # Saving size spectrum plot for non-spatial runs
    ggsave(file.path(results_folder,
                     paste0("size_spectrum_", f, "_searchvol_",
                            fish_param[i,]$search_vol, ".png")),
           size_sp_plot, bg = "white")

    ## Creating growth rate plots -------------------------------------------
    dates_model <- c(min(dbpm_inputs$time)%m-% months(1), dbpm_inputs$time)
    growth_df <- dbpm_output_mat_to_df(init_results, dates_model, "growth") |>
      # Excluding growth value for first time step as it is used for model
      # initialisation only
      filter(time >= min(as_date(dbpm_inputs$time)))


    # Create growth rate plot
    growth_plot <- plot_growth_rate(growth_df, init_results$params, f,
                                    fishing_params = fish_param[i,])

    # Saving size spectrum plot for non-spatial runs
    ggsave(file.path(results_folder,
                     paste0("growth_rates_", f, "_searchvol_",
                            fish_param[i,]$search_vol, ".png")),
           growth_plot, bg = "white")

    #Equilibrium run
    calib_run <- run_model(fish_param[i,], dbpm_inputs, xmin_consumer_u = -3,
                           xmin_consumer_v = -3, include_plankton = T)

    # Save results
    calib_run |>
      write_parquet(file.path(
        results_folder, paste0("dbpm_nonspatial_", f, "_searchvol_",
                                fish_param[i,]$search_vol, "_1801-2000.parquet"))
      )

    # Create plots of biomass (predators, detritivores and detritus)
    biomass_data <- calib_run |>
      select(year, ends_with("biomass"), total_detritus) |>
      group_by(year) |>
      summarise(across(starts_with("total"), ~ mean(.x, na.rm = T)))

    bio_plot <- biomass_data |>
      pivot_longer(!year, names_to = "group", values_to = "values",
                   names_prefix = "total_") |>
      # separate(group, c("group", "type"), sep = "_", fill = "right") |>
      separate_wider_delim(group, delim = "_", names = c("group", "type"),
                           too_few = "align_start") |>
      replace_na(list(type = "detritus")) |>
      ggplot(aes(x = year, y = values, color = group))+
      geom_line()+
      geom_point()+
      facet_grid(type~., scales = "free")+
      labs(title =
             paste0(str_replace_all(str_to_upper(f), "_|-", " "),
                    ": Predator and detritivore biomass, and detritus density"))+
      theme_bw()+
      theme(axis.title.x = element_blank())

    ggsave(file.path(results_folder,
                     paste0("pred-detritus-bio_detritus_", f, "_searchvol_",
                            fish_param[i,]$search_vol, ".png")), bio_plot,
           bg = "white")
  }, .options = furrr_options(seed = T))
}



# Processing outputs per region -------------------------------------------
plan(multisession, workers = availableCores())
for(f in fao_lme){
  results_folder <- file.path(out_folder, f, "equilibrium_run")
  
  init_files <- list.files(results_folder, pattern = "^init_dbpm", 
                           full.names = T)
  
  calib_files <- list.files(results_folder, pattern = "^dbpm_nonspatial_", 
                           full.names = T)
  
  print(paste0("Running DBPM for region: ", f))
  
  # Parallelising output processing within regions of interest
  future_walk(1:length(init_files), ~{
    i <- .x
    init_results <- read_json(init_files[i], simplifyVector = T)
    calib_run <- read_parquet(calib_files[i]) |> 
      select(time, year, ends_with("biomass"), total_detritus)
    min_pred_size <- init_results$params$log10_size_bins[
      init_results$params$ind_min_pred_size]
    min_det_size <- init_results$params$log10_size_bins[
      init_results$params$ind_min_detritivore_size]
    
    # Processing density data
    density_df <- dbpm_output_mat_to_df(init_results, calib_run$time, 
                                        "density") |> 
      select(!detritus) |> 
      mutate(predators = ifelse(size_class < min_pred_size, NA, predators),
             detritivores = ifelse(size_class < min_det_size, NA, 
                                   detritivores)) |> 
      drop_na(detritivores, predators) |> 
      rowwise() |> 
      mutate(total = sum(predators, detritivores, na.rm = T)) |>
      ungroup() |> 
      summarise(across(c(predators, detritivores, total), 
                       ~ mean(.x, na.rm = T)), 
                .by = c(decade, size_class, search_vol)) |> 
      mutate(across(c(predators, detritivores, total), log10)) |>
      pivot_longer(c(predators, detritivores, total), names_to = "group", 
                   values_to = "bio")
    density_df |> 
      write_parquet(file.path(results_folder,
                              paste0("size_spectrum_data_", f, "_searchvol_",
                                     unique(density_df$search_vol), ".parquet")))
   
    # Processing growth rate data
    dates_model <- c(min(calib_run$time)%m-% months(1), calib_run$time)
    growth_df <- dbpm_output_mat_to_df(init_results, dates_model, "growth") |>
      mutate(size_class = 10**size_class) |>
      # Excluding growth value for first time step as it is used for model
      # initialisation only
      filter(time >= min(as_date(calib_run$time)) & size_class >= 0.1 & 
               size_class <= 10**5) |> 
      summarise(across(c(predators, detritivores), ~ mean(.x, na.rm = T)),
                .by = c(decade, size_class, search_vol)) |> 
      pivot_longer(c(predators, detritivores), names_to = "group", 
                   values_to = "growth")
    
    growth_df |> 
      write_parquet(file.path(results_folder,
                              paste0("growth_rates_data", f, "_searchvol_",
                                     unique(growth_df$search_vol), ".parquet")))
  
    # Processing biomass data
    biomass_data <- calib_run |>
      summarise(across(starts_with("total"), ~ mean(.x, na.rm = T)), 
                .by = year) |>
      pivot_longer(!year, names_to = "group", values_to = "values",
                   names_prefix = "total_") |>
      separate_wider_delim(group, delim = "_", names = c("group", "type"),
                           too_few = "align_start") |>
      replace_na(list(type = "detritus")) |> 
      mutate(search_vol = unique(density_df$search_vol))
    
    biomass_data |> 
      write_parquet(file.path(
        results_folder, paste0("plankton-pred-detritus-bio_detritus_", f, 
                               "_searchvol_", unique(density_df$search_vol), 
                               ".parquet")))
    }
  )
}


# Creating diagnostic plots per region per decade -------------------------
for(f in fao_lme){
  results_folder <- file.path(out_folder, f, "equilibrium_run")
  
  print(paste0("Running DBPM for region: ", f))
  
  ## Size spectrum plots --------------------------------------------------
  # Loading data
  density_df <- list.files(results_folder, pattern = "size_spectrum_data_", 
                           full.names = T) |> 
    map(\(x) read_parquet(x)) |> 
    bind_rows() |> 
    # Removing rows with no data
    drop_na(bio) 
  
  # Plotting data for last decade
  density_fig <- density_df |> 
    filter(decade == max(decade)) |>
    ggplot(aes(size_class, bio, colour = factor(search_vol), 
               group = search_vol))+
    geom_line(alpha = 0.5)+
    lims(y = c(-20, NA))+
    facet_grid(~group)+
    labs(colour = "Search volume", 
         y = expression("" *log[10] ~ "abundance density (m"^-3* ")"),
         x = expression("" *log[10] ~ "body mass (g)"),
         title = paste0("Last decade of simulation: Size spectrum - ",
                        str_to_upper(str_replace_all(f, "_|-", " "))))+
    theme_bw()
  
  # Saving plot
  ggsave(file.path(results_folder,
                   paste0("size_spectrum_", f, "_searchvol-check", ".png")),
         density_fig, bg = "white")
  
  ## Growth rate plots ----------------------------------------------------
  # Loading data
  growth_df <- list.files(results_folder, pattern = "growth_rates_data", 
                          full.names = T) |> 
    map(\(x) read_parquet(x)) |> 
    bind_rows() |> 
    # Removing rows with no data
    drop_na(growth) |> 
    # Removing rows where growth rate is equal or less than 0
    filter(growth > 0) 
  
  # Plotting data for last decade
  growth_fig <- growth_df |> 
    filter(decade == max(decade)) |> 
    ggplot(aes(size_class, growth, colour = factor(search_vol), 
               group = search_vol))+
    geom_line(alpha = 0.5)+
    scale_y_continuous(trans = "log10", 
                       name = "Relative growth rate per year")+
    scale_x_continuous(trans = "log10", name = "Body mass (g)")+
    facet_grid(~group)+
    labs(colour = "Search volume", 
         title = paste0("Last decade of simulation: Mean growth rate - ",
                        str_to_upper(str_replace_all(f, "_|-", " "))))+
    theme_bw()
  
  # Saving plot
  ggsave(file.path(results_folder,
                   paste0("growth_rates_", f, "_searchvol-check", ".png")),
         growth_fig, bg = "white")
  
  ## Biomass plots --------------------------------------------------------
  # Loading data
  bio_df <- list.files(results_folder, "plankton-pred-detritus-bio_detritus_", 
                       full.names = T) |> 
    map(\(x) read_parquet(x)) |> 
    bind_rows() |> 
    # Removing rows with no data
    drop_na(values) |> 
    # Removing rows where growth rate is equal or less than 0
    filter(values >= 0) 
  
  good_search_vol <- bio_df |> 
    filter(group != "plankton") |> 
    count(search_vol) |> 
    # A total of 600, which includes 200 timesteps and three groups (predators,
    # detritivores and detritus)
    filter(n == 600) 
  
  # Plotting data for entire simulation
  bio_fig <- bio_df |> 
    filter(search_vol %in% good_search_vol$search_vol) |> 
    ggplot(aes(year, values, colour = factor(search_vol), 
               group = search_vol))+
    geom_line(alpha = 0.5)+
    facet_grid(group~., scales = "free")+
    labs(colour = "Search volume", 
         title = paste0("Last decade of simulation: Estimated biomass per group - ",
                        str_to_upper(str_replace_all(f, "_|-", " "))))+
    theme_bw()+
    theme(axis.title = element_blank())
  
  # Saving plot
  ggsave(file.path(results_folder,
                   paste0("plankton-pred-detritus-bio_detritus_", f, 
                          "_searchvol-check", ".png")),
         bio_fig, bg = "white")
}



# Creating diagnostic plots per region for final timestep -----------------
density_growth_global <- tibble()
for(f in fao_lme){
  results_folder <- file.path(out_folder, f, "equilibrium_run")
  
  print(paste0("Running DBPM for region: ", f))
  

  ## Regional inputs ------------------------------------------------------
  f_inputs <- read_parquet(list.files(
    file.path(base_folder, f, paste0("monthly_weighted", smoothed)), 
    paste0("dbpm_clim-fish-inputs", fn_search, "_", f), full.names = T), 
    col_select = scenario:depth) |> 
    filter(time == min(time))
  
  ## Biomass plots --------------------------------------------------------
  # Loading data
  bio_df <- list.files(results_folder,
                       "plankton-pred-detritus-bio_detritus_.*parquet", 
                       full.names = T) |> 
    map(\(x) read_parquet(x)) |> 
    bind_rows() |> 
    # Removing rows with no data
    drop_na(values) |> 
    # Removing rows where biomass is 0 or more
    filter(values >= 0) 
  
  good_search_vol <- bio_df |> 
    filter(group != "plankton") |> 
    count(search_vol) |> 
    # A total of 600, which includes 200 timesteps and three groups (predators,
    # detritivores and detritus)
    filter(n == 600) 
  
  # Get files for successful runs only
  init_files <- list.files(results_folder, 
                           pattern = paste(
                             paste0("^init_dbpm.*_", good_search_vol$search_vol, 
                                    ".json"), collapse = "|"), full.names = T)
  
  
  density_growth_reg <- tibble()
  
  for(fn in init_files){
    init_results <- read_json(fn, simplifyVector = T)
    
    min_pred_size <- init_results$params$log10_size_bins[
      init_results$params$ind_min_pred_size]
    min_det_size <- init_results$params$log10_size_bins[
      init_results$params$ind_min_detritivore_size]
    
    # Processing density data for size spectrum
    density_growth_df <- data.frame(
      size_class = init_results$params$log10_size_bins, 
      pred_den = init_results$predators[, ncol(init_results$predators)], 
      detrit_den = init_results$detritivores[, ncol(init_results$detritivores)],
      # Growth from time step before end of simulation
      pred_growth = init_results$growth_int_pred[, ncol(init_results$growth_int_pred)-1], 
      detrit_growth = init_results$growth_det[, ncol(init_results$growth_det)-1], 
      search_vol = init_results$params$hr_volume_search) |> 
      mutate(size_class_g = 10**size_class, .before = pred_growth) |>
      mutate(pred_den = ifelse(size_class < min_pred_size, NA, pred_den),
             detrit_den = ifelse(size_class < min_det_size, NA, detrit_den),
             pred_growth = ifelse(size_class < min_pred_size, NA, pred_growth),
             detrit_growth = ifelse(size_class < min_det_size, NA, detrit_growth)) |> 
      drop_na(pred_den:detrit_growth) |> 
      rowwise() |> 
      mutate(total = sum(pred_den, detrit_den, na.rm = T), 
             .after = detrit_den) |>
      ungroup() |> 
      mutate(across(c(pred_den, detrit_den, total), log10)) 
  
    density_growth_reg <- density_growth_reg |> 
      bind_rows(density_growth_df)
  }
  
  # Binding all search volumes together within a single region
  density_growth_reg <- f_inputs |> 
    bind_cols(density_growth_reg)
  
  density_growth_reg |> 
    write_parquet(
      file.path(results_folder, 
                paste0("successful_runs_inputs-outputs_", f, ".parquet")))
}


# Binding all search volumes together within a single region
density_growth_global <- density_growth_global |> 
  bind_rows(density_growth_df)

density_growth_global |> 
  write_parquet(file.path(out_folder, 
                          "successful_runs_inputs-outputs_global.parquet"))


# density_growth_global |> 
#   summarise(min_sv = min(search_vol), max_sv = max(search_vol), .by = region)
  




