# Equilibrium no-fishing runs to find best search volume for all LME-FAO regions

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
library(forcats)


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

# Testing no-fishing and search volume values between 0.64 and 640
fish_param <- data.frame("region" = NA, "fmort_u" = 0, "fmort_v" = 0,
                         "fminx_u" = 0, "fminx_v" = 0, 
                         "search_vol" = 64*seq(0.1, 10, by = 0.1))

# Parallelising runs using parallelly and furrr
plan(multisession, workers = availableCores())

# Initial calibration loop ------------------------------------------------
for(f in fao_lme){
  # Define results folder
  results_folder <- file.path(out_folder, f, "equilibrium_run")
  
  # If the folder does not exist, create a new one
  if(!dir.exists(results_folder)){
    dir.create(results_folder, recursive = T)
  }
  
  # Add region name to variable with fishing parameterss
  fish_param <- fish_param |> 
    mutate(region = f)
  
  # Load inputs - using stable spin-up scenario only
  dbpm_inputs <- read_parquet(
    file.path(base_folder, f, paste0("monthly_weighted", smoothed),
              paste0("dbpm_clim-fish-inputs", fn_search, "_", f, 
                     "_1641-2010.parquet"))) |> 
    filter(str_detect(scenario, "stable")) 
  
  # Print message to console to keep track of region being processed
  print(paste0("Running DBPM for region: ", f))
  
  # Parallelising equilibrium runs per region of interest
  future_walk(1:nrow(fish_param), ~{
    i <- .x
    # Get parameters ready for calibration run
    sizeparam(dbpm_inputs, fish_param[i,], xmin_consumer_u = -3,
              xmin_consumer_v = -3) |>
      #Ensuring up to 10 decimal places are saved in file
      write_json(
        file.path(results_folder,
                  paste0("dbpm_size_params_", f, "_searchvol_",
                         fish_param[i,]$search_vol, ".json")), digits = 10)

    # Run non-spatial DBPM.  This step is necessary to get the initial
    # conditions to be used in the gridded DBPM
    init_results <- run_model(fish_param[i,], dbpm_inputs, withinput = F,
                              xmin_consumer_u = -3, xmin_consumer_v = -3,
                              include_plankton = T)

    # Saving initial results for non-spatial run
    init_results |>
      #Ensuring up to 10 decimal places are saved in file
      write_json(file.path(
        results_folder, paste0("init_dbpm_nonspatial_", f, "_searchvol_",
                               fish_param[i,]$search_vol, ".json")), 
        digits = 10)

    ## Size spectrum plots per group (predators and detritivores) ---------
    # Transform density matrix to data frame to create plots
    density_df <- dbpm_output_mat_to_df(init_results, dbpm_inputs$time, 
                                        "density")
    
    den_data <- plotsizespectrum(density_df, init_results$params, f, 
                                 fishing_params = fish_param[i,], 
                                 mean_decade = T, return_data = T)
    
    write_parquet(den_data$data, file.path(
      results_folder, paste0("size_spectrum_data_", f, "_searchvol_",
                             fish_param[i,]$search_vol, ".parquet")))
    
    
    ## Creating growth rate plots -------------------------------------------
    dates_model <- c(min(dbpm_inputs$time)%m-% months(1), dbpm_inputs$time)
    growth_df <- dbpm_output_mat_to_df(init_results, dates_model, "growth") |>
      # Excluding growth value for first time step as it is used for model
      # initialisation only
      filter(time >= min(as_date(dbpm_inputs$time))) 
    
    plot_growth_rate(growth_df, init_results$params, f,
                     fishing_params = fish_param[i,], return_data = T) |> 
      write_parquet(file.path(
        results_folder, paste0("growth_rates_data", f, "_searchvol_",
                               fish_param[i,]$search_vol, ".parquet")))
    
    #Equilibrium run
    calib_run <- run_model(fish_param[i,], dbpm_inputs, xmin_consumer_u = -3,
                           xmin_consumer_v = -3, include_plankton = T)

    # Save results
    calib_run |> 
      write_parquet(
        file.path(results_folder, 
                  paste0("dbpm_nonspatial_", f, "_searchvol_", 
                         fish_param[i,]$search_vol, "_1641-1840.parquet")))

    # Create plots of biomass (predators, detritivores and detritus)
    calib_run |>
      select(year, ends_with("biomass"), total_detritus) |>
      group_by(year) |>
      summarise(across(starts_with("total"), ~ mean(.x, na.rm = T))) |> 
      pivot_longer(!year, names_to = "group", values_to = "values",
                   names_prefix = "total_") |>
      separate_wider_delim(group, delim = "_", names = c("group", "type"),
                           too_few = "align_start") |>
      replace_na(list(type = "detritus")) |> 
      mutate(search_vol = fish_param[i,]$search_vol) |> 
      write_parquet(file.path(
        results_folder, paste0("plankton-pred-detritus-bio_detritus_", f, 
                               "_searchvol_", fish_param[i,]$search_vol, 
                               ".parquet")))
  }, .options = furrr_options(seed = T))
}


# Creating diagnostic plots per region per decade -------------------------
for(f in fao_lme){
  results_folder <- file.path(out_folder, f, "equilibrium_run")
  
  print(paste0("Running DBPM for region: ", f))
  
  ## Biomass plots --------------------------------------------------------
  # Loading data
  bio_df <- list.files(results_folder, 
                       "plankton-pred-detritus-bio_detritus_.*parquet", 
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
         title = paste0("Last simulated decade: Estimated biomass per group - ",
                        str_to_upper(str_replace_all(f, "_|-", " "))))+
    theme_bw()+
    theme(axis.title = element_blank())
  
  # Saving plot
  ggsave(file.path(results_folder,
                   paste0("plankton-pred-detritus-bio_detritus_", f, 
                          "_searchvol-check", ".png")),
         bio_fig, bg = "white")
  
  ## Size spectrum plots --------------------------------------------------
  # Loading data
  density_df <- list.files(
    results_folder, pattern = paste0("size_spectrum_data.*", 
                                     good_search_vol$search_vol, 
                                     collapse = "|"), full.names = T) |> 
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
  growth_df <- list.files(
    results_folder, pattern = paste0("growth_rates_data.*", 
                                    good_search_vol$search_vol, 
                                    collapse = "|"), full.names = T) |> 
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
  init_files <- list.files(
    results_folder, 
    pattern = paste0("^init_dbpm.*_", good_search_vol$search_vol, ".json", 
                     collapse = "|"), full.names = T)
  
  # Creating empty tibble to store global results
  density_growth_reg <- tibble()
  
  for(fn in init_files){
    init_results <- read_json(fn, simplifyVector = T)
    
    min_pred_size <- init_results$params$min_log10_pred
    min_det_size <- init_results$params$min_log10_detritivore
    
    total_ts <- init_results$params$numb_time_steps
    
    # Processing density data for size spectrum
    density_growth_df <- data.frame(
      size_class = init_results$params$log10_size_bins, 
      pred_den = init_results$predators[,total_ts], 
      detrit_den = init_results$detritivores[,total_ts],
      # Growth from time step before end of simulation
      pred_growth = init_results$growth_int_pred[,total_ts-1], 
      detrit_growth = init_results$growth_det[,total_ts-1], 
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


# Binding all search volume results together for the entire world
global_success <- fao_lme |> 
  map(\(x) list.files(file.path(out_folder, x, "equilibrium_run"), 
                      "successful", full.names = T) |> 
        read_parquet()) |> 
  bind_rows()

# Save all results
global_success |> 
  write_parquet(file.path(
    out_folder, "successful_equilibrium_runs_inputs-outputs_global.parquet"))

# Save inputs only
global_in <- global_success |> 
  select(region, region_name, area_m2:depth) |> 
  distinct() 

global_in |>
  write_parquet(file.path(out_folder, "equilibrium_runs_inputs_global.parquet"))

# Save max and minimum search volume values resulting in successful runs
min_max_sv <- global_success |> 
  summarise(min_sv = min(search_vol), max_sv = max(search_vol), 
            .by = c(region, region_name)) 

min_max_sv |> 
  write_parquet(file.path(
    out_folder, "successful_equilibrium_runs_searchvol_global.parquet"))

plot_global_data <- min_max_sv |> 
  left_join(global_in)
  
# Maximum search volume values leading to successful runs
p1 <- plot_global_data |> 
  ggplot(aes(max_sv, fct_reorder(region, tos), fill = intercept))+
  geom_tile()+
  theme_bw()+
  labs(x = "maximum search volume",
       y = "ordered from warmest at the top to coolest at bottom (tos)",
       title = "")+
  scale_fill_viridis_c() 

ggsave(file.path(out_folder, 
                 "successful_equilibrium_runs_searchvol_global.png"), p1, 
       scale = 1.5, bg = "white")


p2 <- plot_global_data |> 
  ggplot(aes(max_sv, fct_reorder(region, intercept), fill = tos))+
  geom_tile()+
  theme_bw()+
  labs(x = "maximum search volume",
       y = "ordered from top to bottom from most to least productive (intercept)",
       title = "")+
  scale_fill_viridis_c() 

ggsave(file.path(out_folder, 
                 "successful_equilibrium_runs_searchvol_global_int-ordered.png"), 
       p2, scale = 1.5, bg = "white")



# Dynamic equilibrium run for the Arctic ----------------------------------
f <- "fao_lme-64"

dbpm_inputs <- read_parquet(
  file.path(base_folder, f, paste0("monthly_weighted", smoothed),
            paste0("dbpm_dynamic_clim-fish-inputs", fn_search, "_", f, 
                   "_1641-2010.parquet"))) |> 
  filter(scenario == "stable-spin")

results_folder <- file.path(out_folder, f, "dynamic_equilibrium_run")

search_volume <- 12.8

fish_param <- data.frame("region" = str_replace(str_to_upper(f), "-", " "),
                         "fmort_u" = 0, "fmort_v" = 0, "fminx_u" = 0, 
                         "fminx_v" = 0, "search_vol" = search_volume)

# Get parameters ready for calibration run
params <- sizeparam(dbpm_inputs, fish_param, xmin_consumer_u = -3,
                    xmin_consumer_v = -3)

# Saving non-spatial parameters
params |>
  #Ensuring up to 10 decimal places are saved in file
  write_json(file.path(results_folder,
                       paste0("dbpm_size_params_", f, "_searchvol_",
                              fish_param$search_vol, ".json")), digits = 10)

calib_run <- run_model(fish_param, dbpm_inputs, withinput = F, 
                       xmin_consumer_u = -3, xmin_consumer_v = -3,
                       include_plankton = T)

# Saving initial results for non-spatial run
calib_run |>
  #Ensuring up to 10 decimal places are saved in file
  write_json(file.path(results_folder,
                       paste0("init_dbpm_nonspatial_", f, "_searchvol_",
                              fish_param$search_vol, ".json")), digits = 10)


## Size spectrum plots per group (predators and detritivores) ---------
# Transform density matrix to data frame to create plots
density_df <- dbpm_output_mat_to_df(calib_run, dbpm_inputs$time, "density")

# Create size spectrum plots
size_sp_plot <- plotsizespectrum(density_df, calib_run$params, f,
                                 fishing_params = fish_param, mean_decade = T,
                                 nrow = 3)

# Saving size spectrum plot for non-spatial runs
ggsave(file.path(results_folder, paste0("size_spectrum_", f, "_searchvol_", 
                                        fish_param$search_vol, ".png")),
       size_sp_plot, bg = "white")

## Creating growth rate plots -------------------------------------------
dates_model <- c(min(dbpm_inputs$time)%m-% months(1), dbpm_inputs$time)
growth_df <- dbpm_output_mat_to_df(calib_run, dates_model, "growth") |>
  # Excluding growth value for first time step as it is used for model
  # initialisation only
  filter(time >= min(as_date(dbpm_inputs$time)))

# Create growth rate plot
growth_plot <- plot_growth_rate(growth_df, calib_run$params, f,
                                fishing_params = fish_param)

# Saving size spectrum plot for non-spatial runs
ggsave(file.path(results_folder, paste0("growth_rates_", f, "_searchvol_",
                                        fish_param$search_vol, ".png")),
       growth_plot, bg = "white")

#Equilibrium run
equilib_run <- run_model(fish_param, dbpm_inputs, withinput = T,
                         xmin_consumer_u = -3, xmin_consumer_v = -3, 
                         include_plankton = T)

# Save results
equilib_run |>
  write_parquet(file.path(results_folder, 
                          paste0("dbpm_nonspatial_", f, "_searchvol_",
                                 fish_param$search_vol, "_1801-2000.parquet")))

# Create plots of biomass (predators, detritivores and detritus)
biomass_data <- equilib_run |>
  select(year, ends_with("biomass"), total_detritus) |>
  group_by(year) |>
  summarise(across(starts_with("total"), ~ mean(.x, na.rm = T)))

bio_plot <- biomass_data |>
  pivot_longer(!year, names_to = "group", values_to = "values",
               names_prefix = "total_") |>
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
                        fish_param$search_vol, ".png")), bio_plot, bg = "white")

min_pred_size <- calib_run$params$min_log10_pred
min_det_size <- calib_run$params$min_log10_detritivore
total_ts <- calib_run$params$numb_time_steps

# Processing density data for size spectrum
density_growth_df <- data.frame(
  size_class = calib_run$params$log10_size_bins, 
  pred_den = calib_run$predators[, total_ts], 
  detrit_den = calib_run$detritivores[, total_ts],
  # Growth from time step before end of simulation
  pred_growth = calib_run$growth_int_pred[, total_ts-1], 
  detrit_growth = calib_run$growth_det[, total_ts-1], 
  search_vol = calib_run$params$hr_volume_search) |> 
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


size_spec_plot <- density_growth_df |>
  pivot_longer(c(pred_den, detrit_den, total), names_to = "group", 
               values_to = "bio") |> 
  ggplot(aes(size_class, bio, colour = group, linetype = group))+
  geom_line(alpha = 0.5)+
  scale_colour_manual(values = c("#1b9e77", "#d95f02", "#7570b3"))+
  scale_linetype_manual(values = c(2, 6, 1))+
  labs(colour = "Size-structured\ncommunities",
       linetype = "Size-structured\ncommunities")+
  theme(legend.title.position = "left",
        legend.title = element_text(hjust = 0.5, size = 10),
        legend.position = "top", legend.text = element_text(size = 10),
        legend.direction = "horizontal")+
  lims(x = c(calib_run$params$min_log10_detritivore, 
             calib_run$params$max_log10_pred), y = c(-20, NA))+
  labs(y = expression("" *log[10] ~ "abundance density (m"^-3* ")"),
       x = expression("" *log[10] ~ "body mass (g)"),
       title = paste0("Calibration (non-spatial) run - ",
                      str_replace_all(str_to_upper(f), "_|-", " ")),
       caption = paste0("Pelagic preference: ", 
                        round(calib_run$params$pref_pelagic, 3),
                        "\nBenthic preference: ", 
                        round(calib_run$params$pref_benthos, 3)))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        plot.caption = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(results_folder, 
                 paste0("size_spectrum_finaltimestep_", f, "_searchvol_", 
                        fish_param$search_vol, ".png")), size_spec_plot,
       bg = "white")

growth_rates_plot <- density_growth_df |> 
  pivot_longer(ends_with("growth"), names_to = "group", values_to = "growth") |>
  ggplot(aes(size_class_g, growth, colour = group))+
  geom_line()+
  scale_y_continuous(trans = "log10", 
                     name = "Relative growth rate per year")+
  scale_x_continuous(trans = "log10", name = "Body mass (g)")+
  labs(title = paste0(str_replace_all(str_to_upper(f), "_|-", " "),
                      ": Mean growth rate per decade"),
       caption = paste0("Pelagic preference: ", 
                        round(calib_run$params$pref_pelagic, 3),
                        "\nBenthic preference: ", 
                        round(calib_run$params$pref_benthos, 3)))+
  facet_grid(~group, scales = "free")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Saving size spectrum plot for non-spatial runs
ggsave(file.path(results_folder,
                 paste0("growth_rates_finaltimestep_", f, "_searchvol_",
                        fish_param$search_vol, ".png")), growth_rates_plot, 
       bg = "white")



