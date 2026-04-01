# Equilibrium runs

# Choose local R library
# .libPaths("/g/data/vf71/la6889/R_personal_lib/")

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
out_folder <- "/g/data/vf71/fishmip_outputs/ISIMIP3a/fao_lme_outputs/"

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
no_cores <- availableCores()
plan(multisession, workers = no_cores)

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
    density_df <- dbpm_output_mat_to_df(init_results, dbpm_inputs$time, "density")

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

