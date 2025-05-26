source("scripts/useful_functions.R")
library(dplyr)
library(arrow)
library(jsonlite)
library(purrr)
library(tibble)
library(lubridate)
library(zoo)


f <-  "fao-31"
res <- "1deg"
dbpm_inputs <- file.path(base_folder, f, "monthly_weighted", res,
                         paste0("dbpm_clim-fish-inputs_", f, 
                                "_1841-2010.parquet")) |> 
  read_parquet()

dbpm_inputs <- dbpm_inputs |> 
  mutate(time = as_date(time))

new_ts <- data.frame(time = seq(as_date("1841-01-01"), as_date("2010-12-31"), 
                                by = "week")) |> 
  mutate(year = year(time), month = month(time, label = T, abbr = F), 
         day = day(time)) |> 
  group_by(year, month) |> 
  mutate(time = case_when(day == min(day) ~ ymd(paste(year, month, "01", 
                                                      sep = "-")), 
                          T ~ time))

new_ts |> 
  write_csv_arrow("weekly_dates.csv")
  

new_inputs <- dbpm_inputs |> 
  right_join(new_ts, by = c("time", "year", "month")) |> 
  arrange(time) |>
  fill(region, scenario, depth_m, tot_area_m2, depth, area_m2) |> 
  mutate(across(c(expc_bot:tos, nom_active_area_m2_relative), 
                ~ na.approx(.x, na.rm = F))) |> 
  fill(c(expc_bot:tos, nom_active_area_m2_relative))



results_folder <- file.path(base_folder, f, "fishing_params", res,
                            "best_fish_vals_weekly")
  
dir.create(results_folder)

no_iter <- 100
params_calibration_weekly <- LHSsearch(num_iter = no_iter, 
                                forcing_file = new_inputs, 
                                gridded_forcing = NULL, 
                                best_val_folder = results_folder, 
                                best_param = F) |> 
  rowid_to_column("id")



params_corr <- params_calibration_weekly |> 
  split(params_calibration_weekly$id) |>
  map_df(\(x) getError(x, dbpm_inputs, corr = T)) 

params_calibration_weekly <- params_calibration_weekly |> 
  #Removing column to avoid duplication
  select(!region) |> 
  bind_cols(params_corr) |> 
  select(!id) |> 
  #Remove any rows where simulation returned NA values
  filter(catchNA == 0) |> 
  arrange(desc(cor), rmse) |> 
  relocate(region, .before = fmort_u)


params_calibration_weekly |> 
  filter(rmse == min(rmse) & cor >= 0.5) |> 
  corr_calib_plots(dbpm_inputs, results_folder)


fishing_params <- params_calibration_weekly |> 
  # Find parameters with lowest RMSE and correlation of 0.5 or higher
  filter(rmse == min(rmse) & cor >= 0.5)


results_folder <- file.path(base_folder, f, "fishing_params", res,
                            "best_fish_vals_weekly")

params <- sizeparam(new_inputs, fishing_params, xmin_consumer_u = -3, 
                    xmin_consumer_v = -3)

params |> 
  #Ensuring up to 10 decimal places are saved in file
  write_json(file.path(results_folder, 
                       paste0("dbpm_size_params_", f, ".json")), 
             digits = 10)


init_results <- run_model(fishing_params, new_inputs, withinput = F)

pred_initial <- rowMeans(init_results$predators)
detritivore_initial <- rowMeans(init_results$detritivores)
detritus_initial <- mean(init_results$detritus)

gridded_params <- sizeparam(new_inputs, fishing_params, xmin_consumer_u = -3, 
                            xmin_consumer_v = -3, 
                            use_init = T, pred_initial = pred_initial, 
                            detritivore_initial = detritivore_initial, 
                            detritus_initial = detritus_initial,
                            gridded = T)

#Save for use in gridded DBPM (step 05)
gridded_params |> 
  write_json(file.path(results_folder, 
                       paste0("dbpm_gridded_size_params_", f, ".json")),
             digits = 10)        
