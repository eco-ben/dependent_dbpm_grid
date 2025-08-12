source("scripts/useful_functions.R")
library(dplyr)
library(arrow)
library(jsonlite)
library(purrr)
library(tibble)
library(lubridate)
library(zoo)


f <-  "fao-27"
res <- "1deg"
base_folder <- "/g/data/vf71/fishmip_inputs/ISIMIP3a/fao_inputs"
dbpm_inputs <- file.path(base_folder, f, "monthly_weighted", res,
                         paste0("dbpm_clim-fish-inputs_", f, 
                                "_1841-2010.parquet")) |> 
  read_parquet() |> 
  mutate(time = as_date(time))

#Fishing parameters - No fishing
fishing_params <- read_parquet(list.files(file.path(base_folder, f, 
                                                    "fishing_params", res,
                                                    "best_fish_vals"), 
                                          "best-fishing", full.names = T)) |> 
  filter(rmse == min(rmse) & cor >= 0.5) |> 
  mutate(fmort_u = 0, fmort_v = 0)


# Weekly dates - Two columns - dates matching first day of month and 
# new dates not forced to start month at 01
new_ts <- data.frame(new_time = seq(as_date("1741-01-01"), as_date("2010-12-31"), 
                                by = "week")) |> 
  mutate(year = year(new_time), month = month(new_time, label = T, abbr = F), 
         day = day(new_time)) |> 
  group_by(year, month) |> 
  mutate(time = case_when(day == min(day) ~ ymd(paste(year, month, "01", 
                                                      sep = "-")), 
                          T ~ new_time))

new_ts |> 
  write_csv_arrow("scripts/weekly_dates.csv")



# Outlier grid cell -------------------------------------------------------
#Input data from outlier grid cell
dyn_outliers <- read_csv_arrow("scripts/dynamic_inputs_outliers.csv") |> 
  rowid_to_column("id")

#Names of columns where LOESS will be applied
vars_loess <- dyn_outliers |> 
  select(export_ratio:tos) |> 
  names()

#Defining LOESS function
loess_smooth <- function(x, span, df){
  loess(formula = paste(x, "id", sep = "~"), data = df, 
        span = span)$fitted
}  

#Applying LOESS - span 10
dyn_outliers_smooth01 <- as.data.frame(lapply(vars_loess, loess_smooth, 
                                              span = 0.001470588, df = dyn_outliers), 
                                       col.names = vars_loess) |> 
  mutate(time = dyn_outliers$time)

dyn_outliers_smooth02 <- as.data.frame(lapply(vars_loess, loess_smooth, 
                                              span = 0.002, df = dyn_outliers), 
                                       col.names = vars_loess) |> 
  mutate(time = dyn_outliers$time)

dyn_outliers_smooth05 <- as.data.frame(lapply(vars_loess, loess_smooth, 
                                             span = 0.005, df = dyn_outliers), 
                                      col.names = vars_loess) |> 
  mutate(time = dyn_outliers$time)

dyn_outliers_smooth1 <- as.data.frame(lapply(vars_loess, loess_smooth, 
                                             span = 0.01, df = dyn_outliers), 
                                      col.names = vars_loess) |> 
  mutate(time = dyn_outliers$time)

#Applying LOESS - span 20
dyn_outliers_smooth10 <- as.data.frame(lapply(vars_loess, loess_smooth, 
                                              span = 0.1, df = dyn_outliers), 
                                       col.names = vars_loess) |> 
  mutate(time = dyn_outliers$time)

dyn_outliers |>
  select(time, intercept) |>
  mutate(intercept_01 = dyn_outliers_smooth01$intercept,
         intercept_05 = dyn_outliers_smooth05$intercept,
         intercept_1 = dyn_outliers_smooth1$intercept,
         intercept_10 = dyn_outliers_smooth10$intercept) |>
  # filter(year(time) > 1960) |>
  pivot_longer(!time, names_to = "type", values_to = "vals") |>
  ggplot(aes(x = time, y = vals, color = type))+
  geom_line()


#Running model data from outlier grid cell
dbpm_inputs_outlier <- dbpm_inputs |> 
  select(!all_of(vars_loess)) |> 
  left_join(dyn_outliers)

init_results_outlier <- run_model(fishing_params, dbpm_inputs_outlier,
                                   withinput = F, new_detritus_calc = F)

init_results_outlier_new <- run_model(fishing_params, dbpm_inputs_outlier,
                                       withinput = F, new_detritus_calc = T)
#Running model with LOESS - span 0.001
dbpm_inputs_smooth01 <- dbpm_inputs |> 
  select(!all_of(vars_loess)) |> 
  left_join(dyn_outliers_smooth01)

init_results_smooth01 <- run_model(fishing_params, dbpm_inputs_smooth01,
                                   withinput = F, new_detritus_calc = F)

init_results_smooth01_new <- run_model(fishing_params, dbpm_inputs_smooth01,
                                       withinput = F, new_detritus_calc = T)

#Running model with LOESS - span 0.002
dbpm_inputs_smooth02 <- dbpm_inputs |> 
  select(!all_of(vars_loess)) |> 
  left_join(dyn_outliers_smooth02)

init_results_smooth02 <- run_model(fishing_params, dbpm_inputs_smooth02,
                                   withinput = F, new_detritus_calc = F)

init_results_smooth02_new <- run_model(fishing_params, dbpm_inputs_smooth02,
                                       withinput = F, new_detritus_calc = T)

#Running model with LOESS - span 0.05
dbpm_inputs_smooth05 <- dbpm_inputs |> 
  select(!all_of(vars_loess)) |> 
  left_join(dyn_outliers_smooth05)

init_results_smooth05 <- run_model(fishing_params, dbpm_inputs_smooth05,
                                   withinput = F, new_detritus_calc = F)

init_results_smooth05_new <- run_model(fishing_params, dbpm_inputs_smooth05,
                                       withinput = F, new_detritus_calc = T)


#Running model with LOESS - span 10
dbpm_inputs_smooth10 <- dbpm_inputs |> 
  select(!all_of(vars_loess)) |> 
  left_join(dyn_outliers_smooth10)
  
init_results_smooth10 <- run_model(fishing_params, dbpm_inputs_smooth10,
                                   withinput = F, new_detritus_calc = F)

init_results_smooth10_new <- run_model(fishing_params, dbpm_inputs_smooth10,
                                   withinput = F, new_detritus_calc = T)

#Running model with LOESS - span 20
dbpm_inputs_smooth20 <- dbpm_inputs |> 
  select(!all_of(vars_loess)) |> 
  left_join(dyn_outliers_smooth20)

init_results_smooth20 <- run_model(fishing_params, dbpm_inputs_smooth20,
                                   withinput = F, new_detritus_calc = F)
  
init_results_smooth20_new <- run_model(fishing_params, dbpm_inputs_smooth20,
                                       withinput = F, new_detritus_calc = T)


#Plotting detritus
detri_df <- data.frame(time = dbpm_inputs$time, 
                       detritus_original_old = init_results_outlier$detritus,
                       detritus_original_new = init_results_outlier_new$detritus,
                       detritus_smooth01_old = init_results_smooth01$detritus,
                       detritus_smooth01_new = init_results_smooth01_new$detritus,
                       detritus_smooth02_old = init_results_smooth02$detritus,
                       detritus_smooth02_new = init_results_smooth02_new$detritus,
                       detritus_smooth05_old = init_results_smooth05$detritus,
                       detritus_smooth05_new = init_results_smooth05_new$detritus) |> 
  pivot_longer(!time, names_to = "type", values_to = "vals") 


detri_df |> 
  filter(grepl(pattern = "smooth01", type)) |> 
  # filter(year(time) >= 1961) |>
  ggplot(aes(time, vals, color = type, shape = type))+
  geom_line()+
  geom_point()

data.frame(time = dbpm_inputs$time, 
           detritus_original_old = init_results_outlier$detritus,
           detritus_original_new = init_results_outlier_new$detritus) |> 
  pivot_longer(!time, names_to = "type", values_to = "vals") |> 
  # filter(year(time) >= 1961) |> 
  ggplot(aes(time, vals, color = type, shape = type))+
  geom_line()+
  geom_point()



# "Normal" grid cell ------------------------------------------------------
dyn_normal <- read_csv_arrow("scripts/dynamic_inputs_normal.csv") |> 
  rowid_to_column("id")

#Applying LOESS - span 10
dyn_normal_smooth10 <- as.data.frame(lapply(vars_loess, loess_smooth, 
                                              span = 0.1, df = dyn_normal), 
                                       col.names = vars_loess) |> 
  mutate(time = dyn_normal$time)

#Applying LOESS - span 20
dyn_normal_smooth20 <- as.data.frame(lapply(vars_loess, loess_smooth, 
                                              span = 0.2, df = dyn_normal), 
                                       col.names = vars_loess) |> 
  mutate(time = dyn_normal$time)

dyn_normal |>
  select(time, intercept) |>
  mutate(intercept_10 = dyn_normal_smooth10$intercept,
         intercept_20 = dyn_normal_smooth20$intercept) |>
  filter(year(time) > 1960) |>
  pivot_longer(!time, names_to = "type", values_to = "vals") |>
  ggplot(aes(x = time, y = vals, color = type))+
  geom_line()


#Running model data from normal grid cell
dbpm_inputs_normal <- dbpm_inputs |> 
  select(!all_of(vars_loess)) |> 
  left_join(dyn_normal)

init_results_normal <- run_model(fishing_params, dbpm_inputs_normal,
                                  withinput = F, new_detritus_calc = F)

init_results_normal_new <- run_model(fishing_params, dbpm_inputs_normal,
                                      withinput = F, new_detritus_calc = T)

#Running model with LOESS - span 10
dbpm_inputs_smooth10 <- dbpm_inputs |> 
  select(!all_of(vars_loess)) |> 
  left_join(dyn_normal_smooth10)

init_results_smooth10 <- run_model(fishing_params, dbpm_inputs_smooth10,
                                   withinput = F, new_detritus_calc = F)

init_results_smooth10_new <- run_model(fishing_params, dbpm_inputs_smooth10,
                                       withinput = F, new_detritus_calc = T)

#Running model with LOESS - span 20
dbpm_inputs_smooth20 <- dbpm_inputs |> 
  select(!all_of(vars_loess)) |> 
  left_join(dyn_normal_smooth20)

init_results_smooth20 <- run_model(fishing_params, dbpm_inputs_smooth20,
                                   withinput = F, new_detritus_calc = F)

init_results_smooth20_new <- run_model(fishing_params, dbpm_inputs_smooth20,
                                       withinput = F, new_detritus_calc = T)


#Plotting detritus
detri_df_normal <- data.frame(time = dbpm_inputs$time, 
                       detritus_original_old = init_results_normal$detritus,
                       detritus_original_new = init_results_normal_new$detritus,
                       detritus_smooth10_old = init_results_smooth10$detritus,
                       detritus_smooth10_new = init_results_smooth10_new$detritus,
                       detritus_smooth20_old = init_results_smooth20$detritus,
                       detritus_smooth20_new = init_results_smooth20_new$detritus) |> 
  pivot_longer(!time, names_to = "type", values_to = "vals") 


detri_df_normal |> 
  # filter(year(time) >= 1975) |>
  ggplot(aes(time, vals, color = type, shape = type))+
  geom_line()+
  geom_point()






# 
# new_ts <- data.frame(time = seq(as_date("1841-01-01"), as_date("2010-12-31"), 
#                                 by = "week")) |> 
#   mutate(year = year(time), month = month(time, label = T, abbr = F), 
#          day = day(time)) |> 
#   group_by(year, month) |> 
#   mutate(time = case_when(day == min(day) ~ ymd(paste(year, month, "01", 
#                                                       sep = "-")), 
#                           T ~ time))
# 
# new_ts |> 
#   write_csv_arrow("weekly_dates.csv")

new_ts <- read_csv_arrow("weekly_dates.csv")

new_inputs <- dbpm_inputs |> 
  right_join(new_ts, by = c("time", "year", "month")) |> 
  arrange(time) |>
  fill(region, scenario, depth_m, tot_area_m2, depth, area_m2) |> 
  mutate(across(c(expc_bot:tos, nom_active_area_m2_relative), 
                ~ na.approx(.x, na.rm = F))) |> 
  fill(c(expc_bot:tos, nom_active_area_m2_relative))




  
# dir.create(results_folder)




init_results <- run_model(fishing_params, dbpm_inputs, withinput = F, 
                          new_detritus_calc = F)










results_folder <- file.path(base_folder, f, "fishing_params", res,
                            "best_fish_vals_weekly_new_detritus")

params <- sizeparam(new_inputs, fishing_params, xmin_consumer_u = -3, 
                    xmin_consumer_v = -3)

params |> 
  #Ensuring up to 10 decimal places are saved in file
  write_json(file.path(results_folder, 
                       paste0("dbpm_size_params_", f, ".json")), 
             digits = 10)


init_results <- run_model(fishing_params, new_inputs, withinput = F, 
                          new_detritus_calc = T)

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



#Plotting detritus
detri_df <- data.frame(time = new_ts$time, year = new_ts$year, 
                       detritus = init_results$detritus) 

detri_df|> 
  ggplot(aes(time, detritus))+
  geom_line()

ggsave(file.path(results_folder, 
          paste0("detritus_", f, "_1841-2010.png")), device = "png")


detri_df |> 
  filter(year >= 1960) |> 
  ggplot(aes(time, detritus))+
  geom_line()

ggsave(file.path(results_folder, 
                 paste0("detritus_", f, "_1961-2010.png")), device = "png")
