# Processing fishing and effort inputs for Dynamic Benthic-Pelagic Size 
# Spectrum Model (DBPM)

# Loading libraries -------------------------------------------------------
library(dplyr)
library(tidyr)
library(arrow)
library(purrr)
library(janitor)
library(ggplot2)


# Define base variables ---------------------------------------------------
# FAO region number
fao_region <- 58
# Name of region
name_region <- "east_antarctica"
# Resolution
res <- "025deg"

#Creating path to ocean inputs
forcing_folder <- file.path("/g/data/vf71/la6889/dbpm_inputs", name_region, 
                            "monthly_weighted", res)

# Parent folder containing fishing effort and catches
fishing_folder <- "/g/data/vf71/fishmip_inputs/ISIMIP3a"

# Loading DBPM climate inputs ---------------------------------------------
# We will load climate inputs to merge with catch and effort data before saving
# results
clim_forcing_file <- list.files(forcing_folder, pattern = "obsclim|spinup",
                                full.names = T) |>
  map(\(x) read_parquet(x)) |> 
  bind_rows() |> 
  arrange(time) |> 
  clean_names()

# Getting the mean depth and area of the region of interest
depth_area <- clim_forcing_file |> 
  select(depth_m, tot_area_m2) |> 
  distinct()


## Loading effort data ----------------------------------------------------
effort_data <- file.path(fishing_folder, "DKRZ_EffortFiles",
                         "effort_isimip3a_histsoc_1841_2010.csv") |> 
  #Keep only columns with relevant information
  read_csv_arrow(col_select = c("Year", "fao_area", "NomActive")) |> 
  #Selecting data for area of interest
  filter(fao_area == fao_region) |> 
  # calculate sum of effort by area
  group_by(Year, fao_area) |> 
  summarise(total_nom_active = sum(NomActive, na.rm = T), 
            .groups = "drop") |> 
  ungroup() |>
  rename(year = Year, region = fao_area) |> 
  # Adding depth and area information for the area of interest
  mutate(depth = depth_area$depth_m, 
         area_m2 = depth_area$tot_area_m2,
         total_nom_active_area_m2 = total_nom_active/area_m2,
         nom_active_relative = total_nom_active/max(total_nom_active),
         nom_active_area_m2_relative = total_nom_active_area_m2/
         max(total_nom_active_area_m2))


# Loading catches data ----------------------------------------------------
#From Watson et al 2018
catch_watson <- file.path(fishing_folder, "DKRZ_EffortFiles",
                          "catch-validation_isimip3a_histsoc_1850_2004.csv") |> 
  read_csv_arrow(col_select = c("Year", "fao_area", "Reported", "IUU")) |>
  #Selecting area of interest
  filter(fao_area == fao_region) |> 
  # catch is in tonnes. This was checked in "FishingEffort" project
  mutate(catch_tonnes = Reported+IUU) |> 
  group_by(Year, fao_area) |> 
  summarise(catch_tonnes = sum(catch_tonnes), .groups = "drop") |> 
  # also Reg advise to exclude discards 
  ungroup() |> 
  rename(year = Year, region = fao_area) |> 
  # Adding depth and area information for the area of interest
  mutate(depth = depth_area$depth_m, 
         area_m2 = depth_area$tot_area_m2,
         catch_tonnes_area_m2 = catch_tonnes/area_m2)

#From Pauly et al 2020
catch_pauly <- read.csv(list.files(file.path(fishing_folder, "SAU_catch_data"), 
                                   pattern = as.character(fao_region), 
                                   full.names = T)) |> 
  #Keep data up to 2010
  filter(year <= 2010) |> 
  group_by(year) |> 
  #Calculate total tonnes caught per year
  summarise(catch_tonnes_pauly = sum(tonnes, na.rm = T)) 


#From CCAMLR
#Load catch data
catch_ccamlr <- file.path(fishing_folder, 
                          "CCAMLR_Statistical_Bulletin_V36/catch_and_effort",
                          "Catch.csv") |> 
  read_csv_arrow(col_select = c("year", "asd_code", 
                                "greenweight_caught_tonne")) |> 
  #Merge with codes for CCAMLR regions
  left_join(read_csv_arrow(
    file.path(fishing_folder, "CCAMLR_Statistical_Bulletin_V36",
              "catch_and_effort/Area.csv")), by = "asd_code") |> 
  #Keep only data for region and time period of interest
  filter(asd_area_code == fao_region & year <= 2010) |> 
  #Calculate total tonnes caught per year within boundaries of FAO region
  #(including sub regions)
  group_by(year) |> 
  summarise(catch_tonnes_ccamlr = sum(greenweight_caught_tonne, na.rm = T)) |> 
  full_join(catch_pauly, by = "year") |> 
  arrange(year) |> 
  mutate(area_m2 = depth_area$tot_area_m2) |>
  mutate(catch_pauly = catch_tonnes_pauly/area_m2, 
         .after = catch_tonnes_pauly) |> 
  mutate(catch_ccamlr = catch_tonnes_ccamlr/area_m2, 
         .after = catch_tonnes_ccamlr) |>
  select(!area_m2) |> 
  full_join(catch_watson, by = "year") |> 
  select(!region:area_m2) |> 
  arrange(year) |> 
  filter(!if_all(c(catch_tonnes_area_m2, catch_pauly, catch_ccamlr), is.na)) |>  
  rowwise() |>
  mutate(min_catch_density = min(catch_tonnes_area_m2, catch_pauly, 
                                 catch_ccamlr, na.rm = T),
         max_catch_density = max(catch_tonnes_area_m2, catch_pauly,
                                 catch_ccamlr, na.rm = T)) |> 
  select(!catch_tonnes_area_m2)

#Merge catches datasets
catch_data <- catch_watson |>
  full_join(catch_ccamlr, by = "year") |> 
  arrange(year) |> 
  relocate(catch_tonnes, .after = area_m2) |> 
  mutate(region = case_when(is.na(region) ~ fao_region, T ~ region), 
         depth = case_when(is.na(depth) ~ depth_area$depth_m, T ~ depth),
         area_m2 = case_when(is.na(area_m2) ~ depth_area$tot_area_m2, 
                             T ~ area_m2))

rm(catch_pauly, catch_watson, catch_ccamlr)

# Merging catch and effort data -------------------------------------------
DBPM_effort_catch_input <- effort_data |> 
  full_join(catch_data) |> 
  mutate(region = paste0("FAO ", region))

#Saving summarised catch and effort data
DBPM_effort_catch_input |> 
  write_parquet(file.path(forcing_folder, 
                          paste0("dbpm_effort-catch-inputs_fao-",
                                 fao_region, ".parquet")))

#Removing individual data frames
rm(effort_data, catch_data)

#Joining with climate inputs
forcing_file <- clim_forcing_file |> 
  full_join(DBPM_effort_catch_input, by = c("region", "year")) 


## Plotting fish and catch data -------------------------------------------
forcing_file |> 
  ggplot(aes(year, total_nom_active))+
  annotate("rect", xmin = 1841, xmax = 1960, ymin = 0, ymax = Inf, 
           fill = "#b2e2e2", alpha = 0.4)+ 
  annotate("rect", xmin = 1961, xmax = 2010, ymin = 0, ymax = Inf, 
           fill = "#238b45", alpha = 0.4)+ 
  geom_point(size = 1)+
  geom_line()+
  scale_x_continuous(expand = c(.01, 0), breaks = seq(1850, 2010, 20))+
  scale_y_continuous(expand = c(.02, 0))+
  theme_bw()+
  labs(y = "Total nom active", 
       title = paste0("FAO Major Fishing Area # ", fao_region))+
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        axis.title.x = element_blank(), axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) 

#Saving result that matches previous work
ggsave(paste0("new_workflow/outputs/effort_fao-", fao_region, ".pdf"), 
       device = "pdf", dpi = 300)


## Saving catch and effort, and inputs data -------------------------------
forcing_file |> 
  write_parquet(file.path(forcing_folder, 
                          paste0("dbpm_clim-fish-inputs_fao-", fao_region, 
                                 "_1841-2010.parquet")))

