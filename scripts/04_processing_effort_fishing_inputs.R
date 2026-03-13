# Processing fishing and effort inputs for Dynamic Benthic-Pelagic Size 
# Spectrum Model (DBPM)

# Loading libraries -------------------------------------------------------
library(dplyr)
library(tidyr)
library(arrow)
library(purrr)
library(janitor)
library(ggplot2)
library(stringr)
library(lubridate)


# Defining base folder
fishing_folder <- "/g/data/vf71/fishmip_inputs/ISIMIP3a"

# Processing global effort data -------------------------------------------
# This step needs to be completed only once at a global scale
# Creating summaries of effort per year and region of interest
effort_data_global <- file.path(fishing_folder, "DKRZ_EffortFiles",
                         "effort_isimip3a_histsoc_1841_2010.csv") |> 
  read_csv_arrow(col_select = c("Year", "fao_area", "LME", "NomActive")) |> 
  clean_names() |>
  mutate(region = case_when(lme == 0 ~ fao_area+100, .default = lme)) |> 
  # calculate sum of effort by area
  group_by(year, region) |> 
  summarise(total_nom_active = sum(nom_active, na.rm = T)) |> 
  ungroup()

# Saving summarised data
effort_data_global |> 
  write_csv_arrow(
    file.path(fishing_folder, "DKRZ_EffortFiles",
              "yearly_effort_fao-lme_isimip3a_histsoc_1841_2010.csv")) 


# Processing global catch data (Watson) -----------------------------------
# This step needs to be completed only once at a global scale
catch_watson <- file.path(fishing_folder, "DKRZ_EffortFiles",
                          "catch_histsoc_1869_2017_EEZ_addFAO.csv") |> 
  read_csv_arrow(col_select = c("Year", "fao_area", "LME", "Reported", "IUU")) |>
  clean_names() |> 
  mutate(region = case_when(lme == 0 ~ fao_area+100, .default = lme)) |> 
  group_by(year, region) |> 
  summarise(tot_reported = sum(reported, na.rm = T),
            tot_iuu = sum(iuu, na.rm = T)) |> 
  rowwise() |> 
  # catch is in tonnes. This was checked in "FishingEffort" project
  mutate(catch_tonnes = sum(tot_reported, tot_iuu, na.rm = T)) |>
  # also Reg advise to exclude discards 
  ungroup() |> 
  select(!starts_with("tot_")) 

# Saving summarised data
catch_watson |> 
  write_csv_arrow(
    file.path(fishing_folder, "DKRZ_EffortFiles",
              "yearly_catch_fao-lme_isimip3a_histsoc_1869_2017.csv")) 


# Processing size spectrum - catches --------------------------------------
# Size spectra data from Reg (original file name: SizesinLMET2.csv")
ss_catches <- read_csv_arrow(file.path(fishing_folder, "effort_catch_data", 
                                       "size_spectrum_catches_fao-lme.csv")) |> 
  clean_names()
# Loading names of LMEs and FAO regions
lme_names <- read_csv_arrow(file.path("/g/data/vf71/shared_resources", 
                                      "fao_lme_masks/fao-major_lme_keys.csv"),
                      col_select = c("fao_lme", "corrected_name"))

# Finding maximum and minimum weight classes of catches (per year and AOI)
# The 'Log10MidWt' column was selected based on the 'Plot Size Data 8.R' script 
# from Reg that produces the figure of size spectrum plots for each region
ss_catches_summ <- ss_catches |> 
  left_join(lme_names, by = c("area"="fao_lme")) |> 
  group_by(year, area, corrected_name) |> 
  summarise(min_weight_class = min(log10mid_wt, na.rm = T),
            max_weight_class = max(log10mid_wt, na.rm = T), .groups = "drop") |> 
  mutate(region = case_when(area < 100 ~ paste0("LME ", area),
                            .default = paste0("FAO ", area)), .after = year) |> 
  rename(region_name = corrected_name) |> 
  select(!area)

# Saving summarised data
ss_catches_summ |> 
  write_csv_arrow(file.path(fishing_folder, "effort_catch_data", 
                            "summary_size_spectrum_catches_fao-lme.csv"))


# Define base variables ---------------------------------------------------
base_folder <- "/g/data/vf71/fishmip_inputs/ISIMIP3a/fao_lme_inputs"
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

# Applying workflow to all regions
for(f in fao_lme){
  fao_lme_id <- as.numeric(str_extract(f, "[0-9]+"))
  #Creating path to ocean inputs
  forcing_folder <- file.path(base_folder, f, 
                              paste0("monthly_weighted", smoothed))
  
  # Loading DBPM climate inputs ---------------------------------------------
  # We will load climate inputs to merge with catch and effort data before 
  # saving results
  clim_forcing_file <- list.files(forcing_folder, 
                                  pattern = "obsclim|spinup|stable-spin",
                                  full.names = T) |>
    str_subset(paste0("inputs", fn_search, "_fao_lme")) |> 
    map(\(x) read_parquet(x)) |> 
    bind_rows() |> 
    clean_names() |> 
    arrange(time) |> 
    mutate(time = as_date(time)) |> 
    rename(depth = depth_m, area_m2 = tot_area_m2)
  
  # Getting the mean depth and area of the region of interest
  depth_area <- clim_forcing_file |> 
    distinct(depth, area_m2)

  ## Loading effort data ----------------------------------------------------
  effort_data <- read_csv_arrow(
    file.path(fishing_folder, "DKRZ_EffortFiles",
              "yearly_effort_fao-lme_isimip3a_histsoc_1841_2010.csv")) |> 
    #Selecting data for area of interest
    filter(region == fao_lme_id) 
  
  # Extend fishing effort (nom_active_area_m2_relative) starting in 1741
  # Repeat 1841 value for entire stable spin up period
  effort_stable_spin <- effort_data |> 
    filter(year == 1841) |> 
    pull(total_nom_active)
  
  effort_data <- tibble(year = seq(1741, 1840), region = fao_lme_id,
                        total_nom_active = effort_stable_spin) |> 
    bind_rows(effort_data) |> 
    # Adding depth and area information for the area of interest
    mutate(depth = depth_area$depth, 
           area_m2 = depth_area$area_m2,
           total_nom_active_area_m2 = total_nom_active/area_m2,
           nom_active_relative = total_nom_active/max(total_nom_active),
           nom_active_area_m2_relative = total_nom_active_area_m2/
             max(total_nom_active_area_m2))
  
  rm(effort_stable_spin)
  
  # Loading catches data ----------------------------------------------------
  #From Watson et al 2018
  catch_watson <- read_csv_arrow(
    file.path(fishing_folder, "DKRZ_EffortFiles",
              "yearly_catch_fao-lme_isimip3a_histsoc_1869_2017.csv")) |> 
    #Selecting area of interest
    filter(region == fao_lme_id & year <= 2010) |> 
    mutate(depth = depth_area$depth, 
           area_m2 = depth_area$area_m2,
           catch_tonnes_area_m2 = catch_tonnes/area_m2) |> 
    relocate(catch_tonnes, .before = catch_tonnes_area_m2)

  #From Pauly et al 2020
  if(fao_lme_id < 100){
    pat_look <- str_c("LME ", fao_lme_id, " v50-1.csv")
  }else{
    pat_look <- str_c("FAO ", (fao_lme_id-100), " v50-1.csv")
  }
  catch_pauly <- read.csv(
    list.files(file.path(fishing_folder, "SAU_catch_data"), 
               pattern = pat_look, full.names = T)) |> 
    # Keep data up to 2010 and removing discards to match processing of Watson
    # data
    filter(year <= 2010 & catch_type != "Discards") |> 
    group_by(year) |> 
    #Calculate total tonnes caught per year
    summarise(catch_tonnes_pauly = sum(tonnes, na.rm = T))
  
  catch_data <- catch_watson |> 
    full_join(catch_pauly, by = "year") |> 
    mutate(catch_pauly_tonnes_area_m2 = catch_tonnes_pauly/area_m2) |> 
    arrange(year) |> 
    filter(!if_all(c(catch_tonnes_area_m2, catch_pauly_tonnes_area_m2), 
                   is.na)) |> 
    rowwise() |>
    mutate(min_catch_density = min(catch_tonnes_area_m2, 
                                   catch_pauly_tonnes_area_m2, na.rm = T),
           max_catch_density = max(catch_tonnes_area_m2, 
                                   catch_pauly_tonnes_area_m2, na.rm = T)) |> 
    select(!c(region, depth, area_m2))

  rm(catch_pauly, catch_watson)

  # Merging catch and effort data -------------------------------------------
  DBPM_effort_catch_input <- effort_data |> 
    full_join(catch_data, by = "year") |> 
    mutate(region = case_when(fao_lme_id < 100 ~ paste0("LME ", region),
                              T ~ paste0("FAO ", region-100)))
  
  #Saving summarised catch and effort data
  DBPM_effort_catch_input |> 
    write_parquet(file.path(forcing_folder, 
                            paste0("dbpm_effort-catch-inputs", fn_search, "_",
                                   f, ".parquet")))
  
  #Removing individual data frames
  rm(effort_data, catch_data)
  
  #Joining with climate inputs
  forcing_file <- clim_forcing_file |> 
    select(!region) |> 
    full_join(DBPM_effort_catch_input) 


  ## Plotting fish and catch data -------------------------------------------
  forcing_file |> 
    filter(year >= 1841) |> 
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
         title = unique(DBPM_effort_catch_input$region))+
    theme(plot.title = element_text(size = 12, hjust = 0.5),
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 12),
          axis.text = element_text(size = 10), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank()) 
  
  folder_out <- file.path("/g/data/vf71/fishmip_outputs/ISIMIP3a",
                          "fao_lme_outputs", f)
  if(!dir.exists(folder_out)){
    dir.create(folder_out, recursive = T)
  }
  
  #Saving results - only save once per FAO region
  fout <- file.path(folder_out, paste0("effort_", f, ".pdf"))
  if(!file.exists(fout)){
    ggsave(fout, device = "pdf", dpi = 300)
  }

  ## Saving catch and effort, and inputs data -------------------------------
  fout_forcing <- paste0("dbpm_clim-fish-inputs", fn_search, "_", f, "_", 
                         min(forcing_file$year), "-", max(forcing_file$year),
                         ".parquet")
  
  forcing_file |> 
    write_parquet(file.path(forcing_folder, fout_forcing))
}

