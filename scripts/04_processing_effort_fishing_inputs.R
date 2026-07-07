# Processing fishing and effort inputs for Dynamic Benthic-Pelagic Size 
# Spectrum Model (DBPM)

# Activate local R library
.libPaths("/g/data/vf71/la6889/R_personal_lib/")

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
  summarise(min_fished_weight_class = min(log10mid_wt, na.rm = T),
            max_fished_weight_class = max(log10mid_wt, na.rm = T), 
            .groups = "drop") |> 
  mutate(region = case_when(area < 100 ~ paste0("LME ", area),
                            .default = paste0("FAO ", area)), .after = year) |>
  rename(region_name = corrected_name)

# Per-SPECTRUM (U pelagic / V benthic) fished-size window --------------------
# Reg's size_spectrum_catches file is fish-only and taxon-less, so the min/max above
# collapse ALL taxa into ONE window and miss krill/invertebrates. For DBPM's two-spectrum
# fishery we derive a per-spectrum window from the FGroup catch (catch_histsoc, already used
# above), mapping each FGroup to a gram size range and to U vs V:
#   fish FGroups: W = 0.01 * L^3 (g), L = cm size-class bounds
#     <30cm -> L[4,30]  30-90cm -> L[30,90]  >=90cm -> L[90,200]  <90cm -> L[4,90]
#   inverts (fixed): krill[1,2] shrimp[3,60] lobsterscrab[100,4000] cephalopods[20,6000]
#     demersalmollusc[2,500]
#   U = all fish + krill + cephalopods (water-column predators, incl. demersal fish)
#   V = shrimp, lobsterscrab, demersalmollusc (benthos)
# Each FGroup maps to a gram size RANGE [lo,hi]: fish W=0.01*L^3 with the cm-class bounds, where
# the smallest classes' LOWER bound is the realistic minimum FISHED length 10 cm (~10 g, e.g.
# anchovy/sardine gear onset), NOT 4 cm (larvae). Window = [min lo, max hi] over FGroups holding
# >=0.5% of that spectrum-year-region catch (drops trace bycatch). min = lower edge = smallest
# fished size (region-specific: small-pelagic LMEs ~10 g, toothfish LMEs stay high, krill 1 g).
# NA where a spectrum was not fished that year. (U max later overridden by the real WtMax below.)
V_groups <- c("shrimp", "lobsterscrab", "demersalmollusc")
uv_summ <- file.path(fishing_folder, "DKRZ_EffortFiles",
                     "catch_histsoc_1869_2017_EEZ_addFAO.csv") |>
  read_csv_arrow(col_select = c("Year", "fao_area", "LME",
                                "Reported", "IUU", "FGroup")) |>
  clean_names() |>
  mutate(area  = case_when(lme == 0 ~ fao_area + 100, .default = lme),
         catch = reported + iuu,
         lo = case_when(fgroup == "krill" ~ 1, fgroup == "shrimp" ~ 5,     # shrimp 5 g, krill 1 g
                        fgroup == "lobsterscrab" ~ 100, fgroup == "cephalopods" ~ 20,
                        fgroup == "demersalmollusc" ~ 10,                    # krill 1 g = fine mesh
                        str_detect(fgroup, "<30cm")   ~ 0.01 * 5^3,    # 5 cm ~ 1.25 g (graded pelagics)
                        str_detect(fgroup, "30-90cm") ~ 0.01 * 30^3,
                        str_detect(fgroup, ">=90cm")  ~ 0.01 * 90^3,
                        str_detect(fgroup, "<90cm")   ~ 0.01 * 10^3, .default = NA_real_),
         hi = case_when(fgroup == "krill" ~ 2, fgroup == "shrimp" ~ 60,
                        fgroup == "lobsterscrab" ~ 4000, fgroup == "cephalopods" ~ 6000,
                        fgroup == "demersalmollusc" ~ 500,
                        str_detect(fgroup, "<30cm")   ~ 0.01 * 30^3,
                        str_detect(fgroup, "30-90cm") ~ 0.01 * 90^3,
                        str_detect(fgroup, ">=90cm")  ~ 0.01 * 200^3,
                        str_detect(fgroup, "<90cm")   ~ 0.01 * 90^3, .default = NA_real_),
         spectrum = if_else(fgroup %in% V_groups, "V", "U")) |>
  filter(catch > 0, !is.na(lo)) |>
  group_by(area, year, spectrum) |>
  mutate(frac = catch / sum(catch)) |>
  filter(frac >= 0.005) |>
  summarise(mn = log10(min(lo)), mx = log10(max(hi)), .groups = "drop") |>   # lower/upper edges
  select(area, year, spectrum, mn, mx) |>
  pivot_wider(names_from = spectrum, values_from = c(mn, mx)) |>
  rename(min_fished_U = mn_U, max_fished_U = mx_U,
         min_fished_V = mn_V, max_fished_V = mx_V)

ss_catches_summ <- ss_catches_summ |>
  left_join(uv_summ, by = c("area" = "area", "year" = "year")) |>
  # HYBRID max: the pelagic U max uses Reg's real (WtMax-based) max_fished_weight_class rather
  # than the coarse FGroup cm-class max; the benthic V max keeps the invert-group max (the
  # fish-dominated max_fished_weight_class over-extends benthos). min (U,V) stays FGroup-derived.
  mutate(max_fished_U = max_fished_weight_class)

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
  # We will double the timesteps for the stable spinup period
  stable_spin <- list.files(forcing_folder, pattern = "^stable-spin_dbpm",
                            full.names = T) |> 
    read_parquet() |> 
    replicate(2, expr = _, simplify = F) |> 
    bind_rows() |> 
    mutate(time = seq(as_date("1641-01-01"), as_date("1840-12-31"), 
                      by = "month"), 
           year = year(time), month = month(time, label = T, abbr = F))
    
  
  # We will load climate inputs to merge with catch and effort data before 
  # saving results
  clim_forcing_file <- list.files(forcing_folder, 
                                  pattern = "obsclim|spinup",
                                  full.names = T) |>
    str_subset(paste0("inputs", fn_search, "_fao_lme")) |> 
    map(\(x) read_parquet(x)) |> 
    bind_rows(stable_spin) |>
    arrange(time) |> 
    clean_names() |> 
    mutate(time = as_date(time),
           depth = ifelse(is.na(depth_m), mean(depth_m, na.rm = T), depth_m)) |> 
    rename(area_m2 = tot_area_m2) |> 
    select(!depth_m)
  
  ## Dynamic stable spinup period for the Arctic only -----------------------
  if(fao_lme_id == 64){
    spinup <- list.files(forcing_folder, pattern = "spinup", full.names = T) |> 
      read_parquet()
    
    dyn_spinup <- spinup |> 
      group_by(month) |> 
      summarise(across(where(is.double) & !c(year, time), 
                       ~ mean(.x, na.rm = T))) |> 
      mutate(month = factor(month, levels = month.name, ordered = T)) |> 
      arrange(month) |> 
      replicate(200, expr = _, simplify = F) |> 
      bind_rows() |> 
      mutate(region = str_replace(str_to_upper(f), "-", " "), 
             scenario = "stable-spin", 
             time = seq(as_date("1641-01-01"), as_date("1840-12-31"), 
                        by = "month"), year = year(time), .before = month)
    
    clim_forcing_file <- list.files(forcing_folder, pattern = "obsclim",
                                    full.names = T) |> 
      read_parquet() |> 
      bind_rows(dyn_spinup, spinup) |> 
      arrange(time) |> 
      clean_names() |> 
      mutate(time = as_date(time), 
             depth = ifelse(is.na(depth_m), mean(depth_m, na.rm = T),
                            depth_m)) |> 
      rename(area_m2 = tot_area_m2) |> 
      select(!depth_m)
  }
    
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
  
  effort_data <- tibble(year = seq(1641, 1840), region = fao_lme_id,
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
  
  # Load minimum and maximum fish sizes harvested 
  ss_catches_summ <- read_csv_arrow(
    file.path(fishing_folder, "effort_catch_data", 
              "summary_size_spectrum_catches_fao-lme.csv")) |> 
    filter(area == fao_lme_id & year <= 2010) |> 
    select(!c(region, area))
  
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
    select(!c(region, depth, area_m2)) |> 
    full_join(ss_catches_summ, by = "year")
    
  rm(catch_pauly, catch_watson)

  # Merging catch and effort data -------------------------------------------
  DBPM_effort_catch_input <- effort_data |> 
    full_join(catch_data, by = "year") |> 
    mutate(region = case_when(fao_lme_id < 100 ~ paste0("LME ", region),
                              T ~ paste0("FAO ", region-100)),
           region_name = unique(ss_catches_summ$region_name)) |> 
    relocate(region_name, .after = region) |> 
    filter(year <= 2010)
  
  #Saving summarised catch and effort data
  DBPM_effort_catch_input |> 
    write_parquet(file.path(forcing_folder, 
                            paste0("dbpm_effort-catch-inputs", fn_search, "_",
                                   f, ".parquet")))
  
  #Removing individual data frames
  rm(effort_data, catch_data, ss_catches_summ)
  
  #Joining with climate inputs
  forcing_file <- clim_forcing_file |> 
    select(!region) |> 
    full_join(DBPM_effort_catch_input) |> 
    relocate(region, region_name, .after = scenario)


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
  if(fao_lme_id == 64){
    fout_forcing <- paste0("dbpm_dynamic_clim-fish-inputs", fn_search, "_", f, 
                           "_", min(forcing_file$year), "-", 
                           max(forcing_file$year), ".parquet")
  }
  
  forcing_file |> 
    write_parquet(file.path(forcing_folder, fout_forcing))
}

