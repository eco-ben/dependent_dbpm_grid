
# CN 20/08/2023 - create a plotting file for DBPM based on plots prepared by 
# Romain Forestier and JB
library(ggplot2)
library(dplyr)
library(reshape2)
library(data.table)
library(sf)
source("LME_calibration.R")

plotgridbyLME <- function(LMEnumber){
  
  # # trial 
  # LMEnumber <- 2
  
  # load outputs 
  LME_path <- file.path("/g/data/vf71/fishmip_outputs/ISIMIP3a/DBPM/obsclim",
                        paste0("LME_", LMEnumber, "_output"))
  load(list.files(LME_path, "gridded_model_results", full.names = T))
  
  # load inputs and param object 
  load(list.files(LME_path, "grid_inputs", full.names = T))
  lme_inputs_grid <- fread(list.files(LME_path, "LME_inputs_gridded.csv",
                                      full.names = T))
  lme_input_init <- fread(list.files(LME_path, "LME_inputs_non-gridded.csv",
                                     full.names = T))
  
  # Would this be the original time steps ?
  # gridded_params$Neq <- 2040 
  # If so, this should be replaced by:
  gridded_params$Neq <- length(unique(lme_inputs_grid$t))

  if(sum(gridded_params$effort) > 0){
    f.effort = TRUE
  }else{
    f.effort = FALSE
  }
  
  ### WARNING need to check depth integration and Neq - do they match what used
  #for inputs and for runLMEcalibration? 
  out <- getGriddedOutputs(input = lme_inputs_grid, results = grid_results,
                           params = gridded_params)

  # World map
  world <- read_sf("data/world_map.shp")

  ####### Maps of biomass ----
  # averaged across decades for U and V
  biom_df <- out |> 
    select(lat, lon, t, year, TotalVbiomass, TotalUbiomass) |> 
    mutate(totalB = TotalVbiomass + TotalUbiomass,
           decade = (year %/% 10)*10)

  # calculate the mean biomass for each decade  
  biom_decade_avg <- biom_df |>
    filter(decade >= 1950 & decade <= max(decade)) |> 
    group_by(decade, lon, lat) |>
    summarize(avg_Ubiomass = mean(TotalUbiomass, na.rm = T),
              avg_Vbiomass = mean(TotalVbiomass, na.rm = T))

  # Spatial extent
  lat_ext <- c(min(biom_decade_avg$lat), max(biom_decade_avg$lat))
  lon_ext <- c(min(biom_decade_avg$lon), max(biom_decade_avg$lon))
  
  # facet plot U
  p1 <- ggplot()+
    geom_tile(data = biom_decade_avg, aes(lon, lat, fill = avg_Ubiomass))+
    geom_sf(data = world)+
    coord_sf(xlim = lon_ext, ylim = lat_ext)+
    scale_fill_gradient2(low = "white", high = "red", mid = "grey",
                         name = bquote("Biomass (g*"*m^-2*")"))+
    facet_wrap(~decade, nrow = 2)+
    theme_bw()+
    ggtitle("Mean U biomass across decades")+
    theme(legend.position = "bottom", axis.title = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 0.75, hjust = 0.75), 
          plot.title = element_text(hjust = 0.5), legend.title.position = "top")
    
  # facet plot V
  p2 <- ggplot()+
    geom_tile(data = biom_decade_avg, aes(lon, lat, fill = avg_Vbiomass))+
    geom_sf(data = world)+
    coord_sf(xlim = lon_ext, ylim = lat_ext)+
    scale_fill_gradient2(low = "white", high = "red", mid = "grey",
                         name = bquote("Biomass (g*"*m^-2*")"))+
    facet_wrap(~decade, nrow = 2)+
    theme_bw()+
    ggtitle("Mean V biomass across decades")+
    theme(legend.position = "bottom", axis.title = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 0.75, hjust = 0.75), 
          plot.title = element_text(hjust = 0.5), legend.title.position = "top")
  
  ####### Biomass trends ---- 
  # calculate the mean biomass across grid cells
  biom_grid_avg <- biom_df |>
    group_by(year) |>  
    summarize(avg_biomass = mean(totalB)) |> 
    ungroup()

  # time series of biomass
  p3 <- ggplot(biom_grid_avg, aes(x = year, y = avg_biomass))+
    geom_line()+
    labs(y = bquote("Biomass (g*"*m^-2*")"), x = "Year")+
    ggtitle("Mean biomass through time")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))

  ####### Size spectra ---- 
  # U + V averaged per decade per longitude
  # Using grid_results$U and V. Dim are gridcell*size*time (monthly)
  
  totBiom <- grid_results$U + grid_results$V
  
  
  # average per longitude
  uniq_coord <- biom_df |> 
    distinct(lat, lon) 
  
  uniq_lat <- uniq_coord |> 
    distinct(lat)
  
  lat_list <- vector("list", length = nrow(uniq_lat))
  names(lat_list) <- uniq_lat$lat
  
  for(iCell in seq_along(uniq_lat$lat)) {
    lat_list[[iCell]] <- which(uniq_coord$lat == uniq_lat$lat[iCell])
  }
  
  spectra_grid_avg <- array(NA, 
                            dim = c(length(lat_list), ncol(spectra_decade_avg),
                                    length(uniq_dec)),
                            dimnames = list("latitude" = uniq_lat$lat,
                                            "size" = grid_results$params$x,
                                            "decade" = uniq_dec))
  for (iCell in names(lat_list)) {
    tempBiom <- spectra_decade_avg[lat_list[[iCell]], , ]
    if(length(dim(tempBiom)) == 3){
      # if = 2 means only one lat value
      avgBiom <- apply(tempBiom, 2:3, mean) 
    }else{
      avgBiom <- tempBiom
    }
    spectra_grid_avg[iCell, , ] <- avgBiom
  }

  # show sizes only between $ref and $Nx
  ind <- grid_results$params$ref:grid_results$params$Nx
  spectra_grid_avg <- spectra_grid_avg[, ind, ]

  plot_dat <- reshape2::melt(spectra_grid_avg)

  p4 <- ggplot()+
    geom_line(data = plot_dat, aes(x = size, y = value, color = latitude), 
              alpha = .5)+
    facet_wrap(~decade)+
    scale_y_continuous(trans = "log10", name = "Biomass (g)")+
    xlab(bquote("Size ("*log[10]*"g)"))+
    ggtitle("Size spectra averaged per decade and longitude (U+V)")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))

  ####### Growth rate ---- 
  # GG.u + GG.v per decade per gridcell
  totGrowth <- grid_results$GG_u+grid_results$GG_v

  # averaging per decade. First decade is 108 month then the rest is 120
  growth_decade_avg <- array(NA, 
                             dim = c(nrow(totBiom), ncol(totBiom), 
                                     length(uniq_dec)),
                             dimnames = list("grid_cell" = 1:nrow(totBiom),
                                             "size" = grid_results$params$x,
                                             "decade" = uniq_dec))
  
  for(i in seq_along(uniq_dec)){
    tempGrowth <- totGrowth[, , which(decades$decade == uniq_dec[i])]
    avgGrowth <- apply(tempGrowth, 1:2, mean)
    growth_decade_avg[, , i] <- avgGrowth
  }
  
  # average per longitude - using lat_list from previous plot
  growth_grid_avg <- array(NA, 
                           dim = c(length(lat_list), ncol(spectra_decade_avg),
                                   length(uniq_dec)),
                           dimnames = list("latitude" = uniq_lat$lat,
                                           "size" = grid_results$params$x,
                                           "decade" = uniq_dec))

  for(iCell in names(lat_list)){
    tempBiom <- growth_decade_avg[lat_list[[iCell]], , ]
    if(length(dim(tempBiom)) == 3){
      # if = 2 means only one lat value
      avgBiom <- apply(tempBiom, 2:3, mean)
    }else{
      avgBiom <- tempBiom
    }
    growth_grid_avg[iCell, , ] <- avgBiom
  }

  # show sizes only between $ref and #Nx
  growth_grid_avg <- growth_grid_avg[, ind,]

  plot_dat <- reshape2::melt(growth_grid_avg)

  p5 <- ggplot()+
    geom_line(data = plot_dat, aes(x = size, y = value, color = latitude))+
    facet_wrap(~decade)+
    scale_y_continuous(trans = "log10", name = "Relative growth rate per year")+
    xlab(bquote("Size ("*log[10]*"g)"))+
    ggtitle("Mean growth rate per decade and longitude (U+V)")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))

  

  ####### Save plots ---- 
  name <- ifelse(f.effort == FALSE, "plots_no_fishing.pdf", "plots_fishing.pdf")
  pdf(file.path(LME_path, name), height = 10, width = 10, onefile = T)
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
  print(p6)
  print(p7)
  print(p8)
  dev.off()
}

