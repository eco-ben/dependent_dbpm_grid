###############################################################################
# Supporting DBPM functions
# Functions have been adapted from previous DBPM work by JLB, CN, JB and others
# 
# Edited by: Denisse Fierro Arcos
# Date of update: 2024-12-19


# Loading libraries -------------------------------------------------------
library(tidyr)
library(dplyr)
library(stringr)
library(arrow)
library(parallel)
library(ggplot2)
library(lhs)
library(patchwork)

# Getting DBPM model parameters ready -------------------------------------
sizeparam <- function(dbpm_inputs, fishing_params, dx = 0.1, xmin = -12, 
                      xmin_consumer_u = -7, xmin_consumer_v = -7, xmax = 6, 
                      tstepspryr = 48, Ngrid = NA, use_init = FALSE, 
                      pred_initial = NA, detritivore_initial = NA, 
                      detrititus_initial = NA, equilibrium = FALSE,
                      gridded = FALSE){
  
  #Inputs:
  # - dbpm_inputs (data frame) Containing climate and fishing data needed as 
  #   inputs for DBPM. This is produced by script 
  #   "03_processing_effort_fishing_inputs.R"
  # - fishing_params (data frame) Containing fishing parameters: "f_u", "f_v",
  #   "f_minu", "f_minv", and "search_vol"
  # - dx (numeric) Default value is 0.1. Size increment after discretization for 
  #   integration (log body weight)
  # - xmin (numeric) Default value is -12. Minimum log10 body size of plankton
  # - xmin_consumer_u (numeric) Default value is -7. Minimum log10 body size in 
  #   dynamics predators
  # - xmin_consumer_v (numeric) Default value is -7. Minimum log10 body size in
  #   dynamics benthic detritivores
  # - xmax (numeric). Default value is 6. Maximum log10 body size of predators
  # - tstepspryr (numeric) Default value is 48. Number of time steps to include 
  #   within a year
  # - Ngrid (numeric) Optional. Number of grid cells.
  # - use_init (boolean). Default value is FALSE. If set to TRUE, the 
  #   initialisation values for predators (U_initial), detritivores (V.initial) 
  #   and detritus (W.initial) will be used
  # - pred_initial (numeric). Optional. Default is NA. Initialisation value for
  #   predators. If provided, 'use_init' must be set to TRUE
  # - detritivore_initial (numeric). Optional. Default is NA. Initialisation 
  #   value for detritivores. If provided, 'use_init' must be set to TRUE
  # - detrititus_initial (numeric). Optional. Default is NA. Initialisation 
  #   value for detritus If provided, 'use_init' must be set to TRUE
  # - equilibrium (boolean). Default value is FALSE.
  # - gridded(boolean). Default value is FALSE. If set to TRUE, it will provide
  #   params for non-gridded model run.
  #
  # Outputs:
  # er (numeric vector) Export ratio
  
  # Creating empty list to store DBPM parameters:
  param <- list()
  
  # If gridded is selected, the parameters below from "depth" to 
  # "sea_floor_temp" are not included because there are zarr files available 
  # with this information. "pref_benthos" and "pref_pelagic" need to be 
  # calculated from gridded files.
  
  if(!gridded){
    # depth
    param$depth <- mean(dbpm_inputs$depth)  
    # export ratio (er) replaced by fraction of sinking detritus reaching the 
    # seafloor (from export ratio input) (sinking.rate)
    param$sinking_rate <- dbpm_inputs$export_ratio 
    # plankton parameters
    # intercept of phyto-zooplankton spectrum (pp)
    param$int_phy_zoo <- dbpm_inputs$intercept
    # slope of phyto-zooplankton spectrum (r.plank)
    param$slope_phy_zoo <- dbpm_inputs$slope
    # temperature parameters 
    # sea-surface temperature - degrees Celsius (sst)
    param$sea_surf_temp <- dbpm_inputs$tos
    # near sea-floor temperature - degrees Celsius (sft)
    param$sea_floor_temp <- dbpm_inputs$tob
    
    # set predator coupling to benthos, depth dependent - 0.75 above 500 m, 0.5
    # between 500-1800 and 0 below 	1800m (suggestions of values from Clive
    # Trueman based on stable isotope work, and proportion of biomass, 	Rockall 
    # Trough studies) (pref.ben)
    param$pref_benthos <- 0.8*exp(-1/250*param$depth)
    
    # preference for pelagic prey (pref.pel)
    param$pref_pelagic <- 1-param$pref_benthos 
  }
  
  # number of years to run model (tmax)
  param$n_years <- length(unique(dbpm_inputs$year))
  
  # discretisation of year(delta.t)
  param$timesteps_years <- (1/tstepspryr)
  
  # number of time bins (Neq)
  param$numb_time_steps <- param$n_years*tstepspryr
  
  # get rescaled effort
  param$effort <- dbpm_inputs$nom_active_area_m2_relative
  
  # fishing parameters 
  # fishing mortality rate for predators (Fmort.u)
  param$fish_mort_pred <- fishing_params$fmort_u
  # fishing mortality rate for detritivores (Fmort.v)
  param$fish_mort_detritivore <- fishing_params$fmort_v
  # minimum log10 body size fished for predators (min.fishing.size.u)
  param$min_fishing_size_pred <- fishing_params$fminx_u 
  # minimum log10 body size fished for detritivores (min.fishing.size.v)
  param$min_fishing_size_detritivore <- fishing_params$fminx_v 
  
  # Benthic-pelagic coupling parameters
  
  # originally 640, but if 64 then using Quest-fish default of 64 hourly rate
  # volume searched constant m3.yr-1 for fish. need to check this value, its
  # quite large. (A.u)
  param$hr_volume_search <- fishing_params$search_vol
  
  # detritus coupling on? 1 = yes, 0 = no (det.coupling)
  param$detritus_coupling <- TRUE
  
  # feeding and energy budget parameters
  # Mean log10 predator prey mass ratio 100:1. (q0)
  param$log10_pred_prey_ratio <- 2.0
  
  # 0.6 to 1.0 for log normal prey preference function. (sd.q)
  param$log_prey_pref <- 1.0
  
  # hourly rate volume filtered constant m3*yr-1 for benthos. this value yields 
  # believable growth curve. Approximately 10 times less than A.u (A.v)
  param$hr_vol_filter_benthos <- param$hr_volume_search*0.1    
  
  # exponent for metabolic requirements plus swimming for predators (Ware et al 
  # 1978) (alpha.u)
  param$metabolic_req_pred <- 0.82  
  
  # exponent for whole organism basal (sedentary) metabolic rate (see growth.v) 
  # from Peters (1983) and Brown et al. (2004) for detritivores (alpha.v)
  param$metabolic_req_detritivore <- 0.75
  
  # fraction of ingested food that is defecated (Peters,1983) (def.high)
  param$defecate_prop <- 0.3  
  
  # low <- low quality (K) food, high <- high quality (K) food (def.low)
  param$def_low <- 0.5
  
  # net growth conversion efficiency for organisms in the "predator" spectrum
  # from Ware (1978) (K.u)
  param$growth_pred <- 0.3
  
  # net growth conversion efficiency for organisms in the "detritivore" 
  # spectrum (K.v)
  param$growth_detritivore <- 0.2
  
  # net growth conversion efficiency for detritus (K.d)
  param$growth_detritus <- param$growth_detritivore
  
  #fraction of energy required for maintenance & activity etc. (AM.u and AM.v)
  param$energy_pred <- 0.5
  param$energy_detritivore <- 0.7
  
  # if handling time is > 0 Type II functional response, if = 0 linear (no 
  # predator satiation) - 5.7e-7
  param$handling <- 0
  
  # dynamic reproduction on = 1, off = 0 (repro.on)
  param$dynamic_reproduction <- TRUE
  
  # constant used in Jennings et al. 2008 Proc B to standardize 
  # metabolism-temperature effects for Boltzmann equation. Derived from Simon's
  # fit to Andy Clarke's data
  param$c1 <- 25.22 
  
  # activation energy, eV (E)
  param$activation_energy <- 0.63
  
  # Boltzmann's constant (Boltzmann)
  param$boltzmann <- 8.62*10^-5
  
  # "other" mortality parameters
  
  # residual natural mortality (mu0)
  param$natural_mort <- 0.2
  
  # size at senescence (xs)
  param$size_senescence <- 3
  
  # exponent for senescence mortality (p.s)
  param$exp_senescence_mort <- 0.3
  
  # constant for senescence mortality (k.sm)
  param$const_senescence_mort <- 0.2
  
  # Parameters for numerical integration (size & time discretisation)
  # size increment after discretization for integration (log body weight) (dx)
  param$log_size_increase <- dx
  
  # minimum log10 body size of plankton (xmin)
  param$min_log10_plankton <- xmin
  
  # minimum log10 body size in dynamics predators (x1)
  param$min_log10_pred <- xmin_consumer_u
  
  # minimum log10 body size in dynamics benthic detritivores (x1.det)
  param$min_log10_detritivore <- xmin_consumer_v
  
  # maximum log10 body size of predators (xmax)
  param$max_log10_pred <- xmax
  
  # maximum log10 body size before senescence kicks in (departure form 
  # linearity) (xmax2)
  param$max_log10_senescense <- 4.0
  
  # Log10 size bins (x)
  param$log10_size_bins <- seq(xmin, xmax, dx)
  
  # Number of log10 size bins (Nx)
  param$numb_size_bins <- length(param$log10_size_bins)
  
  # index for minimum log10 predator size (ref)
  param$ind_min_pred_size <-
    which(param$log10_size_bins == param$min_log10_pred)
  
  # index for minimum log10 detritivore size (ref.det)
  param$ind_min_detritivore_size <- 
    which(param$log10_size_bins == param$min_log10_detritivore)
  
  # index in F vector corresponding to smallest size fished in U (Fref.u)
  param$ind_min_fish_pred <- ((param$min_fishing_size_pred-
                                 param$min_log10_plankton)/dx)+1
  
  # index in F vector corresponding to smallest size fished in V (Fref.v)
  param$ind_min_fish_det <- ((param$min_fishing_size_detritivore-
                                param$min_log10_plankton)/dx)+1
  
  #short hand for matrix indexing (idx)
  param$idx <- 2:param$numb_size_bins
  
  if(use_init){
    #(U.init)
    param$init_pred <- pred_initial
    #(V.init)
    param$init_detritivores <- detritivore_initial
    #(W.init)
    param$init_detritus <- detrititus_initial
  }else{
    # arbitrary initial value for detritus (W.init)
    param$init_detritus <- 0.00001
    # If gridded is selected and no initial values for predators or detritivores
    # area provided, then the parameters below are not included in the output.
    # These will need to calculated from gridded files
    if(!gridded){
      # set initial detritivore spectrum (V.init)
      param$init_detritivores <- param$sinking_rate[1]*10^param$int_phy_zoo[1]* 
        10^(param$slope_phy_zoo[1]*param$log10_size_bins)
      
      # (phyto+zoo)plankton + pelagic predator size spectrum (U.init)
      param$init_pred <- 10^param$int_phy_zoo[1]*
        10^(param$slope_phy_zoo[1]*param$log10_size_bins)
    }
  }
  
  param$equilibrium <- equilibrium
  
  return(param)
}

# Build a lookup table for diet preference ------
# Looks at all combinations of predator and prey body size: diet preference 
# (in the predator spectrum only)
phi_f <- function(q, log10_pred_prey_ratio, log_prey_pref){
  phi <- ifelse(q > 0, 
                exp(-(q-log10_pred_prey_ratio)*(q-log10_pred_prey_ratio)/
                      (2*log_prey_pref*log_prey_pref))/
                  (log_prey_pref*sqrt(2.0*pi)),
                0) 
  return(phi)
}

# Build lookup tables for (constant) growth ------
# Considers components which remain constant
gphi_f <- function(pred_prey_matrix, log10_pred_prey_ratio, log_prey_pref){
  return(10^(-pred_prey_matrix)*phi_f(pred_prey_matrix, log10_pred_prey_ratio, 
                                      log_prey_pref))
}	

# Build lookup tables for (constant) mortality ------
mphi_f <- function(rev_pred_prey_matrix, log10_pred_prey_ratio, log_prey_pref,
                   metabolic_req_pred){
  return(10^(metabolic_req_pred*rev_pred_prey_matrix)*
           phi_f(rev_pred_prey_matrix, log10_pred_prey_ratio, log_prey_pref))
}	

# Build lookup table for components of 10^(alpha*x) ------
expax_f <- function(log10_size_bins, metabolic_req_pred){
  return(10^(metabolic_req_pred*log10_size_bins)) 
}


# Gravity model -----
gravitymodel <- function(effort, prop_b, depth, iter){
  # redistribute total effort across grid cells according to proportion of
  # biomass in that grid cell using graivity model, Walters & Bonfil, 1999, 
  # Gelchu & Pauly 2007 ideal free distribution - Blanchard et al 2008
  
  eff <- as.vector(effort)
  d <- unlist(as.vector(depth))
  
  for(j in 1:iter){
    suit <- prop_b*(1-d/max(d))
    # rescale:
    rel_suit <- suit/sum(suit)
    neweffort <- eff+rel_suit*eff
    mult <- sum(eff)/sum(neweffort)
    #gradient drives individuals to best locations at equilibrium (assumed to 
    #reached after 10000 iterations)
    eff <- neweffort*mult
  }
  
  return(eff)
}

  
# Run model per grid cell or averaged over an area ------
sizemodel <- function(params, ERSEM_det_input = F, temp_effect = T,
                      use_init = F){
  with(params,{
    # Model for a dynamical ecosystem comprised of: two functionally distinct 
    # size spectra (predators and detritivores), size structured primary 
    # producers and an unstructured detritus resource pool. 
    # time implicit upwind difference discretization algorithm (from Ken 
    # Andersen)
    # (Numerical Recipes and Richards notes, see "C://...standardised test 
    # implicit.txt" for basic example)
    # U represents abundance density (m-2) of "fish" and V abundance density 
    # (m-2) of "detritivores"
    # W is the pool of detritus (expressed in biomass density, g.m-2) - not 
    # size-based
    # Fish feed on other fish and benthic detritivores with size preference 
    # function 
    # Benthic detritivores feed on detritus produced from pelagic spectrum 
    # Senescence Mortality also included.  Options to include dynamic 
    # reproduction and predator handling time (but not currently used).
    #
    # Code modified for global fishing mortality rate application. 
    # JLB 17/02/2014
    # Code modified to include temperature scaling on senescence and detrital 
    # flux. RFH 18/06/2020
    
    # Input parameters to vary:
    #shorthand for matrix referencing
    idx_new <- (ind_min_detritivore_size+1):numb_size_bins
    #Size bin index
    size_bin_index <- 1:numb_size_bins
    #Size bins
    size_bins_vals <- 10^log10_size_bins
    
    #Create matrices for slope_phy_zoo and log10_size_bins so they can be
    #easily multiplied
    slope_phy_zoo_mat <- matrix(rep(slope_phy_zoo, numb_size_bins), 
                                nrow = numb_size_bins, byrow = T)
    log10_size_bins_mat <- matrix(rep(log10_size_bins, numb_time_steps), 
                                  nrow = numb_size_bins)
    
    # Initialising matrices
    
    #q1 is a square matrix holding the log(predatorsize/preysize) for all 
    #combinations of sizes (q1)
    pred_prey_matrix <- matrix(NA, numb_size_bins, numb_size_bins)
    for(i in size_bin_index){
      pred_prey_matrix[, i] <- log10_size_bins[i] - log10_size_bins
    }
    
    #q2 is the reverse matrix holding the log(preysize/predatorsize) for all 
    #combinations of sizes (q2) - Removed, using -pred_prey_matrix (-q1)
    
    #matrix for recording the two size spectra (V - det, U - pred)
    detritivores <- predators <- array(0, c(numb_size_bins, numb_time_steps+1))
    
    #vector to hold detritus biomass density (g.m-2) (W)
    detritus <- array(0, numb_time_steps+1)
    
    #matrix for keeping track of growth (GG_v, GG_u) and reproduction (R_v, R_u) 
    #from ingested food:
    reprod_det <- reprod_pred <- array(0, c(numb_size_bins, numb_time_steps+1))
    growth_det <- growth_int_pred <- array(0, c(numb_size_bins, 
                                                numb_time_steps+1))
    
    #matrix for keeping track of predation mortality (PM_v, PM_u)
    pred_mort_det <- array(0, c(numb_size_bins, numb_time_steps+1))
    pred_mort_pred <- array(0, c(numb_size_bins, numb_time_steps+1))
    
    #matrix for keeping track of catches (Y_v, Y_u)
    catch_det <- catch_pred <- array(0, c(numb_size_bins, numb_time_steps+1))  
    
    #matrix for keeping track of total mortality (Z_v, Z_u)
    tot_mort_det <- tot_mort_pred <- array(0, c(numb_size_bins,
                                                numb_time_steps+1))
    
    #empty vector to hold fishing mortality rates at each size class at time 
    #(Fvec_v, Fvec_u)
    fishing_mort_det <- fishing_mort_pred <- array(0, c(numb_size_bins, 
                                                        numb_time_steps))
    
    #lookup tables for terms in the integrals which remain constant over time
    #(gphi, mphi)
    constant_growth <- gphi_f(pred_prey_matrix, log10_pred_prey_ratio, 
                              log_prey_pref)
    constant_mortality <- mphi_f(-pred_prey_matrix, log10_pred_prey_ratio, 
                                 log_prey_pref, metabolic_req_pred)
    
    rm(pred_prey_matrix)
    
    #lookup table for components of 10^(metabolic_req_pred*log10_size_bins) 
    #(expax)
    met_req_log10_size_bins <- expax_f(log10_size_bins, metabolic_req_pred)
    
    # Numerical integration
    
    # set up with the initial values from param
    #(phyto+zoo)plankton size spectrum - set initial consumer size spectrum (U)
    
    # time series of intercept of plankton size spectrum (estimated from GCM, 
    # biogeophysical model output or satellite data).		
    ui0 <- matrix(rep((10^int_phy_zoo), numb_size_bins), nrow = numb_size_bins, 
                  byrow = T)
    
    #recruitment at smallest consumer mass
    #continuation of plankton hold constant
    predators[1:(ind_min_pred_size-1), 2:(numb_time_steps+1)] <- 
      (ui0*10^(slope_phy_zoo_mat*log10_size_bins_mat))[1:(ind_min_pred_size-1),]
    predators[,1] <- predators[,2]
    predators[ind_min_pred_size:120, 1] <- 
      init_pred[ind_min_pred_size:120]
    
    #remove time series of intercept of plankton size spectrum not in use
    rm(ui0)
    
    # set initial detritivore spectrum (V)
    detritivores[ind_min_detritivore_size:120, 1] <- 
      init_detritivores[ind_min_detritivore_size:120]
    # set initial detritus biomass density (g.m^-3) (W)
    detritus[1] <- init_detritus
    
    if(use_init){
      # set initial consumer size spectrum from previous run
      predators[ind_min_pred_size:numb_size_bins, 1] <- 
        init_pred[ind_min_pred_size:numb_size_bins]
      # set initial detritivore spectrum from previous run
      detritivores[ind_min_detritivore_size:numb_size_bins, 1] <- 
        init_detritivores[ind_min_detritivore_size:numb_size_bins] 
    }
    
    #other (intrinsic) natural mortality (OM.u, OM.v)
    other_mort_det <- other_mort_pred <- 
      natural_mort*10^(-0.25*log10_size_bins)
    
    #senescence mortality rate to limit large fish from building up in the 
    #system
    #same function as in Law et al 2008, with chosen parameters gives similar 
    #M2 values as in Hall et al. 2006 (SM.u, SM.v)
    senes_mort_det <- senes_mort_pred <- const_senescence_mort*
      10^(exp_senescence_mort*(log10_size_bins-size_senescence))
    
    #Fishing mortality (THESE PARAMETERS NEED TO BE ESTIMATED!) (FVec.u, FVec.v)
    # from Benoit & Rochet 2004 
    # here fish_mort_pred and fish_mort_pred = fixed catchability term for 
    # predators and detritivores to be estimated along with ind_min_det and 
    # ind_min_fish_pred
    fishing_mort_pred[ind_min_fish_pred:numb_size_bins,] <- 
      matrix(rep((fish_mort_pred*effort), 
                 length(ind_min_fish_pred:numb_size_bins)), 
             nrow = length(ind_min_fish_pred:numb_size_bins), byrow = T)
    
    fishing_mort_det[ind_min_fish_det:numb_size_bins,] <- 
      matrix(rep((fish_mort_detritivore*effort), 
                 length(ind_min_fish_det:numb_size_bins)),
             nrow = length(ind_min_fish_det:numb_size_bins), byrow = T)
    
    #output fisheries catches per yr at size (Y_u, Y_v)
    catch_pred[ind_min_fish_pred:numb_size_bins, 1] <-
      fishing_mort_pred[ind_min_fish_pred:numb_size_bins, 1]*
      predators[ind_min_fish_pred:numb_size_bins, 1]*
      size_bins_vals[ind_min_fish_pred:numb_size_bins]
    
    #output fisheries catches per yr at size
    catch_det[ind_min_fish_det:numb_size_bins, 1] <- 
      fishing_mort_det[ind_min_fish_det:numb_size_bins, 1]*
      detritivores[ind_min_fish_det:numb_size_bins, 1]*
      size_bins_vals[ind_min_fish_det:numb_size_bins] 
    
    #Temperature effects taken out of loop because they do not need to be 
    #recalculated at every time step
    if(temp_effect){
      pel_tempeffect <- exp(c1-activation_energy/
                              (boltzmann*(sea_surf_temp+273)))
      ben_tempeffect <- exp(c1-activation_energy/
                              (boltzmann*(sea_floor_temp+273)))
    }else{
      pel_tempeffect <- 1
      ben_tempeffect <- 1
    }
    
    #To be applied to feeding rates for pelagics and benthic groups
    feed_mult_pel <- (hr_volume_search*
                        10^(log10_size_bins*metabolic_req_pred)*pref_pelagic)
    
    feed_mult_ben <- (hr_volume_search*
                        10^(log10_size_bins*metabolic_req_pred)*pref_benthos)
    
    #Ingested food
    growth_prop <- 1-defecate_prop
    
    # High quality food
    high_prop <- 1-def_low
    
    #iteration over time, N [days]
    for(i in 1:numb_time_steps){
      # Calculate Growth and Mortality
      # feeding rates yr-1 (f_pel)
      pred_growth <- (predators[,i]*log_size_increase)%*%(constant_growth)
      
      feed_rate_pel <- pel_tempeffect[i]*as.vector(
        (feed_mult_pel*pred_growth)/(1+handling*feed_mult_pel*pred_growth)) 
      
      # yr-1 (f_ben)
      detrit_growth <- (detritivores[,i]*log_size_increase)%*%(constant_growth)
      
      feed_rate_bent <- pel_tempeffect[i]* 
        as.vector((feed_mult_ben*detrit_growth)/ 
                    (1+handling*feed_mult_ben*detrit_growth))
      
      # yr-1 (f_det)
      detritus_multiplier <- (1/size_bins_vals)*hr_vol_filter_benthos*
        10^(log10_size_bins*metabolic_req_detritivore)*detritus[i]
      
      feed_rate_det <- ben_tempeffect[i]*(detritus_multiplier)/
        (1+handling*detritus_multiplier)
      
      # Predator growth integral (GG_u) yr-1 
      growth_int_pred[, i] <- growth_prop*growth_pred*
        feed_rate_pel+high_prop*growth_detritivore*(feed_rate_bent)
      
      # Reproduction (R_u) yr-1
      if(dynamic_reproduction){
        reprod_pred[, i] <- growth_prop*(1-(growth_pred+energy_pred))*
          (feed_rate_pel)+growth_prop*
          (1-(growth_detritivore+energy_detritivore))*feed_rate_bent
      }
      
      # Predator death integrals 
      #Satiation level of predator for pelagic prey
      sat_pel <- ifelse(feed_rate_pel > 0,
                        feed_rate_pel/(feed_mult_pel*pred_growth),
                        0)
      # yr-1 (PM_u)
      pred_mort_pred[, i] <- as.vector(
        (pref_pelagic*hr_volume_search*met_req_log10_size_bins)*
          (predators[, i]*sat_pel*log_size_increase)%*%(constant_mortality))
      
      # yr-1 (Z_u)
      tot_mort_pred[, i] <- pred_mort_pred[, i]+
        pel_tempeffect[i]*other_mort_pred+senes_mort_pred+fishing_mort_pred[,i]
      
      # Benthos growth integral
      # yr-1 (GG_v)
      growth_det[, i] <- high_prop*growth_detritus*feed_rate_det
      
      #reproduction
      # yr-1 (R_v)
      if(dynamic_reproduction){
        reprod_det[, i] <- high_prop*(1-(growth_detritus+energy_detritivore))*
          (feed_rate_det)
      }
      
      # Benthos death integral
      #Satiation level of predator for benthic prey 
      sat_ben <- ifelse(feed_rate_bent > 0,
                        feed_rate_bent/
                          ((hr_volume_search*
                              10^(log10_size_bins*metabolic_req_detritivore)*
                              pref_benthos)* 
                             detrit_growth),
                        0)
      # yr-1 (PM_v)
      pred_mort_det[, i] <- ifelse(sat_ben > 0,
                                   as.vector(
                                     (pref_benthos*hr_volume_search*
                                        met_req_log10_size_bins)*
                                       (predators[, i]*sat_ben*
                                          log_size_increase)%*%
                                       (constant_mortality)),
                                   0)
      
      # yr-1 (Z_v)
      tot_mort_det[, i] <- pred_mort_det[, i]+
        ben_tempeffect[i]*other_mort_det+senes_mort_det+fishing_mort_det[, i]
      
      #detritus output (g.m-2.yr-1)
      # losses from detritivore scavenging/filtering only:
      output_w <- sum(size_bins_vals*feed_rate_det*detritivores[, i]*
                        log_size_increase)   
      
      #total biomass density defecated by pred (g.m-2.yr-1)
      defbypred <- defecate_prop*(feed_rate_pel)*size_bins_vals*
        predators[, i]+def_low*(feed_rate_bent)*size_bins_vals*predators[, i]
      
      # Increment values of detritus, predators & detritivores for next 
      # timestep  
      
      #Detritus Biomass Density Pool - fluxes in and out (g.m-2.yr-1) of 
      #detritus pool and solve for detritus biomass density in next time step 
      if(!ERSEM_det_input){
        #considering pelagic faeces as input as well as dead bodies from both 
        #pelagic and benthic communities and phytodetritus (dying sinking
        #phytoplankton)
        if(detritus_coupling){
          # pelagic spectrum inputs (sinking dead bodies and faeces) - export 
          # ratio used for "sinking rate" + benthic spectrum inputs (dead stuff
          # already on/in seafloor)
          input_w <- (sinking_rate[i]* 
                        (sum(defbypred[ind_min_pred_size:numb_size_bins]*
                               log_size_increase)+
                           sum(pel_tempeffect[i]*other_mort_pred*predators[, i]*
                                 size_bins_vals*log_size_increase) + 
                           sum(pel_tempeffect[i]*senes_mort_pred*predators[, i]*
                                 size_bins_vals*log_size_increase)) +
                        (sum(ben_tempeffect[i]*other_mort_det*detritivores[, i]*
                               size_bins_vals*log_size_increase) + 
                           sum(ben_tempeffect[i]*senes_mort_det*
                                 detritivores[, i]*size_bins_vals*
                                 log_size_increase)))
        }else{
          input_w <- sum(ben_tempeffect[i]*other_mort_det*detritivores[, i]*
                           size_bins_vals*log_size_increase)+ 
            sum(ben_tempeffect[i]*senes_mort_det*detritivores[, i]*
                  size_bins_vals*log_size_increase)
        }
        
        # get burial rate from Dunne et al. 2007 equation 3
        burial <- input_w*(0.013+0.53*input_w^2/(7+input_w)^2)
        
        # losses from detritivory + burial rate (not including remineralisation
        # bc that goes to p.p. after sediment, we are using realised p.p. as
        # inputs to the model) 
        dW <- input_w-(output_w+burial) 
        #biomass density of detritus g.m-2
        detritus[i+1] <- detritus[i]+dW*timesteps_years
      }
      if(ERSEM_det_input){
        detritus[i+1] <- detritus[i]
      }
      
      # Pelagic Predator Density (nos.m-2)- solve for time + timesteps_years 
      # using implicit time Euler upwind finite difference (help from Ken 
      # Andersen and Richard Law)
      
      # Matrix setup for implicit differencing 
      Ai_u <- Bi_u <- Si_u <- rep(0, numb_size_bins)
      
      Ai_u[idx] <- (1/log(10))*-growth_int_pred[idx-1, i]*
        timesteps_years/log_size_increase
      Bi_u[idx] <- 1+(1/log(10))*growth_int_pred[idx, i]*
        timesteps_years/log_size_increase+tot_mort_pred[idx, i]*timesteps_years
      Si_u[idx] <- predators[idx, i]
      
      # Boundary condition at upstream end 
      Ai_u[ind_min_pred_size] <- 0
      Bi_u[ind_min_pred_size] <- 1
      Si_u[ind_min_pred_size] <- predators[ind_min_pred_size, i]
      
      # apply transfer efficiency of 10% *plankton density at same size
      # reproduction from energy allocation
      if(dynamic_reproduction){
        predators[ind_min_pred_size, i+1] <- predators[ind_min_pred_size, i]+
          (sum(reprod_pred[(ind_min_pred_size+1):numb_size_bins, i]*
                 size_bins_vals[(ind_min_pred_size+1):numb_size_bins]*
                 predators[(ind_min_pred_size+1):numb_size_bins, i]*
                 log_size_increase)*timesteps_years)/
          (log_size_increase*size_bins_vals[ind_min_pred_size])-
          (timesteps_years/log_size_increase)*(1/log(10))*
          (growth_int_pred[ind_min_pred_size, i])*
          predators[ind_min_pred_size, i]-timesteps_years*
          tot_mort_pred[ind_min_pred_size, i]*predators[ind_min_pred_size, i]
      }
      
      #main loop calculation
      for(j in (ind_min_pred_size+1):numb_size_bins){
        predators[j, i+1] <- (Si_u[j]-Ai_u[j]*predators[j-1, i+1])/Bi_u[j]
      }
      
      # Benthic Detritivore Density (nos.m-2)
      Ai_v <- Bi_v <- Si_v <- rep(0, numb_size_bins)
      
      Ai_v[idx_new] <- (1/log(10))*-growth_det[idx_new-1,i]*timesteps_years/
        log_size_increase
      Bi_v[idx_new] <- 1+(1/log(10))*growth_det[idx_new, i]*timesteps_years/
        log_size_increase+tot_mort_det[idx_new, i]*timesteps_years
      Si_v[idx_new] <- detritivores[idx_new, i]
      
      #boundary condition at upstream end
      Ai_v[ind_min_detritivore_size] <- 0
      Bi_v[ind_min_detritivore_size] <- 1
      Si_v[ind_min_detritivore_size] <-
        detritivores[ind_min_detritivore_size, i]
      
      #invert matrix
      #recruitment at smallest detritivore mass
      #hold constant continuation of plankton with sinking rate multiplier
      detritivores[1:ind_min_detritivore_size, i+1] <-
        detritivores[1:ind_min_detritivore_size, i]
      
      # apply a very low of transfer efficiency 1%* total biomass of detritus
      #divided by minimum size
      if(dynamic_reproduction){
        detritivores[ind_min_detritivore_size, i+1] <-
          detritivores[ind_min_detritivore_size, i]+
          sum(reprod_det[idx_new, i]*size_bins_vals[idx_new]*
                detritivores[idx_new, i]*log_size_increase)*timesteps_years/
          (log_size_increase*size_bins_vals[ind_min_detritivore_size])-
          (timesteps_years/log_size_increase)*(1/log(10))*
          (growth_det[ind_min_detritivore_size, i])*
          detritivores[ind_min_detritivore_size, i]-
          timesteps_years*tot_mort_det[ind_min_detritivore_size, i]*
          detritivores[ind_min_detritivore_size, i]
      }
      
      #loop calculation
      for(j in idx_new){
        detritivores[j, i+1] <- ((Si_v[j]-Ai_v[j]*detritivores[j-1, i+1])/
                                   Bi_v[j])
      }
      rm(j)
    }
    #end time iteration
    
    #output fisheries catches per yr at size
    catch_pred[ind_min_fish_pred:numb_size_bins, 2:(numb_time_steps+1)] <-
      fishing_mort_pred[ind_min_fish_pred:numb_size_bins,]*
      predators[ind_min_fish_pred:numb_size_bins, 2:(numb_time_steps+1)]*
      size_bins_vals[ind_min_fish_pred:numb_size_bins]
    
    #output fisheries catches per yr at size
    catch_det[ind_min_fish_det:numb_size_bins, 2:(numb_time_steps+1)] <-
      fishing_mort_det[ind_min_fish_det:numb_size_bins,]*
      detritivores[ind_min_fish_det:numb_size_bins, 2:(numb_time_steps+1)]*
      size_bins_vals[ind_min_fish_det:numb_size_bins]
    
    # Subsetting predator, detritivore and detritus results to exclude 
    # initialisation values (i.e., first timestep)
    predators <- predators[,2:(numb_time_steps+1)]
    detritivores <- detritivores[,2:(numb_time_steps+1)]
    detritus <- detritus[2:(numb_time_steps+1)]
    
    return(list(predators = predators[,],
                growth_int_pred = growth_int_pred[,],
                pred_mort_pred = pred_mort_pred[,],
                detritivores = detritivores[,],
                growth_det = growth_det[,],
                pred_mort_det = pred_mort_det[,],
                catch_pred = catch_pred[,],
                catch_det = catch_det[,],
                detritus = detritus[],
                params = params))
  })
  # end with(params)
}
#end size-based model function


# Running model with time series ----
run_model <- function(fishing_params, dbpm_inputs, withinput = T){
  #Inputs:
  # fishing_params (list) - Fishing parameters produced by the `sizeparam` 
  # function
  # dbpm_inputs (data frame) - Climate and fishing forcing data produced in
  # script 03_processing_effort_fishing_inputs.R
  # withinput (boolean) - Default is TRUE. ????
  #
  #Output:
  # If withinput set to TRUE:
  # input (data frame) - ????
  # If withinput set to FALSE:
  # result_set (data frame) - ???
  
  params <- sizeparam(dbpm_inputs, fishing_params, xmin_consumer_u = -3, 
                      xmin_consumer_v = -3, tstepspryr = 12)
  
  # run model through time
  # TO DO IN SIZEMODEL CODE: make fishing function like one in model template
  result_set <- sizemodel(params)
  size_bins <- 10^params$log10_size_bins
  
  if(withinput){
    lims_pred_bio <- params$ind_min_pred_size:params$numb_size_bins
    # JB:  changed inputs to m2 so no need to divide by depth here
    # Timesteps start from index 2 because the first time step contains 
    # initialisation values
    time_steps <- 2:(params$numb_time_steps+1)
    
    dbpm_inputs$total_pred_biomass <- 
      apply(result_set$predators[lims_pred_bio,]*params$log_size_increase*
              size_bins[lims_pred_bio], 2, sum)
    
    lims_det_bio <- params$ind_min_detritivore_size:params$numb_size_bins
    dbpm_inputs$total_detritivore_biomass <- 
      apply(result_set$detritivores[lims_det_bio,]*params$log_size_increase*
              size_bins[lims_det_bio], 2, sum)
    
    dbpm_inputs$total_detritus <- result_set$detritus
    
    #sum catches (currently in grams per m3 per year, across size classes) 
    # keep as grams per m2, then be sure to convert observed from tonnes per m2
    # per year to g.^-m2.^-yr (for each month)
    dbpm_inputs$total_pred_catch <- 
      apply(result_set$catch_pred[, time_steps]*params$log_size_increase,
            2, sum)
    dbpm_inputs$total_detritivore_catch <-
      apply(result_set$catch_det[, time_steps]*params$log_size_increase, 
            2, sum)
    
    #Calculate total catch
    dbpm_inputs <- dbpm_inputs |>
      mutate(total_catch = total_pred_catch + total_detritivore_catch)
    
    return(dbpm_inputs)
    
  }else{
    
    return(result_set)
  }
}

# Comparing observed and predicted fish biomass ----
getError <- function(fishing_params, dbpm_inputs, year_int = 1950, corr = F, 
                     figure_folder = NULL){
  #Inputs:
  # fishing_params (data frame) - Contains fishing parameters
  # dbpm_inputs (data frame) - Climate and fishing forcing data
  # year_int (numeric) - Default is 1949. Used to subset results from this
  # year onwards
  # corr (boolean) - Default is FALSE. If set to TRUE, the correlation between 
  # predicted and observed values is calculated
  # figure_folder (character) - Optional. Full path to the folder where figures
  # comparing observed and predicted data will be stored
  #
  #Output:
  # If corr set to FALSE:
  # rmse (numeric) - RMSE value between observed and predicted catch
  # If corr set to TRUE:
  # corr_nas (numeric) - Correlation between observed and predicted catch
  
  #If a path to save figures is provided, ensure correlation is calculated
  if(!is.null(figure_folder)){
    if(!corr){
      corr <- T
    }
  }
  
  #Getting name of region from dbpm inputs
  region_name <- unique(dbpm_inputs$region)
  
  #Running model
  result <- run_model(fishing_params, dbpm_inputs)
  
  #Aggregate data by year (mean to conserve units)
  error_calc <- result |> 
    filter(year >= year_int) |> 
    group_by(year) |> 
    summarise(mean_total_catch_yr = mean(total_catch),
              mean_obs_catch_yr = mean(catch_tonnes_area_m2, na.rm = T)) |> 
    #Converting units from tonnes to g
    mutate(mean_obs_catch_yr = mean_obs_catch_yr*1e6,
           #calculate and output error 
           #convert from tonnes to grams (m^-2*yr^-1)
           squared_error = (mean_obs_catch_yr-mean_total_catch_yr)^2)
  
  if(corr){
    corr_nas <- tryCatch({
      #Calculate correlation between observed and predicted catches
      corr_nas <- data.frame(cor = cor(error_calc$mean_obs_catch_yr, 
                                       error_calc$mean_total_catch_yr, 
                                       use = "complete.obs"),
                             #Get number of rows with NA values
                             catchNA = sum(is.na(error_calc$mean_total_catch_yr)),
                             region = region_name)
    },
    error = function(e){
      message(paste0("No mean catch estimates available in this simulation ",
                     "from ", year_int))
      message(conditionMessage(e))
      corr_nas <- data.frame(cor = NA,
                             #Get number of rows with NA values
                             catchNA = sum(is.na(error_calc$mean_total_catch_yr)),
                             region = region_name)
      return(corr_nas)
    },
    warning = function(w){
      message(paste0("No mean catch estimates available in this simulation ",
                     "from ", year_int))
      message(conditionMessage(w))
      corr_nas <- data.frame(cor = NA,
                             #Get number of rows with NA values
                             catchNA = sum(is.na(error_calc$mean_total_catch_yr)),
                             region = region_name)
      return(corr_nas)
    })
    
    #If a path to save figures is provided, create figures and save 
    if(!is.null(figure_folder)){
      if(!dir.exists(figure_folder)){
        dir.create(figure_folder)
      }
      #Plotting predicted and observed catches over time
      p1 <- ggplot()+
        geom_line(data = error_calc, aes(x = year, y = mean_total_catch_yr))+
        geom_point(data = error_calc, aes(x = year, y = mean_obs_catch_yr))+
        scale_x_continuous(breaks = seq(min(error_calc$year), 
                                        max(error_calc$year), by = 10))+
        theme_classic()+ 
        theme(axis.text = element_text(colour = "grey20", size = 12),
              text = element_text(size = 15))+
        labs(x = "Year", y = bquote("Mean catch (g*"~yr^-1*"*"*m^-2*")"))
      
      #Plotting predicted vs observed catches
      p2 <- ggplot()+
        geom_point(data = error_calc, 
                   aes(x = mean_total_catch_yr, y = mean_obs_catch_yr))+
        geom_abline(slope = 1, intercept = 0)+
        theme_classic()+
        theme(axis.text = element_text(colour = "grey20", size = 12),
              text = element_text(size = 15))+
        labs(x = "Predicted", y = "Observed")
      
      #Creating a single plot
      p3 <- p1+p2+
        plot_annotation(title = region_name,
                        theme = theme(plot.title = element_text(size = 16, 
                                                                hjust = 0.5)))
      
      #Creating path to save figure
      f_out <- file.path(figure_folder, 
                         paste0("dbpm_pred_obs_catches_", 
                                str_replace(str_to_lower(region_name), " ",
                                            "-"), ".png"))
      
      #Saving composite figure
      ggsave(filename = f_out, plot = p3, width = 15, height = 10)
    }
    
    #Return correlation values
    return(corr_nas)
    
  }else{
    #Calculate RMSE 
    sum_se <- sum(error_calc$squared_error, na.rm = T)
    count <- sum(!is.na(error_calc$squared_error))
    rmse <- sqrt(sum_se/count)
    
    #Return RMSE
    return(rmse)
  }
}

#Carry out LHS param search ----
LHSsearch <- function(num_iter = 1, search_volume = "estimated", seed = 1234,
                      forcing_file = NULL, gridded_forcing = NULL, 
                      best_param = T, best_val_folder = NULL){
  #Inputs:
  # num_iter (integer) - Number of individual runs. Default is 1.
  # search_volume (character or numeric) - Default is "estimated". It also takes
  # a number that will be used as the area searched by predators for prey
  # seed (positive integer) - Default is 1234. Value used to initialise random-
  # number generation
  # forcing_file (character) - Full path to forcing file. This must be 
  # non-gridded data
  # gridded_forcing (character) - Full path to folder containing gridded forcing
  # files
  # best_param (boolean) - Default is True. If True, then only best fishing
  # parameters (i.e., lowest RME) is returned. If False, all fishing parameters
  # tested are returned
  # best_val_folder (character) - Optional. If provided, it must be the full
  # path to the folder where LHS search results will be saved
  #
  #Output:
  # bestvals (data frame) - Contains the values for LHS parameters that resulted
  # in the best performing model based on RMSE values
  
  #Making function reproducible
  set.seed(seed)
  
  #Construct a hypercube with random numbers. Columns represent five specific 
  #parameters needed to run DBPM
  fishing_params <- data.frame(randomLHS(num_iter, 5))
  #Renaming columns 
  colnames(fishing_params) <- c("fmort_u", "fmort_v", "fminx_u", "fminx_v", 
                                "search_vol")
  
  #Adjust range of mi size params, others go from 0-1
  fishing_params <- fishing_params |> 
    mutate(fminx_u = fminx_u*2, 
           fminx_v = fminx_v*2)
  
  if(is.numeric(search_volume)){
    fishing_params <- fishing_params |> 
      mutate(search_vol = search_volume)
  }else{
    # adjust range of search vol, others go from 0-1
    fishing_params <- fishing_params |> 
      mutate(search_vol = search_vol+0.001)
  }
  
  # use below to select a constant value for search.vol
  if(!is.null(forcing_file)){
    dbpm_inputs <- forcing_file
  }
  if(!is.null(gridded_forcing)){
    dbpm_inputs <- gridded_forcing
  }
  
  # parallelise using 75% of cores available using mclapply
  no_cores <- round((detectCores()*.75), 0)
  fishing_params$rmse <- mclapply(1:nrow(fishing_params), 
                                  FUN = function(i) getError(fishing_params[i,],
                                                             dbpm_inputs), 
                                  mc.cores = no_cores) |> 
    unlist()
  
  #Getting name of region from dbpm inputs
  region_name <- unique(dbpm_inputs$region)
  
  # check param set with lowest error
  if(best_param){
    bestvals <- fishing_params |> 
      filter(rmse == min(rmse, na.rm = T)) |> 
      mutate(region = region_name)
  }else{
    bestvals <- fishing_params |> 
      arrange(rmse) |> 
      mutate(region = region_name)
  }
  
  #Print row with lowest RMSE
  print(bestvals)
  
  #If folder to save values is provided - Save results
  if(!is.null(best_val_folder)){
    #Ensure folder exists
    if(!dir.exists(best_val_folder)){
      dir.create(best_val_folder, recursive = T)
    }
    
    #File path to save output
    fout <- file.path(best_val_folder, 
                      paste0("best-fishing-parameters_", 
                             str_replace(str_to_lower(region_name), " ", "-"),
                             "_searchvol_", search_volume, "_numb-iter_", 
                             num_iter, ".parquet"))
    #Save output
    bestvals |> 
      write_parquet(fout)
  }
  
  return(bestvals)
}


# Correlation and calibration plots ----
corr_calib_plots <- function(fishing_params, dbpm_inputs, 
                             figure_folder = NULL){
  #Inputs:
  # fishing_params (named numeric vector) - Single column with named rows 
  # containing LHS parameters
  # forcing_file (character) - Full path to forcing file. This must be 
  # non-gridded data
  # figure_folder (character) - Optional. Full path to the folder where figures 
  # comparing observed and predicted data will be stored
  #
  #Output:
  # corr_nas (data.frame) - Contains the correlation between predicted and
  # observed values
  
  #Calculate correlations with tuned fishing parameters and save plots
  corr_nas <- getError(fishing_params, dbpm_inputs, year_int = 1950,
                       corr = T, figure_folder)
  
  return(corr_nas)
}


# Size spectrum plots ----
plotsizespectrum <- function(modeloutput, params, timeaveraged = F){
  #Inputs:
  # modeloutput (named list) - Output from `sizemodel` function
  # params (data frame) - Single row containing fishing parameters
  # timeaveraged (boolean) - Default is False. If False, the predator and 
  # detritivore values from the last time step are plotted. If True, the average
  # predator and detritivore values are plotted
  #
  #Output:
  # Size spectrum plot
  #
  with(params, {
    pred_size <- ind_min_pred_size:numb_size_bins
    det_size <- ind_min_detritivore_size:numb_size_bins
    # plot changes in the two size spectra over time
    if(timeaveraged){
      predators <- rowMeans(modeloutput$predators[pred_size], na.rm = T)
      detritivores <- rowMeans(modeloutput$detritivores[det_size], na.rm = T)
    }else{
      predators <- modeloutput$predators[pred_size, numb_time_steps]
      detritivores <- modeloutput$detritivores[det_size, numb_time_steps]
    }
    
    maxy <- max(log10(predators), na.rm = T)
    miny <- -20
    
    plot(log10_size_bins[pred_size], log10(predators), type = "l", 
         col = "#004488", cex = 1.6, ylab = "", xlab = "",
         xlim = c(min_log10_detritivore, max_log10_pred), 
         ylim = c(miny, maxy))
    points(log10_size_bins[det_size], log10(detritivores), lty = 2,
           type = "l", col = "#cc3311", cex = 1.6, ylab = "", xlab = "")
    title(ylab = expression("log abundance density [m"^-3* "]"), line = 2.5,
          xlab = "log body mass [g]")
    text(min_log10_detritivore+abs(0.05*min_log10_detritivore),
         miny+abs(0.3*miny), paste("pel.pref = ", round(pref_pelagic, 2),
                                   sep = ""), pos = 4)
    text(min_log10_detritivore+abs(0.05*min_log10_detritivore),
         miny+abs(0.2*miny), paste("ben.pref = ", round(pref_benthos, 2),
                                   sep = ""), pos = 4)
    legend(max_log10_pred-0.3*(max_log10_pred-min_log10_detritivore),
           maxy-0.04*maxy, c("Predators", "Detritivores"),
           col = c("#004488", "#cc3311"), lwd = 1.5, lty = c(1, 2))
  })
}	# end plot function
