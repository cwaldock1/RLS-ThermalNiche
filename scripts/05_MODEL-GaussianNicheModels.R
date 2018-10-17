# Fourth script in thermal-niche RLS analysis for net thermal niche shape analysis.
# This section models the new thermal niche parameters across standardised species ecological performance and thermal ranges . 

# This includes 
# 1. Bayesian models fit in JAGS estimating ecological performance across thermal gradients for 
#      a. P-occ-abun
#      b. P-abun
#      b. P-occ
# 2. The above models fit seperately for tropical and temperate species. 
# 3. The above models fit between 0.8-0.99 quantiles 

# Load libraries and packages ---- 
library(dplyr)
library(tidyr)
library(R2jags)
library(MCMCvis)
library(coda)
library(gridExtra)
library(qgam)
library(ggplot2)

# Load RLS_All from script 2 ----
ThermalNicheData <- readRDS(file = 'data_upload/ThermalNicheData.rds')
RLS_All          <- readRDS(file = 'data_derived/RLS_20-For-Qgams-2019-09-28.rds')
AbundanceDistribution <- read.csv(file = 'data_derived/AbundanceDistributions.csv')
RLS_All <- left_join(RLS_All, AbundanceDistribution)
All_Gams_Preds <- readRDS(file = 'data_derived/All_Gams_Preds_2018-10-08.rds')
# ----------------------------------------------------------------------------

# Figures of Qgams outputs across different abundance centre patterns for panel insets later ----

# Organise data
All_Gams_Preds <- left_join(All_Gams_Preds, AbundanceDistribution)
All_Gams_Preds$ThermalGuild <- ifelse(All_Gams_Preds$Topt < 23, 'temperate', 'tropical')
All_Gams_Preds <- left_join(All_Gams_Preds, ThermalNicheData %>% dplyr::select(SpeciesName, ConfidenceCombined, T_Gam_edf, T_Lower_0.1, T_Upper_0.1, T_Lower_0.25, 
                                                                               T_Upper_0.25, T_Lower_0.5, T_Upper_0.5, T_SD_Upper, T_SD_Lower))
All_Gams_Preds <- All_Gams_Preds 

All_Gams_Preds2 <- All_Gams_Preds %>% group_by(SpeciesName) %>% nest() %>% mutate(ScaledAbundance = purrr::map(data, ~as.numeric(.$Abundance / max(.$Abundance))), 
                                                                                  CentredTemp     = purrr::map(data, ~ 
                                                                                                                 (.$MeanTemp_CoralWatch - .$Topt[1]) / 
                                                                                                                 (mean(c(.$T_SD_Upper[1], .$T_SD_Lower[1]))))
                                                                                  ) %>% 
  unnest()

All_Gams_Preds2$AbundanceDistributionCategory <- revalue(All_Gams_Preds2$AbundanceDistributionCategory, c('Abundance centre' = 'abundance-centre', 
                                                                                                          'Ramped cool-edge' = 'cool-skew',
                                                                                                          'Ramped warm-edge' = 'warm-skew'))

pdf('figures_new/ThermalDistributionsFigure_allConf.pdf', height = 4, width = 8)
ggplot(All_Gams_Preds2 %>% filter(AbundanceDistributionCategory != 'No trend')) + 
  geom_line(aes(x = CentredTemp, y = ScaledAbundance, col = ThermalGuild, group = SpeciesName, alpha = as.factor(as.numeric(ConfidenceCombined)))) + 
  facet_grid(ThermalGuild~AbundanceDistributionCategory, scales = 'free_x') + 
  ylim(c(0,1)) + xlim(c(-4, 4)) + 
  scale_colour_manual(values = c('dark blue', 'dark orange')) + 
  theme_bw() + theme(panel.grid = element_blank(), legend.position = 'none', aspect.ratio = 1, 
                     axis.title.y = element_text(size = 13),
                     axis.title.x = element_text(size = 13), 
                     axis.text = element_text(size = 13),
                     text = element_text(size = 13), 
                     #strip.text = element_blank(), 
                     axis.line.y = element_line(), 
                     axis.line.x = element_line(), 
                     legend.background = element_blank()) + 
  xlab(NULL) + ylab(NULL) + 
  scale_alpha_manual(values = c(0.1,0.1,0.1,0.2,0.3,0.5,1)) + 
  scale_y_continuous(breaks = seq(0, 1, 0.5), limits = c(0,1)) + 
  scale_x_continuous(breaks = seq(-3, 3, 3), limits = c(-4,4))
dev.off()

pdf('figures_new/ThermalDistributionsFigure_allConf_V2.pdf', height = 4, width = 8)
ggplot(All_Gams_Preds2 %>% filter(AbundanceDistributionCategory != 'No trend')) + 
  geom_line(aes(x = CentredTemp, y = ScaledAbundance, col = ThermalGuild, group = SpeciesName, alpha = as.factor(as.numeric(ConfidenceCombined)))) + 
  facet_grid(ThermalGuild~AbundanceDistributionCategory, scales = 'free_x') + 
  ylim(c(0,1)) + xlim(c(-4, 4)) + 
  scale_colour_manual(values = c('dark blue', 'dark orange')) + 
  theme_bw() + theme(panel.grid = element_blank(), legend.position = 'none', aspect.ratio = 1, 
                     axis.title.y = element_text(size = 13),
                     axis.title.x = element_text(size = 13), 
                     axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     axis.line = element_blank(),
                     text = element_text(size = 13), 
                     strip.text = element_blank(), 
                     axis.line.y = element_line(), 
                     axis.line.x = element_line(), 
                     legend.background = element_blank()) + 
  xlab(NULL) + ylab(NULL) + 
  scale_alpha_manual(values = c(0.1,0.1,0.1,0.2,0.3,0.5,1)) + 
  scale_y_continuous(breaks = seq(0, 1, 0.5), limits = c(0,1)) + 
  scale_x_continuous(breaks = seq(-3, 3, 3), limits = c(-4,4))
dev.off()

ggplot(All_Gams_Preds2 %>% filter(AbundanceDistributionCategory == 'No trend')) + 
  geom_line(aes(x = CentredTemp, y = ScaledAbundance, col = ThermalGuild, group = SpeciesName, alpha = as.factor(as.numeric(ConfidenceCombined)))) + 
  facet_grid(ThermalGuild~AbundanceDistributionCategory, scales = 'free_x') + 
  ylim(c(0,1)) + xlim(c(-4, 4)) + 
  scale_colour_manual(values = c('dark blue', 'dark orange')) + 
  xlab(NULL) + ylab(NULL) + 
  scale_alpha_manual(values = c(0.1,0.1,0.1,0.2,0.3,0.5,1)) + 
  scale_y_continuous(breaks = seq(0, 1, 0.5), limits = c(0.5,1)) + 
  scale_x_continuous(breaks = seq(-3, 3, 3), limits = c(-4,4))

# ----

# DEFINE JAGS MODELS ----
# JAGS TPC Model for abundance and occupancy + abundance data ----
# Step 3: Formulate JAGS modelling code
sink("TPC.txt")
cat("
    model{
    
    # Priors for TPC parameters
    T_Opt ~ dnorm(0, 0.001)
    T_SD  ~ dunif(0, 10)
    T_SD_2  ~ dunif(0, 10)
    # T_Upper ~ dnorm(0, 0.001)I(0,10)
    MaxAbun ~ dnorm(1, 0.001)I(0,1)
    Var ~ dnorm(1, 0.001)
    
    # TPC Model
    for(i in 1:nTemp){ 	## loop through sites
    Abundance[i] ~ dnorm(y[i], Var)
    
    # Temp.hi[i]   <-  step(Temp[i]-T_Opt) # Define your breakpoint
    #  y[i]         <-  (1-Temp.hi[i]) *   MaxAbun*(exp(-(((Temp[i] - T_Opt) / (T_SD)))^2)) +      # Model below Topt
    #                   (Temp.hi[i])   *   MaxAbun*(1-((Temp[i] - T_Opt) / (T_Opt - T_Upper))^2)   # Model above Topt
    
    Temp.hi[i]   <-  step(Temp[i]-T_Opt) # Define your breakpoint
    y[i]         <-  (1-Temp.hi[i])   *     MaxAbun*(exp(-(((Temp[i] - T_Opt) / (T_SD)))^2)) +      # Model below Topt
    (Temp.hi[i])     *     MaxAbun*(exp(-(((Temp[i] - T_Opt) / (T_SD_2)))^2))      # Model above Topt
    
    }
    
    # Solve for Topt
    #    YMax <- MaxAbun*(exp(-(((T_Opt - T_Opt) / (T_SD)))^2))
    
    # Make performance estimates from parameters 
    for(j in 1:nPerformance){
    Temp.hi2[j]   <-  step(Temp_pred[j]-T_Opt) # Define your breakpoint
    Performance[j] <- (1-Temp.hi2[j])   *     MaxAbun*(exp(-(((Temp_pred[j] - T_Opt) / (T_SD)))^2)) +      # Model below Topt
    (Temp.hi2[j])     *     MaxAbun*(exp(-(((Temp_pred[j] - T_Opt) / (T_SD_2)))^2))      # Model above Topt
    
    }
    
    # Estimate skew 
    Skew <- T_SD_2 - T_SD
    
    # Estimate niche breadth in sd 
    T_Breadth <- (((T_SD + T_SD_2))/2) 
    
    #3. Discrepancy measures: 
    #   Pearson residuals PRes and ordinary residuals E
    y_hat <- mean(Abundance)
    
    for (i in 1:nTemp) {
    # VarY[i] <- y[i] * (1 - y[i])
    # PRes[i] <- (Abundance[i] - y[i]) / sqrt(VarY[i])   
    S_res[i]    <- (Abundance[i] - y[i])^2
    S_tot[i]    <- (Abundance[i] - y_hat)^2
    }
    
    R_squared <- 1 - (sum(S_res) / sum(S_tot))
    
    }
    ",fill = TRUE)
sink()


# JAGS TPC Model with logit link for occupancy only model ----
sink("TPC_Occupancy.txt")
cat("
    model{
    
    # Priors for TPC parameters
    T_Opt ~ dnorm(0, 0.001)
    T_SD  ~ dunif(0, 10)
    T_SD_2  ~ dunif(0, 10)
    MaxAbun ~ dnorm(1, 0.001)I(0,1)
    # Var ~ dnorm(1, 0.001)
    
    # TPC Model
    for(i in 1:nTemp){ 	# loop through sites
    Occupancy[i] ~ dbern(y[i])
    Temp.hi[i]  <-  step(Temp[i]-T_Opt) # Define your breakpoint
    logit(y[i]) <-  (1-Temp.hi[i])   *     MaxAbun*(exp(-(((Temp[i] - T_Opt) / (T_SD)))^2)) +      # Model below Topt
    (Temp.hi[i])     *     MaxAbun*(exp(-(((Temp[i] - T_Opt) / (T_SD_2)))^2))      # Model above Topt
    }
    
    # Make performance estimates from parameters 
    for(j in 1:nPerformance){
    Temp.hi2[j]   <-  step(Temp_pred[j]-T_Opt) # Define your breakpoint
    Performance[j] <- (1-Temp.hi2[j])   *     MaxAbun*(exp(-(((Temp_pred[j] - T_Opt) / (T_SD)))^2)) +      # Model below Topt
    (Temp.hi2[j])     *     MaxAbun*(exp(-(((Temp_pred[j] - T_Opt) / (T_SD_2)))^2))      # Model above Topt
    }
    
    # Estimate skew 
    Skew <- T_SD_2 - T_SD
    
    # Estimate niche breadth in sd 
    T_Breadth <- (((T_SD + T_SD_2))/2)
    
    #3. Discrepancy measures: 
    #   Pearson residuals PRes and ordinary residuals E
    for (i in 1:nTemp) {
    E[i]    <- Occupancy[i]  - y[i]
    }
    }
    ",fill = TRUE)
sink()




# ----------------------------------------------------------------------------

# DEFINE FUNCTIONS TO RUN AND HANDLE OUTPUTS OF JAGS MODELS ----
# Function to fit JAGS models to quantiles ----
# Function to fit jags models to defined quantiles on certain data and then save with predefined colour scale for plotting. 
RunJAGSModels <- function(modeldata,
                          params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var"),
                          Quantile = c(0.99),
                          Subset = NULL, # Define which subset we are talking about. LABEL ONLY
                          ThermalGuild = 'temperate',# Define which subset we are talking about
                          Scale = 'Species', # Produces different aggregations of data depending on intput here. 
                          n.thin     = 10, 
                          n.chains   = 3,
                          n.burnin   = 400,
                          n.iter     = 500
){
  
  # Define inputs to for loop.
  QuantileDataList <- list()
  JagsDataList <- list()
  TPC_Model <- list()
  PlotDataJags <- list()
  inits  <- function(){list(T_Opt = 0, T_SD = 1, T_SD_2 = 1, T_Upper = 1)}
  
  # Run for loop to estimate model over each level of quantile. 
  for(i in 1:length(Quantile)){
    
    # Create data. 
    #QuantileDataList[[i]] <- RLS_All_HighConf_Temp %>% group_by(Species_MeanTemp_CoralWatch_V2) %>% do(ScaledLogAbundanceAdult40 = quantile(.$ScaledLogAbundanceAdult40, Quantile[i])) %>% unnest(ScaledLogAbundanceAdult40) %>% na.omit()
    
    
    if(Scale == 'Species'){
      # Aggragte by species and temperature groups (species have equal contributions)
      QuantileDataList[[i]] <- modeldata %>% 
        
        # Take quantile for each species and temperature group
        group_by(Species_MeanTemp_CoralWatch_V2, SpeciesName, ThermalGuild) %>% 
        do(ScaledLogAbundanceAdult40 = quantile(.$ScaledLogAbundanceAdult40, Quantile[i])) %>% 
        unnest(ScaledLogAbundanceAdult40) %>% 
        group_by(Species_MeanTemp_CoralWatch_V2, ThermalGuild) #%>% 
      
      # Take average quantile across all species
      # do(ScaledLogAbundanceAdult40    = mean(.$ScaledLogAbundanceAdult40, na.rm = T), 
      #   ScaledLogAbundanceAdult40_SD = sd(.$ScaledLogAbundanceAdult40, na.rm = T)) %>% 
      #unnest(ScaledLogAbundanceAdult40, ScaledLogAbundanceAdult40_SD) 
      
    }else{
      
      # Aggragte over all samples (species have different contributions) 
      QuantileDataList[[i]] <- modeldata %>% 
        
        # Take quantile for each species and temperature group
        group_by(Species_MeanTemp_CoralWatch_V2) %>% 
        do(ScaledLogAbundanceAdult40 = quantile(.$ScaledLogAbundanceAdult40, Quantile[i])) %>% 
        unnest(ScaledLogAbundanceAdult40)
      
    }
    
    # Create jags data list. 
    JagsDataList[[i]]  <- list(Abundance     = QuantileDataList[[i]]$ScaledLogAbundanceAdult40,
                               Temp          = QuantileDataList[[i]]$Species_MeanTemp_CoralWatch_V2,
                               nTemp         = nrow(QuantileDataList[[i]]), 
                               Temp_pred     = seq(min(QuantileDataList[[i]]$Species_MeanTemp_CoralWatch_V2), max(QuantileDataList[[i]]$Species_MeanTemp_CoralWatch_V2), length.out = 500),
                               nPerformance  = 500)
    
    # Run the jags models
    TPC_Model[[i]]   <- jags(data       = JagsDataList[[i]],
                             inits      = inits,
                             parameters = params,
                             model      = "TPC.txt",
                             n.thin     = n.thin, 
                             n.chains   = n.chains,
                             n.burnin   = n.burnin,
                             n.iter     = n.iter)
    
    # Store the outputs plotlines
    PlotDataJags[[i]] <- data.frame(
      Performance       = TPC_Model[[i]]$BUGSoutput$mean$y,
      Performance_Upper = apply(TPC_Model[[i]]$BUGSoutput$sims.list$y, 2, function(x) quantile(x, 0.975)),
      Performance_Lower = apply(TPC_Model[[i]]$BUGSoutput$sims.list$y, 2, function(x) quantile(x, 0.025)), 
      Quantile = as.factor(Quantile[i]), 
      Species_MeanTemp_CoralWatch_V2 = QuantileDataList[[i]]$Species_MeanTemp_CoralWatch_V2, 
      ScaledLogAbundanceAdult40   = QuantileDataList[[i]]$ScaledLogAbundanceAdult40,
      ThermalGuild = ThermalGuild)
    if(is.null(Subset)){NULL}else{PlotDataJags[[i]]$Subset = Subset[i]}
    
  }
  
  ## Aggregate data output for plotting
  # Integrate with data to form ggplot. 
  PlotDataJags <- do.call(rbind, PlotDataJags)
  
  return(list(data = PlotDataJags, models = TPC_Model))
  
}

# Function to fit JAGS mdoel to occupancy patterns using a logit link function for same mathmatical function describing TPC. 
RunJAGSModels_OCC <- function(modeldata,
                              params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var"),
                              Quantile = c(0.99),
                              Subset = NULL, # Define which subset we are talking about. LABEL ONLY
                              ThermalGuild = 'temperate',# Define which subset we are talking about
                              Scale = 'Species', # Produces different aggregations of data depending on intput here. 
                              n.thin     = 10, 
                              n.chains   = 3,
                              n.burnin   = 400,
                              n.iter     = 500
){
  
  # Define inputs to for loop.
  i = 1
  PresentModelData <- data.frame()
  JagsDataList <- list()
  TPC_Model <- data.frame()
  PlotDataJags <- list()
  inits  <- function(){list(T_Opt = 0, T_SD = 1, T_SD_2 = 1, T_Upper = 1)}
  
  # Aggregate by species and temperature groups (species have equal contributions) and take mean of occupancies 
  
  #DataTest <-  modeldata %>% filter(SpeciesName == 'Abudefduf troschelii')
  
  if(Scale == 'Species'){
    PresentModelData <- modeldata %>% 
      # Take mean for each species occupancy rates
      group_by(Species_MeanTemp_CoralWatch_V2, SpeciesName, ThermalGuild) %>% 
      do(Mean_SpeciesPresence = mean(.$Presence, na.rm = T)) %>% 
      unnest(Mean_SpeciesPresence) }else{
        
        PresentModelData <- modeldata %>% 
          # Take mean for each species occupancy rates
          group_by(Species_MeanTemp_CoralWatch_V2, ThermalGuild) %>% 
          do(Mean_SpeciesPresence = mean(.$Presence, na.rm = T)) %>% 
          unnest(Mean_SpeciesPresence) 
      }
  
  
  # Create jags data list. 
  JagsDataList[[i]]  <- list(Abundance = PresentModelData$Mean_SpeciesPresence,
                             Temp          = PresentModelData$Species_MeanTemp_CoralWatch_V2,
                             nTemp         = nrow(PresentModelData), 
                             Temp_pred     = seq(min(PresentModelData$Species_MeanTemp_CoralWatch_V2), max(PresentModelData$Species_MeanTemp_CoralWatch_V2), length.out = 500),
                             nPerformance  = 500)
  
  # Run the jags models
  TPC_Model   <- jags(data       = JagsDataList[[i]],
                      inits      = inits,
                      parameters = params,
                      model      = "TPC.txt",
                      n.thin     = n.thin, 
                      n.chains   = n.chains,
                      n.burnin   = n.burnin,
                      n.iter     = n.iter)
  
  # Store the outputs plotlines
  PlotDataJags[[i]] <- data.frame(
    Performance       = TPC_Model$BUGSoutput$mean$y,
    Performance_Upper = apply(TPC_Model$BUGSoutput$sims.list$y, 2, function(x) quantile(x, 0.975)),
    Performance_Lower = apply(TPC_Model$BUGSoutput$sims.list$y, 2, function(x) quantile(x, 0.025)), 
    Species_MeanTemp_CoralWatch_V2 = PresentModelData$Species_MeanTemp_CoralWatch_V2, 
    Mean_SpeciesPresence   = PresentModelData$Mean_SpeciesPresence,
    ThermalGuild = ThermalGuild)
  
  if(is.null(Subset)){NULL}else{PlotDataJags[[i]]$Subset = Subset}
  
  
  
  ## Aggregate data output for plotting
  # Integrate with data to form ggplot. 
  PlotDataJags <- do.call(rbind, PlotDataJags)
  
  return(list(data = PlotDataJags, models = TPC_Model))
  
}


# Function to obtain the posterior distributions of defined parameters. 
ObtainPosterior <- function(BugsList = list(), ModelData = list(), Parameter = character(), Subset = character(), colour = character()){
  
  OutputPosteriorData  <- list()
  for(i in 1:length(BugsList)){
    
    OutputPosteriorData[[i]]  <- data.frame(Value = BugsList[[i]]$sims.list[[Parameter]], Parameter = Parameter, ThermalGuild = unique(ModelData[[i]][['ThermalGuild']]), Subset = Subset[i], colour = colour[i])
    
  }
  
  OutputPosteriorData <- do.call('rbind', OutputPosteriorData)
  return(OutputPosteriorData)
  
}




# ----------------------------------------------------------------------------

# STANDARDISE DATA ACROSS ECOLOGICAL PERFORMANCE AND TEMPERATURE GRADIENTS WITHIN EACH SPECIES ---- 
# Run data standardisations ----
# Bind in thermal niche parameters to standardised species temperature and abundance ranges. 
RLS_All_JAGS <- RLS_All %>% left_join(., ThermalNicheData %>% dplyr::select(SpeciesName, Topt, T_SD_Upper, T_SD_Lower, ThermalGuild))


# Centre temperature data within species
RLS_All_JAGS_scales <- RLS_All_JAGS %>% 
  group_by(SpeciesName, ThermalGuild) %>% 
  nest() %>% 
  
  # Scale and transpose data to be comparable across
  mutate(Species_MeanTemp_CoralWatch  = purrr::map(data, ~ (.$MeanTemp_CoralWatch - unique(.$Topt))/mean(c(.$T_SD_Upper , .$T_SD_Lower))), 
         ScaledLogAbundanceAdult40 = purrr::map(data, ~ log(.$AbundanceAdult40+1)/max(log(.$AbundanceAdult40+1)))) %>% 
  unnest(Species_MeanTemp_CoralWatch, ScaledLogAbundanceAdult40)

# Assign columns to data. 
RLS_All_JAGS$Species_MeanTemp_CoralWatch     <- RLS_All_JAGS_scales$Species_MeanTemp_CoralWatch
RLS_All_JAGS$ScaledLogAbundanceAdult40       <- RLS_All_JAGS_scales$ScaledLogAbundanceAdult40
RLS_All_JAGS$Species_MeanTemp_CoralWatch_V2  <- plyr::round_any(RLS_All_JAGS_scales$Species_MeanTemp_CoralWatch, 0.1)
RLS_All_JAGS$T_Opt                           <- RLS_All_JAGS$Topt
RLS_All_JAGS$T_SD_2                          <- RLS_All_JAGS$T_SD_Upper
RLS_All_JAGS$T_SD                            <- RLS_All_JAGS$T_SD_Lower

# Fit one model to each thermal guild. 
TempData_JAGS <- RLS_All_JAGS %>% filter(ThermalGuild == 'temperate') #%>% filter(AbundanceAdult40 > 0)
TropData_JAGS <- RLS_All_JAGS %>% filter(ThermalGuild == 'tropical')  #%>% filter(AbundanceAdult40 > 0)
AllData_JAGS <- rbind(TropData_JAGS, TempData_JAGS)

# What is the mean and sd of abundances when standardised in this way? 
mean(AllData_JAGS$ScaledLogAbundanceAdult40[AllData_JAGS$ScaledLogAbundanceAdult40 != 0])
sd(AllData_JAGS$ScaledLogAbundanceAdult40[AllData_JAGS$ScaledLogAbundanceAdult40 != 0])
# ----------------------------------------------------------------------------

# SET UP MODEL DATA ----
# Data subsetting ----
AllData_JAGS$AbundanceAdult40_log <- log(AllData_JAGS$AbundanceAdult40 + 1)

AllData_JAGS_No0s <- AllData_JAGS[AllData_JAGS$ScaledLogAbundanceAdult40 > 0, ]
AllData_JAGS_No0s_Trop <- AllData_JAGS_No0s[which(AllData_JAGS_No0s$ThermalGuild == 'tropical'), ]
AllData_JAGS_No0s_Temp <- AllData_JAGS_No0s[which(AllData_JAGS_No0s$ThermalGuild == 'temperate'), ]

LessThan50 <- AllData_JAGS_No0s %>%  group_by(Species_MeanTemp_CoralWatch_V2) %>% do(N_Spp = length(unique(.$SpeciesName))) %>% unnest() #%>% filter(N_Spp > 100)
LessThan50_Trop <- AllData_JAGS_No0s_Trop %>%  group_by(Species_MeanTemp_CoralWatch_V2) %>% do(N_Spp = length(unique(.$SpeciesName))) %>% unnest() #%>% filter(N_Spp > 100)
LessThan50_Temp <- AllData_JAGS_No0s_Temp %>%  group_by(Species_MeanTemp_CoralWatch_V2) %>% do(N_Spp = length(unique(.$SpeciesName))) %>% unnest() #%>% filter(N_Spp > 100)

AllData_JAGS_No0s <- AllData_JAGS_No0s %>% filter(Species_MeanTemp_CoralWatch_V2 %in% LessThan50$Species_MeanTemp_CoralWatch_V2)
AllData_JAGS_No0s_Trop <- AllData_JAGS_No0s_Trop %>% filter(Species_MeanTemp_CoralWatch_V2 %in% LessThan50_Trop$Species_MeanTemp_CoralWatch_V2)
AllData_JAGS_No0s_Temp <- AllData_JAGS_No0s_Temp %>% filter(Species_MeanTemp_CoralWatch_V2 %in% LessThan50_Temp$Species_MeanTemp_CoralWatch_V2)

# ----

# FIT MODELS FOR P-Abun ---- 
# Create directory for output of JAGS files ----
dir.create('data_derived/jagsmodels/')

# TropJAGS_0.99_ABUN ----
TropJAGS_0.99_ABUN <- RunJAGSModels(modeldata = TropData_JAGS %>% filter(AbundanceAdult40 > 0),
                                    params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var", "R_squared"),
                                    Quantile = c(0.99),
                                    Subset = 'Quantile 0.99', # Define which subset we are talking about
                                    ThermalGuild = 'tropical', 
                                    Scale = 'Species', 
                                    n.thin     = 5, 
                                    n.chains   = 4,
                                    n.burnin   = 2500,
                                    n.iter     = 10000)
saveRDS(TropJAGS_0.99_ABUN, 'data_derived/jagsmodels/TropJAGS_0.99_ABUN.rds')
TropJAGS_0.99_ABUN <- readRDS('data_derived/jagsmodels/TropJAGS_0.99_ABUN.rds')

round(TropJAGS_0.99_ABUN$models[[1]]$BUGSoutput$summary['T_Opt',], 3)
round(TropJAGS_0.99_ABUN$models[[1]]$BUGSoutput$summary['T_Breadth',], 3)
round(TropJAGS_0.99_ABUN$models[[1]]$BUGSoutput$summary['Skew',], 3)
round(TropJAGS_0.99_ABUN$models[[1]]$BUGSoutput$summary['MaxAbun',], 3)
round(TropJAGS_0.99_ABUN$models[[1]]$BUGSoutput$summary['R_squared',], 3)
TropJAGS_0.99_ABUN$models[[1]]$BUGSoutput$DIC
TropJAGS_0.99_ABUN$models[[1]]$BUGSoutput$pD

# TempJAGS_0.99_ABUN ----
TempJAGS_0.99_ABUN <- RunJAGSModels(modeldata = TempData_JAGS %>% filter(AbundanceAdult40 > 0),
                                    params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var", "R_squared"),
                                    Quantile = c(0.99),
                                    Subset = 'Quantile 0.99', # Define which subset we are talking about
                                    ThermalGuild = 'temperate', 
                                    Scale = 'Species', 
                                    n.thin     = 5, 
                                    n.chains   = 4,
                                    n.burnin   = 2500,
                                    n.iter     = 10000)
saveRDS(TempJAGS_0.99_ABUN, 'data_derived/jagsmodels/TempJAGS_0.99_ABUN.rds')
TempJAGS_0.99_ABUN <- readRDS('data_derived/jagsmodels/TempJAGS_0.99_ABUN.rds')

round(TempJAGS_0.99_ABUN$models[[1]]$BUGSoutput$summary['T_Opt',], 3)
round(TempJAGS_0.99_ABUN$models[[1]]$BUGSoutput$summary['T_Breadth',], 3)
round(TempJAGS_0.99_ABUN$models[[1]]$BUGSoutput$summary['Skew',], 3)
round(TempJAGS_0.99_ABUN$models[[1]]$BUGSoutput$summary['MaxAbun',], 3)
round(TempJAGS_0.99_ABUN$models[[1]]$BUGSoutput$summary['R_squared',], 3)
TempJAGS_0.99_ABUN$models[[1]]$BUGSoutput$DIC
TempJAGS_0.99_ABUN$models[[1]]$BUGSoutput$pD

# TropJAGS_0.99_ABUN_AGG ----
TropJAGS_0.99_ABUN_AGG <- RunJAGSModels(modeldata = TropData_JAGS %>% filter(AbundanceAdult40 > 0),
                                        params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var", "R_squared"),
                                        Quantile = c(0.99),
                                        Subset = 'Quantile 0.99', # Define which subset we are talking about
                                        ThermalGuild = 'tropical', 
                                        Scale = 'Aggregated', 
                                        n.thin     = 5, 
                                        n.chains   = 4,
                                        n.burnin   = 2500,
                                        n.iter     = 10000)
saveRDS(TropJAGS_0.99_ABUN_AGG, 'data_derived/jagsmodels/TropJAGS_0.99_ABUN_AGG.rds')
TropJAGS_0.99_ABUN_AGG <- readRDS('data_derived/jagsmodels/TropJAGS_0.99_ABUN_AGG.rds')

round(TropJAGS_0.99_ABUN_AGG$models[[1]]$BUGSoutput$summary['R_squared',], 3)


# TempJAGS_0.99_ABUN_AGG ----
TempJAGS_0.99_ABUN_AGG <- RunJAGSModels(modeldata = TempData_JAGS %>% filter(AbundanceAdult40 > 0),
                                        params = c("T_Opt", "T_SD", "T_SD_2", "T_Upper", "E", "y", 'YMax', "MaxAbun","Skew","T_Breadth", "Var", "R_squared"),
                                        Quantile = c(0.99),
                                        Subset = 'Quantile 0.99', # Define which subset we are talking about
                                        ThermalGuild = 'temperate', 
                                        Scale = 'Aggregated', 
                                        n.thin     = 5, 
                                        n.chains   = 4,
                                        n.burnin   = 2500,
                                        n.iter     = 10000)
saveRDS(TempJAGS_0.99_ABUN_AGG, 'data_derived/jagsmodels/TempJAGS_0.99_ABUN_AGG.rds')
TempJAGS_0.99_ABUN_AGG <- readRDS('data_derived/jagsmodels/TempJAGS_0.99_ABUN_AGG.rds')

round(TempJAGS_0.99_ABUN_AGG$models[[1]]$BUGSoutput$summary['R_squared',], 3)

# Figure integrating the above models ----

# Extract data from plots for plotting
TropJAGS_0.99_ABUN_Data    <- TropJAGS_0.99_ABUN[['data']]
TempJAGS_0.99_ABUN_Data    <- TempJAGS_0.99_ABUN[['data']]

TropJAGS_0.99_ABUN_Data$Subset <- 'Abundance'
TempJAGS_0.99_ABUN_Data$Subset <- 'Abundance'

# Bind in mean abundances for plots 
TropJAGS_0.99_ABUN_Data_V2 <- TropJAGS_0.99_ABUN_Data %>% group_by(Species_MeanTemp_CoralWatch_V2) %>% do(ScaledLogAbundanceAdult40_V2 = mean(.$ScaledLogAbundanceAdult40, na.rm = T), 
                                                                                                          N_species = nrow(.)) %>% unnest() %>% left_join(., TropJAGS_0.99_ABUN_Data)
TempJAGS_0.99_ABUN_Data_V2 <- TempJAGS_0.99_ABUN_Data %>% group_by(Species_MeanTemp_CoralWatch_V2) %>% do(ScaledLogAbundanceAdult40_V2 = mean(.$ScaledLogAbundanceAdult40, na.rm = T), 
                                                                                                          N_species = nrow(.)) %>% unnest() %>% left_join(., TempJAGS_0.99_ABUN_Data)

# Combine guilds for aggregate plot
GuildsJAGS_ABUN <- rbind(TropJAGS_0.99_ABUN_Data, TempJAGS_0.99_ABUN_Data)
GuildsJAGS_ABUN$ThermalGuild <- factor(GuildsJAGS_ABUN$ThermalGuild, levels = c('temperate', 'tropical'))
GuildsJAGS_ABUN_V2 <- rbind(TropJAGS_0.99_ABUN_Data_V2, TempJAGS_0.99_ABUN_Data_V2)
GuildsJAGS_ABUN_V2$ThermalGuild <- factor(GuildsJAGS_ABUN_V2$ThermalGuild, levels = c('temperate', 'tropical'))

# Remove single outlier with 1 species at very high temperature. 
GuildsJAGS_ABUN_V2 <- GuildsJAGS_ABUN_V2[-which(GuildsJAGS_ABUN_V2$ThermalGuild == 'temperate' & GuildsJAGS_ABUN_V2$Species_MeanTemp_CoralWatch_V2 > 5),]
GuildsJAGS_ABUN <- GuildsJAGS_ABUN[-which(GuildsJAGS_ABUN$ThermalGuild == 'temperate' & GuildsJAGS_ABUN$Species_MeanTemp_CoralWatch_V2 > 5),]

pdf('figures_new/Abundance_TPC_ThermalGuild.pdf', width = 8, height = 4, useDingbats = FALSE, bg = 'transparent')
ggplot() + 
  
  # Minor data points
  geom_point(data = GuildsJAGS_ABUN_V2, 
             aes(x = Species_MeanTemp_CoralWatch_V2, y = ScaledLogAbundanceAdult40), pch = 19, size = 0.2, alpha = 0.1) +
  
  # Major data points
  geom_point(data = GuildsJAGS_ABUN_V2 %>% select(Species_MeanTemp_CoralWatch_V2, ScaledLogAbundanceAdult40_V2, N_species, ThermalGuild) %>% unique(), 
             aes(x = Species_MeanTemp_CoralWatch_V2, y = ScaledLogAbundanceAdult40_V2, size = N_species, fill = ThermalGuild),  pch = 21,  alpha = 0.9) +
  
  # Plot confidence intervals.
  geom_ribbon(data  = GuildsJAGS_ABUN , aes(x = Species_MeanTemp_CoralWatch_V2, ymin = Performance_Lower, ymax = Performance_Upper, fill = ThermalGuild), alpha = 0.5) +
  
  # Fitted trend lines
  geom_line(data  = GuildsJAGS_ABUN, aes(x = Species_MeanTemp_CoralWatch_V2, y = Performance, col = ThermalGuild), size = 0.5, lty = 1) +
  
  facet_wrap(~ThermalGuild, scales = 'free_x') + 
  # c('#45BBD1', '#EFF165')
  scale_colour_manual(values = c('dark blue', 'dark orange')) + 
  scale_fill_manual(values = c('dark blue', 'dark orange')) + 
  
  # Hand theme
  theme_bw() +
  theme(axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15), 
        axis.text = element_text(size = 15),
        aspect.ratio = 1, 
        text = element_text(size = 15), 
        #strip.text = element_blank(), 
        axis.line.y = element_line(), 
        axis.line.x = element_line(), 
        legend.background = element_blank(), 
        panel.grid = element_blank()) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) + 
  xlab(bquote(NULL)) +
  ylab(NULL) + 
  scale_size_continuous(breaks = c(100, 200, 300, 400, 500, 600, 700), range = c(0.1, 5)) 
dev.off()




pdf('figures_new/Abundance_TPC_ThermalGuild_temperate.pdf', width = 4, height = 4, useDingbats = FALSE, bg = 'transparent')
ggplot() + 
  # Minor data points
  geom_point(data = GuildsJAGS_ABUN_V2 %>% filter(ThermalGuild == 'temperate'), 
             aes(x = Species_MeanTemp_CoralWatch_V2, y = ScaledLogAbundanceAdult40), pch = 19, size = 0.2, alpha = 0.1) +
  # Major data points
  geom_point(data = GuildsJAGS_ABUN_V2 %>% filter(ThermalGuild == 'temperate') %>% select(Species_MeanTemp_CoralWatch_V2, ScaledLogAbundanceAdult40_V2, N_species, ThermalGuild) %>% unique(), 
             aes(x = Species_MeanTemp_CoralWatch_V2, y = ScaledLogAbundanceAdult40_V2, size = N_species), fill = '#545674',  pch = 21,  alpha = 0.9) +
  # Plot confidence intervals.
  geom_ribbon(data  = GuildsJAGS_ABUN %>% filter(ThermalGuild == 'temperate'), aes(x = Species_MeanTemp_CoralWatch_V2, ymin = Performance_Lower, ymax = Performance_Upper), fill = '#545674', alpha = 0.5) +
  # Fitted trend lines
  geom_line(data  = GuildsJAGS_ABUN %>% filter(ThermalGuild == 'temperate'), aes(x = Species_MeanTemp_CoralWatch_V2, y = Performance), col = '#545674', size = 0.5, lty = 1) +
  # Hand theme
  theme_classic() +
  theme(axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15), 
        axis.text = element_text(size = 15),
        aspect.ratio = 1, 
        text = element_text(size = 15), 
        #strip.text = element_blank(), 
        axis.line.y = element_line(), 
        axis.line.x = element_line(), 
        legend.background = element_blank()) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) + xlim(c(NA,4.9)) + 
  xlab(bquote(NULL)) +
  ylab('Species relative performance') + 
  scale_size_continuous(breaks = c(100, 200, 300, 400, 500, 600, 700), range = c(0.1, 5))
dev.off()



pdf('figures_new/Abundance_TPC_ThermalGuild_tropical.pdf', width = 4, height = 4, useDingbats = FALSE, bg = 'transparent')
ggplot() + 
  # Minor data points
  geom_point(data = GuildsJAGS_ABUN_V2 %>% filter(ThermalGuild == 'tropical'), 
             aes(x = Species_MeanTemp_CoralWatch_V2, y = ScaledLogAbundanceAdult40), pch = 19, size = 0.2, alpha = 0.1) +
  # Major data points
  geom_point(data = GuildsJAGS_ABUN_V2 %>% filter(ThermalGuild == 'tropical') %>% select(Species_MeanTemp_CoralWatch_V2, ScaledLogAbundanceAdult40_V2, N_species, ThermalGuild) %>% unique(), 
             aes(x = Species_MeanTemp_CoralWatch_V2, y = ScaledLogAbundanceAdult40_V2, size = N_species), fill = '#D17F3E',  pch = 21,  alpha = 0.9) +
  # Plot confidence intervals.
  geom_ribbon(data  = GuildsJAGS_ABUN %>% filter(ThermalGuild == 'tropical'), aes(x = Species_MeanTemp_CoralWatch_V2, ymin = Performance_Lower, ymax = Performance_Upper), fill = '#D17F3E', alpha = 0.5) +
  # Fitted trend lines
  geom_line(data  = GuildsJAGS_ABUN %>% filter(ThermalGuild == 'tropical'), aes(x = Species_MeanTemp_CoralWatch_V2, y = Performance), col = '#D17F3E', size = 0.5, lty = 1) +
  # Hand theme
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.line.y.left = element_blank(),
        axis.title.x = element_text(size = 10), 
        axis.text.x = element_text(size = 10),
        aspect.ratio = 1, 
        text = element_text(size = 10), 
        #strip.text = element_blank(), 
        axis.line.y = element_line(), 
        axis.line.x = element_line(), 
        legend.background = element_blank()) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) + 
  xlab(bquote(NULL)) +
  ylab('Species relative performance') + 
  scale_size_continuous(breaks = c(100, 200, 300, 400, 500, 600, 700), range = c(0.1, 5)) 
dev.off()

# ----------------------------------------------------------------------------


# END OF SCRIPT ------