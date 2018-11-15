# Script combines output parameters of the species distribution models and the quantile gams. 
# Produces object used in further modelling 

# 1. Reads in all relevant data sets.
# 2. Combines data.
# 3. Produces quality controls.

# Initiated 07/11/2018
# Author: Conor Waldock

# Load libraries ----
library(lme4)
library(MuMIn)
library(remef)
library(stargazer)
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(viridis)
library(mgcv)



# READ IN ALL RELEVANT DATA --------------------------------------------------------------------------------------------------------------
# OBTAIN DATA FROM MODELS ---- 

# Thermal optima from qgams: 
#load(file = 'data_derived2/AllQgams_2018-09-28.RData')
Quantile_Parameters <- readRDS(file = 'data_derived/qGamModelOutputs-2019-10-03.rds')
head(Quantile_Parameters)
Quantile_Parameters$Topt_SD[is.na(Quantile_Parameters$Topt_SD)] <- 0

# Thermal limits from SDMs: 
SDM_NicheParameters <- readRDS(file = 'data_derived/sdmModelOutputs_2018-10-08.rds')
SDM_NicheParameters[,-1] <- apply(SDM_NicheParameters[,-1], 2, as.numeric)

# DEFINE OBSEVED AND SURVEY THERMAL LIMTIS ----

# Read in RLS data 
RLS_All  <- readRDS(file = 'data_derived/RLS_20-For-Qgams-2019-09-28.rds')

# Create object of thermal niche from observations only (for later merge, some redundancy here)
Observed_Limits <- RLS_All %>% 
  filter(AbundanceAdult40 > 0) %>% 
  group_by(SpeciesName) %>% 
  nest() %>% 
  mutate(T_Upper_Obs     = purrr::map(data, ~quantile(.$MeanTemp_CoralWatch, 0.975)),
         T_Lower_Obs     = purrr::map(data, ~quantile(.$MeanTemp_CoralWatch, 0.025)),
         T_Midpoint_Obs  = purrr::map(data, ~    mean(.$MeanTemp_CoralWatch, na.rm = T)), 
         T_Upper_Obs_MAX = purrr::map(data, ~quantile(.$MaxTemp_CoralWatch, 0.975)),
         T_Lower_Obs_MIN = purrr::map(data, ~quantile(.$MinTemp_CoralWatch, 0.025))) %>% 
  unnest(T_Upper_Obs, T_Lower_Obs, T_Midpoint_Obs, T_Upper_Obs_MAX, T_Lower_Obs_MIN) %>% 
  dplyr::select(-data)

N_Samples <- RLS_All %>% 
  group_by(SpeciesName) %>% 
  nest() %>% 
  mutate(N_Presence = purrr::map(data, ~sum(.$Presence == 1)), 
         N_Absence = purrr::map(data, ~sum(.$Presence == 0)), 
         N_Samples = purrr::map(data, ~nrow(.)))%>% 
  unnest(N_Presence, N_Absence, N_Samples) %>% 
  dplyr::select(-data)

N_Samples$Prevalence <- with(N_Samples, N_Presence / N_Samples)

NicheLimitSampling <- RLS_All %>% dplyr::select(SpeciesName, N_Absences_T_Lower, N_Absences_T_Upper) %>% unique()

# Create sampling uppers and sampling lowers
AbsenceLimits <- RLS_All %>% 
  dplyr::select(SpeciesName, MeanTemp_CoralWatch, Presence) %>% group_by(SpeciesName) %>% 
  do(Absence_TUpper = max(.$MeanTemp_CoralWatch), 
     Absence_TLower = min(.$MeanTemp_CoralWatch)) %>% 
  unnest(Absence_TUpper, Absence_TLower)


# COMBINE ThermalNicheData_Review ----
ThermalNicheData_Review <- left_join(left_join(left_join(left_join(Quantile_Parameters, Observed_Limits), SDM_NicheParameters), N_Samples), NicheLimitSampling)
apply(ThermalNicheData_Review, 2, function(x) sum(is.na(x))) # There are 13 species' for which SDMs did not fit. 

# Combine with thermal niche data new
ThermalNicheData_Review <- left_join(ThermalNicheData_Review, AbsenceLimits)


# Create thermal guilds ---- 
ThermalNicheData_Review$ThermalGuild <- ifelse(ThermalNicheData_Review$Topt >= 23, 'tropical', 'temperate')
# ----------------------------------------------------------------------------------------------------------------------------------------









# PRODUDE QUALTIY CONTROL V2 ----
# Quality control objects ----
Observed_limits_QC <- data.frame(SpeciesName                 = ThermalNicheData_Review$SpeciesName, 
                                 Conf_absence_limits          = NA, 
                                 Conf_qgam_deviance          = NA, 
                                 Conf_qgam_bootstrap         = NA, 
                                 Conf_qgam_samplinglimit     = NA, 
                                 ConfidenceCombined_observed = NA)

SDM_limits_QC <- data.frame(SpeciesName               = ThermalNicheData_Review$SpeciesName, 
                            ConfSDM_T_limits          = NA, 
                            ConfSDM_sdm_measure       = NA, 
                            ConfSDM_thermal_rangesize = NA, 
                            ConfidenceCombined_SDM    = NA)
#
# Observed thermal limit quality controls ----

# 1. Are there more than 10 samples above and below the species observed thermal niche limit? 
Observed_limits_QC$Conf_absence_limits <- ifelse(ThermalNicheData_Review$N_Absences_T_Lower < 10 | ThermalNicheData_Review$N_Absences_T_Upper < 10, -1, 0)

# 2. Deviance explained by quantile gam is > 75th quantile of deviance explained across all species
Observed_limits_QC$Conf_qgam_deviance  <- ifelse(ThermalNicheData_Review$T_Gam_deviance.exp < quantile(ThermalNicheData_Review$T_Gam_deviance.exp, na.rm = T, 0.75), 0, -1)

# 3. Standard deviation in estimates of Topt, due to bootstrapping, is not > 0.5
Observed_limits_QC$Conf_qgam_bootstrap <- ifelse(ThermalNicheData_Review$Topt_SD <= 0.5, 0, -1)

# 4. Estimated Topt is not at the limits of sampled temperatures within a species' range.
Observed_limits_QC$Conf_qgam_samplinglimit <- ifelse(ThermalNicheData_Review$Absence_TUpper - ThermalNicheData_Review$Topt < 0 | 
                                                     ThermalNicheData_Review$Topt - ThermalNicheData_Review$Absence_TLower < 0, 
                                                     -1, 0)


# Sum together confidence estimates
Observed_limits_QC$ConfidenceCombined_Observed <- 3 + apply(Observed_limits_QC[,2:5], 1, sum)

table(Observed_limits_QC$ConfidenceCombined_Observed)
# 0   1   2   3 
# 13 178 330 181 

# SDM quality controls ----

# 1. Tmax or Tmin are 3°C above or below sampling limits
SDM_limits_QC$ConfSDM_T_limits <- ifelse(ThermalNicheData_Review$T_Upper_0.1 - ThermalNicheData_Review$Absence_TUpper >  3 | 
                                         ThermalNicheData_Review$T_Lower_0.1 - ThermalNicheData_Review$Absence_TLower < -3, 
                                         -1, 0)

# 2. Specificity or sensitivity from SDMs are < 0.7
SDM_limits_QC$ConfSDM_sdm_measure <- ifelse(ThermalNicheData_Review$specificity < 0.7 | ThermalNicheData_Review$sensitivity < 0.7 | is.na(ThermalNicheData_Review$sensitivity), -1, 0)

# 3. Range of temperatures between Tmin and Tmax is not 0.5 – 2x the sampling range
T_Range_RATIO <- (ThermalNicheData_Review$T_Upper_0.1 - ThermalNicheData_Review$T_Lower_0.1) / (ThermalNicheData_Review$T_Upper_Obs - ThermalNicheData_Review$T_Lower_Obs)
SDM_limits_QC$ConfSDM_thermal_rangesize <- ifelse(T_Range_RATIO <= 2, ifelse(T_Range_RATIO <= 2 & T_Range_RATIO > 0.5, 0, -1), -1)

# Sum together confidence estimates
SDM_limits_QC$ConfidenceCombined_SDM <- apply(SDM_limits_QC[,2:4], 1, sum)

table(SDM_limits_QC$ConfidenceCombined_SDM)
# -2  -1   0 
# 16 321 365 


# Combine together confidence scoring criteria onto ThermalNicheData_Review ----
ThermalNicheData_Review$ConfidenceCombined     <- Observed_limits_QC$ConfidenceCombined_Observed
ThermalNicheData_Review$ConfidenceCombined_SDM <- Observed_limits_QC$ConfidenceCombined_Observed + SDM_limits_QC$ConfidenceCombined_SDM

table(ThermalNicheData_Review$ConfidenceCombined)
# 0   1   2   3 
# 13 178 330 181 

table(ThermalNicheData_Review$ConfidenceCombined_SDM)
# -1   0   1   2   3 
# 12 103 233 255  99 


# Figure for SDM summary statistics for supporting materials ----
pdf('figures_new/SDM_statistic_comparison.pdf', width = 7, height = 5)
gridExtra::grid.arrange(
  ggplot() + geom_histogram(data = ThermalNicheData_Review[SDM_limits_QC$ConfidenceCombined_SDM != 0,], aes(x = specificity)) +
    geom_histogram(data = ThermalNicheData_Review[SDM_limits_QC$ConfidenceCombined_SDM == 0,], aes(x = specificity), fill = 'red3')+ theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank()),
  ggplot() + geom_histogram(data = ThermalNicheData_Review[SDM_limits_QC$ConfidenceCombined_SDM != 0,], aes(x = sensitivity)) + 
    geom_histogram(data = ThermalNicheData_Review[SDM_limits_QC$ConfidenceCombined_SDM == 0,], aes(x = sensitivity), fill = 'red3')+ theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank()),
  ggplot() + geom_histogram(data = ThermalNicheData_Review[SDM_limits_QC$ConfidenceCombined_SDM != 0,], aes(x = AUC)) + 
    geom_histogram(data = ThermalNicheData_Review[SDM_limits_QC$ConfidenceCombined_SDM == 0,], aes(x = AUC), fill = 'red3')+ theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank()),
  ggplot() + geom_histogram(data = ThermalNicheData_Review[SDM_limits_QC$ConfidenceCombined_SDM != 0,], aes(x = TSS)) + 
    geom_histogram(data = ThermalNicheData_Review[SDM_limits_QC$ConfidenceCombined_SDM == 0,], aes(x = TSS), fill = 'red3')+ theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank()),
  ncol = 2, nrow = 2) + theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank())
dev.off()





# Summaries ----

mean(ThermalNicheData_Review$AUC[which(ThermalNicheData_Review$ConfidenceCombined_SDM == 3)], na.rm = T) # 0.82
mean(ThermalNicheData_Review$TSS[which(ThermalNicheData_Review$ConfidenceCombined_SDM == 3)], na.rm = T) # 0.60
summary(ThermalNicheData_Review$AUC); sd(ThermalNicheData_Review$AUC)
summary(ThermalNicheData_Review$TSS); sd(ThermalNicheData_Review$TSS)

# Estimate amount of skew for 1. observed limits, 2. observed seasonal limits, 3. SDM based limits ---- 
# Re-estimate SDs based on reestimation of Topt and criteria for thermal limits. 
ThermalNicheData_Review$T_SD_Upper <- abs(ThermalNicheData_Review$T_Upper_0.1 - ThermalNicheData_Review$Topt) / 1.96
ThermalNicheData_Review$T_SD_Lower <- abs(ThermalNicheData_Review$T_Lower_0.1 - ThermalNicheData_Review$Topt) / 1.96

# Re-estimate SDs based on reestimation of Topt and criteria for thermal limits. 
ThermalNicheData_Review$T_SD_Upper_OBS <- abs(ThermalNicheData_Review$T_Upper_Obs - ThermalNicheData_Review$Topt) / 1.96
ThermalNicheData_Review$T_SD_Lower_OBS <- abs(ThermalNicheData_Review$T_Lower_Obs - ThermalNicheData_Review$Topt) / 1.96

# Re-estimate SDs based on reestimation of Topt and criteria for thermal limits including seasonality
ThermalNicheData_Review$T_SD_Upper_OBS_season <- abs(ThermalNicheData_Review$T_Upper_Obs_MAX - ThermalNicheData_Review$Topt) / 1.96
ThermalNicheData_Review$T_SD_Lower_OBS_season <- abs(ThermalNicheData_Review$T_Lower_Obs_MIN - ThermalNicheData_Review$Topt) / 1.96


# Create column estimating skew
ThermalNicheData_Review$T_Skew            <- ThermalNicheData_Review$T_SD_Upper - ThermalNicheData_Review$T_SD_Lower
ThermalNicheData_Review$T_Skew_OBS        <- ThermalNicheData_Review$T_SD_Upper_OBS - ThermalNicheData_Review$T_SD_Lower_OBS
ThermalNicheData_Review$T_Skew_OBS_season <- ThermalNicheData_Review$T_SD_Upper_OBS_season - ThermalNicheData_Review$T_SD_Lower_OBS_season

ThermalNicheData_Review$T_SD_Upper <- NULL
ThermalNicheData_Review$T_SD_Lower <- NULL
ThermalNicheData_Review$T_SD_Upper_OBS <- NULL
ThermalNicheData_Review$T_SD_Lower_OBS <- NULL
ThermalNicheData_Review$T_SD_Upper_OBS_season <- NULL
ThermalNicheData_Review$T_SD_Lower_OBS_season <- NULL

# Save thermal niche data object for modelling ----
saveRDS(ThermalNicheData_Review, file = 'data_derived/ThermalNicheData_Review.rds')
saveRDS(ThermalNicheData_Review, file = 'data_upload/ThermalNicheData.rds') # Save here for ease of running the next script if not wanting to run this one. 

# ----


##### END OF SCRIPT ######




# ----
# OLD QUALITY CONTROLS ---------------------------------------------------------------------------------------------------------------
# These match the workflow in the supporting materials 
# 1. Confidence in estimation of upper and lower limit ----

# Create confidence limits (first are they modelled, second are they within 3°C of absence records limits
# So this says, if our estimate of niche limit is above or below our sampling limits by too much, how can we be sure that we have sampled the limits
# of species distributions well enough? (this isn't to do with the presence limits of the species, but the limits of the sampled range). 
ThermalNicheData_Review$Conf_T_Upper <- ifelse(ThermalNicheData_Review$T_Upper_0.1 - ThermalNicheData_Review$Absence_TUpper > 3, -1, 0)
ThermalNicheData_Review$Conf_T_Lower <- ifelse(ThermalNicheData_Review$T_Lower_0.1 - ThermalNicheData_Review$Absence_TLower < -3, -1, 0)

# New confidence scores are based on metrics of sensitivity and specificity
mean(ThermalNicheData_Review$specificity, na.rm = T) # A false absence (i.e., missing presences) is not as critical for our niche parameter estimates. 
mean(ThermalNicheData_Review$sensitivity, na.rm = T) # A false presence is worse for our parameters than a false absence. (leads to overestimating range size)

# 219 species have high specificity and sensitivity
ThermalNicheData_Review$Conf_SDM <- ifelse(ThermalNicheData_Review$specificity < 0.7 | ThermalNicheData_Review$sensitivity < 0.7 | is.na(ThermalNicheData_Review$sensitivity), -2, 0)
sum(ThermalNicheData_Review$Conf_SDM==0)
mean(ThermalNicheData_Review$AUC[which(ThermalNicheData_Review$Conf_SDM == 0)], na.rm = T) # 
mean(ThermalNicheData_Review$TSS[which(ThermalNicheData_Review$Conf_SDM == 0)], na.rm = T) # 

# Plot of histograms of sensitivity and specificity and consequences of threshold for AUC and TSS selection. 
pdf('figures_new/SDM_statistic_comparison.pdf', width = 7, height = 5)
gridExtra::grid.arrange(
  ggplot() + geom_histogram(data = ThermalNicheData_Review[ThermalNicheData_Review$Conf_SDM != 0,], aes(x = specificity)) +
    geom_histogram(data = ThermalNicheData_Review[ThermalNicheData_Review$Conf_SDM == 0,], aes(x = specificity), fill = 'red3')+ theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank()),
  ggplot() + geom_histogram(data = ThermalNicheData_Review[ThermalNicheData_Review$Conf_SDM != 0,], aes(x = sensitivity)) + 
    geom_histogram(data = ThermalNicheData_Review[ThermalNicheData_Review$Conf_SDM == 0,], aes(x = sensitivity), fill = 'red3')+ theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank()),
  ggplot() + geom_histogram(data = ThermalNicheData_Review[ThermalNicheData_Review$Conf_SDM != 0,], aes(x = AUC)) + 
    geom_histogram(data = ThermalNicheData_Review[ThermalNicheData_Review$Conf_SDM == 0,], aes(x = AUC), fill = 'red3')+ theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank()),
  ggplot() + geom_histogram(data = ThermalNicheData_Review[ThermalNicheData_Review$Conf_SDM != 0,], aes(x = TSS)) + 
    geom_histogram(data = ThermalNicheData_Review[ThermalNicheData_Review$Conf_SDM == 0,], aes(x = TSS), fill = 'red3')+ theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank()),
  ncol = 2, nrow = 2) + theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank())
dev.off()

# Test for correlations amost SDM parameters and number of presences, absences and prevalence (relative presence)
par(mfrow = c(2,3))
plot(ThermalNicheData_Review$specificity ~ log10(ThermalNicheData_Review$N_Presence))
plot(ThermalNicheData_Review$specificity ~ log10(ThermalNicheData_Review$N_Absence))
plot(ThermalNicheData_Review$specificity ~ log10(ThermalNicheData_Review$Prevalence))
plot(ThermalNicheData_Review$sensitivity ~ log10(ThermalNicheData_Review$N_Presence))
plot(ThermalNicheData_Review$sensitivity ~ log10(ThermalNicheData_Review$N_Absence))
plot(ThermalNicheData_Review$sensitivity ~ log10(ThermalNicheData_Review$Prevalence))
dev.off()

mean(ThermalNicheData_Review$AUC); sd(ThermalNicheData_Review$AUC); range(ThermalNicheData_Review$AUC)
mean(ThermalNicheData_Review$TSS); sd(ThermalNicheData_Review$TSS); range(ThermalNicheData_Review$TSS)

# Plot correlations between 
CorThermalEdges <- ThermalNicheData_Review[,c('T_Lower_0.1', 'T_Lower_0.25', 'T_Lower_0.5', 'T_Upper_0.1', 'T_Upper_0.25', 'T_Upper_0.5')]
names(CorThermalEdges) <- c('Tmin - 0.1', 'Tmin - 0.25', 'Tmin - 0.5', 'Tmax - 0.1', 'Tmax - 0.25', 'Tmax - 0.5')
pdf('figures_new/ThermalEdges-Correlations.pdf', width = 5, height = 5)
psych::pairs.panels(CorThermalEdges)
dev.off()

# 2. Confidence in sampling of upper and lower limit ---- 
ThermalNicheData_Review$Conf_Absence_Sampling <- ifelse(ThermalNicheData_Review$N_Absences_T_Lower < 10 | ThermalNicheData_Review$N_Absences_T_Upper < 10, -1, 0)

# 2. Confidence in T_Breadth ----
ThermalNicheData_Review$T_Range <- ThermalNicheData_Review$T_Upper_0.1 - ThermalNicheData_Review$T_Lower_0.1

# Estimate ratio between observed range and estimated range. 
ThermalNicheData_Review$T_Range_RATIO <- ThermalNicheData_Review$T_Range / (ThermalNicheData_Review$T_Upper_Obs - ThermalNicheData_Review$T_Lower_Obs)

# We see that the observed T_Range is now rarely much bigger than the modelled T_Range.
hist(ThermalNicheData_Review$T_Range_RATIO, breaks = 50)

# Assign confidence scores. 
ThermalNicheData_Review$Conf_T_Breadth <- ifelse(ThermalNicheData_Review$T_Range_RATIO <= 2, ifelse(ThermalNicheData_Review$T_Range_RATIO <= 2 & ThermalNicheData_Review$T_Range_RATIO > 0.5, 0, -1), -1)
# ThermalNicheData_Review %>% select(ThermalGuild, Conf_T_Breadth, SpeciesName) %>% unique(.) %>% filter(Conf_T_Breadth != 0) %>% .$ThermalGuild %>% table(.) # 

# Use realized (sampling) thermal limits to define Tmin or Tmax with previous addition of difference between confidence 3 species sampling limits and modelled thermal limits. 
#ThermalNicheData_Review$T_Upper[which((ThermalNicheData_Review$Conf_T_Upper == -1 | ThermalNicheData_Review$Conf_T_Breadth == -1) & ThermalNicheData_Review$ThermalGuild == 'tropical' )] <- 
#  ThermalNicheData_Review$T_Upper_Obs[which((ThermalNicheData_Review$Conf_T_Upper == -1 | ThermalNicheData_Review$Conf_T_Breadth == -1) & ThermalNicheData_Review$ThermalGuild == 'tropical' )] + median(Tropical_Difference_Upper)
#ThermalNicheData_Review$T_Upper[which((ThermalNicheData_Review$Conf_T_Upper == -1 | ThermalNicheData_Review$Conf_T_Breadth == -1) & ThermalNicheData_Review$ThermalGuild == 'temperate' )] <- 
#  ThermalNicheData_Review$T_Upper_Obs[which((ThermalNicheData_Review$Conf_T_Upper == -1 | ThermalNicheData_Review$Conf_T_Breadth == -1) & ThermalNicheData_Review$ThermalGuild == 'temperate' )] + median(Temperate_Difference_Upper)

#ThermalNicheData_Review$T_Lower[which((ThermalNicheData_Review$Conf_T_Lower == -1 | ThermalNicheData_Review$Conf_T_Breadth == -1) & ThermalNicheData_Review$ThermalGuild == 'tropical' )] <- 
#  ThermalNicheData_Review$T_Lower_Obs[which((ThermalNicheData_Review$Conf_T_Lower == -1 | ThermalNicheData_Review$Conf_T_Breadth == -1) & ThermalNicheData_Review$ThermalGuild == 'tropical' )] + median(Tropical_Difference_Lower)
#ThermalNicheData_Review$T_Lower[which((ThermalNicheData_Review$Conf_T_Lower == -1 | ThermalNicheData_Review$Conf_T_Breadth == -1) & ThermalNicheData_Review$ThermalGuild == 'temperate' )] <- 
#  ThermalNicheData_Review$T_Lower_Obs[which((ThermalNicheData_Review$Conf_T_Lower == -1 | ThermalNicheData_Review$Conf_T_Breadth == -1) & ThermalNicheData_Review$ThermalGuild == 'temperate' )] + median(Temperate_Difference_Lower)

# 3. Confidence in qgam estimate of Topt ---- 

# Confidence_Qgam
table(ThermalNicheData_Review$T_Gam_pvalue < 0.05)
# FALSE  TRUE # Significant vs. non-significant. 
# 248    442  # Most are significant. 

hist(ThermalNicheData_Review$T_Gam_deviance.exp, 50)
summary(ThermalNicheData_Review$T_Gam_deviance.exp)

ThermalNicheData_Review$Conf_Qgam <- ifelse(ThermalNicheData_Review$T_Gam_deviance.exp < quantile(ThermalNicheData_Review$T_Gam_deviance.exp, na.rm = T, 0.75), 0, -1)
table(ThermalNicheData_Review$Conf_Qgam)
ThermalNicheData_Review$Conf_Qgam2 <- ifelse(ThermalNicheData_Review$Topt_SD <= 0.5, 0, -1)
table(ThermalNicheData_Review$Conf_Qgam2)


#ThermalNicheData_Review %>% select(ThermalGuild, Conf_Qgam, SpeciesName) %>% unique(.) %>% filter(Conf_Qgam != 0) %>% .$ThermalGuild %>% table(.) #
# temperate  tropical 
# 26        66  

# Redefine Topt to midpoint when non-significant or no confidence in upper/lower. 
# ThermalNicheData_Review$Topt[which(ThermalNicheData_Review$Conf_Qgam == -1)] <- ThermalNicheData_Review$T_Midpoint_Obs[which(ThermalNicheData_Review$Conf_Qgam == -1)]

# Re-estimate SDs based on reestimation of Topt and criteria for thermal limits. 
ThermalNicheData_Review$T_SD_Upper <- abs(ThermalNicheData_Review$T_Upper_0.1 - ThermalNicheData_Review$Topt) / 1.96
ThermalNicheData_Review$T_SD_Lower <- abs(ThermalNicheData_Review$T_Lower_0.1 - ThermalNicheData_Review$Topt) / 1.96

# Re-estimate SDs based on reestimation of Topt and criteria for thermal limits. 
ThermalNicheData_Review$T_SD_Upper_OBS <- abs(ThermalNicheData_Review$T_Upper_Obs - ThermalNicheData_Review$Topt) / 1.96
ThermalNicheData_Review$T_SD_Lower_OBS <- abs(ThermalNicheData_Review$T_Lower_Obs - ThermalNicheData_Review$Topt) / 1.96

# Re-estimate SDs based on reestimation of Topt and criteria for thermal limits including seasonality
ThermalNicheData_Review$T_SD_Upper_OBS_season <- abs(ThermalNicheData_Review$T_Upper_Obs_MAX - ThermalNicheData_Review$Topt) / 1.96
ThermalNicheData_Review$T_SD_Lower_OBS_season <- abs(ThermalNicheData_Review$T_Lower_Obs_MIN - ThermalNicheData_Review$Topt) / 1.96


# Create column estimating skew
ThermalNicheData_Review$T_Skew <- ThermalNicheData_Review$T_SD_Upper - ThermalNicheData_Review$T_SD_Lower
ThermalNicheData_Review$T_Skew_OBS <- ThermalNicheData_Review$T_SD_Upper_OBS - ThermalNicheData_Review$T_SD_Lower_OBS
ThermalNicheData_Review$T_Skew_OBS_season <- ThermalNicheData_Review$T_SD_Upper_OBS_season - ThermalNicheData_Review$T_SD_Lower_OBS_season

# 4. Confidence in Topt (should not buffer exactly the minimum temperature sampled) ----
ThermalNicheData_Review$T_Opt_Difference_Upper_Sampling<- ThermalNicheData_Review$Absence_TUpper - ThermalNicheData_Review$Topt
ThermalNicheData_Review$T_Opt_Difference_Lower_Sampling <- ThermalNicheData_Review$Topt    - ThermalNicheData_Review$Absence_TLower

# Confidence_T_Opt_Difference_Upper + Confidence_T_Opt_Difference_Lower
ThermalNicheData_Review$Confidence_T_Opt_Difference_Upper <- ifelse(ThermalNicheData_Review$T_Opt_Difference_Upper_Sampling < 0, -1, 0)
ThermalNicheData_Review$Confidence_T_Opt_Difference_Lower <- ifelse(ThermalNicheData_Review$T_Opt_Difference_Lower_Sampling < 0, -1, 0)

# Topt shouldn't exceed SDM estimate of niche limit. Suggests the projection of Topt combined with the SDM is an inaccurate representation of the niche shape.  
ThermalNicheData_Review$T_Opt_Difference_Upper1<- ThermalNicheData_Review$T_Upper_0.1 - ThermalNicheData_Review$Topt
ThermalNicheData_Review$T_Opt_Difference_Lower1 <- ThermalNicheData_Review$Topt - ThermalNicheData_Review$T_Lower_0.1
ThermalNicheData_Review$Confidence_T_Opt_Difference_Upper1 <- ifelse(ThermalNicheData_Review$T_Opt_Difference_Upper1 < 0, -1, 0)
ThermalNicheData_Review$Confidence_T_Opt_Difference_Lower1 <- ifelse(ThermalNicheData_Review$T_Opt_Difference_Lower1 < 0, -1, 0)

# Combine confidence scores ----

ThermalNicheData_Review$ConfidenceCombined <- 3 + 
  ThermalNicheData_Review$Conf_SDM + 
  ThermalNicheData_Review$Conf_T_Upper + 
  ThermalNicheData_Review$Conf_T_Lower +
  ThermalNicheData_Review$Conf_T_Breadth + 
  ThermalNicheData_Review$Conf_Qgam + ThermalNicheData_Review$Conf_Qgam2 + 
  ThermalNicheData_Review$Confidence_T_Opt_Difference_Upper + ThermalNicheData_Review$Confidence_T_Opt_Difference_Upper1 + 
  ThermalNicheData_Review$Confidence_T_Opt_Difference_Lower + ThermalNicheData_Review$Confidence_T_Opt_Difference_Lower1 + 
  ThermalNicheData_Review$Conf_Absence_Sampling

table(ThermalNicheData_Review$ConfidenceCombined)

# SUMMARY OF CONFIDENCE CRITERIA ----
table(ThermalNicheData_Review$ConfidenceCombined)
# -3  -2  -1   0   1   2   3 
#  3  16 115 129 198 145  96 

# ----------------------------------------------------------------------------------------------------------------------------------------




# ----

# Relationship between Confidence and parameters of models ----
ggplot(ThermalNicheData_Review) + 
  geom_boxplot(aes(x = ConfidenceCombined, y = AUC, group = ConfidenceCombined))

ggplot(ThermalNicheData_Review) + 
  geom_boxplot(aes(x = ConfidenceCombined, y = T_Gam_deviance.exp, group = ConfidenceCombined))

ggplot(ThermalNicheData_Review) + 
  geom_boxplot(aes(x = ConfidenceCombined, y = sqrt(sqrt(T_Gam_pvalue)), group = ConfidenceCombined))

ggplot(ThermalNicheData_Review) + 
  geom_boxplot(aes(x = ConfidenceCombined, y = log10(N_Presence), group = ConfidenceCombined))

gridExtra::grid.arrange(
  ggplot(ThermalNicheData_Review) + 
    geom_point(aes(y = T_Upper_Obs_MAX, x = Topt), col = 'red') + 
    geom_point(aes(y = T_Lower_Obs_MIN, x = Topt), col = 'blue') + 
    geom_abline(),
  
  ggplot(ThermalNicheData_Review) + 
    geom_point(aes(y = T_Upper_Obs, x = Topt), col = 'red') + 
    geom_point(aes(y = T_Lower_Obs, x = Topt), col = 'blue') + 
    geom_abline(),
  
  ggplot(ThermalNicheData_Review) + 
    geom_point(aes(y = T_Upper_0.1, x = Topt), col = 'orange3') + 
    geom_point(aes(y = T_Lower_0.1, x = Topt), col = 'darkblue') + 
    geom_abline()
)

# ---- 



# SAVE THERMAL NICHE DATA OBJECT ----
saveRDS(ThermalNicheData_Review, file = 'data_derived/ThermalNicheData_Review.rds')
saveRDS(ThermalNicheData_Review, file = 'data_upload/ThermalNicheData.rds') # Save here for ease of running the next script if not wanting to run this one. 
# ----



