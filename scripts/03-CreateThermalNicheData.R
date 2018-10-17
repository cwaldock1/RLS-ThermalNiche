# Script combines output parameters of the species distribution models and the quantile gams. 
# Produces object used in further modelling 

# 1. Reads in all relevant data sets.
# 2. Combines data.
# 3. Produces quality controls.

# Initiated 08/10/2018
# Author: Conor Waldock

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
  mutate(T_Upper_Obs = purrr::map(data, ~max(.$MeanTemp_CoralWatch)),
         T_Lower_Obs = purrr::map(data, ~min(.$MeanTemp_CoralWatch)),
         T_Midpoint_Obs = purrr::map(data, ~mean(.$MeanTemp_CoralWatch, na.rm = T))) %>% 
  unnest(T_Upper_Obs, T_Lower_Obs, T_Midpoint_Obs) %>% 
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


# COMBINE ThermalNicheData_Review ----
ThermalNicheData_Review <- left_join(left_join(left_join(Quantile_Parameters, Observed_Limits), SDM_NicheParameters), N_Samples)
apply(ThermalNicheData_Review, 2, function(x) sum(is.na(x))) # There are 13 species' for which SDMs did not fit. 

# Create therma guilds ---- 
ThermalNicheData_Review$ThermalGuild <- ifelse(ThermalNicheData_Review$Topt >= 23, 'tropical', 'temperate')
# ----------------------------------------------------------------------------------------------------------------------------------------


# PRODUCE QUALITY CONTROLS ---------------------------------------------------------------------------------------------------------------
# These match the workflow in the supporting materials 
# 1. Confidence in upper and lower limit ----

# Create sampling uppers and sampling lowers
AbsenceLimits <- RLS_All %>% dplyr::select(SpeciesName, MeanTemp_CoralWatch, Presence) %>% group_by(SpeciesName) %>% 
  do(Absence_TUpper = max(.$MeanTemp_CoralWatch), 
     Absence_TLower = min(.$MeanTemp_CoralWatch)) %>% 
  unnest(Absence_TUpper, Absence_TLower)

# Combine with thermal niche data new
ThermalNicheData_Review <- left_join(ThermalNicheData_Review, AbsenceLimits)

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
mean(ThermalNicheData_Review$AUC[which(ThermalNicheData_Review$Conf_SDM == 0)], na.rm = T) # 0.89
mean(ThermalNicheData_Review$TSS[which(ThermalNicheData_Review$Conf_SDM == 0)], na.rm = T) # 0.74

# Plot of histograms of sensitivity and specificity and consequences of threshold for AUC and TSS selection. 
pdf('figures_new/SDM_statistic_comparison.pdf', width = 7, height = 5)
gridExtra::grid.arrange(
  ggplot() + geom_histogram(data = ThermalNicheData_Review, aes(x = specificity)) + geom_histogram(data = ThermalNicheData_Review[ThermalNicheData_Review$Conf_SDM == 0,], aes(x = specificity), col = 'red3')+ theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank()),
  ggplot() + geom_histogram(data = ThermalNicheData_Review, aes(x = sensitivity)) + geom_histogram(data = ThermalNicheData_Review[ThermalNicheData_Review$Conf_SDM == 0,], aes(x = sensitivity), col = 'red3')+ theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank()),
  ggplot() + geom_histogram(data = ThermalNicheData_Review, aes(x = AUC)) + geom_histogram(data = ThermalNicheData_Review[ThermalNicheData_Review$Conf_SDM == 0,], aes(x = AUC), col = 'red3')+ theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank()),
  ggplot() + geom_histogram(data = ThermalNicheData_Review, aes(x = TSS)) + geom_histogram(data = ThermalNicheData_Review[ThermalNicheData_Review$Conf_SDM == 0,], aes(x = TSS), col = 'red3')+ theme_bw() + theme(aspect.ratio = 0.75, panel.grid = element_blank()),
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

# Create column estimating skew
ThermalNicheData_Review$T_Skew <- ThermalNicheData_Review$T_SD_Upper - ThermalNicheData_Review$T_SD_Lower

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
  ThermalNicheData_Review$Confidence_T_Opt_Difference_Lower + ThermalNicheData_Review$Confidence_T_Opt_Difference_Lower1

table(ThermalNicheData_Review$ConfidenceCombined)

# SUMMARY OF CONFIDENCE CRITERIA ----
table(ThermalNicheData_Review$ConfidenceCombined)
# -3  -2  -1   0   1   2   3 
# 1   8  59 143 159 205 127 



# ----------------------------------------------------------------------------------------------------------------------------------------




# ----
# Relationship between Confidence and parameters of models ----
ggplot(ThermalNicheData) + 
  geom_boxplot(aes(x = ConfidenceCombined, y = AUC, group = ConfidenceCombined))

ggplot(ThermalNicheData) + 
  geom_boxplot(aes(x = ConfidenceCombined, y = T_Gam_deviance.exp, group = ConfidenceCombined))

ggplot(ThermalNicheData) + 
  geom_boxplot(aes(x = ConfidenceCombined, y = sqrt(sqrt(T_Gam_pvalue)), group = ConfidenceCombined))

ggplot(ThermalNicheData) + 
  geom_boxplot(aes(x = ConfidenceCombined, y = log10(N_Presence), group = ConfidenceCombined))

# ---- 
# SAVE THERMAL NICHE DATA OBJECT ----
saveRDS(ThermalNicheData_Review, file = 'data_derived/ThermalNicheData_Review.rds')
saveRDS(ThermalNicheData_Review, file = 'data_upload/ThermalNicheData.rds') # Save here for ease of running the next script if not wanting to run this one. 
# ----











### PRELIMINARY COMPARISON OF SKEWS -----
ggplot(ThermalNicheData_Review %>% filter(ConfidenceCombined == 3)) +
  geom_linerange(aes(ymax = T_Upper_0.1, ymin = T_Lower_0.1, x = Topt), col = 'black') + 
  geom_point(aes(y = T_Lower_0.1, x = Topt, pch = ThermalGuild), col = 'blue', size = 5) + 
  geom_point(aes(y = T_Upper_0.1, x = Topt, pch = ThermalGuild), col = 'red', size = 5) + 
  stat_smooth(aes(y = T_Lower_0.1, x = Topt, group = ThermalGuild), col = 'black', size = 1, method = 'lm') + 
  stat_smooth(aes(y = T_Upper_0.1, x = Topt, group = ThermalGuild), col = 'black', size = 1, method = 'lm') + 
  #geom_point(aes(y = Topt, x = Topt), col = 'black', size = 0.5) + 
  theme_bw() + theme(aspect.ratio = 0.75) + xlab('Topt') +
  ylab('SDM thermal niche limits') + 
  geom_abline()

ggplot(ThermalNicheData_Review %>% filter(ConfidenceCombined == 3)) +
  geom_point(aes(y = T_Skew, x = Topt, col = ThermalGuild), size = 2) + 
  stat_smooth(aes(y = T_Skew, x = Topt, col = ThermalGuild), size = 2, method = 'lm') + 
  #geom_point(aes(y = T_Upper_0.1, x = Topt), col = 'red', size = 2) + 
  #geom_linerange(aes(ymax = T_Lower_0.1, ymin = Topt, x = Topt), col = 'blue') + 
  #geom_linerange(aes(ymax = T_Upper_0.1, ymin = Topt, x = Topt), col = 'red') + 
  #geom_point(aes(y = Topt, x = Topt), col = 'black', size = 0.5) + 
  theme_bw() + theme(aspect.ratio = 0.75) + xlab('Topt') +
  ylab('SDM thermal niche limits')


gridExtra::grid.arrange(ggplot(ThermalNicheData_Review %>% filter(N_Presence > 100)) +
                          geom_point(aes(y = T_Lower_0.5, x = Topt), col = 'blue', size = 2) + 
                          geom_point(aes(y = T_Upper_0.5, x = Topt), col = 'red', size = 2) + 
                          geom_linerange(aes(ymax = T_Lower_0.5, ymin = Topt, x = Topt), col = 'blue') + 
                          geom_linerange(aes(ymax = T_Upper_0.5, ymin = Topt, x = Topt), col = 'red') + 
                          geom_point(aes(y = Topt, x = Topt), col = 'black', size = 0.5) + 
                          theme_bw() + theme(aspect.ratio = 0.75) + xlab('Topt') +
                          ylab('SDM thermal niche limits'), 
                        
                        ggplot(ThermalNicheData_Review) + #%>% filter(T_Lower_NEW >= Topt_NEW | T_Upper_NEW <= Topt_NEW)) + 
                          geom_point(aes(y = T_Lower, x = Topt_NEW, pch = ThermalGuild), col = 'blue', size = 2) + 
                          geom_point(aes(y = T_Upper, x = Topt_NEW, pch = ThermalGuild), col = 'red', size = 2) + 
                          geom_linerange(aes(ymax = T_Lower, ymin = Topt_NEW, x = Topt_NEW), col = 'blue') + 
                          geom_linerange(aes(ymax = T_Upper, ymin = Topt_NEW, x = Topt_NEW), col = 'red') + 
                          geom_point(aes(y = Topt_NEW, x = Topt_NEW), col = 'black', size = 0.5) + 
                          theme_bw() + theme(aspect.ratio = 0.75) + xlab('Topt') +
                          ylab('Occupancy mixed model niche limit'),
                        
                        ggplot(ThermalNicheData_Review) + #%>% filter(T_Lower_NEW >= Topt_NEW | T_Upper_NEW <= Topt_NEW)) + 
                          geom_point(aes(y = T_Lower_Obs, x = Topt_NEW, pch = ThermalGuild), col = 'blue', size = 2) + 
                          geom_point(aes(y = T_Upper_Obs, x = Topt_NEW, pch = ThermalGuild), col = 'red', size = 2) + 
                          geom_linerange(aes(ymax = T_Lower_Obs, ymin = Topt_NEW, x = Topt_NEW), col = 'blue') + 
                          geom_linerange(aes(ymax = T_Upper_Obs, ymin = Topt_NEW, x = Topt_NEW), col = 'red') + 
                          geom_point(aes(y = Topt_NEW, x = Topt_NEW), col = 'black', size = 0.5) + 
                          ylab('Sampling thermal limit') + xlab('Topt') +
                          theme_bw() + theme(aspect.ratio = 0.75), 
                        nrow = 3)






