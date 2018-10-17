# Second script in thermal-niche RLS analysis for species level analysis.
# This section models all the data in qgams and occupancy models and produces the dataset of realized thermal niche estimates for later use. 

# This includes 
# 1. Qgam analysis of abundance
# 2. Occupancy models for estimating niche limits
# 3. Derivation of thermal niche dataset
# 4. Plot of TPCs. 

# Initiated 27/03/2018
# Author: Conor Waldock

dir.create('figures_new')
dir.create('data_derived')

# PREAMBLE ----
# Load libraries and packages ----

# Detach packages from script 1. 
detach("package:geosphere", unload=TRUE)
detach("package:rgeos", unload=TRUE)

# Managing and manipulating data. 
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(summarytools)

# Fitting PCAs 
library(pcaMethods)
library(factoextra)

# Fitting quantile gams. 
library(qgam)

# Fitting occupancy models. 
library(glmmTMB)

# Doing parallel stuff
library(doParallel)

# Initial data management performed in external scripts ----
load(file = 'data_upload/01-organsiing-data24092018.RData')
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------




# CREATE COVARIATE PCA SCORES ----
# Explore the structure of covariates across all of RLS sites ----

# Read in covariate data create in '00' script. 
Covariates <- as_data_frame(read.csv(file =  'data_upload/RLS_Site_Covariates_V2_2018-09-26.csv', row.names = NULL))

# View dataframe in summary-tools package
view(dfSummary(Covariates))

# Filter sites to exclude those that won't be used in modelling  
Covariates <- Covariates %>% filter(SiteCode %in% unique(RLS_19$SiteCode))

# Transform variables with distributions that are far from normal in distributions
Covariates$HumPop50km2 <- log10(Covariates$HumPop50km2 + 1)
Covariates$ReefAreaIn15km <- log10(Covariates$ReefAreaIn15km + 1)
Covariates$npp_mean <- log(Covariates$npp_mean + 1)
Covariates$npp_max <- log(Covariates$npp_max + 1)
Covariates$Silicate.Mean <- log(Covariates$Silicate.Mean + 1)
Covariates$Iron.Mean <- log10(Covariates$Iron.Mean*1000 + 1)
Covariates$Current.Velocity.Mean <- log10(Covariates$Current.Velocity.Mean + 1)
Covariates$Current.Velocity.Lt.max <- log10(Covariates$Current.Velocity.Lt.max + 1)

# NOTE:
# Salinity mean vs. max is highly correalted 
# NPP mean vs. max is highly correlated 
# Oxygen mean and max are highly correlated
# Nitrate shows high number of extreme points that will swamp variation in PCA. It is highly correlated with phosphate anyway. 
# If highly correlated values are input then this will swamp the PCA without and combined information getting through. 

# Remove highly correlated variables. 
Covariates_V2 <- Covariates %>% dplyr::select(SiteCode, HumPop50km2, ReefAreaIn15km, 
                                              npp_mean, Silicate.Mean, Salinity.Mean, 
                                              Phosphate.Mean, Iron.Mean,
                                              Dissolved.oxygen.Mean, Current.Velocity.Mean, Current.Velocity.Lt.max, 
                                              pH)

# Use the psych package to produce pairs plots. 
pdf('figures_new/CovariatePairsPlot.pdf', width = 15, height = 15)
psych::pairs.panels(Covariates_V2[, -1])
dev.off()

# Perform PCA analysis across all sites in RLS ----

# Set up covariates as a matrix. 
PCA_data.input <- as.matrix(Covariates_V2[, -1])

# Compute PCA 
PCA_results <- prcomp(na.omit(PCA_data.input), scale = TRUE) # Just testing a subset. 

# Visualise variance explained by each component using the factoextra package. 
pdf('figures_new/PCA-variance.pdf', width = 4, height = 3)
fviz_eig(PCA_results)
dev.off()


# Plotting PCA axis and saving plots ----

# How much is explained by first 3 components? 
round(summary(PCA_results)$importance, 2)# First 3 components explain 62% variation. 

# Extract results for variables
res.var <- get_pca_var(PCA_results)

# Combine results for the variables. 
VariablePCA_Coords <- cbind(melt(res.var$coord), melt(res.var$contrib)[,3])
levels(VariablePCA_Coords$Var1) <- gsub('.Mean','', VariablePCA_Coords$Var1)
names(VariablePCA_Coords)[c(3,4)] <- c('coord', 'contrib')

# Get the PCA scores for individuals sites. 
res.ind <- get_pca_ind(PCA_results)

# Coordinates: this is what I want to model with. 
res.var1 <- res.var$coord %>% data.frame
res.var1$vars <- rownames(res.var1)


# Produce boxplots of coordinates from indepedent axis. 
pdf('figures_new/PCA-Coords-All-Axis-Boxplot.pdf', width= 10, height = 10)
ggplot(VariablePCA_Coords) + 
  geom_bar(aes(x = Var1, y = coord, fill = contrib), stat = 'identity') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', aspect.ratio = 1)+
  scale_fill_gradientn(colours = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  ylab('Dimension 1') + xlab(NULL) +
  facet_wrap(~Var2, scale = 'free_y')
dev.off() 

# Next take the top three axis and plot up boxplots of contribtuons against individual axis scores for each site
BarAxis1 <- ggplot(res.var1) + 
  geom_bar(aes(x = vars, y = Dim.1, fill = Dim.1), stat = 'identity') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', aspect.ratio = 1)+
  scale_fill_gradientn(colours = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  ylab('Dimension 1') + xlab(NULL) 

BarAxis2 <- ggplot(res.var1) + 
  geom_bar(aes(x = vars, y = Dim.2, fill = Dim.2), stat = 'identity') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', aspect.ratio = 1)+
  scale_fill_gradientn(colours = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  ylab('Dimension 2') + xlab(NULL)

BarAxis3 <- ggplot(res.var1) + 
  geom_bar(aes(x = vars, y = Dim.3, fill = Dim.3), stat = 'identity') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', aspect.ratio = 1)+
  scale_fill_gradientn(colours = c("#00AFBB", "#E7B800", "#FC4E07")) + 
  ylab('Dimension 3') + xlab(NULL)

# This shows the correlation between variables. Positively correlated variables point to the same sides of the plots. 
Arrow1 <- fviz_pca_var(PCA_results, axes = c(1, 2), col.var = "contrib", 
                       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, 
                       title = 'Community Niche PCA axis 1-2') + theme(aspect.ratio = 1, legend.position = 'none')     # Avoid text overlapping)

Arrow2 <- fviz_pca_var(PCA_results, axes = c(2, 3), col.var = "contrib", 
                       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, 
                       title = 'Community Niche PCA axis 2-3') + theme(aspect.ratio = 1, legend.position = 'none')     # Avoid text overlapping)

Arrow3 <- fviz_pca_var(PCA_results, axes = c(1, 3), col.var = "contrib", 
                       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, 
                       title = 'Community Niche PCA axis 1-3') + theme(aspect.ratio = 1, legend.position = 'none')     # Avoid text overlapping)

Inds1 <- fviz_pca_ind(PCA_results, axes = c(1, 2), col.ind = "contrib", geom = c("point"),
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, 
                      title = 'Community Niche PCA axis 1-2') + theme(aspect.ratio = 1, legend.position = 'none')     # Avoid text overlapping)

Inds2 <- fviz_pca_ind(PCA_results, axes = c(2, 3), col.ind  = "contrib", geom = c("point"),
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, 
                      title = 'Community Niche PCA axis 2-3') + theme(aspect.ratio = 1, legend.position = 'none')     # Avoid text overlapping)

Inds3 <- fviz_pca_ind(PCA_results, axes = c(1, 3), col.ind  = "contrib", geom = c("point"),
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, 
                      title = 'Community Niche PCA axis 1-3') + theme(aspect.ratio = 1, legend.position = 'none')     # Avoid text overlapping)

# Produce a multipanel plot summarising the first three PCA axis. 
pdf('figures_new/PCA-multiplot.pdf', width = 15, height = 15)
gridExtra::grid.arrange(
  BarAxis1 , Arrow1  , Arrow3, 
  Inds1   , BarAxis2, Arrow2, 
  Inds3   , Inds2  , BarAxis3, 
  nrow = 3, ncol = 3)
dev.off()

# Combind in site level coordinates
res.ind <- get_pca_ind(PCA_results)
Covariates_V2$PCA_1 <- res.ind$coord[,1]
Covariates_V2$PCA_2 <- res.ind$coord[,2]
Covariates_V2$PCA_3 <- res.ind$coord[,3]
Covariates_V2$PCA_4 <- res.ind$coord[,4]


# What is the relationship between PCA scores and temperature across sites? ----

# Is it likley that these axis just capture information that is already in temperature as a covariate? 
Covariates_WithTemp <- RLS_19 %>% dplyr::select(SiteCode, MeanTemp_CoralWatch) %>% unique() %>% left_join(., Covariates_V2 %>% dplyr::select(SiteCode, PCA_1, PCA_2, PCA_3, PCA_4))

# At a global scale PCA1 is correlated with temperature. But, there is a huge about of variation in temperature not
# explained in PC1. And, for each species there will be less of a signal. 
# This suggests we can model with PC1 and temperature together. 
pdf(file = 'figures_new/PairsWithTemp.pdf', width = 8, height = 8)
psych::pairs.panels(Covariates_WithTemp[-1], scale = T)
dev.off()

# Testing code for PCA on single species before running across all species ----
#LabroidesSites <- RLS_19 %>%  filter(SpeciesName == unique(.$SpeciesName)[1]) %>% dplyr::select(SiteCode) %>% unique() %>% .$SiteCode %>% as.character()
#filter(SpeciesName == 'Labroides dimidiatus') %>% dplyr::select(SiteCode) %>% unique() %>% .$SiteCode %>% as.character()
#Covariates_V2_Labroides <- Covariates_V2 %>% filter(SiteCode %in% LabroidesSites) %>% dplyr::select(-PCA_1, -PCA_2, -PCA_3, -PCA_4)
#PCA_Labroides <- prcomp(na.omit(as.matrix(Covariates_V2_Labroides[,-1])), scale = TRUE) # Just testing a subset. 
#fviz_eig(PCA_Labroides)
#round(summary(PCA_Labroides)$importance, 2)# First 3 components explain 83% variation. 
#res.var_Lab <- get_pca_var(PCA_Labroides)
#res.var_Lab <- res.var_Lab$coord %>% data.frame; 
#res.var_Lab$vars <- rownames(res.var_Lab)
#ggplot(res.var_Lab) + geom_bar(aes(x = vars, y = Dim.1, fill = Dim.1), stat = 'identity') + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', aspect.ratio = 1)
#ggplot(res.var_Lab) + geom_bar(aes(x = vars, y = Dim.2, fill = Dim.2), stat = 'identity') + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', aspect.ratio = 1)
#ggplot(res.var_Lab) + geom_bar(aes(x = vars, y = Dim.3, fill = Dim.3), stat = 'identity') + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', aspect.ratio = 1)
#ggplot(res.var_Lab) + geom_bar(aes(x = vars, y = Dim.4, fill = Dim.4), stat = 'identity') + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', aspect.ratio = 1)
#ggplot(res.var_Lab) + geom_bar(aes(x = vars, y = Dim.5, fill = Dim.5), stat = 'identity') + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', aspect.ratio = 1)
#ggplot(res.var_Lab) + geom_bar(aes(x = vars, y = Dim.6, fill = Dim.6), stat = 'identity') + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', aspect.ratio = 1)
#ggplot(res.var_Lab) + geom_bar(aes(x = vars, y = Dim.7, fill = Dim.7), stat = 'identity') + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none', aspect.ratio = 1)
#res.ind_Lab <- get_pca_ind(PCA_Labroides)
#Covariates_V2_Labroides$PCA_1 <- res.ind_Lab$coord[,1]
#Covariates_V2_Labroides$PCA_2 <- res.ind_Lab$coord[,2]
#Covariates_V2_Labroides$PCA_3 <- res.ind_Lab$coord[,3]
#Covariates_V2_Labroides$PCA_4 <- res.ind_Lab$coord[,4]
#Covariates_WithTemp_SPP <- RLS_19 %>% filter(SiteCode %in% LabroidesSites) %>% dplyr::select(SiteCode, MeanTemp_CoralWatch) %>% unique() %>% left_join(., Covariates_V2_Labroides %>% dplyr::select(SiteCode, PCA_1, PCA_2, PCA_3, PCA_4))
#psych::pairs.panels(Covariates_WithTemp_SPP[-1], scale = T)

# Estimate PCA scores for all species individually across sites within their range ----

# Example input data use testing function.
#Species_Sites_Temp <- RLS_19 %>% filter(SpeciesName == unique(.$SpeciesName)[1]) %>% dplyr::select(SpeciesName, SiteCode, MeanTemp_CoralWatch) %>% unique()
#Species_Sites_Temp <- TestSpecies %>% dplyr::select(SpeciesName, SiteCode, MeanTemp_CoralWatch) %>% unique()
#SpeciesPCAScore(Species_Sites_Temp = Species_Sites_Temp, covariates = Covariates_V2)

# Function to estimate PCA scores per species. 
SpeciesPCAScore <- function(Species_Sites_Temp,          # Input a dataframe of temperature by sitecodes
                            covariates = Covariates_V2){
  
  # Convert species subset to a matrix 
  Covariates_Species <- covariates %>% filter(SiteCode %in% Species_Sites_Temp$SiteCode)
  
  if(sum(apply(Covariates_Species[,-1], 2, sd) == 0) != 0){
    #If a column is singluar it breaks the pc. Remove all singluar columns
    Covariates_Species <- Covariates_Species[, -which(apply(Covariates_Species, 2, sd) == 0)]
  }else{NULL}
  PCA_Species <- prcomp(na.omit(as.matrix(Covariates_Species[,-1])), scale = TRUE) # Just testing a subset. 
  
  # Join up temperature with site codes to ensure consistent sitecode order in single object used for PCA
  Species_Sites_Temp <- left_join(Covariates_Species %>% dplyr::select(SiteCode), Species_Sites_Temp, by = 'SiteCode')
  
  # Create output of PCA for variables to decide how many. 
  PCA_Variances <- round(summary(PCA_Species)$importance, 2)
  PCAs <- colnames( PCA_Variances[, -which(PCA_Variances[2,] < 0.1) ] )
  
  # Create columns in temperature by site code based on PCAs. 
  IndividualPCA <- get_pca_ind(PCA_Species)
  
  ToFill <- matrix(NA, nrow(Covariates_Species), length(PCAs)+3)
  ToFill[,1] <- Species_Sites_Temp$SpeciesName
  ToFill[,2] <- Species_Sites_Temp$SiteCode
  ToFill[,3] <- Species_Sites_Temp$MeanTemp_CoralWatch
  for(i in 1:length(PCAs)){
    ToFill[,i+3] <- IndividualPCA$coord[,i]
  }
  
  # Create dataframe to link with species covariates. 
  PCA_Covarites <- data.frame(ToFill, stringsAsFactors = F)
  names(PCA_Covarites) <- c('SpeciesName', 'SiteCode','MeanTemp_CoralWatch', PCAs) 
  PCA_Covarites[,-c(1:2)] <- apply(PCA_Covarites[,-c(1:2)], 2, as.numeric)
  # psych::pairs.panels(PCA_Covarites[,-1])
  # str(PCA_Covarites)
  
  # Determine correlation between PCAs and temperature 
  TempPCAcor <- matrix(NA, 1, length(PCAs) + 1)
  TempPCAcor[,1] <- unique(Species_Sites_Temp$SpeciesName)
  for(i in 1:length(PCAs)){
    TempPCAcor[,i+1] <- cor(PCA_Covarites$MeanTemp_CoralWatch, 
                            PCA_Covarites[,i+3])
  }
  TempPCAcor <- data.frame(TempPCAcor)
  names(TempPCAcor) <- c('SpeciesName', paste0(PCAs, 'Temp_cor'))
  
  # Bundle up for output
  PCA_Covarites_DF <- as_data_frame(PCA_Covarites) %>% dplyr::select(-MeanTemp_CoralWatch) %>% nest(.key = "PCA_Covarites")
  PCA_Covarites_DF$PCA_Temp_Cors <- list(TempPCAcor) #%>% nest(.key = "PCA_Temp_Cors")
  return(PCA_Covarites_DF) 
}

# Run function in loop
PCA_Outputs <- list()
for(i in 1:length(unique(RLS_19$SpeciesName))){
  PCA_Outputs[[i]] <- tryCatch(SpeciesPCAScore(Species_Sites_Temp = RLS_19 %>% filter(SpeciesName == unique(.$SpeciesName)[i]) %>% dplyr::select(SiteCode, MeanTemp_CoralWatch, SpeciesName) %>% unique(), 
                                               covariates = Covariates_V2), error = function(e) NA)
}

# Obtain species list which produce errors. 
PCA_Fails <- which(cbind(lapply(PCA_Outputs, function(x) length(class(x)) == 1)) == T)

# Species which run PCAs
PCA_Outputs_Run <- do.call(rbind, PCA_Outputs)

# Combine the outputs of the runctions. 
PCA_Covariates_Final  <- PCA_Outputs_Run %>% unnest(PCA_Covarites)
PCA_Correlations_Final <- PCA_Outputs_Run %>% unnest(PCA_Temp_Cors)
hist(as.numeric(PCA_Correlations_Final$PC1Temp_cor))
hist(as.numeric(PCA_Correlations_Final$PC2Temp_cor))
hist(as.numeric(PCA_Correlations_Final$PC3Temp_cor))

# These are the species which produce NAs for PCAs. 
MissingSpecies <- unique(RLS_19$SpeciesName)[!unique(RLS_19$SpeciesName) %in% unique(PCA_Correlations_Final$SpeciesName)]
TestSpecies <- RLS_19 %>% filter(SpeciesName == MissingSpecies[1])

# Combine species x site PCA scores with RLS data 
RLS_20 <- left_join(RLS_19, PCA_Covariates_Final, by = c('SpeciesName', 'SiteCode'))

# Which species are removed? These bugs are fixed in the funciton now. 
#sum(is.na(RLS_20$PC1)) # These are the species' removed. 
#identical(round(RLS_20$MeanTemp_CoralWatch.x) , round(RLS_20$MeanTemp_CoralWatch.y))
#hist(RLS_20$MeanTemp_CoralWatch.x - RLS_20$MeanTemp_CoralWatch.y)

# Save RLS_20 object using in q-gams and SDM script ----
saveRDS(RLS_20, file = 'data_derived/RLS_20-For-Qgams-2019-09-28.rds')
# ----




# NEW MODELLING QGAMS SECTION ----
# Function to constrain absences ----
ConstrainAbsences <- function(Data){
  
  NrowPresences <- Data %>% filter(AbundanceAdult40 > 0) %>% nrow(.)
  NrowAbsences <-  Data %>% filter(AbundanceAdult40 == 0) %>% nrow(.)
  
  if(NrowPresences <= NrowAbsences){
    Data <- rbind(Data %>% filter(AbundanceAdult40 > 0), Data %>% filter(AbundanceAdult40 == 0) %>% .[sample(1:nrow(.), NrowPresences),])
  }else{NULL}
  
  return(Data)}

# For a subset of species (choose by N) fit the qgam model and explore the inflight of PCA on temperature effect
Top10 <- RLS_20 %>% filter(Presence == 1) %>%  group_by(SpeciesName) %>% do(N_Sites = nrow(.)) %>% unnest(N_Sites) %>% .[order(.$N_Sites, decreasing = T),] %>% .$SpeciesName %>% .[1:10]
Bottom10 <- RLS_20 %>% filter(Presence == 1) %>%  group_by(SpeciesName) %>% do(N_Sites = nrow(.)) %>% unnest(N_Sites) %>% .[order(.$N_Sites, decreasing = F),] %>% .$SpeciesName %>% .[1:10]

Species1 <- RLS_20 %>% filter(SpeciesName == Top10[9])# %>% ConstrainAbsences(.)
Species1 <- RLS_20 %>% filter(SpeciesName == Bottom10[9])# %>% ConstrainAbsences(.)
sum(is.na(Species1$PC1)) # These are the species' removed. 

# transform, scale and centre columns for modelling
Species1$AbundanceAdult40_log           <- as.numeric(log(Species1$AbundanceAdult40 + 1))
Species1$SamplingIntensity_Scaled       <- as.numeric(scale(Species1$SamplingIntensity))
Species1$MeanTemp_CoralWatch_Scaled     <- as.numeric(scale(Species1$MeanTemp_CoralWatch))
Species1$Depth_Site                     <- as.numeric(scale(Species1$Depth_Site))
Species1$NEOLI[is.na(Species1$NEOLI)]   <- 0
Species1$PC1                            <- as.numeric(Species1$PC1)
Species1$PC2                            <- as.numeric(Species1$PC2)
Species1$PC3                            <- as.numeric(Species1$PC3)
Species1 <- Species1[,-which(apply(apply(Species1, 2, is.na), 2, sum) > 1)]

# Fit qgam model 
summary(glmmTMB(AbundanceAdult40_log ~ poly(Depth_Site,2) + poly(PC1,2) + poly(PC2,2) + poly(PC3,2) + SamplingIntensity_Scaled + NEOLI,
                data = Species1, family = gaussian, zi = ~ 1))
GLMM_Resid <- resid(glmmTMB(AbundanceAdult40_log ~ poly(Depth_Site,2) + poly(PC1,2) + poly(PC2,2) + poly(PC3,2) + SamplingIntensity_Scaled + NEOLI,
                            data = Species1, family = gaussian, zi = ~ 1))
hist(GLMM_Resid)

Species1$GLMM_Resid <- GLMM_Resid
Species1_GAM <- qgam(GLMM_Resid ~ s(MeanTemp_CoralWatch_Scaled, k = 4), qu = 0.8, data = Species1)
Species1_GAM_V2 <- qgam(AbundanceAdult40_log ~ 
                          s(MeanTemp_CoralWatch_Scaled, k = 4) + 
                          s(Depth_Site, k = 4) + 
                          s(PC1, k = 4) + 
                          s(PC2, k = 4) + 
                          s(PC3, k = 4) + 
                          SamplingIntensity_Scaled + 
                          NEOLI
                        , qu = 0.8, data = Species1)

summary(Species1_GAM)
summary(Species1_GAM_V2)

plot(Species1_GAM, pages = 1)
plot(Species1_GAM_V2, pages = 1)

#ResidQgam <- resid(Species1_GAM)

par(mfrow = c(1,2))
plot(AbundanceAdult40_log ~ MeanTemp_CoralWatch_Scaled, data = Species1)
#plot(ResidQgam ~ MeanTemp_CoralWatch_Scaled, data = Species1)
plot(GLMM_Resid ~ MeanTemp_CoralWatch_Scaled, data = Species1)



# Function to fit quantile gams accounting for covariates ----

# This function fits all quantile gam models. Bootstrapping available where 
# subsetting of occupancy is necessary using the NRUNS term. 
FitModels_27092015 <- function(Species1, k = 4, q = 0.8, NRUNS = 1, Folder = NULL){
  
  # Check before with a particular dataset - but at present there are 2 NAs which cause problems with 1 species. 
  Species1 <- Species1[!is.na(Species1$MeanTemp_CoralWatch),]
  
  # Ensure balanced design of positive abundance values and absences. 
  NrowPresences <- Species1 %>% filter(AbundanceAdult40 > 0) %>% nrow(.)
  NrowAbsences  <- Species1 %>% filter(AbundanceAdult40 == 0) %>% nrow(.)
  
  # Scale and transform variables as appropriate
  Species1$AbundanceAdult40_log           <- as.numeric(log(Species1$AbundanceAdult40 + 1))
  Species1$SamplingIntensity_Scaled       <- as.numeric(scale(Species1$SamplingIntensity))
  Species1$MeanTemp_CoralWatch_Scaled     <- as.numeric(scale(Species1$MeanTemp_CoralWatch))
  Species1$Depth_Site                     <- as.numeric(scale(Species1$Depth_Site))
  Species1$NEOLI[is.na(Species1$NEOLI)]   <- 0
  
  # Remove and columsn with NA values as this malfunctions qgam. 
  Species1 <- Species1[,-which(apply(apply(Species1, 2, is.na), 2, sum) > 1)]
  
  # Convert varying number of PCs to numerics  
  PC_Columns <- Species1[,grepl('PC', colnames(Species1))]
  PC_Columns <- apply(PC_Columns, 2, as.numeric)
  for(i in 1:ncol(PC_Columns)){
    Species1[,colnames(Species1[,grepl('PC', colnames(Species1))])[i]] <- PC_Columns[,i]
  }
  
  # Create a model structure that can be flexible to number of PCs. 
  BaselineFormula <- 'AbundanceAdult40_log ~ poly(Depth_Site, 2) + SamplingIntensity_Scaled + NEOLI'
  PC_Columns <- Species1[,grepl('PC', colnames(Species1))]
  for(i in 1:ncol(PC_Columns)){
    BaselineFormula <- paste0(BaselineFormula, ' + ', 'poly(', colnames(Species1[,grepl('PC', colnames(Species1))])[i], ',2)')
  }
  BaselineFormula <- formula(BaselineFormula)
  
  # Columns for backtransformation of temperature data
  Species1$MeanTemp_CoralWatch_scaleBT  <- attr(scale(Species1$MeanTemp_CoralWatch), "scaled:scale")
  Species1$MeanTemp_CoralWatch_centerBT <- attr(scale(Species1$MeanTemp_CoralWatch), "scaled:center")
  
  # Bootstrapping procedure for estimation of Topt and its SD.
  Temps <- Species1 %>% filter(AbundanceAdult40 > 0) %>% do(MinTemp = min(.$MeanTemp_CoralWatch_Scaled), MaxTemp = max(.$MeanTemp_CoralWatch_Scaled)) %>% unnest()
  TempData <- data.frame(MeanTemp_CoralWatch_Scaled = seq(Temps$MinTemp, Temps$MaxTemp, length.out = 1000))
  
  # Objects to assign inside bootstrap 
  Topt <- c()
  T_Gam_pvalue <- c()
  T_Gam_edf <- c()
  T_Gam_deviance.exp <- c()
  
  
  # If there is subsampling then need to bootstrap to get more reliable value for Topt. 
  if(NrowPresences < NrowAbsences){
    for(i in 1:NRUNS){
      
      Species1_Subset <- rbind(Species1 %>% filter(AbundanceAdult40 > 0), Species1 %>% filter(AbundanceAdult40 == 0) %>% .[sample(1:nrow(.), NrowPresences),])
      
      # Extract residuals from fitted values. 
      GLMM_Resid <- resid(glmmTMB(BaselineFormula,
                                  data = Species1_Subset, family = gaussian, zi = ~ 1))
      # Fit model to residuals
      Species1_GAM <- qgam(GLMM_Resid ~ s(MeanTemp_CoralWatch_Scaled, k = k), qu = q, data = Species1_Subset)
      TempData$Abundance <- predict(Species1_GAM, TempData, se = T)$fit
      Topt[i]         <- (TempData[which(TempData$Abundance == max(TempData$Abundance)), 'MeanTemp_CoralWatch_Scaled'] * Species1$MeanTemp_CoralWatch_scaleBT[1]) + Species1$MeanTemp_CoralWatch_centerBT[1]
      T_Gam_pvalue[i] <- summary(Species1_GAM)$s.table[4]
      T_Gam_edf[i]   <- summary(Species1_GAM)$s.table[1]
      T_Gam_deviance.exp[i] <- summary(Species1_GAM)$dev.expl
    }
    
    
  }else{ # If not then use single model. 
    
    # Extract residuals from fitted values. 
    GLMM_Resid <- resid(glmmTMB(AbundanceAdult40_log ~ poly(Depth_Site,2) + poly(PC1,2) + poly(PC2,2) + poly(PC3,2) + SamplingIntensity_Scaled + NEOLI,
                                data = Species1, family = gaussian, zi = ~ 1))
    # Fit model to residuals
    Species1_GAM <- qgam(GLMM_Resid ~ s(MeanTemp_CoralWatch_Scaled, k = k), qu = q, data = Species1)
    TempData$Abundance <- predict(Species1_GAM, TempData, se = T)$fit
    Topt         <- (TempData[which(TempData$Abundance == max(TempData$Abundance)), 'MeanTemp_CoralWatch_Scaled'] * Species1$MeanTemp_CoralWatch_scaleBT[1]) + Species1$MeanTemp_CoralWatch_centerBT[1]
    T_Gam_pvalue <- summary(Species1_GAM)$s.table[4]
    T_Gam_edf    <- summary(Species1_GAM)$s.table[1]
    T_Gam_deviance.exp <- summary(Species1_GAM)$dev.expl
    
  }
  rm(TempData)
  
  # Extract temperture at maximum abundance, and maximum abundance at optimum temperature. 
  #Temps <- Species1 %>% filter(AbundanceAdult40 > 0) %>% do(MinTemp = min(.$MeanTemp_CoralWatch_Scaled), MaxTemp = max(.$MeanTemp_CoralWatch_Scaled)) %>% unnest()
  #TempData <- data.frame(MeanTemp_CoralWatch_Scaled = seq(Temps$MinTemp, Temps$MaxTemp, length.out = 1000))
  #AbunPredictions <- predict(Species1_GAM, TempData, se = T)
  #TempData$Abundance <- AbunPredictions$fit
  #TempData$SE <- AbunPredictions$se
  #TempData$Upr <- AbunPredictions$fit + AbunPredictions$se*1.96
  #TempData$Lwr <- AbunPredictions$fit - AbunPredictions$se*1.96
  
  # Estimate optimum temperature based on maximum abundance predicted 
  #Topt_Raw <- TempData[which(TempData$Abundance == max(TempData$Abundance)), 'MeanTemp_CoralWatch_Scaled']# * Species1$MeanTemp_CoralWatch_scaleBT[1]) + Species1$MeanTemp_CoralWatch_centerBT[1]
  #Topt <- (TempData[which(TempData$Abundance == max(TempData$Abundance)), 'MeanTemp_CoralWatch_Scaled'] * Species1$MeanTemp_CoralWatch_scaleBT[1]) + Species1$MeanTemp_CoralWatch_centerBT[1]
  
  # Take SE at Topt
  #SE <- TempData$SE[which.min(abs((TempData$MeanTemp_CoralWatch_Scaled*Species1$MeanTemp_CoralWatch_scaleBT[1] + Species1$MeanTemp_CoralWatch_centerBT[1]) - Topt))]
  #(quantile(rnorm(mean = Topt_Raw, sd = SE, n = 10000), c(0.025, 0.975)) * Species1$MeanTemp_CoralWatch_scaleBT[1]) + Species1$MeanTemp_CoralWatch_centerBT[1]
  
  # Estimate variation in temperature below thermal optimum. OLD METHOD. 
  # Tsd <- Species1 %>% filter(MeanTemp_CoralWatch < Topt) %>% filter(AbundanceAdult40 > 0) %>% .$MeanTemp_CoralWatch %>% sd(., na.rm = T)
  # Tsd_V2 <- Species1 %>% filter(AbundanceAdult40 > 0) %>% .$MeanTemp_CoralWatch %>% sd(., na.rm = T)/2
  
  ### --- Estimate max abundance
  MaxAbundance <- Species1 %>% .$MaxAbundance %>% unique
  
  ### --- Estimate confidence score criteria 4. 
  T_Opt_Difference_Upper <- max(Species1$MeanTemp_CoralWatch) - mean(Topt)
  T_Opt_Difference_Lower <- mean(Topt) - min(Species1$MeanTemp_CoralWatch)
  
  
  ### --- Estimate confidence score criteria 5. 
  # Example extraction
  #T_Gam_pvalue <- summary(Species1_GAM)$s.table[4] # < 0.05
  #T_Gam_edf    <- summary(Species1_GAM)$s.table[1] # > 1
  
  ### --- Create output dataframe of predictions from model to save with data (and plot up later)
  PredData <- Species1 %>% filter(AbundanceAdult40 > 0) %>% 
    do(SpeciesName = unique(.$SpeciesName), 
       MeanTemp_CoralWatch        = seq(min(.$MeanTemp_CoralWatch), max(.$MeanTemp_CoralWatch), length.out = 100), 
       MeanTemp_CoralWatch_Scaled = seq(min(.$MeanTemp_CoralWatch_Scaled), max(.$MeanTemp_CoralWatch_Scaled), length.out = 100)) %>% 
    unnest(SpeciesName) %>% unnest(MeanTemp_CoralWatch, MeanTemp_CoralWatch_Scaled)
  
  AbunPredictions <- predict(Species1_GAM, PredData, se = T)
  PredData$Abundance <- AbunPredictions$fit
  PredData$SE  <- AbunPredictions$se
  PredData$Upr <- AbunPredictions$fit + AbunPredictions$se*1.96
  PredData$Lwr <- AbunPredictions$fit - AbunPredictions$se*1.96
  PredData$Topt <- mean(Topt)
  PredData$Topt_SD <- sd(Topt)
  PredData$MaxAbundance <- MaxAbundance
  PredData$q <- q
  PredData$k <- k
  
  OutputData <- list(model = Species1_GAM, predictions = PredData)
  
  ### --- Save qgam outputs for if needed later
  dir.create(paste0(Folder))
  saveRDS(OutputData, file = paste(paste0(Folder, '/'), gsub(' ' , '_', unique(Species1$SpeciesName)), gsub('-','_',Sys.Date()), 'q_is', gsub('0.','',q),'k_is',k,'.rds', sep = '_'))
  rm(Species1_GAM)
  
  Topt2 <- mean(Topt)
  Topt_SD <- sd(Topt)
  return(data_frame(SpeciesName = unique(Species1$SpeciesName), 
                    Topt = Topt2, 
                    Topt_SD = Topt_SD,
                    MaxAbundance = MaxAbundance, 
                    T_Opt_Difference_Upper = T_Opt_Difference_Upper,
                    T_Opt_Difference_Lower = T_Opt_Difference_Lower,
                    T_Gam_pvalue = mean(T_Gam_pvalue), 
                    T_Gam_edf = mean(T_Gam_edf), 
                    T_Gam_deviance.exp = mean(T_Gam_deviance.exp)))
}

# Testing function 
TestSpecies <- RLS_20 %>% filter(SpeciesName == unique(.$SpeciesName)[349]) # %>% ConstrainAbsences(.)
TestSpecies <- RLS_20 %>% filter(SpeciesName == Bottom10[9])# %>% ConstrainAbsences(.)
FitModels_27092015(Species1 = TestSpecies, Folder = 'TEST')

# save.image('data_derived/SaveObjectsScript2-2018-27-09.RData')
# load('data_derived/SaveObjectsScript2-2018-27-09.RData')

# Fitting models for all species - ESTIMATING TOPT ----

# Read in RLS data to model with. 
RLS_20 <- readRDS(file = 'data_derived/RLS_20-For-Qgams-2019-09-28.rds')

view(dfSummary(RLS_20))

# Runs function in parallel and saves qgams into a new folder. 
a <- Sys.time()
cl <- makeCluster(4)
registerDoParallel(cl)
ModelOutputs <- foreach(i=1:length(unique(RLS_20$SpeciesName)), 
                        .packages=c('tidyr', 'dplyr', 'glmmTMB', 'qgam')) %dopar% {
                          tryCatch(FitModels_27092015(RLS_20[which(RLS_20$SpeciesName == unique(RLS_20$SpeciesName)[i]),], 
                                                      Folder = 'data_derived/AllQgamModels_2018-10-03', NRUNS = 25), error = function(e) NA)
                        }
stopCluster(cl)
b <- Sys.time() - a # Takes ~5 hours to run. 
b

# Combine output and save as a file. This is the parameters for the Topt. 
Quantile_Parameters <- do.call(rbind, ModelOutputs)
saveRDS(Quantile_Parameters, file = 'data_derived/qGamModelOutputs-2019-10-03.rds')
Quantile_Parameters <- readRDS(file = 'data_derived/qGamModelOutputs-2019-10-03.rds')
# ---- 



# PLOT OLD TOPT AND NEW TOPT VALUES FOR REVIEW ----
ThermalNicheData_Old <- readRDS(file = 'data_upload/ThermalNicheData_Old.rds')
ThermalNicheData_Old$Topt_OLD <- ThermalNicheData_Old$Topt

load(file = 'data_derived/qGamModelOutputs-2019-10-03.rds')

ToptComparison <- left_join(ThermalNicheData_Old[,c('SpeciesName', 'Topt_OLD', 'ConfidenceCombined')], Quantile_Parameters[,c('SpeciesName', 'Topt')])

# Plot models and see relationship. 
pdf(file = 'figures_new/ComparisonTopt-NewVsOld.pdf', height = 5, width = 5)
ggplot() + 
  geom_point(data = ToptComparison, aes(x = Topt_OLD, y = Topt), alpha = 0.5, pch = 1) + 
  geom_point(data = ToptComparison %>% filter(ConfidenceCombined==3), aes(x = Topt_OLD, y = Topt), col = 'red2', size = 3) + 
  theme_bw() + theme(panel.grid = element_blank(), aspect.ratio = 0.75) + 
  xlab('Topt - no covariates') + ylab('Topt - with covariates') + 
  geom_abline()
dev.off()

summary(lm(ToptComparison$Topt ~ ToptComparison$Topt_OLD))

#Residuals:
#   Min      1Q  Median      3Q     Max 
# -6.5402 -0.5012 -0.0630  0.6040  4.4961 
#
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)             -1.10976    0.33421  -3.321 0.000945 ***
#   ToptComparison$Topt_OLD  1.03942    0.01296  80.205  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

# Residual standard error: 1.389 on 701 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.9017,	Adjusted R-squared:  0.9016 
# F-statistic:  6433 on 1 and 701 DF,  p-value: < 2.2e-16


# END OF SCRIPT ----
