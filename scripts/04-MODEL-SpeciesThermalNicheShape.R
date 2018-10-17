# Third script in thermal-niche RLS analysis for species level analysis.
# This section models the thermal niche parameters derived in script 2. 

# This includes 
# 1. Tests of abundance declines at thermal niche edges (X^2 tests)
# 2. Analysis of thermal niche skew vs. Topt (figure 2) 
# 3. Supporting analyses with different data subsets (SOM)

# This code runs the analyses described in section 2.2 and 2.4 of the main manuscript. 

# Initiated 09/10/2018
# Author: Conor Waldock

# ----------------------------------------------------------------------------
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

# Custom functions ----
myvis.gam <- function (x, view = NULL, cond = list(), n.grid = 30, too.far = 0, 
                       col = NA, color = "heat", contour.col = NULL, se = -1, type = "link", 
                       plot.type = "persp", zlim = NULL, nCol = 50, ..., ManualCols = c('dark orange', 'gray80', 'dark blue'))
{
  fac.seq <- function(fac, n.grid) {
    fn <- length(levels(fac))
    gn <- n.grid
    if (fn > gn) 
      mf <- factor(levels(fac))[1:gn]
    else {
      ln <- floor(gn/fn)
      mf <- rep(levels(fac)[fn], gn)
      mf[1:(ln * fn)] <- rep(levels(fac), rep(ln, fn))
      mf <- factor(mf, levels = levels(fac))
    }
    mf
  }
  dnm <- names(list(...))
  v.names <- names(x$var.summary)
  if (is.null(view)) {
    k <- 0
    view <- rep("", 2)
    for (i in 1:length(v.names)) {
      ok <- TRUE
      if (is.matrix(x$var.summary[[i]])) 
        ok <- FALSE
      else if (is.factor(x$var.summary[[i]])) {
        if (length(levels(x$var.summary[[i]])) <= 1) 
          ok <- FALSE
      }
      else {
        if (length(unique(x$var.summary[[i]])) == 1) 
          ok <- FALSE
      }
      if (ok) {
        k <- k + 1
        view[k] <- v.names[i]
      }
      if (k == 2) 
        break
    }
    if (k < 2) 
      stop("Model does not seem to have enough terms to do anything useful")
  }
  else {
    if (sum(view %in% v.names) != 2) 
      stop(paste(c("view variables must be one of", v.names), 
                 collapse = ", "))
    for (i in 1:2) if (!inherits(x$var.summary[[view[i]]], 
                                 c("numeric", "factor"))) 
      stop("Don't know what to do with parametric terms that are not simple numeric or factor variables")
  }
  ok <- TRUE
  for (i in 1:2) if (is.factor(x$var.summary[[view[i]]])) {
    if (length(levels(x$var.summary[[view[i]]])) <= 1) 
      ok <- FALSE
  }
  else {
    if (length(unique(x$var.summary[[view[i]]])) <= 1) 
      ok <- FALSE
  }
  if (!ok) 
    stop(paste("View variables must contain more than one value. view = c(", 
               view[1], ",", view[2], ").", sep = ""))
  if (is.factor(x$var.summary[[view[1]]])) 
    m1 <- fac.seq(x$var.summary[[view[1]]], n.grid)
  else {
    r1 <- range(x$var.summary[[view[1]]])
    m1 <- seq(r1[1], r1[2], length = n.grid)
  }
  if (is.factor(x$var.summary[[view[2]]])) 
    m2 <- fac.seq(x$var.summary[[view[2]]], n.grid)
  else {
    r2 <- range(x$var.summary[[view[2]]])
    m2 <- seq(r2[1], r2[2], length = n.grid)
  }
  v1 <- rep(m1, n.grid)
  v2 <- rep(m2, rep(n.grid, n.grid))
  newd <- data.frame(matrix(0, n.grid * n.grid, 0))
  for (i in 1:length(x$var.summary)) {
    ma <- cond[[v.names[i]]]
    if (is.null(ma)) {
      ma <- x$var.summary[[i]]
      if (is.numeric(ma)) 
        ma <- ma[2]
    }
    if (is.matrix(x$var.summary[[i]])) 
      newd[[i]] <- matrix(ma, n.grid * n.grid, ncol(x$var.summary[[i]]), 
                          byrow = TRUE)
    else newd[[i]] <- rep(ma, n.grid * n.grid)
  }
  names(newd) <- v.names
  newd[[view[1]]] <- v1
  newd[[view[2]]] <- v2
  if (type == "link") 
    zlab <- paste("linear predictor")
  else if (type == "response") 
    zlab <- type
  else stop("type must be \"link\" or \"response\"")
  fv <- predict.gam(x, newdata = newd, se.fit = TRUE, type = type)
  z <- fv$fit
  if (too.far > 0) {
    ex.tf <- exclude.too.far(v1, v2, x$model[, view[1]], 
                             x$model[, view[2]], dist = too.far)
    fv$se.fit[ex.tf] <- fv$fit[ex.tf] <- NA
  }
  if (is.factor(m1)) {
    m1 <- as.numeric(m1)
    m1 <- seq(min(m1) - 0.5, max(m1) + 0.5, length = n.grid)
  }
  if (is.factor(m2)) {
    m2 <- as.numeric(m2)
    m2 <- seq(min(m1) - 0.5, max(m2) + 0.5, length = n.grid)
  }
  if (se <= 0) {
    old.warn <- options(warn = -1)
    av <- matrix(c(0.5, 0.5, rep(0, n.grid - 1)), n.grid, 
                 n.grid - 1)
    options(old.warn)
    max.z <- max(z, na.rm = TRUE)
    z[is.na(z)] <- max.z * 10000
    z <- matrix(z, n.grid, n.grid)
    surf.col <- t(av) %*% z %*% av
    surf.col[surf.col > max.z * 2] <- NA
    if (!is.null(zlim)) {
      if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
        stop("Something wrong with zlim")
      min.z <- zlim[1]
      max.z <- zlim[2]
    }
    else {
      min.z <- min(fv$fit, na.rm = TRUE)
      max.z <- max(fv$fit, na.rm = TRUE)
    }
    surf.col <- surf.col - min.z
    surf.col <- surf.col/(max.z - min.z)
    surf.col <- round(surf.col * nCol)
    con.col <- 1
    if (color == "heat") {
      pal <- heat.colors(nCol)
      con.col <- 3
    }
    else if (color == "topo") {
      pal <- topo.colors(nCol)
      con.col <- 2
    }
    else if (color == "cm") {
      pal <- cm.colors(nCol)
      con.col <- 1
    }
    else if (color == "terrain") {
      pal <- terrain.colors(nCol)
      con.col <- 2
    }
    else if (color == "gray" || color == "bw") {
      pal <- gray(seq(0.1, 0.9, length = nCol))
      con.col <- 1
    }
    ### customized here
    else if (color == 'manual') {
      pal <- colorRampPalette(ManualCols, space = "Lab")(nCol)
      con.col = 1
    }
    ####
    else stop("color scheme not recognised")
    if (is.null(contour.col)) 
      contour.col <- con.col
    surf.col[surf.col < 1] <- 1
    surf.col[surf.col > nCol] <- nCol
    if (is.na(col)) 
      col <- pal[as.array(surf.col)]
    z <- matrix(fv$fit, n.grid, n.grid)
    if (plot.type == "contour") {
      stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                    ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), 
                    ifelse("main" %in% dnm, "", ",main=zlab"), ",...)", 
                    sep = "")
      if (color != "bw") {
        txt <- paste("image(m1,m2,z,col=pal,zlim=c(min.z,max.z)", 
                     stub, sep = "")
        eval(parse(text = txt))
        txt <- paste("contour(m1,m2,z,col=contour.col,zlim=c(min.z,max.z)", 
                     ifelse("add" %in% dnm, "", ",add=TRUE"), ",...)", 
                     sep = "")
        eval(parse(text = txt))
      }
      else {
        txt <- paste("contour(m1,m2,z,col=1,zlim=c(min.z,max.z)", 
                     stub, sep = "")
        eval(parse(text = txt))
      }
    }
    else {
      stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                    ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), 
                    ifelse("main" %in% dnm, "", ",zlab=zlab"), ",...)", 
                    sep = "")
      if (color == "bw") {
        op <- par(bg = "white")
        txt <- paste("persp(m1,m2,z,col=\"white\",zlim=c(min.z,max.z) ", 
                     stub, sep = "")
        eval(parse(text = txt))
        par(op)
      }
      else {
        txt <- paste("persp(m1,m2,z,col=col,zlim=c(min.z,max.z)", 
                     stub, sep = "")
        eval(parse(text = txt))
      }
    }
  }
  else {
    if (color == "bw" || color == "gray") {
      subs <- paste("grey are +/-", se, "s.e.")
      lo.col <- "gray"
      hi.col <- "gray"
    }
    else {
      subs <- paste("red/green are +/-", se, "s.e.")
      lo.col <- "green"
      hi.col <- "red"
    }
    if (!is.null(zlim)) {
      if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
        stop("Something wrong with zlim")
      min.z <- zlim[1]
      max.z <- zlim[2]
    }
    else {
      z.max <- max(fv$fit + fv$se.fit * se, na.rm = TRUE)
      z.min <- min(fv$fit - fv$se.fit * se, na.rm = TRUE)
    }
    zlim <- c(z.min, z.max)
    z <- fv$fit - fv$se.fit * se
    z <- matrix(z, n.grid, n.grid)
    if (plot.type == "contour") 
      warning("sorry no option for contouring with errors: try plot.gam")
    stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                  ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), ifelse("zlab" %in% 
                                                                         dnm, "", ",zlab=zlab"), ifelse("sub" %in% dnm, 
                                                                                                        "", ",sub=subs"), ",...)", sep = "")
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% 
                                                             dnm, "", ",border=lo.col"), stub, sep = "")
    eval(parse(text = txt))
    par(new = TRUE)
    z <- fv$fit
    z <- matrix(z, n.grid, n.grid)
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% 
                                                             dnm, "", ",border=\"black\""), stub, sep = "")
    eval(parse(text = txt))
    par(new = TRUE)
    z <- fv$fit + se * fv$se.fit
    z <- matrix(z, n.grid, n.grid)
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% 
                                                             dnm, "", ",border=hi.col"), stub, sep = "")
    eval(parse(text = txt))
  }
}

GaussianFunction <- function(Topt = Topt, TUpper = TUpper, Tsd = Tsd, Tsd_2 = Tsd_2, MaxPerformance = MaxPerformance){
  
  Temp = seq(0, TUpper+5, length.out = 1000)
  
  Performance <- data.frame(Temp = Temp, 
                            Performance = NA, 
                            Topt = Topt, 
                            TUpper = TUpper, 
                            Tsd = Tsd, 
                            Tsd_2 = Tsd_2,
                            MaxPerformance = MaxPerformance)
  
  for(i in 1:length(Temp)){
    if(Temp[i] <= Topt){
      Performance$Performance[i] <- MaxPerformance*(exp(-(((Temp[i] - Topt) / (Tsd)))^2)) # Remove the 2* TSD as we estimate the whole SD 
    }else{
      # Performance$Performance[i] <- MaxPerformance*(1-((Temp[i] - Topt) / (Topt - if(Confidence_Upper != 0){Topt+UpperDiff}else{TUpper}))^2)
      Performance$Performance[i] <- MaxPerformance*(exp(-(((Temp[i] - Topt) / (Tsd_2)))^2)) # Remove the 2* TSD as we estimate the whole SD 
    }
  }
  
  Performance$Performance[which(Performance$Performance < 0)] <- 0
  
  return(Performance)
}

# Load objects from previous scripts to run below script ----
ThermalNicheData <- readRDS(file = 'data_upload/ThermalNicheData.rds')
RLS_All          <- readRDS(file = 'data_derived/RLS_20-For-Qgams-2019-09-28.rds')
# ----------------------------------------------------------------------------



# ANALYSIS OF SHAPE PATTERNS IN QGAMS AND 'ECOLOGICAL PERFORMANCE CURVES' ----
# Extract model runs and predict curves from models for conf 3 ----

# THE BELOW SCRIPT REQUIRES RUNNING ALL THE QGAMS AND SAVING LOCALLY. 
# ALTERNATIVELY SEE LINES XXX FOR RUN OF OBJECT THAT IS THE QGAM PREDICTIONS BOUND TOGETHER

# Finds the most recent model runs. 
#QuantileModels <- list.files('data_derived/AllQgamModels_2018-10-03') # For on CW local computer. 
QuantileModels <- list.files('###INSERT FOLDER WITH SPECIES LEVEL QGAMS HERE###')

# Which species' have confidence scores of > 3. 
HighConfidenceSpecies <- ThermalNicheData$SpeciesName[which(ThermalNicheData$ConfidenceCombined == 3)]
HighConfidenceSpecies2 <- gsub(' ', '_', ThermalNicheData$SpeciesName[which(ThermalNicheData$ConfidenceCombined == 3)])

# Extract the location within the vector for each model of a high-confidence species. 
QuantileModels_Species <- do.call(rbind, lapply(HighConfidenceSpecies2, function(x) grep(x, QuantileModels)))

# Subsets list of all models and returns file names of all relevant models. 
qgamFiles_Conf3 <- QuantileModels[QuantileModels_Species]
qgamFiles_Conf3 <- paste('data_derived2/AllQgamModels_2018-10-03/', qgamFiles_Conf3, sep = '')

# Create list of model outputs for all species. 
Conf3_Gams <- list()
for(i in 1:length(qgamFiles_Conf3)){ 
  Conf3_Gams[[i]] <- readRDS(file = qgamFiles_Conf3[i])
}

# Extract predicitons from models and combine together (and deviance explained)
Conf3_Gams_Preds <- do.call(rbind, lapply(Conf3_Gams, function(x) x[2][[1]]))
#saveRDS(Conf3_Gams_Preds, file = 'data_derived2/Conf3_Gams_Preds_2018-10-08.rds')
Conf3_Gams_Preds <- readRDS(file = 'data_derived2/Conf3_Gams_Preds_2018-10-08.rds')



# Extract model runs and predict curves from models for all ----
# Extract for all species
#AllSpecies <- gsub(' ', '_', ThermalNicheData$SpeciesName)

# Extract the location within the vector for each model of a high-confidence species. 
#QuantileModels_AllSpecies <- do.call(rbind, lapply(AllSpecies, function(x) grep(x, QuantileModels)))

# Subsets list of all models and returns file names of all relevant models. 
#qgamFiles_ConfAll <- QuantileModels[QuantileModels_AllSpecies]
qgamFiles_ConfAll <- paste('data_derived2/AllQgamModels_2018-10-03/', QuantileModels, sep = '')

# Create list of model outputs for all species. 
All_Gams <- list()
for(i in 1:length(qgamFiles_ConfAll)){
  print(i)
  All_Gams[[i]] <- readRDS(file = qgamFiles_ConfAll[i])$predictions
}

# Extract predicitons from models and combine together (and deviance explained)
#All_Gams_Preds <- do.call(rbind, lapply(All_Gams, function(x) x[2][[1]]))
All_Gams_Preds <- do.call(rbind, All_Gams)
#saveRDS(All_Gams_Preds, file = 'data_derived2/All_Gams_Preds_2018-10-08.rds')

# Set up data to obtain model runs for all species ----

All_Gams_Preds <- readRDS(file = 'data_derived/All_Gams_Preds_2018-10-08.rds')

# IF WANT RESULTS FOR ALL  SPECIES ONLY RUN THIS PART
HighConfidenceSpecies <- ThermalNicheData$SpeciesName
Conf3_Gams_Preds <- All_Gams_Preds

# IF WANT RESULTS FOR ONLY HIGH CONFIDENCE SPECIES RUN THIS PART
HighConfidenceSpecies <- ThermalNicheData %>% filter(ConfidenceCombined == 3) %>% .$SpeciesName
Conf3_Gams_Preds <- All_Gams_Preds %>% filter(SpeciesName %in% HighConfidenceSpecies)



# How many species have non-significant qgam? ---- 
print(paste(round(-sum(ThermalNicheData$Conf_Qgam)/704*100,3), 'percent'))
sum(ThermalNicheData$Conf_Qgam == 0)

# Test 
HighConfidence <- sum(ThermalNicheData$ConfidenceCombined == 3)
print(paste(round(-sum(ThermalNicheData$T_Gam_pvalue[which(ThermalNicheData$ConfidenceCombined == 3)])/HighConfidence*100,3)+100, 'percent'))

sum(ThermalNicheData$T_Gam_pvalue[which(ThermalNicheData$ConfidenceCombined == 3)] < 0.05) #


# Taken predictions and define abundance groups for both 25% and 50% thresholds ----

# Define high confidence species
#HighConfidenceSpecies <- ThermalNicheData$SpeciesName[which(ThermalNicheData$ConfidenceCombined == 3)]

# Define categories
Category_25 <- c()
Abundance_percent_cool <- c()
Abundance_percent_warm <- c()
TG <- c() # Thermal guild for a subset. 

for(i in 1:length(unique(HighConfidenceSpecies))){
  
  QGAMfunction <- Conf3_Gams_Preds %>% filter(SpeciesName == HighConfidenceSpecies[i])
  RLS_Species <- RLS_All %>% filter(SpeciesName == HighConfidenceSpecies[i])
  ThermalNicheData_Spp <- ThermalNicheData %>% filter(SpeciesName == HighConfidenceSpecies[i])
  
  # Extract abundance and scale
  Abundance <- QGAMfunction$Abundance[which(QGAMfunction$Abundance > 0)]
  Temp <- QGAMfunction$MeanSiteSST_NOAA[which(QGAMfunction$Abundance > 0)]
  scaleAbundance <- Abundance/max(Abundance)
  
  TG[i] <- ThermalNicheData_Spp$ThermalGuild
  
  # No trend: mean percentage decline from peak model performance from 
  MaxRaw <- log(ThermalNicheData_Spp$MaxAbundance + 1)
  Warm_Drop <-1 - ( Abundance[length(Abundance)] / max(Abundance) )
  Cool_Drop <-1 - ( Abundance[1] / max(Abundance) )
  NoTrend_Percentage <- max(Cool_Drop, Warm_Drop)
  NoTrend <- ifelse(NoTrend_Percentage > 0.25, FALSE, TRUE)
  
  # Ramped edges
  WarmEdge <- 1 - scaleAbundance[length(Abundance)] # Amount of decline from max abun
  ColdEdge <- 1 - scaleAbundance[1]   # Amount of decline from max abun
  
  Abundance_percent_warm[i] <- scaleAbundance[length(Abundance)] # Amount of decline from max abun
  Abundance_percent_cool[i] <- scaleAbundance[1]   # Amount of decline from max abun
  
  # Test is edge falls by at least 25% of maximum. 
  ColdEdgeRamp <- ifelse(ColdEdge > 0.25, FALSE, TRUE)
  WarmEdgeRamp <- ifelse(WarmEdge > 0.25, FALSE, TRUE)
  
  # Test if both edges are ramped and therefore are at the centre. 
  AbundanceCentre <- ifelse(ColdEdgeRamp == F & WarmEdgeRamp == F, TRUE, FALSE)
  
  AbundanceEdge   <- ifelse(AbundanceCentre == F & NoTrend == F, TRUE, FALSE)
  
  Categories <- c('No trend', 'Ramped cool-edge', 'Ramped warm-edge', 'Abundance centre', 'Abundance edge')
  
  if(NoTrend == T){ 
    Category_25[i] <- 'No trend'
  }else{
    Category_25[i] <- Categories[which(c(NoTrend, ColdEdgeRamp, WarmEdgeRamp, AbundanceCentre, AbundanceEdge))]
  }
}

# Define categories
Category_50 <- c()
TG_50 <- c() # Thermal guild for a subset. 
for(i in 1:length(unique(HighConfidenceSpecies))){
  
  QGAMfunction <- Conf3_Gams_Preds %>% filter(SpeciesName == HighConfidenceSpecies[i])
  RLS_Species <- RLS_All %>% filter(SpeciesName == HighConfidenceSpecies[i])
  ThermalNicheData_Spp <- ThermalNicheData %>% filter(SpeciesName == HighConfidenceSpecies[i])
  
  # Extract abundance and scale
  Abundance <- QGAMfunction$Abundance
  Temp <- QGAMfunction$MeanSiteSST_NOAA
  scaleAbundance <- Abundance/max(Abundance)
  
  TG_50[i] <- ThermalNicheData_Spp$ThermalGuild
  
  # No trend: mean percentage decline from peak model performance from 
  MaxRaw <- log(ThermalNicheData_Spp$MaxAbundance + 1)
  Warm_Drop <-1 - ( Abundance[length(Abundance)] / max(Abundance) )
  Cool_Drop <-1 - ( Abundance[1] / max(Abundance) )
  NoTrend_Percentage <- max(Cool_Drop, Warm_Drop)
  NoTrend <- ifelse(NoTrend_Percentage > 0.5, FALSE, TRUE)
  
  # Ramped edges
  WarmEdge <- 1 - scaleAbundance[length(Abundance)] # Amount of decline from max abun
  ColdEdge <- 1 - scaleAbundance[1]   # Amount of decline from max abun
  
  # Test is edge falls by at least 25% of maximum. 
  ColdEdgeRamp <- ifelse(ColdEdge > 0.25, FALSE, TRUE)
  WarmEdgeRamp <- ifelse(WarmEdge > 0.25, FALSE, TRUE)
  
  # Test if both edges are ramped and therefore are at the centre. 
  AbundanceCentre <- ifelse(ColdEdgeRamp == F & WarmEdgeRamp == F, TRUE, FALSE)
  
  AbundanceEdge   <- ifelse(AbundanceCentre == F & NoTrend == F, TRUE, FALSE)
  
  Categories <- c('No trend', 'Ramped cool-edge', 'Ramped warm-edge', 'Abundance centre', 'Abundance edge')
  
  if(NoTrend == T){ 
    Category_50[i] <- 'No trend'
  }else{
    Category_50[i] <- Categories[which(c(NoTrend, ColdEdgeRamp, WarmEdgeRamp, AbundanceCentre, AbundanceEdge))]
  }
}

# Tests and summaries for 25% thresholds ----

# Create dataframe of categories for other analyses. 
AbundanceDistribution <- data.frame(SpeciesName = HighConfidenceSpecies, AbundanceDistributionCategory = Category_25)
write.csv(AbundanceDistribution, file = 'data_derived/AbundanceDistributions.csv', row.names = F)

# Proportions of each grouping.
round(table(Category_25) / length(Category_25) * 100, 1)

# Number of species in in group
table(Category_25)

# chisq.test without additional group of no trend. 
chisq.test(x = c(table(Category_25)))
chisq.test(x = c(table(Category_25[Category_25 != 'No trend'])))

# What are the abundance reductions at the edges of species ranges? 
# Estimate abundances at the edges of species ranges as a scaled percentage of maximum from model.

# OVERALL:
round(mean(Abundance_percent_cool), 3)*100
round(mean(Abundance_percent_warm), 3)*100

# ABUNDANCE CENTRE:
round(mean(Abundance_percent_cool[Category_25=='Abundance centre']), 3)*100
round(mean(Abundance_percent_warm[Category_25=='Abundance centre']), 3)*100

# COOL RAMPED:
round(mean(Abundance_percent_cool[Category_25=='Ramped cool-edge']), 3)*100
round(mean(Abundance_percent_warm[Category_25=='Ramped cool-edge']), 3)*100

# WARM RAMPED:
round(mean(Abundance_percent_cool[Category_25=='Ramped warm-edge']), 3)*100
round(mean(Abundance_percent_warm[Category_25=='Ramped warm-edge']), 3)*100

# Number with peak at extreme of distribution
sum(Abundance_percent_cool == 1)
sum(Abundance_percent_warm == 1)


# How do these patterns vary between thermal guilds? 
Category_25_Temperate <- Category_25[which(ThermalNicheData %>% filter(SpeciesName %in% HighConfidenceSpecies) %>% .$ThermalGuild == 'temperate')]
Category_25_Tropical <- Category_25[which(ThermalNicheData %>% filter(SpeciesName %in% HighConfidenceSpecies) %>% .$ThermalGuild == 'tropical')]

table(Category_25_Temperate) / length(Category_25_Temperate)
table(Category_25_Tropical) / length(Category_25_Tropical) 

chisq.test(x = c(table(Category_25_Temperate)))
chisq.test(x = c(table(Category_25_Tropical)))

round(mean(Abundance_percent_cool[TG == 'tropical']), 3)*100
round(mean(Abundance_percent_warm[TG == 'tropical']), 3)*100
round(mean(Abundance_percent_cool[TG == 'temperate']), 3)*100
round(mean(Abundance_percent_warm[TG == 'temperate']), 3)*100

round(mean(Abundance_percent_cool[Category_25=='Abundance centre' & TG == 'tropical']), 3)*100
round(mean(Abundance_percent_warm[Category_25=='Abundance centre' & TG == 'tropical']), 3)*100
round(mean(Abundance_percent_cool[Category_25=='Abundance centre' & TG == 'temperate']), 3)*100
round(mean(Abundance_percent_warm[Category_25=='Abundance centre' & TG == 'temperate']), 3)*100

round(mean(Abundance_percent_cool[Category_25=='Ramped cool-edge' & TG == 'tropical']), 3)*100
round(mean(Abundance_percent_warm[Category_25=='Ramped cool-edge' & TG == 'tropical']), 3)*100
round(mean(Abundance_percent_cool[Category_25=='Ramped cool-edge' & TG == 'temperate']), 3)*100
round(mean(Abundance_percent_warm[Category_25=='Ramped cool-edge' & TG == 'temperate']), 3)*100

round(mean(Abundance_percent_cool[Category_25=='Ramped warm-edge'& TG == 'tropical']), 3)*100
round(mean(Abundance_percent_warm[Category_25=='Ramped warm-edge'& TG == 'tropical']), 3)*100
round(mean(Abundance_percent_cool[Category_25=='Ramped warm-edge'& TG == 'temperate']), 3)*100
round(mean(Abundance_percent_warm[Category_25=='Ramped warm-edge'& TG == 'temperate']), 3)*100

sum(Abundance_percent_warm == 1) / length(Abundance_percent_warm)
sum(Abundance_percent_warm[TG == 'temperate'] == 1) / length(Abundance_percent_warm)
sum(Abundance_percent_warm[TG == 'tropical'] == 1) / length(Abundance_percent_warm)

sum(Abundance_percent_cool == 1) / length(Abundance_percent_cool)
sum(Abundance_percent_cool[TG == 'temperate'] == 1) / length(Abundance_percent_cool)
sum(Abundance_percent_cool[TG == 'tropical'] == 1) / length(Abundance_percent_cool)

# Where niche edges are not contrained by geographic space
sum(c(
sum(Abundance_percent_warm[TG == 'temperate'] == 1) / length(Abundance_percent_warm),
sum(Abundance_percent_cool[TG == 'tropical'] == 1) / length(Abundance_percent_cool)
)
)


# Barplot of %s ---- 
Temperate_cat <- data.frame(table(Category_25_Temperate) / length(Category_25_Temperate)); Temperate_cat$ThermalGuild = 'temperate'; names(Temperate_cat)[1] <- 'DistributionType'
Temperate_cat<- rbind(Temperate_cat, data.frame(DistributionType = 'No trend', Freq = 0, ThermalGuild = 'temperate'))
Tropical_cat  <- data.frame(table(Category_25_Tropical) / length(Category_25_Tropical));   Tropical_cat$ThermalGuild = 'tropical';  names(Tropical_cat)[1] <- 'DistributionType'

Barplot_25_data <- rbind(Temperate_cat, Tropical_cat)
Barplot_25_data$DistributionType <- revalue(Barplot_25_data$DistributionType, c('Ramped cool-edge' = 'Cool skew', 'Ramped warm-edge' = 'Warm skew', 'Abundance centre' = 'Centre'))
Barplot_25_data$DistributionType <- factor(Barplot_25_data$DistributionType, unique(Barplot_25_data$DistributionType)[c(1,3,4,2)])

pdf('figures_new/DistributionType-barplots.pdf', width = 5, height = 4, useDingbats = FALSE, bg = 'transparent')
ggplot(Barplot_25_data) + 
  geom_bar(aes(y = Freq*100, x = DistributionType, fill = ThermalGuild), stat = 'identity', position = 'dodge') + 
  scale_fill_manual(values = c('dark blue', 'dark orange')) + 
  theme_classic() + 
  ylab('% distribution type') + 
  xlab(NULL) + 
  theme(aspect.ratio = 0.75, axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.y = element_text(size = 15),axis.text.y = element_text(size = 12),
        text = element_text(size = 10))#+ 
  #facet_wrap(~ThermalGuild, nrow = 2)
dev.off()

# Boxplots of thermal distribution edges and t-tests ---- 
BoxplotData <- data.frame(Category_25, Abundance_percent_cool, Abundance_percent_warm, TG, SpeciesName = HighConfidenceSpecies)

BoxplotData <- reshape2::melt(BoxplotData, measure.vars = c('Abundance_percent_cool', 'Abundance_percent_warm'))

BoxplotData$variable <- revalue(BoxplotData$variable, c('Abundance_percent_warm' = 'Warm edge', 'Abundance_percent_cool' = 'Cool edge'))
BoxplotData$Category_25 <- revalue(BoxplotData$Category_25, c('Abundance centre' = 'Centre', 
                                                           'Ramped warm-edge' = 'Warm skew', 
                                                           'Ramped cool-edge' = 'Cool skew', 
                                                           'No trend' = 'No trend'))

#BoxplotData <- rbind(BoxplotData, data.frame(Category_25 = 'No trend', TG = 'temperate', SpeciesName = 'filler', variable = c('Warm edge', 'Cool edge'), value = 2))

pdf('figures_new/Edge-abundance-boxplots.pdf', width = 5, height = 4, useDingbats = FALSE, bg = 'transparent')
ggplot(BoxplotData) + 
  geom_boxplot(data = BoxplotData, aes(x = variable, y = value, fill = TG)) + 
  facet_wrap(~Category_25) + 
  scale_colour_manual(values = c('dark blue', 'dark orange')) + 
  scale_fill_manual(values = c('dark blue', 'dark orange')) + 
  theme_classic() + 
  ylab('Proportion of maximum \n modelled abundance') + 
  xlab(NULL) + 
  theme(aspect.ratio = 0.75, axis.text.x = element_text(angle = 45, hjust = 1, size = 13), axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 10),
        text = element_text(size = 10), legend.position = 'none', strip.text = element_text(size = 13)) + 
  coord_cartesian(ylim = c(0,1))
dev.off()

t.test_WR_W <- BoxplotData %>% filter(Category_25 == 'Warm skew', variable == 'Warm edge')
t.test_CR_W <- BoxplotData %>% filter(Category_25 == 'Cool skew', variable == 'Warm edge')
t.test_NT_W <- BoxplotData %>% filter(Category_25 == 'No trend',         variable == 'Warm edge')
t.test_AC_W <- BoxplotData %>% filter(Category_25 == 'Centre', variable == 'Warm edge')

t.test_WR_C <- BoxplotData %>% filter(Category_25 == 'Warm skew', variable == 'Cool edge')
t.test_CR_C <- BoxplotData %>% filter(Category_25 == 'Cool skew', variable == 'Cool edge')
t.test_NT_C <- BoxplotData %>% filter(Category_25 == 'No trend',         variable == 'Cool edge')
t.test_AC_C <- BoxplotData %>% filter(Category_25 == 'Centre', variable == 'Cool edge')

# Warm edges
t.test(t.test_WR_W$value ~ t.test_WR_W$TG)# 0.0002062 ***
t.test(t.test_CR_W$value ~ t.test_CR_W$TG)# 0.03623   * 
t.test(t.test_NT_W$value ~ t.test_NT_W$TG)# 0.01241   * 
t.test(t.test_AC_W$value ~ t.test_AC_W$TG)# 0.3303    NS

# Cool edges
t.test(t.test_WR_C$value ~ t.test_WR_C$TG)# 0.01047   ** 
t.test(t.test_CR_C$value ~ t.test_CR_C$TG)# 0.003201  **
t.test(t.test_NT_C$value ~ t.test_NT_C$TG)# 0.003643  ** 
t.test(t.test_AC_C$value ~ t.test_AC_C$TG)# 0.004205  **


# Tests and summaries for 50% thresholds ----

# Proportions of each grouping.
round(table(Category_50) / length(Category_50) * 100, 1)

# Number of species in in group
table(Category_50)

# chisq.test without additional group of no trend. 
chisq.test(x = c(table(Category_50)))
chisq.test(x = c(table(Category_50[Category_50 != 'No trend'])))

round(mean(Abundance_percent_cool[Category_25=='Abundance centre']), 3)*100
round(mean(Abundance_percent_warm[Category_25=='Abundance centre']), 3)*100

# COOL RAMPED:
round(mean(Abundance_percent_cool[Category_25=='No trend']), 3)*100
round(mean(Abundance_percent_warm[Category_25=='No trend']), 3)*100


# ----------------------------------------------------------------------------



# ANALYSIS OF THERMAL NICHE SHAPE FOR EACH SPECIES. SKEW VS. TOPT ---- 
# Define species habitat associations ----

# Summaries covariates per species 
Species_SitePreferences <- RLS_All %>% filter(AbundanceAdult40 > 0) %>% 
  group_by(SpeciesName) %>% 
  nest() %>% 
  mutate(#SSTmean_Spp = purrr::map(data, ~mean(.$SSTmean_Site, na.rm = T)),
         #SSTrange_Spp = purrr::map(data, ~mean(.$SSTrange_Site, na.rm = T)),
         #HumanPop_Spp = purrr::map(data, ~mean(.$HumanPop_Site, na.rm = T)),
         #Exposure_Spp = purrr::map(data, ~mean(.$Exposure_Site, na.rm = T)), 
         #PredatorB_Spp = purrr::map(data, ~mean(.$PredatorB_Site, na.rm = T)), 
         LiveCoralCover_Spp = purrr::map(data, ~mean(.$LiveCoralCover_Site, na.rm = T)), 
         #ComplexCoralCover_Spp = purrr::map(data, ~mean(.$ComplexCoralCover_Site, na.rm = T)), 
         AlgalCover_Spp = purrr::map(data, ~mean(.$AlgalCover_Site, na.rm = T))) %>% 
  unnest(#SSTmean_Spp, 
         #SSTrange_Spp, 
         #HumanPop_Spp, 
         #Exposure_Spp, 
         #PredatorB_Spp, 
         LiveCoralCover_Spp, 
         #ComplexCoralCover_Spp, 
         AlgalCover_Spp) %>% 
  dplyr::select(-data)

ThermalNicheData$LiveCoralCover_Spp <- NULL
ThermalNicheData$AlgalCover_Spp <- NULL

ThermalNicheData <- left_join(ThermalNicheData, Species_SitePreferences)

# SOM figure: Plot correlations between species habitat associations and thermal traits ---- 
pdf(file = 'figures_new/SOM_Habitat associations vs. Topt.pdf', width = 4, height = 6.5)
gridExtra::grid.arrange(
  ggplot() + 
    geom_point(data = ThermalNicheData , aes(x = Topt, y = LiveCoralCover_Spp)) +
    stat_smooth(data = ThermalNicheData , aes(x = Topt, y = LiveCoralCover_Spp), method = lm , colour = 'gray80', se = F) + 
    theme_classic() + theme(text = element_text(size = 15)) + 
    ylab('Species mean live coral cover') + 
    xlab('Species thermal optima')+ ylim(0,NA),
  
  ggplot() + 
    geom_point(data = ThermalNicheData , aes(x = Topt, y = AlgalCover_Spp)) + 
    stat_smooth(data = ThermalNicheData , aes(x = Topt, y = AlgalCover_Spp), method = lm , colour = 'gray80', se = F) + 
    theme_classic() + theme(text = element_text(size = 15)) + 
    ylab('Species mean algae cover') + 
    xlab('Species thermal optima') + ylim(0,NA)
)
dev.off()
# SOM figure: Plot correlations between site temperature and site habitat ----
# What are the habitat differences across all sites at a larger scale? 
pdf(file = 'figures_new/SOM_GlobalHabitat_ThermalGuild.pdf', width = 5, height = 4)
ggplot(data  = RLS_All[which(RLS_All$AlgalCover_Site<101),] %>% dplyr::select(MeanTemp_CoralWatch, SiteCode, LiveCoralCover_Site, AlgalCover_Site) %>% unique(), aes(x = MeanTemp_CoralWatch)) + 
  geom_point(aes(y = LiveCoralCover_Site, colour = "Coral"), alpha = 0.3) + 
  geom_point(aes(y = AlgalCover_Site,  colour = "Algae"), alpha = 0.3) + 
  stat_smooth(aes(y = AlgalCover_Site, colour = "Algae"), method = glm, se = F) +
  stat_smooth(aes(y = LiveCoralCover_Site, colour = "Coral"), method = glm, se = F) + 
  ylab('Live coral cover') + 
  scale_y_continuous(sec.axis = sec_axis(~.*1, name = "Algae cover"), limits = c(0,100)) + 
  scale_colour_manual(values = c("dark blue", "dark orange")) + 
  labs(colour = "Habitat cover")  + 
  theme_classic() + 
  theme(legend.position = 'bottom') + 
  xlab('Site sea surface temperature')
dev.off()

# Set up data for further analyses below ----

# Filter to high confidence
ThermalNicheData_conf3 <- ThermalNicheData %>% filter(ConfidenceCombined == 3)

# Combine in taxonomic information
SpeciesTaxonomy <- read.csv('data_upload/RLS_SpeciesTaxonomy.csv')
ThermalNicheData_conf3 <- left_join(ThermalNicheData_conf3, SpeciesTaxonomy)
ThermalNicheData <- left_join(ThermalNicheData, SpeciesTaxonomy)

mean(ThermalNicheData_conf3$T_Skew)
sd(ThermalNicheData_conf3$T_Skew)


# ----




# GLOBAL ANALYSIS PRESENTED IN THE MAIN PAPER: Thermal niche shapes for all species (no filtering by confidence scores) ----
# Average skew value ---- 
summary(ThermalNicheData$T_Skew)
median(ThermalNicheData$T_Skew)
IQR(ThermalNicheData$T_Skew)
# Model t-skew with global model (across all species conf = all) TROPICAL ----

ThermalNicheData_tropical_allconf <- ThermalNicheData %>% filter(ThermalGuild == 'tropical')

# Tropical models
TropicalModel_allconf_coral <- lmer(T_Skew ~  
                                      Topt + LiveCoralCover_Spp + AlgalCover_Spp + 
                                      (1|Order / Family / Genus ), data = ThermalNicheData_tropical_allconf) # 

# Test effect of algae. 
TropicalModel_allconf_coral2 <- lmer(T_Skew ~  
                                       Topt + LiveCoralCover_Spp + 
                                       (1|Order / Family / Genus), data = ThermalNicheData_tropical_allconf) # 
anova(TropicalModel_allconf_coral, TropicalModel_allconf_coral2)
# Df    AIC    BIC   logLik deviance  Chisq Chi Df Pr(>Chisq)    
# TropicalModel_allconf_coral2  7 1398.9 1428.8 -692.47   1384.9                           
# TropicalModel_allconf_coral   8 1395.2 1429.2 -689.57   1379.2 5.7977      1    0.01605 *

# Test for effect of live coral cover
TropicalModel_allconf_coral3 <- lmer(T_Skew ~  
                                       Topt + AlgalCover_Spp + 
                                       (1|Order / Family / Genus), data = ThermalNicheData_tropical_allconf) # 

anova(TropicalModel_allconf_coral, TropicalModel_allconf_coral3) # Significant coral effect.  
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# TropicalModel_allconf_coral3  7 1402.8 1432.7 -694.41   1388.8                            
# TropicalModel_allconf_coral   8 1395.2 1429.2 -689.57   1379.2 9.6817      1   0.001861 **

# Test effect of Topt.  
TropicalModel_allconf_coral4 <- lmer(T_Skew ~  
                                       LiveCoralCover_Spp + AlgalCover_Spp + 
                                       (1|Order / Family / Genus ), data = ThermalNicheData_tropical_allconf) # 
anova(TropicalModel_allconf_coral, TropicalModel_allconf_coral4) # Significant topt effect.  
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# TropicalModel_allconf_coral4  7 1887.5 1917.4 -936.76   1873.5                             
# TropicalModel_allconf_coral   8 1395.2 1429.2 -689.57   1379.2 494.38      1  < 2.2e-16 ***


# Remove phylogeny 
TropicalModel_allconf_coral_ML_nophylo <- lm(T_Skew ~  
                                               Topt + AlgalCover_Spp + LiveCoralCover_Spp, data = ThermalNicheData_tropical_allconf, REML = 'ML') # 
TropicalModel_allconf_coral_ML <- lmer(T_Skew ~  
                                         Topt +
                                         LiveCoralCover_Spp +
                                         (1|Order / Family / Genus), data = ThermalNicheData_tropical_allconf, REML = F) # 
AICc(TropicalModel_allconf_coral_ML_nophylo, TropicalModel_allconf_coral_ML)
# df     AICc
# TropicalModel_allconf_coral_ML_nophylo  5 1389.699
# TropicalModel_allconf_coral_ML          7 1399.160

# Define final model
TropicalModel_allconf_coral

# Keep all remaining terms. 
drop1(TropicalModel_allconf_coral)
#                    Df    AIC
# <none>                1395.2
# Topt                1 1887.5
# LiveCoralCover_Spp  1 1402.8
# AlgalCover_Spp      1 1398.9

summary(TropicalModel_allconf_coral)
# Fixed effects:
#                     Estimate Std. Error t value
# (Intercept)        17.745258   0.817763  21.700
# Topt               -0.674938   0.023627 -28.566
# LiveCoralCover_Spp  0.029509   0.009496   3.107
# AlgalCover_Spp     -0.018994   0.007785  -2.440

r.squaredGLMM(TropicalModel_allconf_coral)
# R2m       R2c
# 0.6266017 0.631729
# Check residuals
plot(TropicalModel_allconf_coral)
plot(fitted(TropicalModel_allconf_coral), resid(TropicalModel_allconf_coral, type = 'pearson'))
qqnorm(y = resid(TropicalModel_allconf_coral)); qqline(y = resid(TropicalModel_allconf_coral)) 

# Do we get the same effect when fitted to < 26°C
TropicalModel_allconf_coral_Final_UpperTruncated <- lmer(formula(TropicalModel_allconf_coral), data = ThermalNicheData_tropical_allconf %>% filter(Topt < median(Topt)))
# Fixed Effects:
#  (Intercept)                Topt  LiveCoralCover_Spp      AlgalCover_Spp  
# 12.50596            -0.44447             0.01176            -0.02469
qqnorm(y = resid(TropicalModel_allconf_coral)); qqline(y = resid(TropicalModel_allconf_coral)) 



# Model t-skew with global model (across all species conf = all) TEMPERATE ----

# Select temperate species only 
ThermalNicheData_temp_allconf <- ThermalNicheData %>% filter(ThermalGuild == 'temperate') %>% na.omit()

# Temperate models for algae cover. 
TemperateModel_allconf_algae <- lmer(T_Skew ~  
                                       Topt + AlgalCover_Spp + LiveCoralCover_Spp + 
                                       (1|Order / Family / Genus), data = ThermalNicheData_temp_allconf) # 

# Significant coral effect. 
TemperateModel_allconf_algae2 <- lmer(T_Skew ~  
                                        Topt + AlgalCover_Spp +
                                        (1|Order / Family / Genus), data = ThermalNicheData_temp_allconf) # 
anova(TemperateModel_allconf_algae, TemperateModel_allconf_algae2) # Siginificant temperature effect. 
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# TemperateModel_allconf_algae2  7 649.02 671.25 -317.51   635.02                             
# TemperateModel_allconf_algae   8 635.96 661.37 -309.98   619.96 15.057      1  0.0001043 ***

# Test algae effect 
TemperateModel_allconf_algae3 <- lmer(T_Skew ~  
                                        Topt + LiveCoralCover_Spp + 
                                        (1|Order / Family / Genus), data = ThermalNicheData_temp_allconf) # 
anova(TemperateModel_allconf_algae, TemperateModel_allconf_algae3) # Significant algae effect  
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
# TemperateModel_allconf_algae3  7 638.18 660.41 -312.09   624.18                           
# TemperateModel_allconf_algae   8 635.96 661.37 -309.98   619.96 4.2136      1     0.0401 *


# Test Topt effect 
TemperateModel_allconf_algae4 <- lmer(T_Skew ~  
                                        AlgalCover_Spp + LiveCoralCover_Spp + 
                                        (1|Order / Family / Genus), data = ThermalNicheData_temp_allconf) # 
anova(TemperateModel_allconf_algae, TemperateModel_allconf_algae4) # Highly significant algae effect  
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# TemperateModel_allconf_algae4  7 724.65 746.88 -355.32   710.65                             
# TemperateModel_allconf_algae   8 635.96 661.37 -309.98   619.96 90.682      1  < 2.2e-16 ***

# Remove phylogeny 
TemperateModel_allconf_algae_ML_nophylo <- lm(T_Skew ~  
                                                Topt +
                                                AlgalCover_Spp, data = ThermalNicheData_temp_allconf, REML = F) # 
TemperateModel_allconf_algae_ML <- lmer(T_Skew ~  
                                          Topt +
                                          AlgalCover_Spp +
                                          (1|Order / Family / Genus), data = ThermalNicheData_temp_allconf[!is.na(ThermalNicheData_temp_allconf$AlgalCover_Spp),], REML = F) # 
AICc(TemperateModel_allconf_algae_ML_nophylo, TemperateModel_allconf_algae_ML)
# df     AICc
# TemperateModel_allconf_algae_ML_nophylo  4 647.9228
# TemperateModel_allconf_algae_ML          7 649.6834

# Define final model
TemperateModel_allconf_algae

summary(TemperateModel_allconf_algae)
#                     Estimate Std. Error t-value
# (Intercept)        10.979567   1.003580  10.940
# Topt               -0.496109   0.045559 -10.889
# AlgalCover_Spp     -0.017104   0.008891  -1.924
# LiveCoralCover_Spp  0.090580   0.022813   3.970

r.squaredGLMM(TemperateModel_allconf_algae)
# R2m      R2c
# 0.40    0.58

# Check residuals
plot(TemperateModel_allconf_algae)
plot(fitted(TemperateModel_allconf_algae), resid(TemperateModel_allconf_algae, type = 'pearson'))
qqplot(y = resid(TemperateModel_allconf_algae), x = rnorm(1000)) # Errors are normally distributed
qqnorm(y = resid(TemperateModel_allconf_algae)); qqline(y = resid(TemperateModel_allconf_algae)) # Errors are normally distributed

# Model t-skew with global model (across all species conf = all) TROPICAL < 26°C Topt ----

# Subset data for < 26°C global data 
ThermalNicheData_allconf_Trop26 <- ThermalNicheData %>% filter(ThermalGuild == 'tropical') %>% filter(Topt < median(Topt)) %>% na.omit()

# Global models for algal cover. 
GlobalModel_allconf_conservative <- lmer(T_Skew ~  
                                           Topt + AlgalCover_Spp + LiveCoralCover_Spp + 
                                           (1|Order / Family / Genus), data =ThermalNicheData_allconf_Trop26)

# Test influence of algae association 
GlobalModel_allconf_conservative2 <- lmer(T_Skew ~  
                                            Topt + LiveCoralCover_Spp + 
                                            (1|Order / Family / Genus), data = ThermalNicheData_allconf_Trop26)

anova(GlobalModel_allconf_conservative, GlobalModel_allconf_conservative2) # 
#Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
#GlobalModel_allconf_conservative2  7 807.31 832.29 -396.66   793.31                           
#GlobalModel_allconf_conservative   8 805.51 834.06 -394.76   789.51 3.7975      1    0.05133 .

# Test influence of coral cover. 
GlobalModel_allconf_conservative3 <- lmer(T_Skew ~  
                                            Topt + AlgalCover_Spp + 
                                            (1|Order / Family / Genus), data = ThermalNicheData_allconf_Trop26)
anova(GlobalModel_allconf_conservative, GlobalModel_allconf_conservative3, test="LRT") # Significant interaction between T-opt and thermal guild 
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# GlobalModel_allconf_conservative3  7 804.17 829.15 -395.08   790.17                         
# GlobalModel_allconf_conservative   8 805.51 834.06 -394.76   789.51 0.6535      1     0.4189

GlobalModel_allconf_conservative4 <- lmer(T_Skew ~  
                                            AlgalCover_Spp + LiveCoralCover_Spp + 
                                            (1|Order / Family / Genus), data = ThermalNicheData_allconf_Trop26)
anova(GlobalModel_allconf_conservative, GlobalModel_allconf_conservative4, test="LRT") # Significant effect of algae (dissimilar across realms) 
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# GlobalModel_allconf_conservative4  7 852.05 877.03 -419.03   838.05                             
# GlobalModel_allconf_conservative   8 805.51 834.06 -394.76   789.51 48.539      1  3.239e-12 ***

# Fit final model. 
GlobalModel_allconf_conservative_final <- lmer(T_Skew ~  
                                                 Topt + 
                                                 (1|Order / Family / Genus), data =ThermalNicheData_allconf_Trop26)


drop1(GlobalModel_allconf_conservative)               # All terms are important. 
# T_Skew ~ Topt + AlgalCover_Spp + LiveCoralCover_Spp + (1 | Order/Family/Genus)
# Df    AIC
# <none>                805.51
# Topt                1 852.05
# AlgalCover_Spp      1 807.31
# LiveCoralCover_Spp  1 804.17

anova(GlobalModel_allconf_conservative, test = 'LRT') 
# Analysis of Variance Table
#                    Df Sum Sq Mean Sq F value
# Topt                1 46.386  46.386 39.5875
# AlgalCover_Spp      1 13.780  13.780 11.7599
# LiveCoralCover_Spp  1  0.750   0.750  0.6403

summary(GlobalModel_allconf_conservative)
#                    Estimate Std. Error t value
# (Intercept)        12.50596    1.79201   6.979
# Topt               -0.44447    0.06172  -7.201
# AlgalCover_Spp     -0.02469    0.01277  -1.934
# LiveCoralCover_Spp  0.01176    0.01470   0.800

r.squaredGLMM(GlobalModel_allconf_conservative)
# R2m   R2c
# 0.16 0.21

# Check residuals
plot(GlobalModel_allconf_conservative)
plot(fitted(GlobalModel_allconf_conservative), resid(GlobalModel_allconf_conservative, type = 'pearson'))
qqplot(y = resid(GlobalModel_allconf_conservative), x = rnorm(1000)) # Errors are normally distributed
qqnorm(resid(GlobalModel_allconf_conservative)); qqline(resid(GlobalModel_allconf_conservative)) # Errors are normally distributed

# Explore shape of interaction between covariates in tropical and temperate realms ----  
GAM_Interaction_Coral_Topt <- gamm(T_Skew ~  
                                                 te(Topt, LiveCoralCover_Spp),
                                               random = list(Order = ~1/Family/Genus), data = ThermalNicheData)
GAM_Interaction_Coral_Topt2 <- gamm(T_Skew ~  
                                                  s(Topt) + s(LiveCoralCover_Spp),
                                                random = list(Order = ~1/Family/Genus), data = ThermalNicheData)

GAM_Interaction_Coral_Topt3 <- gamm(T_Skew ~  
                                       LiveCoralCover_Spp + s(Topt),
                                    random = list(Order = ~1/Family/Genus), data = ThermalNicheData)

GAM_Interaction_Coral_Topt4 <- gamm(T_Skew ~  
                                      Topt + s(LiveCoralCover_Spp),
                                    random = list(Order = ~1/Family/Genus), data = ThermalNicheData)

# This interaction is not significant across all species. 
AICc(GAM_Interaction_Coral_Topt, GAM_Interaction_Coral_Topt2, GAM_Interaction_Coral_Topt3, GAM_Interaction_Coral_Topt4)

# Explore in all species. 
GAM_Interaction_Algae_Topt <- gamm(T_Skew ~  
                                                  te(Topt, AlgalCover_Spp), data = ThermalNicheData)
GAM_Interaction_Algae_Topt2 <- gamm(T_Skew ~  
                                                   s(Topt) + s(AlgalCover_Spp), data = ThermalNicheData)

GAM_Interaction_Algae_Topt3 <- gamm(T_Skew ~  
                                      s(Topt) + AlgalCover_Spp, data = ThermalNicheData)

GAM_Interaction_Algae_Topt4 <- gamm(T_Skew ~  
                                      Topt + s(AlgalCover_Spp), data = ThermalNicheData)

# This interaction is not significant across all species. 
AICc(GAM_Interaction_Algae_Topt, GAM_Interaction_Algae_Topt2, GAM_Interaction_Algae_Topt3, GAM_Interaction_Algae_Topt4)

png('figures_new/TropicalSkew_Coral_allconf.png', res = 300, width = 1500, height = 1600)
myvis.gam(GAM_Interaction_Coral_Topt2$gam, color="manual", theta=45, 
          xlab = 'Thermal optima', ylab = 'Coral association', too.far = 0.15, plot.type = 'contour', labcex = 1, contour.col='white', nCol = 10000,
          main = 'T-skew')
dev.off()
png('figures_new/TropicalSkew_Algae_allconf.png', res = 300, width = 1500, height = 1600)
myvis.gam(GAM_Interaction_Algae_Topt2$gam, color="manual", theta=45, 
          xlab = 'Thermal optima', ylab = 'Algae association', too.far = 0.15, plot.type = 'contour', labcex = 1, contour.col='white', nCol = 10000,
          main = 'T-skew')
dev.off()


# Test if slope coefficients are significant in global models ----
# Estimate if difference between parameter values is significant (from https://stats.stackexchange.com/questions/55501/test-a-significant-difference-between-two-slope-values)
z_score = (summary(TropicalModel_allconf_coral)$coef[2,1] - summary(TemperateModel_allconf_algae)$coef[2,1]) / 
  (sqrt((summary(TropicalModel_allconf_coral)$coef[2,2]^2) + (summary(TemperateModel_allconf_algae)$coef[2,2]^2)))
pnorm(-abs(z_score))
hist(rnorm(-abs(z_score), n = 100000), 1000)

# Predict data from tropical models -----

# Create dataframe to predict over. 
ModelTrop_Pred_confall <- plyr::ddply(ThermalNicheData_tropical_allconf, 
                                      .(ThermalGuild), 
                                      plyr::summarize,
                                      Topt = seq(min(Topt), max(Topt), length = 100))

ModelTrop_Pred_confall$LiveCoralCover_Spp <- mean(ThermalNicheData_tropical_allconf$LiveCoralCover_Spp, na.rm = T)
ModelTrop_Pred_confall$AlgalCover_Spp     <- mean(ThermalNicheData_tropical_allconf$AlgalCover_Spp, na.rm = T)

# Predict from model coefficients
se.fit <- predict(TropicalModel_allconf_coral, ModelTrop_Pred_confall, re.form = NA, se = T)$se.fit
preds <- predict(TropicalModel_allconf_coral, ModelTrop_Pred_confall, re.form = NA)
ModelTrop_Pred_confall$y <- preds
ModelTrop_Pred_confall$SeUp  <- preds + se.fit*1.96
ModelTrop_Pred_confall$SeLo  <- preds - se.fit*1.96

# Create partial residuals. 
TropModel_PartialResiduals_confall <- remef(TropicalModel_allconf_coral, fix = c('LiveCoralCover_Spp', 'AlgalCover_Spp'), ran = list("Genus:(Family:Order)" = 1), grouping = T)
TropModel_PartialResiduals_confall = TropModel_PartialResiduals_confall + (mean(ThermalNicheData_tropical_allconf$T_Skew) - mean(TropModel_PartialResiduals_confall))

# Predict data for temperate models ----
ModelTemp_Pred_confall <- plyr::ddply(ThermalNicheData_temp_allconf, 
                                      .(ThermalGuild), 
                                      plyr::summarize,
                                      Topt = seq(min(Topt), max(Topt), length = 100))

ModelTemp_Pred_confall$LiveCoralCover_Spp <- mean(ThermalNicheData_temp_allconf$LiveCoralCover_Spp, na.rm = T)
ModelTemp_Pred_confall$AlgalCover_Spp     <- mean(ThermalNicheData_temp_allconf$AlgalCover_Spp, na.rm = T)

# Predict from model coefficients
se.fit <- predict(TemperateModel_allconf_algae, ModelTemp_Pred_confall, re.form = NA, se = T)$se.fit
preds <- predict(TemperateModel_allconf_algae, ModelTemp_Pred_confall, re.form = NA)
ModelTemp_Pred_confall$y <- preds
ModelTemp_Pred_confall$SeUp  <- preds + se.fit*1.96
ModelTemp_Pred_confall$SeLo  <- preds - se.fit*1.96

# Create partial residuals. 
TempModel_PartialResiduals_confall <- remef(TemperateModel_allconf_algae, fix = c('LiveCoralCover_Spp', 'AlgalCover_Spp'), ran = list("Genus:(Family:Order)" = 1), grouping = T)
TempModel_PartialResiduals_confall = TempModel_PartialResiduals_confall + (mean(ThermalNicheData_temp_allconf$T_Skew) - mean(TempModel_PartialResiduals_confall))

# Create plot of global model ----
Model_Predictions_confall      <- rbind.fill(ModelTrop_Pred_confall, ModelTemp_Pred_confall)

# Select example species near to fitted line
Skew_TempSpp_min <- 'Tetractenos glaber' #<- ThermalNicheData_temp_allconf %>% filter(ThermalGuild == 'temperate') %>% .[order(.$Topt, decreasing = T),] %>% .$SpeciesName %>% .[4]
Skew_TempSpp_max <- 'Pempheris affinis' #ThermalNicheData_temp_allconf %>% filter(ThermalGuild == 'temperate') %>% .[order(.$Topt),] %>% .$SpeciesName %>% .[4]
Skew_TropSpp_min <- 'Chromis flavomaculata'#ThermalNicheData_tropical_allconf %>% filter(ThermalGuild == 'tropical') %>% .[order(.$Topt, decreasing = T),] %>% .$SpeciesName %>% .[1]
Skew_TropSpp_max <- 'Cephalopholis urodeta'#ThermalNicheData_tropical_allconf %>% filter(ThermalGuild == 'tropical') %>% .[order(.$Topt),] %>% .$SpeciesName %>% .[2]

# Select example species near to fitted line
pdf(file = 'figures_new/Topt_GlobalModel_confall.pdf', width = 3.22835, height = 2.421263, useDingbats = F)
ggplot() + 
  geom_point(data = ThermalNicheData_temp_allconf,  
             aes(x = Topt, y = TempModel_PartialResiduals_confall, alpha = as.factor(ConfidenceCombined)), colour = 'darkblue', size = 2) + 
  geom_point(data = ThermalNicheData_tropical_allconf, 
             aes(x = Topt, y = TropModel_PartialResiduals_confall, alpha = as.factor(ConfidenceCombined)), colour = 'dark orange', size = 2) + 
  
  geom_line(data = Model_Predictions_confall %>% filter(ThermalGuild == 'temperate'), aes(y = y, x = Topt), colour = 'black') + 
  geom_line(data = Model_Predictions_confall %>% filter(ThermalGuild == 'tropical'), aes(y = y, x = Topt), colour = 'black') + 
  
  geom_ribbon(data = Model_Predictions_confall %>% filter(ThermalGuild == 'temperate'), aes(x = Topt, ymax = SeUp, ymin = SeLo), alpha = 0.5, fill = 'gray50') + 
  geom_ribbon(data = Model_Predictions_confall %>% filter(ThermalGuild == 'tropical'), aes(x = Topt,  ymax = SeUp, ymin = SeLo), alpha = 0.5, fill = 'gray50') + 
  
  # Plot representative species 
  geom_point(data = ThermalNicheData_temp_allconf %>% filter(SpeciesName == Skew_TempSpp_min),  
             aes(x = Topt, y = TempModel_PartialResiduals_confall[which(ThermalNicheData_temp_allconf$SpeciesName == Skew_TempSpp_min)]), col = 'darkblue' , fill = 'white', pch = 21, size = 3) + 
  geom_point(data = ThermalNicheData_tropical_allconf %>% filter(SpeciesName == Skew_TropSpp_min), 
             aes(x = Topt, y = TropModel_PartialResiduals_confall[which(ThermalNicheData_tropical_allconf$SpeciesName == Skew_TropSpp_min)]),col = 'dark orange', fill = 'white', pch = 21, size = 3) + 
  geom_point(data = ThermalNicheData_temp_allconf %>% filter(SpeciesName == Skew_TempSpp_max),  
             aes(x = Topt, y = TempModel_PartialResiduals_confall[which(ThermalNicheData_temp_allconf$SpeciesName == Skew_TempSpp_max)]), col = 'darkblue' , fill = 'gray75', pch = 21, size = 3) + 
  geom_point(data = ThermalNicheData_tropical_allconf %>% filter(SpeciesName == Skew_TropSpp_max), 
             aes(x = Topt, y = TropModel_PartialResiduals_confall[which(ThermalNicheData_tropical_allconf$SpeciesName == Skew_TropSpp_max)]),col = 'dark orange', fill = 'gray75', pch = 21, size = 3) + 
  
  theme_classic() + theme(text = element_text(size = 10), aspect.ratio = 0.75, legend.position = 'none') + 
  xlab(expression(T["opt"])) + 
  ylab(expression(T["skew"])) + 
  scale_x_continuous(breaks = seq(round(min(Model_Predictions_confall$Topt)), round(max(Model_Predictions_confall$Topt)), by = 5)) + 
  scale_alpha_manual(values = c(0.1,0.1,0.1,0.1,0.1,0.5,1))
dev.off()

# Example species' curves to compare with the above plot ----
# Apply thermal performance curve to all species

# Vector of 4 species that are highlighted in figure 2
Skew_Species <- c(Skew_TempSpp_min, Skew_TempSpp_max, Skew_TropSpp_min, Skew_TropSpp_max)

# Estimate gaussian function for skew species. 
Skew_TPCS <- ThermalNicheData %>% 
  filter(SpeciesName%in%Skew_Species) %>%
  group_by(SpeciesName) %>% 
  do(ThermalPerformanceCurve = GaussianFunction(Topt = .$Topt,
                                                TUpper = .$T_Upper_0.1,
                                                Tsd = .$T_SD_Lower,
                                                Tsd_2 = .$T_SD_Upper,
                                                MaxPerformance = .$MaxAbundance)) %>% 
  ungroup() %>% 
  unnest(ThermalPerformanceCurve)

Trop_MaxSkew       <- Skew_TPCS %>% filter(SpeciesName == Skew_TropSpp_min)
Trop_MaxSkew_warm  <- Skew_TPCS %>% filter(SpeciesName == Skew_TropSpp_max)

Temp_MinSkew       <- Skew_TPCS %>% filter(SpeciesName == Skew_TempSpp_max)
Temp_MinSkew_cool  <- Skew_TPCS %>% filter(SpeciesName == Skew_TempSpp_min)

Trop_MaxSkew_V2_label    <- Trop_MaxSkew[which(Trop_MaxSkew$Performance/Trop_MaxSkew$MaxPerformance == max(Trop_MaxSkew$Performance/Trop_MaxSkew$MaxPerformance)), ]
Trop_MaxSkew_warm_label <- Trop_MaxSkew_warm[which(Trop_MaxSkew_warm$Performance/Trop_MaxSkew_warm$MaxPerformance == max(Trop_MaxSkew_warm$Performance/Trop_MaxSkew_warm$MaxPerformance)), ]
Temp_MinSkew_V2_label    <- Temp_MinSkew[which(Temp_MinSkew$Performance/Temp_MinSkew$MaxPerformance == max(Temp_MinSkew$Performance/Temp_MinSkew$MaxPerformance)), ]
Temp_MinSkew_cool_label <- Temp_MinSkew_cool[which(Temp_MinSkew_cool$Performance/Temp_MinSkew_cool$MaxPerformance == max(Temp_MinSkew_cool$Performance/Temp_MinSkew_cool$MaxPerformance)), ]

Labels <- rbind(Trop_MaxSkew_V2_label, Trop_MaxSkew_warm_label, Temp_MinSkew_V2_label, Temp_MinSkew_cool_label)

pdf('figures_new/Figure2-skew-examples_allconf.pdf', width = 3.22835, height = 2.421263, useDingbats = F)
ggplot() + 
  geom_ribbon(data = Trop_MaxSkew_warm, aes(x = Temp, ymin = 0, ymax = Performance/MaxPerformance), alpha = 0.75, col = 'dark orange', fill = 'gray75', size = 1) + 
  geom_ribbon(data = Temp_MinSkew, aes(x = Temp, ymin = 0, ymax = Performance/MaxPerformance), alpha = 0.75, fill = 'gray75', col = 'darkblue', size = 1) + 
  geom_ribbon(data = Temp_MinSkew_cool, aes(x = Temp, ymin = 0, ymax = Performance/MaxPerformance), alpha = 0.75, col = 'darkblue', fill = 'white', size = 1) + 
  geom_ribbon(data = Trop_MaxSkew, aes(x = Temp, ymin = 0, ymax = Performance/MaxPerformance), alpha = 0.75, fill = 'white', col = 'dark orange', size = 1) +
  geom_text(data = Labels[1,], aes(x = Topt, y=c(1.1), label = SpeciesName), size = 2, colour = 'dark orange', fontface='italic') + 
  geom_text(data = Labels[2,], aes(x = Topt, y=c(1.05), label = SpeciesName), size = 2, colour = 'dark orange', fontface='italic') + 
  geom_text(data = Labels[3,], aes(x = Topt-2, y=c(1.05), label = SpeciesName), size = 2, colour = 'darkblue', fontface='italic') + 
  geom_text(data = Labels[4,], aes(x = Topt-2, y=c(1.1), label = SpeciesName), size = 2, colour = 'darkblue', fontface='italic') + 
  xlim(12, 33) +
  scale_y_continuous(breaks = c(0, 0.2,.4,.6,.8, 1)) + 
  theme_classic() + 
  ylab('Ecological performance') + 
  xlab('Temperature') + 
  theme(text = element_text(size = 10), aspect.ratio = 0.75)
dev.off()

# Create plot of all species distribution limits ----
pdf(file = 'figures_new/Topt_ThermalRange_Comparisons_confall.pdf', width = 3.22835, height = 2.421263, useDingbats = F)
ggplot() + 
  
  # All confidence levels 
  geom_linerange(data = ThermalNicheData %>% filter(ThermalGuild == 'temperate'),  
                 aes(x = Topt, ymax = T_Upper_0.1, ymin = Topt, alpha = as.factor(ConfidenceCombined)), colour = 'darkblue', size = 0.2) + 
  geom_point(data = ThermalNicheData %>% filter(ThermalGuild == 'temperate'),  
             aes(x = Topt, y = T_Upper_0.1, alpha = as.factor(ConfidenceCombined)), colour = 'darkblue', size = 2, pch = 24, fill = 'gray90') + 
  
  geom_linerange(data = ThermalNicheData %>% filter(ThermalGuild == 'tropical'),  
                 aes(x = Topt, ymax = T_Upper_0.1, ymin = Topt, alpha = as.factor(ConfidenceCombined)), colour = 'dark orange', size = 0.2) + 
  geom_point(data = ThermalNicheData %>% filter(ThermalGuild == 'tropical'), 
             aes(x = Topt, y = T_Upper_0.1, alpha = as.factor(ConfidenceCombined)), colour = 'dark orange', size = 2, pch = 24, fill = 'gray90') + 
  
  geom_linerange(data = ThermalNicheData %>% filter(ThermalGuild == 'temperate'),  
                 aes(x = Topt, ymin = T_Lower_0.1, ymax = Topt, alpha = as.factor(ConfidenceCombined)), colour = 'darkblue', size = 0.2) + 
  geom_point(data = ThermalNicheData %>% filter(ThermalGuild == 'temperate'),  
             aes(x = Topt, y = T_Lower_0.1, alpha = as.factor(ConfidenceCombined)), colour = 'darkblue', size = 2,pch = 21, fill = 'gray90') + 
  
  geom_linerange(data = ThermalNicheData %>% filter(ThermalGuild == 'tropical'),  
                 aes(x = Topt, ymin = T_Lower_0.1, ymax = Topt, alpha = as.factor(ConfidenceCombined)), colour = 'dark orange',size = 0.2) + 
  geom_point(data = ThermalNicheData %>% filter(ThermalGuild == 'tropical'), 
             aes(x = Topt, y = T_Lower_0.1, alpha = as.factor(ConfidenceCombined)), colour = 'dark orange', size = 2,pch = 21, fill = 'gray90' ) + 
  
  # High confidence species
  geom_point(data = ThermalNicheData %>% filter(ThermalGuild == 'temperate', ConfidenceCombined == 3),  
             aes(x = Topt, y = T_Upper_0.1), colour = 'darkblue', size = 2, pch = 24, fill = 'gray90', alpha = 1) + 
  geom_point(data = ThermalNicheData %>% filter(ThermalGuild == 'tropical', ConfidenceCombined == 3), 
             aes(x = Topt, y = T_Upper_0.1), colour = 'dark orange', size = 2, pch = 24, fill = 'gray90', alpha = 1) + 
  geom_point(data = ThermalNicheData %>% filter(ThermalGuild == 'temperate', ConfidenceCombined == 3),  
             aes(x = Topt, y = T_Lower_0.1), colour = 'darkblue', size = 2, alpha = 1, pch = 21, fill = 'gray90') + 
  geom_point(data = ThermalNicheData %>% filter(ThermalGuild == 'tropical', ConfidenceCombined == 3), 
             aes(x = Topt, y = T_Lower_0.1), colour = 'dark orange', size = 2, alpha = 1, pch = 21, fill = 'gray90') + 
  
  # Lines
  #stat_smooth(data = ThermalNicheData, 
  #            aes(x = Topt, y = T_Lower_0.1), colour = 'black', alpha = 1, se = F) + 
  #stat_smooth(data = ThermalNicheData, 
  #            aes(x = Topt, y = T_Upper_0.1), colour = 'black', alpha = 1, se = F) + 
  
  #stat_smooth(data = ThermalNicheData %>% filter(ThermalGuild == 'temperate'),  
  #            aes(x = Topt, y = T_Lower_0.1), colour = 'black', alpha = 1, method = lm, se = F) + 
  #stat_smooth(data = ThermalNicheData %>% filter(ThermalGuild == 'tropical'), 
  #            aes(x = Topt, y = T_Upper_0.1), colour = 'black', alpha = 1, method = lm, lty = 2, se = F) + 
#stat_smooth(data = ThermalNicheData %>% filter(ThermalGuild == 'temperate'),  
#            aes(x = Topt, y = T_Upper_0.1), colour = 'black', alpha = 1, method = lm, lty = 2, se = F) + 

geom_abline() + 
  
  theme_classic() + theme(text = element_text(size = 10), aspect.ratio = 0.75, legend.position = 'none') + 
  xlab(expression(T["opt"]))+ 
  ylab('Realized niche edges (°C)') +
  scale_x_continuous(breaks = seq(round(min(ThermalNicheData$Topt)), round(max(ThermalNicheData$Topt)), by = 5)) + 
  scale_alpha_manual(values = c(0.1,0.1,0.1,0.1,0.1,0.5,1))
dev.off()

# Create final models ----
stargazer(TropicalModel_allconf_coral, TemperateModel_allconf_algae, GlobalModel_allconf_conservative,
          type = 'html', out = 'figures_new/Tskew analysis outputs-allconf.htm', 
          dep.var.labels=c("T-skew"), 
          covariate.labels = c('T-opt', 
                               'Coral association', 
                               'Algae association', 
                               'Intercept'),
          column.labels = c('Tropical', 'Temperate', 'Tropical < 26°C'), 
          digits = 2)
# ----

# SUPPORTING ANALYSIS WITH HIGH CONFIDENCE SUBSET ----
# MAIN ANALYSIS 1. Model t-skew with global model (across all species conf = 3) ----

png('figures_new/paris_AllSpecies_conf3.png', res = 300, width = 1500, height = 1500)
psych::pairs.panels(na.omit(ThermalNicheData_conf3) %>% dplyr::select(ThermalGuild,Topt, AlgalCover_Spp, LiveCoralCover_Spp)); dev.off()
png('figures_new/paris_tropical_conf3.png', res = 300, width = 1500, height = 1500)
psych::pairs.panels(na.omit(ThermalNicheData_conf3) %>% filter(ThermalGuild == 'tropical') %>% dplyr::select(ThermalGuild,Topt, AlgalCover_Spp, LiveCoralCover_Spp)); dev.off()
png('figures_new/paris_temperate_conf3.png', res = 300, width = 1500, height = 1500)
psych::pairs.panels(na.omit(ThermalNicheData_conf3) %>% filter(ThermalGuild == 'temperate') %>% dplyr::select(ThermalGuild,Topt, AlgalCover_Spp, LiveCoralCover_Spp)); dev.off()

# Perform model selection for tropical species. ----
ThermalNicheData_conf3_trop <- ThermalNicheData_conf3 %>% filter(ThermalGuild == 'tropical')

# Fit the full model
GlobalModel_tropical <- lmer(T_Skew ~  
                               Topt + 
                               AlgalCover_Spp +
                               LiveCoralCover_Spp +
                               (1|Order / Family / Genus), data = na.omit(ThermalNicheData_conf3_trop))

# Significant of Topt - highly significant. 
GlobalModel_tropical2 <- lmer(T_Skew ~  
                               AlgalCover_Spp +
                               LiveCoralCover_Spp +
                               (1|Order / Family / Genus), data = na.omit(ThermalNicheData_conf3_trop))
anova(GlobalModel_tropical, GlobalModel_tropical2, test = 'LRT')
# Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)    
# GlobalModel_tropical2  7 339.89 357.39 -162.95   325.89                            
# GlobalModel_tropical   8 267.81 287.81 -125.91   251.81 74.08      1  < 2.2e-16 ***
  
# Significance of algae association - non-significant. 
GlobalModel_tropical3 <- lmer(T_Skew ~  
                               Topt + 
                               LiveCoralCover_Spp +
                               (1|Order / Family / Genus), data = na.omit(ThermalNicheData_conf3_trop))
anova(GlobalModel_tropical, GlobalModel_tropical3, test = 'LRT')
# Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
# GlobalModel_tropical3  7 265.81 283.31 -125.91   251.81                        
# GlobalModel_tropical   8 267.81 287.81 -125.91   251.81 1e-04      1     0.9909

# Significance of coral association - highly significant. 
GlobalModel_tropical4 <- lmer(T_Skew ~  
                               Topt + 
                               AlgalCover_Spp +
                               (1|Order / Family / Genus), data = na.omit(ThermalNicheData_conf3_trop))
anova(GlobalModel_tropical, GlobalModel_tropical4, test = 'LRT')
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# GlobalModel_tropical4  7 272.86 290.36 -129.43   258.86                            
# GlobalModel_tropical   8 267.81 287.81 -125.91   251.81 7.0466      1   0.007942 **

# Drop 1 for final model structure. 
drop1(GlobalModel_tropical3, test = 'Chi')
# Df    AIC    LRT   Pr(Chi)    
# <none>                265.81                     
# Topt                1 343.33 79.517 < 2.2e-16 ***
# LiveCoralCover_Spp  1 277.00 13.189 0.0002816 ***
  
# Model diagnostics
plot(GlobalModel_tropical3)
plot(GlobalModel_tropical3)
qqnorm(resid(GlobalModel_tropical3)); qqline(resid(GlobalModel_tropical3))

# Remove taxonomy 
GlobalModel_tropical3_nophylo <- lm(T_Skew ~  
                                      Topt + 
                                      LiveCoralCover_Spp, data = ThermalNicheData_conf3_trop)

GlobalModel_tropical3_ML <- lmer(T_Skew ~  
                                   Topt + 
                                   LiveCoralCover_Spp +
                                   (1|Order / Family / Genus), data = na.omit(ThermalNicheData_conf3_trop), REML = F)

# Taxonomic terms are not essential
AICc(GlobalModel_tropical3_nophylo, GlobalModel_tropical3_ML)

#                               df     AICc
# GlobalModel_tropical3_nophylo  4 260.5299
# GlobalModel_tropical3_ML       7 267.1808


# Perform model selection for temperate species. ----

# Doesn't make sense to fit coral cover here. 

# Select temperate species. 
ThermalNicheData_conf3_temp <- ThermalNicheData_conf3 %>% filter(ThermalGuild == 'temperate')

# Fit the full model
GlobalModel_temperate <- lmer(T_Skew ~  
                               Topt + 
                               AlgalCover_Spp +
                               (1|Order / Family / Genus), data = na.omit(ThermalNicheData_conf3_temp))

# Significant of Topt
GlobalModel_temperate2 <- lmer(T_Skew ~  
                                AlgalCover_Spp +
                                 (1|Order / Family / Genus), data = na.omit(ThermalNicheData_conf3_temp))
anova(GlobalModel_temperate, GlobalModel_temperate2, test = 'LRT')
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# GlobalModel_temperate2  6 136.56 146.22 -62.279  124.557                             
# GlobalModel_temperate   7 101.75 113.03 -43.874   87.748 36.809      1  1.303e-09 ***
  
# Significance of algae association  
GlobalModel_temperate3 <- lmer(T_Skew ~  
                                 Topt +
                                (1|Order / Family / Genus), data = na.omit(ThermalNicheData_conf3_temp))
anova(GlobalModel_temperate, GlobalModel_temperate3, test = 'LRT')
# Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)  
# GlobalModel_temperate3  6 102.46 112.13 -45.231   90.461                          
# GlobalModel_temperate   7 101.75 113.03 -43.874   87.748 2.713      1    0.09953 .

# Drop 1 for final model structure. 
drop1(GlobalModel_temperate, test = 'Chi')
# Df    AIC    LRT   Pr(Chi)    
# <none>            101.75                     
# Topt            1 136.56 36.809 1.303e-09 ***
# AlgalCover_Spp  1 102.46  2.713   0.09953 .  

# Model diagnostics
plot(GlobalModel_temperate_final_GAM$lme)
plot(GlobalModel_temperate_final_GAM$gam)
qqnorm(resid(GlobalModel_temperate_final_GAM$lme)); qqline(resid(GlobalModel_temperate_final_GAM$lme))

# Check for the effect of taxonomic structure 
# Remove taxonomy 
GlobalModel_temperate3_nophylo <- lm(T_Skew ~  
                                      Topt, data = ThermalNicheData_conf3_temp)

GlobalModel_temperate3_ML <- lmer(T_Skew ~  
                                   Topt +
                                   (1|Order / Family / Genus), data = na.omit(ThermalNicheData_conf3_temp), REML = F)

# Taxonomic terms are not essential
AICc(GlobalModel_temperate3_nophylo, GlobalModel_temperate3_ML)

# df      AICc
# GlobalModel_temperate3_nophylo  3  97.41774
# GlobalModel_temperate3_ML       6 105.26145

# Compare model outputs ----
summary(GlobalModel_tropical3)
summary(GlobalModel_temperate3)

# Estimate if difference between parameter values is significant (from https://stats.stackexchange.com/questions/55501/test-a-significant-difference-between-two-slope-values)
z_score = (summary(GlobalModel_tropical3)$coef[2,1] - summary(GlobalModel_temperate3)$coef[2,1]) / 
  (sqrt((summary(GlobalModel_tropical3)$coef[2,2]^2) + (summary(GlobalModel_temperate3)$coef[2,2]^2)))
pnorm(-abs(z_score))
hist(rnorm(-abs(z_score), n = 100000), 1000)

# Predict data from tropical models -----

# Create dataframe to predict over. 
ModelTrop_Pred <- plyr::ddply(ThermalNicheData_conf3_trop, 
                                .(ThermalGuild), 
                                plyr::summarize,
                                Topt = seq(min(Topt), max(Topt), length = 100))

ModelTrop_Pred$LiveCoralCover_Spp <- mean(ThermalNicheData_conf3_trop$LiveCoralCover_Spp, na.rm = T)

# Predict from model coefficients
se.fit <- predict(GlobalModel_tropical3, ModelTrop_Pred, re.form = NA, se = T)$se.fit
preds <- predict(GlobalModel_tropical3, ModelTrop_Pred, re.form = NA)
ModelTrop_Pred$y <- preds
ModelTrop_Pred$SeUp  <- preds + se.fit*1.96
ModelTrop_Pred$SeLo  <- preds - se.fit*1.96

# Create partial residuals. 
TropModel_PartialResiduals <- remef(GlobalModel_tropical3, fix = c('LiveCoralCover_Spp'), ran = list("Genus:(Family:Order)" = 1), grouping = T)
TropModel_PartialResiduals = TropModel_PartialResiduals + (mean(ThermalNicheData_conf3_trop$T_Skew) - mean(TropModel_PartialResiduals))

# Predict data for temperate models ----
ModelTemp_Pred <- plyr::ddply(ThermalNicheData_conf3_temp, 
                              .(ThermalGuild), 
                              plyr::summarize,
                              Topt = seq(min(Topt), max(Topt), length = 100))

# Predict from model coefficients
se.fit <- predict(GlobalModel_temperate3, ModelTemp_Pred, re.form = NA, se = T)$se.fit
preds <- predict(GlobalModel_temperate3, ModelTemp_Pred, re.form = NA)
ModelTemp_Pred$y <- preds
ModelTemp_Pred$SeUp  <- preds + se.fit*1.96
ModelTemp_Pred$SeLo  <- preds - se.fit*1.96

# Create partial residuals. 
TempModel_PartialResiduals <- remef(GlobalModel_temperate3, ran = list("Genus:(Family:Order)" = 1), grouping = T)
TempModel_PartialResiduals = TempModel_PartialResiduals + (mean(ThermalNicheData_conf3_temp$T_Skew) - mean(TempModel_PartialResiduals))

# FIGURE 2. Plot of skews with model fits ----
Model_Predictions      <- rbind.fill(ModelTemp_Pred, ModelTrop_Pred)

# Select example species near to fitted line
Skew_TempSpp_min <- ThermalNicheData_conf3 %>% filter(ThermalGuild == 'temperate') %>% .[order(.$Topt, decreasing = T),] %>% .$SpeciesName %>% .[3]
Skew_TempSpp_max <- ThermalNicheData_conf3 %>% filter(ThermalGuild == 'temperate') %>% .[order(.$Topt),] %>% .$SpeciesName %>% .[1]
Skew_TropSpp_min <- ThermalNicheData_conf3 %>% filter(ThermalGuild == 'tropical') %>% .[order(.$Topt, decreasing = T),] %>% .$SpeciesName %>% .[4]
Skew_TropSpp_max <- ThermalNicheData_conf3 %>% filter(ThermalGuild == 'tropical') %>% .[order(.$Topt),] %>% .$SpeciesName %>% .[4]

pdf(file = 'figures_new/Topt_GlobalModel_V2.pdf', width = 3.22835, height = 2.421263, useDingbats = F)
ggplot() + 
  geom_ribbon(data = Model_Predictions %>% filter(ThermalGuild == 'temperate'), aes(x = Topt, ymax = SeUp, ymin = SeLo), alpha = 0.1, fill = 'darkblue') + 
  geom_ribbon(data = Model_Predictions %>% filter(ThermalGuild == 'tropical'), aes(x = Topt,  ymax = SeUp, ymin = SeLo), alpha = 0.1, fill = 'dark orange') + 
  
  geom_point(data = ThermalNicheData_conf3_temp,  
             aes(x = Topt, y = TempModel_PartialResiduals), colour = 'darkblue', alpha = 1, size = 2) + 
  geom_point(data = ThermalNicheData_conf3_trop, 
             aes(x = Topt, y = TropModel_PartialResiduals), colour = 'dark orange', alpha = 1, size = 2) + 
  
  geom_line(data = Model_Predictions %>% filter(ThermalGuild == 'temperate'), aes(y = y, x = Topt), colour = 'black') + 
  geom_line(data = Model_Predictions %>% filter(ThermalGuild == 'tropical'), aes(y = y, x = Topt), colour = 'black') + 
  
  geom_point(data = ThermalNicheData_conf3_temp %>% filter(SpeciesName == Skew_TempSpp_min),  
             aes(x = Topt, y = TempModel_PartialResiduals[which(ThermalNicheData_conf3_temp$SpeciesName == Skew_TempSpp_min)]), col = 'darkblue' , fill = 'white', pch = 21, size = 3) + 
  geom_point(data = ThermalNicheData_conf3_trop %>% filter(SpeciesName == Skew_TropSpp_min), 
             aes(x = Topt, y = TropModel_PartialResiduals[which(ThermalNicheData_conf3_trop$SpeciesName == Skew_TropSpp_min)]),col = 'dark orange', fill = 'white', pch = 21, size = 3) + 
  
  geom_point(data = ThermalNicheData_conf3_temp %>% filter(SpeciesName == Skew_TempSpp_max),  
             aes(x = Topt, y = TempModel_PartialResiduals[which(ThermalNicheData_conf3_temp$SpeciesName == Skew_TempSpp_max)]), col = 'darkblue' , fill = 'gray75', pch = 21, size = 3) + 
  geom_point(data = ThermalNicheData_conf3_trop %>% filter(SpeciesName == Skew_TropSpp_max), 
             aes(x = Topt, y = TropModel_PartialResiduals[which(ThermalNicheData_conf3_trop$SpeciesName == Skew_TropSpp_max)]),col = 'dark orange', fill = 'gray75', pch = 21, size = 3) + 
  
  theme_classic() + theme(text = element_text(size = 10), aspect.ratio = 0.75) + 
  xlab(expression(T["opt"])) + 
  ylab(expression(T["skew"])) + 
  scale_x_continuous(breaks = seq(round(min(Model_Predictions$Topt)), round(max(Model_Predictions$Topt)), by = 2))
dev.off()




# Plot of distribution limits for 'FIGURE 2' ---- 
pdf(file = 'figures_new/Topt_ThermalRange_Comparisons.pdf', width = 3.22835, height = 2.421263, useDingbats = F)
ggplot() + 
  
  geom_linerange(data = ThermalNicheData_conf3 %>% filter(ThermalGuild == 'temperate'),  
                 aes(x = Topt, ymax = T_Upper_0.25, ymin = Topt), colour = 'darkblue', alpha = 1, size = 0.2) + 
  geom_point(data = ThermalNicheData_conf3 %>% filter(ThermalGuild == 'temperate'),  
             aes(x = Topt, y = T_Upper_0.25), colour = 'darkblue', alpha = 1, size = 2, pch = 21, fill = 'gray90') + 
  stat_smooth(data = ThermalNicheData_conf3 %>% filter(ThermalGuild == 'temperate'),  
              aes(x = Topt, y = T_Upper_0.25), colour = 'black', alpha = 1, method = lm, lty = 2, se = F) + 
  
  geom_linerange(data = ThermalNicheData_conf3 %>% filter(ThermalGuild == 'tropical'),  
                 aes(x = Topt, ymax = T_Upper_0.25, ymin = Topt), colour = 'dark orange', alpha = 1, size = 0.2) + 
  geom_point(data = ThermalNicheData_conf3 %>% filter(ThermalGuild == 'tropical'), 
             aes(x = Topt, y = T_Upper_0.25), colour = 'dark orange', alpha = 1, size = 2, pch = 21, fill = 'gray90') + 
  stat_smooth(data = ThermalNicheData_conf3 %>% filter(ThermalGuild == 'tropical'), 
              aes(x = Topt, y = T_Upper_0.25), colour = 'black', alpha = 1, method = lm, lty = 2, se = F) + 
  
  
  geom_linerange(data = ThermalNicheData_conf3 %>% filter(ThermalGuild == 'temperate'),  
                 aes(x = Topt, ymin = T_Lower_0.25, ymax = Topt), colour = 'darkblue', alpha = 1, size = 0.2) + 
  geom_point(data = ThermalNicheData_conf3 %>% filter(ThermalGuild == 'temperate'),  
             aes(x = Topt, y = T_Lower_0.25), colour = 'darkblue', alpha = 1, size = 2) + 
  stat_smooth(data = ThermalNicheData_conf3 %>% filter(ThermalGuild == 'temperate'),  
              aes(x = Topt, y = T_Lower_0.25), colour = 'black', alpha = 1, method = lm, se = F) + 
  
  geom_linerange(data = ThermalNicheData_conf3 %>% filter(ThermalGuild == 'tropical'),  
                 aes(x = Topt, ymin = T_Lower_0.25, ymax = Topt), colour = 'dark orange', alpha = 1, size = 0.2) + 
  geom_point(data = ThermalNicheData_conf3 %>% filter(ThermalGuild == 'tropical'), 
             aes(x = Topt, y = T_Lower_0.25), colour = 'dark orange', alpha = 1, size = 2) + 
  stat_smooth(data = ThermalNicheData_conf3 %>% filter(ThermalGuild == 'tropical'), 
              aes(x = Topt, y = T_Lower_0.25), colour = 'black', alpha = 1, method = lm, se = F) + 
  
  geom_abline() + 
  
  theme_classic() + theme(text = element_text(size = 10), aspect.ratio = 0.75) + 
  xlab(expression(T["opt"]))+ 
  ylab('Realized niche edges (°C)') +
  scale_x_continuous(breaks = seq(round(min(Model_Predictions$Topt)), round(max(Model_Predictions$Topt)), by = 2))
dev.off()


# Plot skew examples for 'FIGURE 2' ---- 
# Apply thermal performance curve to all species

# Vector of 4 species that are highlighted in figure 2
Skew_Species <- c(Skew_TempSpp_min, Skew_TempSpp_max, Skew_TropSpp_min, Skew_TropSpp_max)

# Estimate gaussian function for skew species. 
Skew_TPCS <- ThermalNicheData %>% 
  filter(SpeciesName%in%Skew_Species) %>%
  group_by(SpeciesName) %>% 
  do(ThermalPerformanceCurve = GaussianFunction(Topt = .$Topt,
                                                TUpper = .$T_Upper_0.1,
                                                Tsd = .$T_SD_Lower,
                                                Tsd_2 = .$T_SD_Upper,
                                                MaxPerformance = .$MaxAbundance)) %>% 
  ungroup() %>% 
  unnest(ThermalPerformanceCurve)

Trop_MaxSkew       <- Skew_TPCS %>% filter(SpeciesName == Skew_TropSpp_min)
Trop_MaxSkew_warm  <- Skew_TPCS %>% filter(SpeciesName == Skew_TropSpp_max)

Temp_MinSkew       <- Skew_TPCS %>% filter(SpeciesName == Skew_TempSpp_max)
Temp_MinSkew_cool  <- Skew_TPCS %>% filter(SpeciesName == Skew_TempSpp_min)

Trop_MaxSkew_V2_label    <- Trop_MaxSkew[which(Trop_MaxSkew$Performance/Trop_MaxSkew$MaxPerformance == max(Trop_MaxSkew$Performance/Trop_MaxSkew$MaxPerformance)), ]
Trop_MaxSkew_warm_label <- Trop_MaxSkew_warm[which(Trop_MaxSkew_warm$Performance/Trop_MaxSkew_warm$MaxPerformance == max(Trop_MaxSkew_warm$Performance/Trop_MaxSkew_warm$MaxPerformance)), ]
Temp_MinSkew_V2_label    <- Temp_MinSkew[which(Temp_MinSkew$Performance/Temp_MinSkew$MaxPerformance == max(Temp_MinSkew$Performance/Temp_MinSkew$MaxPerformance)), ]
Temp_MinSkew_cool_label <- Temp_MinSkew_cool[which(Temp_MinSkew_cool$Performance/Temp_MinSkew_cool$MaxPerformance == max(Temp_MinSkew_cool$Performance/Temp_MinSkew_cool$MaxPerformance)), ]

Labels <- rbind(Trop_MaxSkew_V2_label, Trop_MaxSkew_warm_label, Temp_MinSkew_V2_label, Temp_MinSkew_cool_label)

pdf('figures_new/Figure2-skew-examples.pdf', width = 3.22835, height = 2.421263, useDingbats = F)
ggplot() + 
  geom_ribbon(data = Trop_MaxSkew_warm, aes(x = Temp, ymin = 0, ymax = Performance/MaxPerformance), alpha = 0.75, col = 'dark orange', fill = 'gray75', size = 1) + 
  geom_ribbon(data = Temp_MinSkew, aes(x = Temp, ymin = 0, ymax = Performance/MaxPerformance), alpha = 0.75, fill = 'gray75', col = 'darkblue', size = 1) + 
  geom_ribbon(data = Temp_MinSkew_cool, aes(x = Temp, ymin = 0, ymax = Performance/MaxPerformance), alpha = 0.75, col = 'darkblue', fill = 'white', size = 1) + 
  geom_ribbon(data = Trop_MaxSkew, aes(x = Temp, ymin = 0, ymax = Performance/MaxPerformance), alpha = 0.75, fill = 'white', col = 'dark orange', size = 1) +
  geom_text(data = Labels[1:2,], aes(x = Topt, y=c(1.15, 1.05), label = SpeciesName), size = 2, colour = 'dark orange', fontface='italic') + 
  geom_text(data = Labels[3:4,], aes(x = Topt, y=c(1.05, 1.1), label = SpeciesName), size = 2, colour = 'darkblue', fontface='italic') + 
  xlim(12, 33) +
  scale_y_continuous(breaks = c(0,0.2,.4,.6,.8,1)) + 
  theme_classic() + 
  ylab('Ecological performance') + 
  xlab('Temperature') + 
  theme(text = element_text(size = 10), aspect.ratio = 0.75)
dev.off()

# Fit models excluding species' with Topt < 26 TROPICAL SPECIES ----
ThermalNicheData_conf3_trop_SUBSET <- ThermalNicheData_conf3 %>% filter(ThermalGuild == 'tropical') %>% filter(Topt < median(.$Topt))

# Fit the full model
GlobalModel_tropical_SUBSET <- lmer(T_Skew ~  
                               Topt + LiveCoralCover_Spp + AlgalCover_Spp + 
                               (1|Order / Family / Genus), data = na.omit(ThermalNicheData_conf3_trop_SUBSET))

# Significant of Algae - non-significant
GlobalModel_tropical2_SUBSET <- lmer(T_Skew ~  
                                       Topt + 
                                LiveCoralCover_Spp +
                                (1|Order / Family / Genus), data = na.omit(ThermalNicheData_conf3_trop_SUBSET))
anova(GlobalModel_tropical_SUBSET, GlobalModel_tropical2_SUBSET, test = 'LRT')
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# GlobalModel_tropical2_SUBSET  7 151.21 163.85 -68.603   137.21                         
# GlobalModel_tropical_SUBSET   8 151.59 166.04 -67.794   135.59 1.6178      1     0.2034

# Significance of temperature association - non-significant. 
GlobalModel_tropical3_SUBSET <- lmer(T_Skew ~  
                                       Topt + 
                                       AlgalCover_Spp +
                                       (1|Order / Family / Genus), data = na.omit(ThermalNicheData_conf3_trop_SUBSET))
anova(GlobalModel_tropical_SUBSET, GlobalModel_tropical3_SUBSET, test = 'LRT')
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# GlobalModel_tropical3_SUBSET  7 150.04 162.69 -68.022   136.04                         
# GlobalModel_tropical_SUBSET   8 151.59 166.04 -67.794   135.59 0.4559      1     0.4995

# Significance of coral association - No effect 
GlobalModel_tropical4_SUBSET <- lmer(T_Skew ~  
                                       AlgalCover_Spp +
                                       LiveCoralCover_Spp +
                                       (1|Order / Family / Genus), data = na.omit(ThermalNicheData_conf3_trop_SUBSET))
anova(GlobalModel_tropical_SUBSET, GlobalModel_tropical4_SUBSET, test = 'LRT')
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# GlobalModel_tropical4_SUBSET  7 152.19 164.84 -69.096   138.19                         
# GlobalModel_tropical_SUBSET   8 151.59 166.04 -67.794   135.59 2.6041      1     0.1066

# Produce tables from models ----
summary(GlobalModel_tropical3)
summary(GlobalModel_temperate3)

stargazer(GlobalModel_tropical3, GlobalModel_temperate3, GlobalModel_tropical_SUBSET,
          type = 'html', out = 'figures_new/Tskew analysis outputs.htm', 
          column.labels = c('Tropical', 'Temperate', 'Tropical < 26°C'), digits = 2,  
          dep.var.labels=c("T-skew"), 
          covariate.labels = c('T-opt', 
                               'Coral association', 
                               'Algae association', 
                               'Intercept'))

# ----------------------------------------------------------------------------



# SOM plot of fundamental thermal niche data from globtherm (data not publically available) ----

# Read in Amanda Bates compiled Ctmax and Ctmin data. 
Bates_TND <- read.csv('data_raw/Ctmax and Ctmin.csv')
Bates_TND <- Bates_TND %>% filter(taxon == 'fish')

# Read in GlobTherm compiled data
GlobTherm_TND <- read.csv('data_raw/ctmax_ctmin_globthermo.csv')
unique(GlobTherm_TND$max_metric)
GlobTherm_TND <- GlobTherm_TND %>% filter(max_metric == 'ctmax', max_pretreatment != 'F')
unique(GlobTherm_TND$Tmax)
GlobTherm_TND$max_pretreatment <- as.numeric(as.character(GlobTherm_TND$max_pretreatment))

# Plot together 
pdf('figures_new/SOM_acclimation-temperatures-critical-limits.pdf', width = 5, height = 5)
ggplot() + 
  geom_point(data = Bates_TND %>% filter(tolerance == 'hot'),  aes(x = haccl, y = Critical.Limit.haccl.AEB), col = 'red') + 
  geom_point(data = Bates_TND %>% filter(tolerance == 'cold'), aes(x = laccl, y = Critical.Limit.laccl.AEB), col = 'blue') + 
  #  stat_smooth(data = Bates_TND %>% filter(tolerance == 'hot'),  aes(x = haccl, y = Critical.Limit.haccl.AEB), col = 'red') + 
  #  stat_smooth(data = Bates_TND %>% filter(tolerance == 'cold'), aes(x = haccl, y = Critical.Limit.haccl.AEB), col = 'blue') + 
  geom_point(data = GlobTherm_TND,  aes(x = max_pretreatment, y = Tmax), col = 'red', pch = 22) + 
  geom_point(data = GlobTherm_TND, aes( x = min_pretreatment, y = Tmin), col = 'blue', pch = 22) +  
  #  stat_smooth(data = GlobTherm_TND,  aes(x = max_pretreatment, y = Tmax), col = 'red', pch = 21) + 
  #  stat_smooth(data = GlobTherm_TND, aes( x = max_pretreatment, y = Tmin), col = 'blue', pch = 21) + 
  xlab('Acclimation Temperature (°C)') + 
  ylab('Critical limits') + 
  theme_bw() + 
  theme(aspect.ratio = 0.75) 
dev.off()



# ----------------------------------------------------------------------------

# Data creation 17/10/2018 ----
# Organise excel file and descriptors 

library(dataframes2xls)

ThermalNiche_refined <- ThermalNicheData %>% dplyr::select(SpeciesName, ThermalGuild, 
                                                           Topt, T_Upper_0.1, T_Lower_0.1, 
                                                           T_Upper_Obs, T_Lower_Obs, T_Midpoint_Obs, 
                                                           T_Skew, ConfidenceCombined)

names(ThermalNiche_refined) <- c('SpeciesName', 'ThermalGuild',
                                 'Topt', 'Tmax', 'Tmin', 
                                 'Tmax_obs', 'Tmin_obs', 'Topt_obs',
                                 'Tskew', 'ConfidenceScore')

#ThermalNiche_refined$Topt_conf <- ifelse(3 + ThermalNicheData$Confidence_T_Opt_Difference_Upper + ThermalNicheData$Confidence_T_Opt_Difference_Lower + ThermalNicheData$Conf_Qgam == 3, 1, 0)

ThermalNiche_refined <- as.data.frame(ThermalNiche_refined)


# Create column descriptions
Column_descriptions <- data.frame(Column = names(ThermalNiche_refined)) 

Column_descriptions$Description = c('Species names', 
                                    'Thermal guilds are identified as Topt < 23°C = temperate, Topt > 23°C = tropical', 
                                    'Topt is defined from quantile generalized additive models', 
                                    'Tmax is defined from SDMs',
                                    'Tmin is defined from SDMs',
                                    'Tmax_obs the maximum temperature of species sampled occurrence',
                                    'Tmin_obs the minimum temperature of species sampled occurrence',
                                    'Topt_obs the midpoint of species temperature distribution, the median temperature of species occurence', 
                                    'Tskew is the difference between sigma-Tmax and sigma-Tmin', 
                                    'Confidence score as described in appendix 1')

dataframes2xls::write.xls(c(ThermalNiche_refined, Column_descriptions), "data_upload/RLS-thermal-niche-estimates.xls")


# ----

# END OF SCRIPT ----------------------------------











# FIGURES FOR SUPPORTING ONLINE MATIERIALS ---------------------------------------------------------------------------

# SOM plot of thermal niche limits (results) ----
# Testing variation in limits and edges ----

par(mfrow = c(2,3))
sd(ThermalNicheData_conf3_trop$Topt); hist(ThermalNicheData_conf3_trop$Topt, main = 'Tropical opt; sd = 1.53')
sd(ThermalNicheData_conf3_trop$T_Lower); hist(ThermalNicheData_conf3_trop$T_Lower, main = 'Tropical lower; sd = 0.83')
sd(ThermalNicheData_conf3_trop$T_Upper); hist(ThermalNicheData_conf3_trop$T_Upper, main = 'Tropical upper; sd = 0.71')

sd(ThermalNicheData_conf3_temp$Topt); hist(ThermalNicheData_conf3_temp$Topt, main = 'Temperate opt; sd = 1.67')
sd(ThermalNicheData_conf3_temp$T_Lower); hist(ThermalNicheData_conf3_temp$T_Lower, main = 'Temperate lower; sd = 1.25')
sd(ThermalNicheData_conf3_temp$T_Upper); hist(ThermalNicheData_conf3_temp$T_Upper, main = 'Temperate upper; sd = 1.29')

Hist_data <- ThermalNicheData_conf3 %>% filter(ConfidenceCombined == 3) %>% dplyr::select(SpeciesName, ThermalGuild, Topt, T_Lower, T_Upper)

Hist_data <- tidyr::gather(Hist_data, Parameter, Value, Topt:T_Upper)
Hist_data$Parameter <- factor(Hist_data$Parameter, levels = c('T_Lower', 'Topt', 'T_Upper'), labels = c('Lower', 'Peak', 'Upper'))

Hist_data <- Hist_data %>% group_by(ThermalGuild, Parameter) %>% nest() %>% mutate(SD = purrr::map(data, ~sd(.$Value))) %>% unnest(SD) %>% unnest(data)

pdf('figures_new/SOM_histogram-of-thermal-niche-parameters.pdf', width = 6, height = 4)
ggplot(data = Hist_data) + 
  geom_histogram(aes(Value, fill = ThermalGuild)) + 
  facet_wrap(ThermalGuild ~ Parameter) + 
  geom_text(aes(x = unique(mean(Value)), y = 30, label = paste('SD = ', round(SD, 3), sep = ''))) + 
  scale_fill_manual(values = c('dark blue', 'dark orange')) + 
  theme_light() + 
  theme(aspect.ratio = 1, legend.position = 'none')
dev.off()

# Variation in topt across taxonomic groups show with boxplots ----

# Attach species' taxonomy
ThermalNicheData_v2 <- left_join(ThermalNicheData, SpeciesTaxonomy)
ThermalNicheData_v2$Family <- as.character(ThermalNicheData_v2$Family)

# Extract families with > 20 individuals. 
ThermalNicheData_v2 <- ThermalNicheData_v2 %>% filter(ConfidenceCombined == 3)
Families <- rownames(table(ThermalNicheData_v2$Family)[which(table(ThermalNicheData_v2$Family)>5)])
ThermalNicheData_v3 <- ThermalNicheData_v2 %>% filter(Family %in% Families)
ThermalNicheData_v3_temp <- ThermalNicheData_v3 %>% filter(ThermalGuild == 'temperate')
ThermalNicheData_v3_trop <- ThermalNicheData_v3 %>% filter(ThermalGuild == 'tropical')

Median_Topts_Temp <- aggregate(Topt ~ Family, data = ThermalNicheData_v3_temp, FUN = median)
names(Median_Topts_Temp)[2] <- 'Topt_Median'
ThermalNicheData_v3_temp <- left_join(ThermalNicheData_v3_temp, Median_Topts_Temp)
ThermalNicheData_v3_temp$Family <- factor(ThermalNicheData_v3_temp$Family, levels = unique(ThermalNicheData_v3_temp$Family[order(ThermalNicheData_v3_temp$Topt_Median)]))

Median_Topts_Trop <- aggregate(Topt ~ Family, data = ThermalNicheData_v3_trop, FUN = median)
names(Median_Topts_Trop)[2] <- 'Topt_Median'
ThermalNicheData_v3_trop <- left_join(ThermalNicheData_v3_trop, Median_Topts_Trop)
ThermalNicheData_v3_trop$Family <- factor(ThermalNicheData_v3_trop$Family, levels = unique(ThermalNicheData_v3_trop$Family[order(ThermalNicheData_v3_trop$Topt_Median)]))

pdf('figures_new/SOM-Topt-across-taxonomy.pdf', width = 5, height = 5)
gridExtra::grid.arrange(
  
  ggplot(ThermalNicheData_v3_trop) + 
    geom_boxplot(aes(x = Family, y = Topt)) + 
    facet_wrap(~ThermalGuild) +    
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), aspect.ratio = 0.5) , 
  
  ggplot(ThermalNicheData_v3_temp) + 
    geom_boxplot(aes(x = Family, y = Topt)) + 
    facet_wrap(~ThermalGuild) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), aspect.ratio = 0.5), 
  
  nrow = 2)
dev.off()

# ----------------------------------------------------------------------------






# SOM plot of comparison of thermal range breadths between tropical and temperate species ----
pdf(file = 'figures_new/SOM_thermal-range-comparison.pdf', width = 3, height = 3, useDingbats = F)
ggplot(data = ThermalNicheData %>% filter(ConfidenceCombined == 3), aes(x = ThermalGuild, y = T_Range, fill = ThermalGuild)) + 
  geom_violin(alpha = 0.5, col = NA) + 
  geom_jitter(aes(col = ThermalGuild), width = 0.1, height = 0) + 
  scale_fill_manual(values = c('dark blue', 'dark orange 2')) + 
  scale_colour_manual(values = c('dark blue', 'dark orange 2')) + 
  theme_classic() + 
  ylab('Thermal range') + 
  xlab(NULL) + 
  theme(legend.position = 'none', aspect.ratio = 1)
dev.off()

t.test(ThermalNicheData$T_Range[which(ThermalNicheData$ConfidenceCombined == 3)] ~ ThermalNicheData$ThermalGuild[which(ThermalNicheData$ConfidenceCombined == 3)])
t.test(ThermalNicheData$T_Range ~ ThermalNicheData$ThermalGuild)

# ----------------------------------------------------------------------------





# SOM plot of comparison of rarity between tropical and temperate species ---- 

SppOccupancyRate <- left_join(RLS_All, ThermalNicheData %>% dplyr::select(SpeciesName, Topt)) %>% 
  group_by(SpeciesName) %>% 
  nest() %>%
  mutate(SppOcc = purrr::map(data, ~mean(.$Presence))) %>% 
  unnest(SppOcc) %>% dplyr::select(-data)

ThermalNicheData_Occ <- left_join(SppOccupancyRate, ThermalNicheData)

pdf(file = 'figures_new/SOM-max abundance comparison.pdf', width = 6, height = 3, useDingbats = F)
grid.arrange(
  
  ggplot(data = ThermalNicheData) + 
    geom_boxplot(aes(x = ThermalGuild, y = log10(MaxAbundance), fill = ThermalGuild), draw_quantiles = c(0.25, 0.5, 0.95), col = 'black', alpha = 0.5) + 
    geom_jitter(aes(x = ThermalGuild, y = log10(MaxAbundance), col = ThermalGuild), width = 0.2, height = 0, pch = 21) +
    theme_classic() + 
    theme(legend.position = 'none') + 
    scale_fill_manual(values = c('darkblue', 'dark orange')) + 
    scale_colour_manual(values = c('darkblue', 'dark orange')) + 
    ylab('Log10 maximum abundance') + 
    xlab(NULL),
  
  ggplot(data = ThermalNicheData_Occ) + 
    geom_jitter(aes(x = ThermalGuild, y = SppOcc, col = ThermalGuild), width = 0.2, height = 0, pch = 21) +
    geom_boxplot(aes(x = ThermalGuild, y = SppOcc, fill = ThermalGuild), draw_quantiles = c(0.25, 0.5, 0.95), col = 'black', alpha = 0.5) + 
    theme_classic() + 
    theme(legend.position = 'none') + 
    scale_fill_manual(values = c('darkblue', 'dark orange')) + 
    scale_colour_manual(values = c('darkblue', 'dark orange')) + 
    ylab('Species range occupancy rate') + 
    xlab(NULL), 
  
  ncol = 2
  
)
dev.off()

# No difference in means. 
t.test(SppOcc ~ ThermalGuild, ThermalNicheData_Occ)
t.test(MaxAbundance ~ ThermalGuild, ThermalNicheData_Occ)

# ----------------------------------------------------------------------------





# SOM plot of fundamental thermal niche data from globtherm ----

# Read in Amanda Bates compiled Ctmax and Ctmin data. 
Bates_TND <- read.csv('data_raw/Ctmax and Ctmin.csv')
Bates_TND <- Bates_TND %>% filter(taxon == 'fish')

# Read in GlobTherm compiled data
GlobTherm_TND <- read.csv('data_raw/ctmax_ctmin_globthermo.csv')
unique(GlobTherm_TND$max_metric)
GlobTherm_TND <- GlobTherm_TND %>% filter(max_metric == 'ctmax', max_pretreatment != 'F')
unique(GlobTherm_TND$Tmax)
GlobTherm_TND$max_pretreatment <- as.numeric(as.character(GlobTherm_TND$max_pretreatment))

# Plot together 
pdf('figures_new/SOM_acclimation-temperatures-critical-limits.pdf', width = 5, height = 5)
ggplot() + 
  geom_point(data = Bates_TND %>% filter(tolerance == 'hot'),  aes(x = haccl, y = Critical.Limit.haccl.AEB), col = 'red') + 
  geom_point(data = Bates_TND %>% filter(tolerance == 'cold'), aes(x = laccl, y = Critical.Limit.laccl.AEB), col = 'blue') + 
#  stat_smooth(data = Bates_TND %>% filter(tolerance == 'hot'),  aes(x = haccl, y = Critical.Limit.haccl.AEB), col = 'red') + 
#  stat_smooth(data = Bates_TND %>% filter(tolerance == 'cold'), aes(x = haccl, y = Critical.Limit.haccl.AEB), col = 'blue') + 
  geom_point(data = GlobTherm_TND,  aes(x = max_pretreatment, y = Tmax), col = 'red', pch = 22) + 
  geom_point(data = GlobTherm_TND, aes( x = min_pretreatment, y = Tmin), col = 'blue', pch = 22) +  
#  stat_smooth(data = GlobTherm_TND,  aes(x = max_pretreatment, y = Tmax), col = 'red', pch = 21) + 
#  stat_smooth(data = GlobTherm_TND, aes( x = max_pretreatment, y = Tmin), col = 'blue', pch = 21) + 
  xlab('Acclimation Temperature (°C)') + 
  ylab('Critical limits') + 
  theme_bw() + 
  theme(aspect.ratio = 0.75) 
dev.off()
  


# ----------------------------------------------------------------------------





# Data release 05/06/2018 ----
# Organise excel file and descriptors 

library(dataframes2xls)

ThermalNiche_refined <- ThermalNicheData %>% dplyr::select(SpeciesName, ThermalGuild, 
                                Topt, T_Upper, T_Lower, 
                                T_Upper_Obs, T_Lower_Obs, T_Midpoint_Obs, 
                                T_Upper_Mod,T_Lower_Mod, T_Skew 
                                )
names(ThermalNiche_refined) <- c('SpeciesName', 'ThermalGuild',
                                 'Topt', 'Tmax', 'Tmin', 
                                 'Tmax_obs', 'Tmin_obs', 'Topt_obs', 
                                 'Tmax_model', 'Tmin_model', 'Tskew' )

ThermalNiche_refined$Topt_conf <- ifelse( 3 + ThermalNicheData$Confidence_T_Opt_Difference_Upper + ThermalNicheData$Confidence_T_Opt_Difference_Lower + ThermalNicheData$Conf_Qgam == 3, 1, 0)
ThermalNiche_refined <- as.data.frame(ThermalNiche_refined)


# Create column descriptions
Column_descriptions <- data.frame(Column = names(ThermalNiche_refined)) 

Column_descriptions$Description = c('Species names', 
  'Thermal guilds are identified as Topt < 23°C = temperate, Topt > 23°C = tropical', 
  'Topt is defined from quantile generalized additive models where Topt_conf = 1, and from species distribution midpoints where Topt_conf = 0', 
  'Tmax is defined from occupancy models when tmax_model is != NA, and Tmax_obs otherwise',
  'Tmin is defined from occupancy models when tmin_model is != NA, and Tmin_obs otherwise',
  'Tmax_obs the maximum temperature of species occurrence',
  'Tmin_obs the minimum temperature of species occurrence',
  'Topt_obs the midpoint of species temperature distribution, the median temperature of species occurence', 
  'Tmax_model is defined from occupancy model of species upper thermal distribution limits', 
  'Tmin_model is defined from occupancy model of species lower thermal distribution limits', 
  'Tskew is the difference between sigma-Tmax and sigma-Tmin')


dataframes2xls::write.xls(c(ThermalNiche_refined, Column_descriptions), "data_derived/RLS-thermal-niche-estimates.xls")


# ----




# Save progress ----
#save.image('data_derived/script3_save-image.RData')
load('data_derived/script3_save-image.RData')
# ----





