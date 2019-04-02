#Morphological strategies of corals in the anthropocene: data setup

rm(list = ls())

#Set up ====

#.....Libraries ====

library(tidyverse)
library(modelr)
library(broom)

#.....Functions ====

source("Functions.r")

#' @title PRESS
#' @author Thomas Hopper
#' @description Returns the PRESS statistic (predictive residual sum of squares).
#'              Useful for evaluating predictive power of regression models.
#' @param linear.model A linear regression model (class 'lm'). Required.
#'
PRESS <- function(linear.model) {
  #' calculate the predictive residuals
  pr <- residuals(linear.model)/(1-lm.influence(linear.model)$hat)
  #' calculate the PRESS
  PRESS <- sum(pr^2)

  return(PRESS)
}

### This calculate the Predictive r-squared

#' @title Predictive R-squared
#' @author Thomas Hopper
#' @description returns the predictive r-squared. Requires the function PRESS(), which returns
#'              the PRESS statistic.
#' @param linear.model A linear regression model (class 'lm'). Required.
#'
pred_r_squared <- function(linear.model) {
  #' Use anova() to get the sum of squares for the linear model
  lm.anova <- anova(linear.model)
  #' Calculate the total sum of squares
  tss <- sum(lm.anova$'Sum Sq')
  # Calculate the predictive R^2
  pred.r.squared <- 1-PRESS(linear.model)/(tss)

  return(pred.r.squared)
}


#.....Data ====

LITDat <- read.csv("2.DataSetup/LIT_CoralsWithGrowthForm.csv")
MorphoDat_All <- read.csv("2.DataSetup/MorphologyData_HighQualityWhole5GFs.csv")
#.====
#Build predictive models ====

#.....setup ====

#..........Variables/transforms/formulas ====

ShapeVars <- c("ColonyVol", "OuterSA", "Sphericity", "Convexity", "Packing", "FractalDimension", "SAVol", "ColonyVol2ndMomentVertScaled", "ColonySA2ndMomentVertScaled")
Transforms <- c("log10", "log10", "logit", "logit", "log10", "logit", "log10", "log10", "log10")
Formulas <- paste0(Transforms, "(", ShapeVars, ") ~ log10(PlanarArea) * GrowthForm") %>%
  map(., as.formula)

#..........Load data, select coloumns, gather shape variables, and nest dataframe ====

MorphoDat <- read.csv("2.DataSetup/MorphologyData_HighQualityWhole5GFs.csv") %>%
  select(ScanID, PlanarArea, GrowthForm, ShapeVars) %>%
  filter(GrowthForm %in% "BranchingClosed" == F) %>%
  gather("ShapeVariable", "Value", ShapeVars) %>%
  group_by(ShapeVariable) %>%
  nest()

#..........update value to shape variable name in nested data ====

lapply(1:length(ShapeVars), function(x) names(MorphoDat$data[[x]])[4] <<- ShapeVars[x])

#.....Modelling ====

#..........add formulas transforms and models ====

MorphoDat <- MorphoDat %>%
  mutate(ModelFormula = Formulas,
         Transform = Transforms,
         Model = map2(Formulas, data, lm))

#..........Model summary stats ====

MorphoDat <- MorphoDat %>%
  mutate(glance = map(Model, glance))

#.....Model plots ====

modelPlots <- lapply(1:nrow(MorphoDat), function(x) {
  ggplot(data = MorphoDat$data[[x]], aes(x = log10(PlanarArea),
                                         y = get(MorphoDat$ShapeVariable[x]),
                                         fill = GrowthForm,
                                         colour = GrowthForm,
                                         label = ScanID)) +
    geom_point(size = 3, alpha = 0.7, shape = 21, colour = "black") +
    #geom_label(colour = "black") +
    geom_smooth(method = "lm", alpha = 0.5) +
    geom_smooth(method = "gam", linetype = "dashed") +
    scale_y_continuous(trans = MorphoDat$Transform[x]) +
    scale_fill_brewer(palette = "Set2") +
    scale_colour_brewer(palette = "Set2") +
    labs(y = MorphoDat$ShapeVariable[x],
         subtitle = paste0("adjR2: ", round(MorphoDat$glance[[x]]$adj.r.squared,3))
         ) +
    facet_wrap(~GrowthForm) +
    theme_minimal()
})

#walk(modelPlots, ~plot(.))

MorphoDat$Plots <- modelPlots

#..........Predicted R2 ====

MorphoDat <- MorphoDat %>%
  mutate(PredR2 = map(.$Model, pred_r_squared))


#.....Export model objects ====

saveRDS(MorphoDat, "3.PredictingShapeVariables/ModelObject.rds")

#.....Model tables ====

Vars <- MorphoDat$ShapeVariable
Coefs <- map(MorphoDat$Model, tidy)

coefsClean <- function(CoefTable, ReferenceVarName = "RefVar", SlopeVar = "SlopeVar") {
  RefInt <- CoefTable[1,]
  RefSlope <- CoefTable[2,]
  RefInt$term <- ReferenceVarName
  RefSlope$term <- paste0(SlopeVar, ":", ReferenceVarName)

  OtherInts <- CoefTable[3:(nrow(CoefTable)/2)+1,]
  OtherSlopes <- CoefTable[((nrow(CoefTable)/2)+2):nrow(CoefTable),]

  OtherInts$estimate <- OtherInts$estimate + RefInt$estimate
  OtherSlopes$estimate <- OtherSlopes$estimate + RefSlope$estimate
  NewCoefTable <- rbind(RefInt,OtherInts, RefSlope,OtherSlopes) %>% as.data.frame() %>% select(-p.value, -statistic)
  return(NewCoefTable)
}


CoefsNew <- map(Coefs, coefsClean, "Arborescent", "log10(PlanarArea)")

CoefsNew <- lapply(1:length(Vars), function(x) {
  Coefs <- CoefsNew[[x]]
  Names <- names(Coefs)
  Names[2] <- paste0(Transforms[[x]] ,"(",Vars[x], ")_estimate")
  Names[3] <- paste0(Transforms[[x]] ,"(", Vars[x], ")_std.error")
  names(Coefs) <- Names
  Coefs
})

lapply(1:length(CoefsNew), function(x) {
  CoefsNew[[x]]$term <<- gsub(pattern = "GrowthForm", replacement = "", CoefsNew[[x]]$term)
  CoefsNew[[x]]$term <<- gsub(pattern = "PlanarArea", replacement = "PlanarArea_cm2", CoefsNew[[x]]$term )
})

CoefsAll <- full_join(CoefsNew[[1]], CoefsNew[[4]], by = "term") %>% full_join(., CoefsNew[[6]], by = "term")

#.....Model stats ====

ModelStats <- lapply(1:length(Vars), function(x) {
  Table <- cbind(MorphoDat$ShapeVariable[[x]], glance(MorphoDat$Model[[x]]), pred_r_squared(MorphoDat$Model[[x]]))
  names(Table)[1] <-"variable"
  names(Table)[13] <-"pred.r.squared"
  Table <- Table %>% select(variable, r.squared, adj.r.squared, pred.r.squared, df, df.residual)
})

StatsAll <- do.call(rbind, ModelStats)

#StatsAll <- StatsAll[c(1,4,6),]

#..........Export stats ====

write.csv(StatsAll, file = "3.PredictingShapeVariables/ModelStats.csv", row.names = F)
#.====
#Predict shape variables on LIT data ====

#.....Remove LIT growth forms missing from morphodat ====

LITFiltered <- filter(LITDat, GrowthForm %in% unique(MorphoDat$data[[1]]$GrowthForm))

#.....Predict shape variables ====

LITDatPredictedShapeVars <- lapply(1:nrow(MorphoDat), function(x) {
  predict.lm(MorphoDat$Model[[x]],
             newdata = dplyr::select(LITFiltered, PlanarArea, GrowthForm))
  })

names(LITDatPredictedShapeVars) <- ShapeVars

MorphoDat <- MorphoDat %>%
  mutate(LITData_PlusShapeVars = map(LITDatPredictedShapeVars, ~cbind(LITFiltered,.)))

#.....update value to shape variable name in nested data ====

lapply(1:length(ShapeVars), function(x)
  names(MorphoDat$LITData_PlusShapeVars[[x]])[ncol(MorphoDat$LITData_PlusShapeVars[[x]])] <<- ShapeVars[x]
  )

#.....Calculate weighted averages ====

MorphoDat$LITData_PlusShapeVars <- lapply(1:nrow(MorphoDat), function(x) {
  MorphoDat$LITData_PlusShapeVars[[x]] %>%
    group_by(campaign,site) %>%
    mutate(ShapeWeightedAverage = (get(ShapeVars[x]) * intercept_cm)/sum(intercept_cm))
})

lapply(1:length(ShapeVars), function(x)
  names(MorphoDat$LITData_PlusShapeVars[[x]])
  [ncol(MorphoDat$LITData_PlusShapeVars[[x]])] <<- paste0(ShapeVars[x],"_Weighted")
)

MorphoDat$ShapeWeightedAverages <- lapply(1:nrow(MorphoDat), function(x) {
  MorphoDat$LITData_PlusShapeVars[[x]] %>%
    group_by(campaign,site) %>%
    summarise(ShapeWeightedAverageSum = sum(get(paste0(ShapeVars[x],"_Weighted")))) %>%
    mutate(ShapeWeightedAverageSum_OrigScale = if(MorphoDat$Transform[x] == "logit") { invLogit(ShapeWeightedAverageSum) } else { 10^(ShapeWeightedAverageSum) } )
})

names(MorphoDat$ShapeWeightedAverages) <- ShapeVars

lapply(1:length(ShapeVars), function(x) {
  names(MorphoDat$ShapeWeightedAverages[[x]])[(ncol(MorphoDat$ShapeWeightedAverages[[x]])-1)] <<- ShapeVars[x]
  names(MorphoDat$ShapeWeightedAverages[[x]])[ncol(MorphoDat$ShapeWeightedAverages[[x]])] <<- paste0(ShapeVars[x],"_OrigScale")
})

#.....Export work ====

saveRDS(MorphoDat, "3.PredictingShapeVariables/PredictiveModels_LITShapeData.rds")

#Convexity data for luisa

# ConvexityData <- MorphoDat[4,]
#
# write.csv(ConvexityData$ShapeWeightedAverages, "ConvexityWA_Data_ForLuisa.csv", row.names = F)
# saveRDS(ConvexityData, "ConvexityWA_RModels_ForLuisa.rds")

