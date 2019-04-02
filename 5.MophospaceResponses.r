#Morphological strategies of corals in the anthropocene: shape distribution analysis

rm(list = ls())

#Set up ====

#.....Libraries ====

library(plyr)
library(tidyverse)
library(viridis)
library(ggpubr)
library(grid)
library(gridExtra)
library(hypervolume)
library(factoextra)

#.....Functions ====

source("Functions.r")

#.....Data ====

LITDat <- readRDS("3.PredictingShapeVariables/PredictiveModels_LITShapeData.rds")
BioGeoDat <-  read.csv("2.DataSetup/BioGeoDat_PlusLIBearing.csv")
MorphoDat <- read.csv("2.DataSetup/MorphologyData_HighQualityWhole5GFs.csv")

#PCA Analysis ====

#.....create LIT PCA Dataframe ====

#Only uses sites with conisitent data from 2011 to 2017

LITPCADatList <- lapply(1:nrow(LITDat), function(x) {
  LITDat$LITData_PlusShapeVars[[x]] %>%
    ungroup() %>%
    select(-campaign, -observer, -site, -transect, -Species,
           -intercept_cm, -GrowthForm, -PlanarArea, -BenthicType)
})

LITPCADat <- bind_cols(LITPCADatList) %>%
  bind_cols(select(LITDat$LITData_PlusShapeVars[[1]], campaign, observer, site, transect, Species,
                   intercept_cm, GrowthForm, PlanarArea, BenthicType), .) %>%
  filter(site %in% c("North Reef", "Lizard Head", "Horseshoe", "Lagoon2", "Trimodal", "Lagoon 1")) %>%
  filter(campaign %in% c("2011","2014","2015","2016","2017"))

#.....transform morphology data ====

MorphoDat_OrigScale <- MorphoDat

MorphoDat <- MorphoDat %>% mutate(Sphericity = logit(Sphericity),
                                  FractalDimension = logit(FractalDimension),
                                  ColonyVol = log10(ColonyVol),
                                  SAVol = log10(SAVol),
                                  ColonySA2ndMomentVertScaled = log10(ColonySA2ndMomentVertScaled),
                                  Convexity = logit(Convexity))


#.....PCA Set up ====

PCAVars <- c("Sphericity", "FractalDimension", "SAVol")
PCAVarsLabels <- c("S", "FD", "SA:V", "Vol")

PCAVars <- c("Convexity", "FractalDimension", "ColonySA2ndMomentVertScaled")
PCAVarsLabels <- c("Volume\ncompactness", "Surface\ncomplexity", "Top-heaviness\nsurface area")

generatePCA <- function(Data, Vars, VarLabels = NA, PCALabel = "PCA") {
  #if(is.na(VarLabels)) VarLabels <- Vars
  PCA <- prcomp(Data[,Vars], center = T, scale. = T)
  Loadings <- as.data.frame(cor(Data[,Vars], PCA$x))
  VariableArrows <- Arrows <- data.frame(xS = rep(0,nrow(Loadings)), #loading arrow x center
                                         yS = rep(0,nrow(Loadings)), #loading arrow y center
                                         xE = Loadings$PC1, #loading arrow x coordinate
                                         yE = Loadings$PC2) #loading arrow y coordinate
  VarLabels <- data.frame(Labels = VarLabels,
                          X = VariableArrows$xE *1.1, # x coordinate for variable labels, scaled away from arrow end
                          Y = VariableArrows$yE *1.1) # y coordinate for variable lables, scale away from arrow end
  PCADataframe <- data.frame(PCA$x,
                             Group = PCALabel)
  GeneratedPCA <- list(PCALabel = PCALabel, Data = Data[,Vars], PCA = PCA, PCADataframe =  PCADataframe,
                       Loadings = Loadings, VariableArrows = VariableArrows,VarLabels = VarLabels)
  return(GeneratedPCA)
}

predictPCA <- function(Data, Vars,  PCA, PCALabel = "PCA") {
  PCADataframe <- data.frame(predict(PCA, newdata = Data[,Vars]), Group = PCALabel)
  return(PCADataframe)
}

ArrowLabelScaler <- 5

plotPCA <- function(GeneratedPCA, GroupingVar = F, VarLabels = T,
                    PC1Shift = 0, PC2Shift = 0, Points = F, Bary = F, Spread = T) {
  ggplot(data = GeneratedPCA$PredictedPCAData, aes(x = PC1 + PC1Shift, y = PC2 + PC2Shift, colour = Group, group = Group)) + #ggplot set up
    #PlotTheme + #plotting theme
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey31") + #horizontal reference line
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey31") + #verticalreference line
    {if(Points) geom_jitter(size = 5, alpha = 0.05, show.legend = T, width = 0.25, height = 0.25)} + #Scatter plot of observed PC1 and 2 values
    geom_segment(data = GeneratedPCA$VariableArrows, #Outwards loading arrows
                 aes(x = xS, y = yS, xend = xE * ArrowLabelScaler, yend = yE * ArrowLabelScaler),
                 colour = "grey21",
                 size = 0.85,
                 arrow = arrow(length = unit(0.25,"cm"),
                               type = "closed"),
                 inherit.aes = F) +
    geom_segment(data = GeneratedPCA$VariableArrows, #Inwards loading arrows
                 aes(x = xS, y = yS, xend = -xE * ArrowLabelScaler, yend = -yE * ArrowLabelScaler),
                 linetype = "dashed",
                 colour = "grey21",
                 size = 0.85,
                 inherit.aes = F) +
   {if(Spread) stat_conf_ellipse(alpha = 0, #95% confidence ellipse around the mean
                                 linetype = "solid",
                                 geom = "polygon",
                                 colour = "black",
                                 size = 1.5,
                                 bary = F,
                                 level= 0.95)} +
   {if(Spread)stat_conf_ellipse(alpha = 0, #95% confidence ellipse around the mean
                                linetype = "solid",
                                geom = "polygon",
                                size = 1.25,
                                bary = F,
                                level= 0.95)} +
    {if(Bary)stat_conf_ellipse(alpha = 0, #95% confidence ellipse around the mean
                               linetype = "solid",
                               geom = "polygon",
                               colour = "black",
                               size = 1.5,
                               bary = T,
                               level= 0.95)} +
   {if(Bary)stat_conf_ellipse(alpha = 0, #95% confidence ellipse around the mean
                              linetype = "solid",
                              geom = "polygon",
                              size = 1.25,
                              bary = T,
                              level= 0.95)} +
    {if(VarLabels) geom_text(data = GeneratedPCA$VarLabels, #Loading variable text labels
                           aes(x =  X * ArrowLabelScaler, y = Y * ArrowLabelScaler, label = Labels),
                           size = 5,
                           colour = "black" ,
                           inherit.aes = F,
                           parse = F)} +
    #geom_text(aes(label = Labels), colour = "black") +
    scale_fill_viridis(discrete = T, begin = 0.2, end = 0.95, drop= F,  option = "magma") + #set pallete colours
    scale_color_viridis(discrete = T, begin = 0.2, end = 0.95, drop = F, option = "magma") + #set pallete colours
    # theme(axis.ticks = element_blank(), #remove axis ticks
    #       panel.grid.major = element_blank()) + #remove panel grid lines
    #coord_fixed(xlim = c(-11,11), ylim = c(-11,9)) + #setplot x and y limits
    guides(fill = F, colour = F) +
    theme_classic() +
    theme(axis.text = element_blank(), axis.title = element_blank(), plot.title = element_text(size = 10, hjust = 0.5))
}

#.....Generate PCA ====

MorphoPCA <- generatePCA(MorphoDat, PCAVars, PCALabel = "PossibleSpace", VarLabels = PCAVarsLabels)
LITPCA <- generatePCA(LITPCADat, PCAVars, PCALabel = "PossibleSpace", VarLabels = PCAVarsLabels)


#.....All sites PCA ====

AllSitesPCAData <- MorphoPCA
AllSitesPCAData <- LITPCA

AllSitesPCAData <- LITPCADat %>%
  group_by(campaign) %>%
  nest()

AllSitesPCAData <- AllSitesPCAData %>%
  mutate(PredPCA = map2(data, campaign, predictPCA, Vars = PCAVars, PCA = LITPCA$PCA)) %>% #Predict PCA positions for each dataset from MorphoPCA
  unnest() %>% #unnest the list columns
  full_join(.,expand(., campaign)) %>% #add missing site~campaign combinations
  mutate(Group = campaign) %>%
  select(PC1, PC2, Group) %>%
  mutate(Group = factor(Group))

AllSitesPCA <- LITPCA
AllSitesPCA$PredictedPCAData <- AllSitesPCAData


AllSitesPCAPlot <- plotPCA(AllSitesPCA,T, PC1Shift = 0.65, PC2Shift = -0.25, Bary = T, Spread = T, Points = T)

(AllSitesPCAPlot <- AllSitesPCAPlot +
  # geom_point(data = MorphoPCA$PCADataframe, aes(x =  PC1, y = PC2),
  #            shape = 21,
  #            colour = "grey",
  #            fill = "white",
  #            size = 2) +
  guides(colour = guide_legend(keyheight = 0.7,
                               title.position = "top",
                               label.position = "right",
                               direction = "vertical",
                               title.hjust = 0.5,
                               ncol = 1)
         ) +
  scale_colour_viridis(discrete = T, begin = 0.2, end = 0.95, drop = F, option = "magma",
                       labels = c("2011 - Pre Ita",
                                  "2014 - Post Ita | Pre Nathan",
                                  "2015 - Post Nathan | Pre Bleaching",
                                  "2016 - Post Bleaching",
                                  "2017 - 1 Year Post Bleaching")) +
  coord_fixed(xlim = c(-11,11),ylim = c(-7,7)) +
  theme_classic2(base_size = 12) +
  theme(axis.text = element_text(),
        axis.title = element_text(),
        plot.title = element_blank(),
        legend.position = "bottom") +
  labs(colour = "Year", x = "PC1 (65%)", y = "PC2 (30%)"))

# AllSitesPCAPlot$layers <- c(AllSitesPCAPlot$layers[[11]],AllSitesPCAPlot$layers)
# AllSitesPCAPlot$layers[[12]] <- NULL

#ggsave(AllSitesPCAPlot, filename = "5.MorphospaceResponses/ChangesInMorphospace_AllSites_ShapeVars_Points.png", height = 6, width = 6, scale = 1.5)
#ggsave(AllSitesPCAPlot, filename = "5.MorphospaceResponses/ChangesInMorphospace_AllSites_ShapeVars_Points.pdf", height = 6, width = 6, scale = 1.5)


#.....Quantifying changes in PCAs ====

fviz_contrib(LITPCA$PCA, "var", axes = 1)
fviz_contrib(LITPCA$PCA, "var", axes = 2)
fviz_screeplot(LITPCA$PCA)

#..........Dispersion and PERMANOVA of multivariate spaces ====

library(vegan)

DispersionDat <- cbind(LITPCA$Data, LITPCADat$campaign) %>%
  dplyr::rename(campaign = 'LITPCADat$campaign') %>%
  filter(campaign %in% c("2011","2017"))

AllSitesDist <- DispersionDat %>% dplyr::select(-campaign) %>% dist()

AllSites_Dispersion <- betadisper(AllSitesDist, DispersionDat$campaign)
anova(AllSites_Dispersion)
summary(AllSites_Dispersion)
boxplot(AllSites_Dispersion)

AllSites_PERMANOVA <- adonis(AllSitesDist ~ DispersionDat$campaign)
anova(AllSites_PERMANOVA)
AllSites_PERMANOVA


summary(AllSites_Dispersion)
TukeyHSD(AllSites_Dispersion)

summary(AllSites_Similarity)
TukeyHSD(AllSites_Similarity)

summary(AllSites_PERMANOVA)
TukeyHSD(AllSites_PERMANOVA)

