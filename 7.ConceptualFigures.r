#Conceptual figures plots

#Morphological strategies of corals in the anthropocene: data setup

rm(list = ls())

#Set up ====

#.....Libraries ====

library(tidyverse)
library(ggpubr)
library(viridis)

#.....Functions ====

source("Functions.r")

summariseData <- function(Data, Var) {
  SummDat <- Data %>%
    select(GrowthForm, Value) %>%
    group_by(GrowthForm) %>%
    summarise(Mean = mean(Value),
              Median = median(Value),
              SD = sd(Value))
  names(SummDat)[2:4] <- paste0(names(SummDat)[2:4],"_",Var)
  SummDat$GrowthForm <- factor(SummDat$GrowthForm, levels = SummDat$GrowthForm[
    c(order(as.data.frame(
      SummDat[,paste0("Mean_",Var)])))
    ]
    )
  return(SummDat)
}

reorderVar <- function(MainData, SummaryData) {
  MainData$GrowthForm <- factor(MainData$GrowthForm, levels = SummaryData$GrowthForm[
    c(order(as.data.frame(SummaryData[,2])))
    ]
  )
  return(MainData)
}

densPlotObs <- function(Data, Var) {
  Min <- min(Data$Value)
  Max <- max(Data$Value)
  ggplot(data = Data, aes(x = Value, colour = GrowthForm)) +
    geom_line(stat = "density", size = 2, alpha = 0.7) +
    scale_color_viridis(discrete = T) +
    theme_classic2() +
    theme(panel.background = element_rect(colour = "black")) +
    coord_cartesian(xlim = c(Min,Max), expand = T) +
    ggtitle(paste0(Var," distribution observed")) +
    labs(x = Var)
}

resampleMeanSD <- function(SummaryData) {
  DataResampled <- lapply(1:nrow(SummaryData), function(x) {
    Mean <- as.numeric(SummaryData[x,2])
    SD <- as.numeric(SummaryData[x,4])
    data.frame(GrowthForm = SummaryData$GrowthForm[x],
               Value = rnorm(1000000, Mean, SD)
               )
  }) %>%
  bind_rows()
  return(DataResampled)
}

densPlotResampled <- function(Data, Var) {
  Min <- min(Data$Value)
  Max <- max(Data$Value)
  ggplot(data = Data, aes(x = Value, colour = GrowthForm)) +
    geom_line(stat = "density", size = 2, alpha = 0.7) +
    scale_color_viridis(discrete = T) +
    theme_classic2() +
    theme(panel.background = element_rect(colour = "black")) +
    coord_cartesian(xlim = c(Min,Max), expand = T) +
    ggtitle(paste0(Var," distribution resampled mean/sd")) +
    labs(x = Var)
}

#.....Data ====

MorphoDat <- read.csv("2.DataSetup/MorphologyData_HighQualityWhole5GFs.csv")
Vars <- c("Sphericity", "Convexity", "FractalDimension", "Packing",
          "ColonyVol2ndMomentVertScaled", "ColonySA2ndMomentVertScaled")

MorphoDat <- MorphoDat%>%
  filter(GrowthForm %in% "BranchingClosed" == F)


#Shape variable distribution ====

MorphoDatNested <- MorphoDat %>%
  select(GrowthForm, Vars) %>%
  rename(SecondMomentVolume = ColonyVol2ndMomentVertScaled, SecondMomentArea = ColonySA2ndMomentVertScaled) %>%
  gather(key = "Trait", value = "Value", -GrowthForm) %>%
  group_by(Trait) %>%
  nest(GrowthForm, Value)

MorphoDatNested$DataSummary <- map2(MorphoDatNested$data, MorphoDatNested$Trait, summariseData)

MorphoDatNested$data <- map2(MorphoDatNested$data, MorphoDatNested$DataSummary, reorderVar)

MorphoDatNested$DataResampled <- map(MorphoDatNested$DataSummary, resampleMeanSD)

MorphoDatNested$DensPlotObs <- map2(MorphoDatNested$data, MorphoDatNested$Trait, densPlotObs)

MorphoDatNested$DensPlotResampled <- map2(MorphoDatNested$DataResampled, MorphoDatNested$Trait, densPlotResampled)

MorphoDatSummary <- MorphoDatNested %>%
  select(DataSummary) %>%
  do.call(rbind, .) %>%
  do.call(cbind, .) %>%
  .[,!duplicated(colnames(.), fromLast = T)] %>%
  select(GrowthForm, everything())

#write.csv(MorphoDatSummary, file = "MorphologyByGrowthForm_ForDama.csv", row.names = F)

#.....Convexity ====

ConDat <- MorphoDat %>%
  select(GrowthForm, Convexity)

ConDatMean <- ConDat %>%
  group_by(GrowthForm) %>%
  summarise(ConvexityMean = mean(Convexity),
            ConvexityMedian = median(Convexity),
            ConvexitySD = sd(Convexity))

ConDat2 <- lapply(1:nrow(ConDatMean), function(x) {
  data.frame(GrowthForm = ConDatMean$GrowthForm[x],
             Convexity = rnorm(1000000, ConDatMean$ConvexityMedian[x], ConDatMean$ConvexitySD[x]))
}) %>%
  bind_rows()

ConDat$GrowthForm <- factor(ConDat$GrowthForm , levels = ConDatMean$GrowthForm[c(order(ConDatMean$ConvexityMean))])
ConDat2$GrowthForm <- factor(ConDat2$GrowthForm , levels = ConDatMean$GrowthForm[c(order(ConDatMean$ConvexityMean))])

ConDens <- ggplot(data = ConDat, aes(x = Convexity, colour = GrowthForm)) +
  geom_line(stat = "density", size = 2, alpha = 0.7) +
  scale_color_viridis(discrete = T) +
  theme_classic2() +
  theme(panel.background = element_rect(colour = "black")) +
  coord_cartesian(xlim = c(0,1), ylim =c(0,12), expand = F) +
  ggtitle("Convexity distribution observed")


# ggsave(ConDens, file = "7.ConceptualFigures/ConvexityDist_GrowthForm_Observed.png", width = 6, height = 3, scale = 1)
# ggsave(ConDens, file = "7.ConceptualFigures/ConvexityDist_GrowthForm_Observed.pdf", width = 6, height = 3, scale = 1)


SampleConDens <- ggplot(data = ConDat2, aes(x = Convexity, colour = GrowthForm)) +
  geom_line(stat = "density", size = 2, alpha = 0.7) +
  scale_color_viridis(discrete = T) +
  theme_classic2() +
  theme(panel.background = element_rect(colour = "black")) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,11), expand = F)# +
  #ggtitle("Convexity distribution, resampled median and sd, n = 1000000")

# ggsave(SampleConDens, file = "7.ConceptualFigures/ConvexityDist_GrowthForm_Sampled.png", width = 6, height = 3, scale = 1)
# ggsave(SampleConDens, file = "7.ConceptualFigures/ConvexityDist_GrowthForm_Sampled.pdf", width = 6, height = 3, scale = 1)

#.....FracDim ====

FDDat <- MorphoDat %>%
  select(GrowthForm, FractalDimension) %>%
  mutate(FractalDimension = FractalDimension + 2)

FDDatMean <- FDDat %>%
  group_by(GrowthForm) %>%
  summarise(FractalDimensionMean = mean(FractalDimension),
            FractalDimensionMedian = median(FractalDimension),
            FractalDimensionSD = sd(FractalDimension))

FDDat2 <- lapply(1:nrow(FDDatMean), function(x) {
  data.frame(GrowthForm = FDDatMean$GrowthForm[x],
             FractalDimension = rnorm(1000000, FDDatMean$FractalDimensionMean[x], FDDatMean$FractalDimensionSD[x]))
}) %>%
  bind_rows()

FDDat$GrowthForm <- factor(FDDat$GrowthForm , levels = FDDatMean$GrowthForm[c(order(FDDatMean$FractalDimensionMean))])
FDDat2$GrowthForm <- factor(FDDat2$GrowthForm , levels = FDDatMean$GrowthForm[c(order(FDDatMean$FractalDimensionMean))])

#FDDat2 <- sapply(   FractalDimension = map2(FractalDimensionMean, FractalDimensionSD, rnorm, n = 100))

ggplot(data = FDDat, aes(x = FractalDimension, colour = GrowthForm)) +
  geom_line(stat = "density", size = 2, alpha = 0.7) +
  scale_color_viridis(discrete = T) +
  theme_classic2() +
  theme(panel.background = element_rect(colour = "black")) +
  coord_cartesian(xlim = c(2,3))

SampleFDDens <- ggplot(data = FDDat2, aes(x = FractalDimension, colour = GrowthForm)) +
  geom_line(stat = "density", size = 2, alpha = 0.7) +
  scale_color_viridis(discrete = T) +
  theme_classic2() +
  theme(panel.background = element_rect(colour = "black")) +
  coord_cartesian(xlim = c(2,3), ylim = c(0,7.5), expand = F)

# ggsave(SampleFDDens, file = "7.ConceptualFigures/FractalDimensionDist_GrowthForm_Sampled.png", width = 6, height = 3, scale = 1)
# ggsave(SampleFDDens, file = "7.ConceptualFigures/FractalDimensionDist_GrowthForm_Sampled.pdf", width = 6, height = 3, scale = 1)

#.....Vvol ====

VvolDat <- MorphoDat %>%
  select(GrowthForm, ColonyVol2ndMomentVertScaled)

VvolDatMean <- VvolDat %>%
  group_by(GrowthForm) %>%
  summarise(ColonyVol2ndMomentVertScaledMean = mean(ColonyVol2ndMomentVertScaled),
            ColonyVol2ndMomentVertScaledMedian = median(ColonyVol2ndMomentVertScaled),
            ColonyVol2ndMomentVertScaledSD = sd(ColonyVol2ndMomentVertScaled))

VvolDat2 <- lapply(1:nrow(VvolDatMean), function(x) {
  data.frame(GrowthForm = VvolDatMean$GrowthForm[x],
             ColonyVol2ndMomentVertScaled = rnorm(1000000, VvolDatMean$ColonyVol2ndMomentVertScaledMean[x], VvolDatMean$ColonyVol2ndMomentVertScaledSD[x]))
}) %>%
  bind_rows()

VvolDat$GrowthForm <- factor(VvolDat$GrowthForm , levels = VvolDatMean$GrowthForm[c(order(VvolDatMean$ColonyVol2ndMomentVertScaledMean))])
VvolDat2$GrowthForm <- factor(VvolDat2$GrowthForm , levels = VvolDatMean$GrowthForm[c(order(VvolDatMean$ColonyVol2ndMomentVertScaledMean))])

#VvolDat2 <- sapply(   ColonyVol2ndMomentVertScaled = map2(ColonyVol2ndMomentVertScaledMean, ColonyVol2ndMomentVertScaledSD, rnorm, n = 100))

ggplot(data = VvolDat, aes(x = ColonyVol2ndMomentVertScaled, colour = GrowthForm)) +
  geom_line(stat = "density", size = 2, alpha = 0.7) +
  scale_color_viridis(discrete = T) +
  theme_classic2() +
  theme(panel.background = element_rect(colour = "black")) #+
  #coord_cartesian(xlim = c(0,1))

SampleVvolDens <- ggplot(data = VvolDat2, aes(x = ColonyVol2ndMomentVertScaled, colour = GrowthForm)) +
  geom_line(stat = "density", size = 2, alpha = 0.7) +
  scale_color_viridis(discrete = T) +
  theme_classic2() +
  theme(panel.background = element_rect(colour = "black")) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,7.5), expand = F)

# ggsave(SampleVvolDens, file = "7.ConceptualFigures/ColonyVol2ndMomentVertScaledDist_GrowthForm_Sampled.png", width = 6, height = 3, scale = 1)
# ggsave(SampleVvolDens, file = "7.ConceptualFigures/ColonyVol2ndMomentVertScaledDist_GrowthForm_Sampled.pdf", width = 6, height = 3, scale = 1)

#.....Vsa ====

VsaDat <- MorphoDat %>%
  select(GrowthForm, ColonyVol2ndMomentVertScaled)

VsaDatMean <- VsaDat %>%
  group_by(GrowthForm) %>%
  summarise(ColonyVol2ndMomentVertScaledMean = mean(ColonyVol2ndMomentVertScaled),
            ColonyVol2ndMomentVertScaledSD = sd(ColonyVol2ndMomentVertScaled))

VsaDat2 <- lapply(1:nrow(VsaDatMean), function(x) {
  data.frame(GrowthForm = VsaDatMean$GrowthForm[x],
             ColonyVol2ndMomentVertScaled = rnorm(1000000, VsaDatMean$ColonyVol2ndMomentVertScaledMean[x], VsaDatMean$ColonyVol2ndMomentVertScaledSD[x]))
}) %>%
  bind_rows()

VsaDat$GrowthForm <- factor(VsaDat$GrowthForm , levels = VsaDatMean$GrowthForm[c(order(VsaDatMean$ColonyVol2ndMomentVertScaledMean))])
VsaDat2$GrowthForm <- factor(VsaDat2$GrowthForm , levels = VsaDatMean$GrowthForm[c(order(VsaDatMean$ColonyVol2ndMomentVertScaledMean))])

#VsaDat2 <- sapply(   ColonyVol2ndMomentVertScaled = map2(ColonyVol2ndMomentVertScaledMean, ColonyVol2ndMomentVertScaledSD, rnorm, n = 100))

ggplot(data = VsaDat, aes(x = ColonyVol2ndMomentVertScaled, colour = GrowthForm)) +
  geom_line(stat = "density", size = 2, alpha = 0.7) +
  scale_color_viridis(discrete = T) +
  theme_classic2() +
  theme(panel.background = element_rect(colour = "black")) #+
#coord_cartesian(xlim = c(0,1))

SampleVsaDens <- ggplot(data = VsaDat2, aes(x = ColonyVol2ndMomentVertScaled, colour = GrowthForm)) +
  geom_line(stat = "density", size = 2, alpha = 0.7) +
  scale_color_viridis(discrete = T) +
  theme_classic2() +
  theme(panel.background = element_rect(colour = "black")) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,7.5), expand = F)

# ggsave(SampleVsaDens, file = "7.ConceptualFigures/ColonyVol2ndMomentVertScaledDist_GrowthForm_Sampled.png", width = 6, height = 3, scale = 1)
# ggsave(SampleVsaDens, file = "7.ConceptualFigures/ColonyVol2ndMomentVertScaledDist_GrowthForm_Sampled.pdf", width = 6, height = 3, scale = 1)


#.....Added time multi panel! ====

library(gridExtra)

SampleConDens <- SampleConDens + labs(colour = "Growth form")
SampleFDDens <- SampleFDDens + labs(x = "Fractal dimension", colour = "Growth form")
SampleVvolDens <- SampleVvolDens + labs(x = "2nd moment of volume", colour = "Growth form")


MultiPan <- grid.arrange(SampleConDens, SampleFDDens, SampleVvolDens, ncol = 1)

ggsave(MultiPan, file = "7.ConceptualFigures/ThreeVars_Dists_Sampled.png", width = 6, height = 9, scale = 1)
ggsave(MultiPan, file = "7.ConceptualFigures/ThreeVars_Dists_Sampled.pdf", width = 6, height = 9, scale = 1)

