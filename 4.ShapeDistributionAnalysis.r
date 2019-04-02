#Morphological strategies of corals in the anthropocene: shape distribution analysis

rm(list = ls())

#Set up ====

#.....Libraries ====

library(tidyverse)
library(modelr)
library(broom)
library(ggmap)
library(maptools)
library(viridis)
library(grid)
library(gridExtra)
library(png)
library(ggrepel)
library(ggpubr)

#.....Functions ====

source("Functions.r")

#.....Data ====

LITDat <- readRDS("3.PredictingShapeVariables/PredictiveModels_LITShapeData.rds")
BioGeoDat <-  read.csv("2.DataSetup/BioGeoDat_PlusLIBearing.csv")

#Set up data ====

#.....lit data ====

LITDat$LITData_PlusShapeVars[[1]]$campaign <- factor(LITDat$LITData_PlusShapeVars[[1]]$campaign)

#.....biogeo dat ====

#.....remove sites with no LIT data ====

BioGeoDat <- BioGeoDat %>%
  filter(site != "North 2" & site !=  "North 3") %>%
  droplevels

#.....remove corner beach ====

#Corner beach only has data for 2017

BioGeoDat <- filter(BioGeoDat, site != "Corner Beach")

#Weighted averages and shape distributions ====

#.....By site and year ====
#
# LITDat$ShapeHistograms <- lapply(1:nrow(LITDat), function(x) {
#   ggplot(data = LITDat$LITData_PlusShapeVars[[x]], aes_string(x = LITDat$ShapeVariable[x], colour = "factor(campaign)")) +
#     geom_vline(data = LITDat$ShapeWeightedAverages[[x]], aes_string(xintercept =  LITDat$ShapeVariable[x], colour = "factor(campaign)"),
#                linetype = "dashed", size = 1) +
#     #geom_histogram(fill = "white", colour = "black", binwidth = 1) +
#     # geom_segment(data = plotDat %>% group_by(campaign, site) %>% summarise(mean = mean(ln_area_cm2), sd = sd(ln_area_cm2)),
#     #              aes(y = 0.5, yend = 0.5, x = mean-sd, xend = mean + sd),
#     #              colour = "blue", size = 1) +
#     # geom_vline(data = plotDat %>% group_by(campaign, site) %>% summarise(mean = mean(ln_area_cm2), sd = sd(ln_area_cm2)),
#     #            aes(xintercept =  mean - sd),
#     #            linetype = "dashed", colour = "#006677", size = 0.5) +
#     # geom_vline(data = plotDat %>% group_by(campaign, site) %>% summarise(mean = mean(ln_area_cm2), sd = sd(ln_area_cm2)),
#     #            aes(xintercept =  mean + sd),
#     #            linetype = "dashed", colour = "#006677", size = 0.5) +
#     geom_density(fill = "white", alpha = 0) +
#     facet_wrap(~site) +
#     scale_color_brewer(palette = "RdYlBu", direction = -1) +
#     labs(title = paste( LITDat$ShapeVariable[x], " distributions and shifts in weighted average")) +
#     #coord_cartesian(ylim = c(0,1)) +
#     theme_minimal(base_size = 8)
#   })
#


#.....Changes on a map ====

#..........load map ====

Lizard <- readShapePoly("dataInput/LizardIslandMap/Lizard_polylines.shp")

LizardDF <- fortify(Lizard)

#..........plot lizard island plus sites ====

LizardMap <- ggplot(data = BioGeoDat, aes(x = Long, y = Lat, label = site)) +
  geom_path(data = LizardDF, aes(x = long, y = lat, group = group), inherit.aes = F) +
  geom_point(size = 5, shape = 21, alpha = 0.7) +
  geom_label_repel() +
  coord_fixed(xlim = c(145.42, 145.50)) +
  theme_classic()


#ggsave(LizardMap, filename = "4.ShapeDistributionAnalysis/LizardMap.png", height = 6, width = 6, scale = 2)
#ggsave(LizardMap, filename = "4.ShapeDistributionAnalysis/LizardMap.pdf", height = 6, width = 6, scale = 2)

#..........site-wise mutli shape variable plots ====

#...............Scale all shape data to 0,1 ====

LITDat$ShapeWeightedAveragesScaled <- lapply(1:nrow(LITDat), function(x) {
  Dat <- LITDat$ShapeWeightedAverages[[x]]
  Dat[,paste0(LITDat$ShapeVariable[[x]], "_Scaled")] <- rescale(Dat[,LITDat$ShapeVariable[[x]]])
  Dat[,paste0(LITDat$ShapeVariable[[x]], "_Scaled_OrigScale")] <- rescale(Dat[,paste0(LITDat$ShapeVariable[[x]],"_OrigScale")])
  Dat
})

#...............create colony count dataframe ====

ColonyCounts <- LITDat$LITData_PlusShapeVars[[1]] %>%
  select(site, campaign) %>%
  table() %>%
  as.data.frame() %>%
  rename(Counts = Freq)

TransectCounts <- LITDat$LITData_PlusShapeVars[[1]] %>%
  group_by(site, campaign) %>%
  summarise(NTransects = length(unique(transect)))
TransectCounts$NTransectsAdded <- ifelse(TransectCounts$NTransects >= 6, TransectCounts$NTransects, 6)

ColonyCountsDF <- left_join(TransectCounts, ColonyCounts, by = c("site", "campaign")) %>%
  mutate(ShapeVariable = "MeanColoniesPerTransect_Scaled")
ColonyCountsDF$Value <- ColonyCountsDF$Counts/ColonyCountsDF$NTransectsAdded
ColonyCountsDF$Value <- rescale(ColonyCountsDF$Value)
ColonyCountsDF <- ColonyCountsDF %>% select(site, campaign, ShapeVariable, Value)

CoverDF <- LITDat$LITData_PlusShapeVars[[1]] %>%
  select(site, campaign, cover) %>%
  mutate(ShapeVariable = "Cover") %>%
  rename(Value = cover) %>%
  unique()
CoverDF$Value <- rescale(CoverDF$Value)


#...............Plotting dataframe

MultiVarYearPlotDat <- multiJoin(LITDat$ShapeWeightedAveragesScaled) %>%
  select(site, campaign, ColonyVol_Scaled_OrigScale, Convexity_Scaled_OrigScale, FractalDimension_Scaled_OrigScale, ColonySA2ndMomentVertScaled_Scaled_OrigScale) %>%
  #mutate(Convexity_Scaled_OrigScale = 1-Convexity_Scaled_OrigScale) %>%
  gather(ColonyVol_Scaled_OrigScale, Convexity_Scaled_OrigScale, FractalDimension_Scaled_OrigScale, ColonySA2ndMomentVertScaled_Scaled_OrigScale, key = "ShapeVariable", value = "Value") %>%
  mutate(campaign = factor(campaign)) %>%
  bind_rows(CoverDF) %>%
  filter(site %in% c("North 3") == F) %>%
  filter(ShapeVariable != "ColonyVol_Scaled_OrigScale")

#...............List column approach ====

DeltaShapeBySite <- MultiVarYearPlotDat %>%
  group_by(site) %>%
  full_join(.,expand(., campaign)) %>% #fill in missing year data levels to keep x axes consistent between plots
  mutate(campaignInteger = as.integer(as.character(campaign))) %>% #create integer version of campaign
  nest() #nest to list columns

#..........Plotting ====

VariablePalette <- c("#8D0197", #purple
                "#EE0000", #red
                "lightgrey",
                "#026969", #cyan
                "#299700" #green
                )

ItaCol <- "grey3"
NathanCol <- "grey3"
BleachingCol <- "grey3"

deltaShapeByYearPlot <- function(Site, Data, Points = T) {
  ggplot(data = Data, aes(x = campaign, y = Value, fill = ShapeVariable)) +
    geom_point(size = ifelse(Points,3,0), shape = 21, alpha = 0.7, colour = "black")  +
    geom_vline(data = Data,
               aes(xintercept = 4.5),
               size = 1.5,
               alpha = 0.5,
               colour = ItaCol,
               linetype = "solid") +
    geom_vline(data = Data,
               aes(xintercept = 5.5),
               size = 1.5,
               alpha = 0.5,
               colour = NathanCol,
               linetype = "dotted") +
    geom_vline(data = Data,
               aes(xintercept = 6.5),
               size = 1.5,
               alpha = 0.5,
               colour = BleachingCol,
               linetype = "dashed") +
    geom_line(data = Data,
              aes(x = as.integer(campaign), y = Value, group = ShapeVariable),
              size = 1.7,
              colour = "black",
              alpha = 0.5,
              inherit.aes = F) +
    geom_line(data = Data,
              aes(x = as.integer(campaign), y = Value, colour = ShapeVariable),
              size = 1.5,
              alpha = 1,
              inherit.aes = F) +
    # geom_vline(aes(xintercept = which(levels(campaign) %in% '2014')),
    #            linetype = "dashed") +
    #scale_color_manual(values = VariablePalette) +
    #scale_fill_manual(values = VariablePalette) +
    scale_colour_viridis(discrete = T) +
    scale_fill_viridis(discrete = T) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    #scale_x_discrete(breaks = paste(levels(DeltaShapeBySite$data[[1]]$campaign))) +
    theme_classic() +
    guides(fill = F, colour = F) +
    ggtitle(Site) +
    theme(axis.text = element_blank(), axis.title = element_blank(), plot.title = element_text(size = 10, hjust = 0.5)) +
    coord_cartesian(ylim = c(-0.1,1.1))
}

DeltaShapeBySite <- DeltaShapeBySite %>%
  mutate(Plots = map2(DeltaShapeBySite$site, DeltaShapeBySite$data, deltaShapeByYearPlot, F))

names(DeltaShapeBySite$Plots) <- DeltaShapeBySite$site

BlankPlot <- deltaShapeByYearPlot(DeltaShapeBySite$site[[12]], DeltaShapeBySite$data[[12]])

BlankPlot <- BlankPlot +
  guides(fill = guide_legend(keyheight = 0.7)) +
  theme(axis.text = element_text(),
        axis.title = element_text(),
        plot.title = element_blank(),
        legend.position= c (0.5, 0.5)) +
  #scale_fill_manual(values = VariablePalette, labels = c("Volume", "Fractal Dimension", "Colonies", "SA:Vol", "Sphericity")) +
  #scale_fill_manual(values = VariablePalette, labels = c("Size", "Surface complexity", "Colonies", "SA:Vol", "Volume compactness")) +
  #scale_fill_manual(values = VariablePalette, labels = c("Top heaviness", "Colony Size", "Coral cover", "Volume compactness", "Surface complexity")) +
  scale_fill_viridis(discrete = T, labels = c("Top heaviness",  "Volume compactness", "Coral cover", "Surface complexity")) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "Year",
       y = "Scaled \nvariable value",
       fill = "") +
  ggtitle("Key") +
  theme_classic2(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(size = 10, hjust = 0.5), legend.background = element_blank())


Legend <- get_legend(BlankPlot)

BlankPlot$layers[[1]]$aes_params$alpha <- 0
BlankPlot$layers[[5]] <- NULL
BlankPlot$layers[[5]] <- NULL
# BlankPlot$layers[[2]] <- NULL
# BlankPlot$layers[[2]] <- NULL
# BlankPlot$layers[[2]] <- NULL

BlankPlot <- BlankPlot +
  guides(fill = F) +
  annotation_custom(Legend, xmin = 2, xmax = 5, ymin = 0.6, ymax = 0.8)

DeltaShapeBySite$Plots$`North Reef`

#Grob

DeltaShapeBySite <- DeltaShapeBySite %>%
  mutate(Grobs = map(DeltaShapeBySite$Plots, ggplotGrob))


#...............create left, top and right hand side plot grobs ====

ImageXCoords <- BioGeoDat$Long# + LizardSitesXYShifts$XShift
ImageYCoords <- BioGeoDat$Lat #+ LizardSitesXYShifts$YShift

LizardMap <- ggplot(data = BioGeoDat, aes(x = Long, y = Lat)) +
  geom_path(data = LizardDF, aes(x = long, y = lat, group = group), inherit.aes = F) +
  geom_point(size = 3, alpha = 0.7, colour = "grey3") +
  #lapply(1:length(PNGsLoaded), function(p) annotation_custom(PNGsLoaded[[p]],PNGXMin[p],PNGXMax[p],PNGYMin[p],PNGYMax[p])) +
  geom_label_repel(aes(x = ImageXCoords, y = ImageYCoords, label = BioGeoDat$site)) +
  geom_segment(aes(x = 145.485, xend = 145.47, y = -14.645, yend = -14.65),
               size = 5,
               colour = ItaCol,
               alpha = 0.2,
               arrow = arrow(length = unit(0.5,"cm"), type = "closed")) +
  geom_segment(aes(x = 145.485, xend = 145.47, y = -14.70, yend = -14.695),
               size = 5,
               colour = NathanCol,
               alpha = 0.2,
               arrow = arrow(length = unit(0.5,"cm"), type = "closed")) +
  coord_fixed(xlim = c(145.425, 145.4875), ylim = c(-14.64, -14.71)) +
  theme_classic() +
  theme(axis.title = element_blank(), axis.text = element_blank())

# LizardMap <- LizardMap +
#   annotation_custom(ggplotGrob(BlankPlot), xmin = 145.425, xmax = 145.4425, ymin = -14.7125, ymax = -14.6975)

LizardMap <- ggplotGrob(LizardMap)

#..........Big map plot with panels ====

BigMap <-  grid.arrange(arrangeGrob(DeltaShapeBySite$Grobs$`Turtle Beach`,
                                    DeltaShapeBySite$Grobs$`Cooks Path`, #left column
                                    DeltaShapeBySite$Grobs$Resort,
                                    ggplotGrob(BlankPlot), #rectGrob(gp = gpar(col = "white", fill = "white")),
                                    DeltaShapeBySite$Grobs$Osprey,
                                    DeltaShapeBySite$Grobs$Vickies,
                                    DeltaShapeBySite$Grobs$Horseshoe,
                                    nrow = 7),
                        arrangeGrob(arrangeGrob(DeltaShapeBySite$Grobs$`Mermaid Cove`, #top row
                                                DeltaShapeBySite$Grobs$`North Reef`,
                                                DeltaShapeBySite$Grobs$`Washing Machine`,
                                                ncol = 3),
                        LizardMap,
                        arrangeGrob(DeltaShapeBySite$Grobs$`Lagoon 2`,             #bottom row
                                    DeltaShapeBySite$Grobs$Trimodal,
                                    DeltaShapeBySite$Grobs$`Lagoon 1`,
                                    ncol = 3), #Map
                        nrow = 3,
                        layout_matrix = cbind(c(1,2,2,2,2,2,3))),
                        arrangeGrob(DeltaShapeBySite$Grobs$`Gnarly Tree`,
                                    DeltaShapeBySite$Grobs$`North of Paradise`,  #right column
                                    DeltaShapeBySite$Grobs$`No Mans Land`,
                                    DeltaShapeBySite$Grobs$`Easter Point`,
                                    DeltaShapeBySite$Grobs$`Lizard Head`,
                                    DeltaShapeBySite$Grobs$Southeast,
                                    DeltaShapeBySite$Grobs$`South Island`,
                                    nrow = 7),
                        ncol = 3,
                        layout_matrix = rbind(c(1,2,2,2,3))
                        )

#ggsave(BigMap, filename = "4.ShapeDistributionAnalysis/BigMap_TakeTwo_PlusVsa.png", height = 6, width = 6, scale = 2)
#ggsave(BigMap, filename = "4.ShapeDistributionAnalysis/BigMap_takeTwo_PlusVsa.pdf", height = 6, width = 6, scale = 2)

