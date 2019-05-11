#Morphological strategies of corals in the anthropocene: data setup

rm(list = ls())

#Set up ====

#.....Libraries ====

library(tidyverse)

#.....Functions ====

source("Functions.r")

#.....Data ====

#LITDat <- read.csv("dataInput/lit_merged_20180412.csv")
LITDat <- read.csv("dataInput/lit_lizard_meps.csv")

CTDB <- read.csv("dataInput/ctdb_1.1.1.1_NomenclatureUpdate_data.csv") #DOI Here:
MorphoDat <- read.csv("dataInput/Morphology_WholeAndMedtoHighQualityMeshes_Ready.csv") #DOI here:
BioGeoDat <-  read.csv("dataInput/Site_Biogeography_UpdatedSiteNames.csv") #Available from GitHub page


#Process data ====

#.....Trait data ====

#Extract typical growth forms for species in LIT data and convert to factor

TraitGrowthForms <- filter(CTDB, location_name %in% "Global estimate") %>%
  filter(., trait_name %in% "Growth form typical") %>%
  filter(., specie_name %in% LITDat$species) %>%
  mutate(value = factor(as.character(value)))


#check levels and set to match morphology data naming scheme

levels(TraitGrowthForms$value)

levels(TraitGrowthForms$value) <- c("BranchingClosed",
                                    "Arborescent",
                                    "Columnar",
                                    "Corymbose",
                                    "Digitate",
                                    "Encrusting",
                                    "EncrustingLongUprights",
                                    "Hispidose",
                                    "Laminar",
                                    "Massive",
                                    "Submassive",
                                    "Tabular")

names(TraitGrowthForms)[6] <- "species"

#.....LIT data ====

#Add growth forms to LIT data

LITDat <- full_join(LITDat, select(TraitGrowthForms,"species","value"), "species")

#Check missing species

filter(LITDat, is.na(value)) %>% .$species %>% unique(.) %>%  sort(.)

#Update species names

LITDat[which(LITDat$Species == "Montastrea curta")] <- "Montastraea curta"

#Manually add growth form where appropriate

source("2a.DataSetup_ManualUpdateGrowthForms.r")

#Create planar area from intercept, assuming intercept = diameter, multiply by 100 for planar area in mm

LITDat <- mutate(LITDat, PlanarArea = (pi * (intercept_cm/2)^2)*100)

#rename LIT variables to match morphology data

LITDat <- rename(LITDat, Species = species, GrowthForm = value)

#Add benthic type data

OtherBenthicTypes <- read.csv("dataInput/AdditionalBenthicTypes.csv")

GrowthFormBenthicTypes <- LITDat %>% filter(is.na(GrowthForm) == F) %>% select(Species) %>% unique(.) %>% mutate(BenthicType = "HardCoral")

BenthicTypesTable <- bind_rows(OtherBenthicTypes, GrowthFormBenthicTypes)

LITDat$BenthicType <- sapply(1:nrow(LITDat), function(x) BenthicTypesTable$BenthicType[which(BenthicTypesTable$Species %in% LITDat$Species[x])]) #%>% unlist()


sum(is.na(LITDat$GrowthForm))/nrow(LITDat) #proportion of data with no growth form

#check transects

TransectCheck <- LITDat %>%
  group_by(site, campaign) %>%
  summarise(NTransect = length(unique(transect)))

TransectCheck2 <- LITDat %>%
  group_by(site, campaign, transect) %>%
  summarise(LengthTransect = sum(intercept_cm))

#calculate cover by site and year

LITDat <- LITDat %>%
  group_by(site, campaign) %>%
  mutate(ntrans = length(unique(transect))) %>%
  mutate(ntransadded = ifelse(ntrans >= 6, ntrans, 6)) %>%
  filter(BenthicType == "HardCoral") %>%
  mutate(length = sum(intercept_cm),
         cover = length/(ntransadded*1000))

#Export full lit data

write.csv(LITDat, "2.DataSetup/LIT_All_PlusGrowthFormBenthicType.csv", row.names = F)

#Remove non-coral and no growth form

LITDat <- LITDat %>% filter(is.na(GrowthForm) == F)

#Export setup data

write.csv(LITDat, "2.DataSetup/LIT_CoralsWithGrowthForm.csv", row.names = F)


#.....Morphology data ====

#create SAVol variable

MorphoDat <- mutate(MorphoDat, SAVol = OuterSA/ColonyVol)

#minus 2 from FracDim for future logit transform

MorphoDat$FractalDimension <- MorphoDat$FractalDimension - 2

#Select growth forms with >= 5 observations

FiveGFs <- table(MorphoDat$GrowthForm) >= 5
FiveGFs <- names(FiveGFs)[which(FiveGFs == T)]

MorphoDat <- filter(MorphoDat, GrowthForm %in% FiveGFs)

#export dataset

write.csv(MorphoDat, "2.DataSetup/MorphologyData_HighQualityWhole5GFs.csv", row.names = F)


#.....Biogeography data ====

#LI coords from google maps

LICoords<- data.frame(Lat = -14.668058, Long = 145.463566)

#Calculate bearing of each site from LI center point

BioGeoDat$bearing2LI <- map2_dbl(BioGeoDat$Lat, BioGeoDat$Long, bearingLatLong, Lat2 = LICoords$Lat[1], Long2 = LICoords$Long[1])

#save

write.csv(BioGeoDat, "2.DataSetup/BioGeoDat_PlusLIBearing.csv", row.names = F)

