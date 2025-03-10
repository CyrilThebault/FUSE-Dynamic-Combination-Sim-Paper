#! ---------------------------------------------------------------------------------------
#!
#! Description       :
#!
#! Authors           : Cyril Thebault <cyril.thebault@inrae.fr>
#!
#! Creation date     : 2025-12-02 12:00:57
#! Modification date :
#!
#! Comments          :
#!
#! ---------------------------------------------------------------------------------------

setwd("/Users/cyrilthebault/Postdoc_Ucal/02_DATA/git/FUSE-Dynamic-Combination-Sim-Paper")

#! ----------------------------- package loading

library(sf)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(ggpubr)
library(dplyr)
library(tidyr)
library(abind)

source(file = "Metrics.R")


##-----------------------------------------
##-------------- VARIABLES ----------------
##-----------------------------------------

dir_FUSE = '00_DATA'

su_HQ = loadRData(file.path(dir_FUSE, "SamplingUncertainty", "SU_HQ.Rdata"))
su_HQ <- su_HQ[order(rownames(su_HQ)), ]
su_HQ$range = su_HQ$p95 - su_HQ$p05

su_LQ = loadRData(file.path(dir_FUSE, "SamplingUncertainty", "SU_LQ.Rdata"))
su_LQ <- su_LQ[order(rownames(su_LQ)), ]
su_LQ$range = su_LQ$p95 - su_LQ$p05

su_WA = loadRData(file.path(dir_FUSE, "SamplingUncertainty", "SU_WA.Rdata"))
su_WA <- su_WA[order(rownames(su_WA)), ]
su_WA$range = su_WA$p95 - su_WA$p05

su_Mosa = loadRData(file.path(dir_FUSE, "SamplingUncertainty", "SU_Mosa.Rdata"))
su_Mosa <- su_Mosa[order(rownames(su_Mosa)), ]
su_Mosa$range = su_Mosa$p95 - su_Mosa$p05

su_Comp = loadRData(file.path(dir_FUSE, "SamplingUncertainty", "SU_Comp.Rdata"))
su_Comp <- su_Comp[order(rownames(su_Comp)), ]
su_Comp$range = su_Comp$p95 - su_Comp$p05

# Creating a function to classify Mosa vs WA
classify_performance <- function(p05_M, p95_M, score_M, p05_W, p95_W, score_W) {
  # Check for missing values
  if (is.na(p05_M) | is.na(p95_M) | is.na(p05_W) | is.na(p95_W)) {
    return(NA) # Assign NA if any value is missing
  }
  
  highest_score = max(c(score_M, score_W))
  lowest_score = min(c(score_M, score_W))
  
  p05 = ifelse(score_M == highest_score, p05_M, p05_W)
  p95 = ifelse(score_M == highest_score, p95_M, p95_W)
  
  
  # Check conditions for classification
  if (score_W == highest_score & p05 > lowest_score) {
    return("Better")
  } else if (score_M == highest_score & p05 > lowest_score) {
    return("Lower")
  } else {
    return("Equivalent")
  }
}



######## Code catchments

catchments = as.character(unname(unlist(read.table(file.path(dir_FUSE, "liste_BV_CAMELS_559.txt")))))

######## Shapefiles

dir_DB = 'Shp'

## catchments outlets
csv_outlets_DB <- readLines(file.path(dir_DB,"gauge_information.txt"))
csv_outlets_DB = strsplit(csv_outlets_DB, "\t")
csv_outlets_DB[[1]] = c("HUC_02", "GAGE_ID", "GAGE_NAME", "LAT", "LON", "DRAINAGE AREA (KM^2)")
csv_outlets_DB = data.frame(matrix(unlist(csv_outlets_DB), nrow=length(csv_outlets_DB), byrow=TRUE),stringsAsFactors=FALSE)
colnames(csv_outlets_DB) = csv_outlets_DB[1,]
csv_outlets_DB = csv_outlets_DB[! csv_outlets_DB$HUC_02 == "HUC_02",]
csv_outlets_DB$Station_lon = csv_outlets_DB$LON
csv_outlets_DB$Station_lat = csv_outlets_DB$LAT
csv_outlets_DB$Station_id = csv_outlets_DB$GAGE_ID
csv_outlets_DB$Country = "USA"

shp_outlets_DB <- sf::st_as_sf(csv_outlets_DB, coords = c('Station_lon','Station_lat'), crs = 4326)

modif = nchar(shp_outlets_DB$Station_id) == 7
shp_outlets_DB$Station_id[shp_outlets_DB$Country == 'USA' & modif] = paste0("0", shp_outlets_DB$Station_id[shp_outlets_DB$Country == 'USA' & modif])

shp_outlets_DB$Codes = paste(shp_outlets_DB$Country, shp_outlets_DB$Station_id, sep='_')

## rivers

## north-america
NorthAm = file.path(dir_DB,"boundaries_p_2021_v3.shp")
NorthAm = read_sf(NorthAm)
NorthAm = st_transform(NorthAm, crs = 4326)
NorthAm = NorthAm[NorthAm$COUNTRY == "USA",]

drop_these = c('US-AK', 'US-HI', 'US-PR', 'US-VI')
NorthAm = NorthAm[! NorthAm$STATEABB %in% drop_these,]


##################################
##            MAP 1             ##
##################################


# Apply to vectors
comparison_vector <- mapply(
  classify_performance, 
  su_Comp$p05, 
  su_Comp$p95,
  su_Comp$score,
  su_WA$p05, 
  su_WA$p95,
  su_WA$score
)

shp_outlets = subset(shp_outlets_DB, Codes %in% catchments)

# Add classification results to spatial dataset
shp_outlets$Diff_WA_Comp <- factor(comparison_vector, 
                           levels = c("Better", "Equivalent", "Lower"),
                           labels = c("Dynamic Combination is better", 
                                      "Both approaches are equivalent", 
                                      "FUSE - KGEcomp is better"))

shp_outlets_sub <- shp_outlets[!is.na(shp_outlets$Diff_WA_Comp), ]

# Update the color scheme
gg1 <- ggplot() +
  geom_sf(data = NorthAm, fill = "grey90", color = "grey40") +
  geom_sf(data = shp_outlets_sub, aes(fill = Diff_WA_Comp), size = 3, shape = 21, color = "black") +
  scale_fill_manual(values = c("Dynamic Combination is better" = "firebrick4",
                               "Both approaches are equivalent" = "white",
                               "FUSE - KGEcomp is better" = "darkorange")) +
  guides(fill = guide_legend(title = NULL)) +
  labs(title = "")+
  annotation_scale(location = "bl", style = "bar") +
  annotation_north_arrow(location = "br", style = north_arrow_fancy_orienteering, 
                         height = unit(1.1, "cm"), width = unit(1.1, "cm")) +
  theme_bw()+
  theme(legend.key.height = unit(1, "cm"))  # Adjust the value as needed


ggsave(plot = gg1, filename = "99_Figures/Fig11_MapEqui.png",
       width = 10, height = 4, dpi = 300)
