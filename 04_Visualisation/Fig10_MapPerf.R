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
library(viridis)

source(file = "Metrics.R")

#! ---------------------------- bassins

##-----------------------------------------
##-------------- VARIABLES ----------------
##-----------------------------------------

######## Code catchments

dir_FUSE = '00_DATA'
catchments = as.character(unname(unlist(read.table(file.path(dir_FUSE, "liste_BV_CAMELS_559.txt")))))

######## Shapefiles

dir_DB = paste0("Shp")

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


## north-america
NorthAm = file.path(dir_DB,"boundaries_p_2021_v3.shp")
NorthAm = read_sf(NorthAm)
NorthAm = st_transform(NorthAm, crs = 4326)
NorthAm = NorthAm[NorthAm$COUNTRY == "USA",]

drop_these = c('US-AK', 'US-HI', 'US-PR', 'US-VI')
NorthAm = NorthAm[! NorthAm$STATEABB %in% drop_these,]

######## Subset

shp_outlets = subset(shp_outlets_DB, Codes %in% catchments)


###### Add performance to the subset
Eval_FUSE_long_WA = loadRData(file.path(dir_FUSE, "WA", "Eval_FUSE.Rdata"))%>%
  mutate(`Cal : KGEcomp` = (`Cal : KGE :  1` + `Cal : KGE : -1`) / 2,
         `Eval : KGEcomp` = (`Eval : KGE :  1` + `Eval : KGE : -1`) / 2)

decision_medians_WA <- Eval_FUSE_long_WA %>%
  group_by(ModelDecisions) %>%
  summarise(median_Eval_KGEcomp = median(`Eval : KGEcomp`, na.rm = TRUE))

best_decision_WA <- decision_medians_WA %>%
  filter(median_Eval_KGEcomp == max(median_Eval_KGEcomp)) %>%
  pull(ModelDecisions)

Eval_FUSE_best_WA <- Eval_FUSE_long_WA %>%
  filter(ModelDecisions == best_decision_WA) %>%
  dplyr::select(Codes,`Eval : KGE :  1`, `Eval : KGE : -1`) %>%
  rename(`KGE(Q)` = `Eval : KGE :  1`, `KGE(1/Q)` = `Eval : KGE : -1`) %>%
  mutate(Type = "Dynamic Combination")


shp_outlets <- shp_outlets %>%
  left_join(Eval_FUSE_best_WA, by = "Codes")

###### Plot

shp_outlets = shp_outlets[! (is.na(shp_outlets$`KGE(Q)`)|is.na(shp_outlets$`KGE(1/Q)`)),]

gg1 <- ggplot() +
  geom_sf(data = NorthAm, fill = "grey90", color = "grey40") +
  geom_sf(data = shp_outlets, aes(fill = `KGE(Q)`), size = 3, shape = 21, color = "black") +
  scale_fill_gradientn(
    limits = c(-0.42, 1),
    colors = c(rep("purple",2), viridis(80)),
    breaks = seq(-0.4,1,0.1),
    labels = c("< -0.4", seq(-3,10,1)/10),
    oob = scales::squish,
    guide = guide_colorbar(title = "KGE (Q)", barwidth = 1.5, barheight = 10)
  ) +
  labs(title = "a)")+
  annotation_scale(location = "bl", style = "bar") +
  annotation_north_arrow(location = "br", style = north_arrow_fancy_orienteering, height = unit(1.1, "cm"), width = unit(1.1, "cm"))+
  theme_bw()

gg2 <- ggplot() +
  geom_sf(data = NorthAm, fill = "grey90", color = "grey40") +
  geom_sf(data = shp_outlets, aes(fill = `KGE(1/Q)`), size = 3, shape = 21, color = "black") +
  scale_fill_gradientn(
    limits = c(-0.42, 1),
    colors = c(rep("purple",2), viridis(80)),
    breaks = seq(-0.4,1,0.1),
    labels = c("< -0.4", seq(-3,10,1)/10),
    oob = scales::squish,
    guide = guide_colorbar(title = "KGE (1/Q)", barwidth = 1.5, barheight = 10)
  ) +
  labs(title = "b)")+
  annotation_scale(location = "bl", style = "bar") +
  annotation_north_arrow(location = "br", style = north_arrow_fancy_orienteering, height = unit(1.1, "cm"), width = unit(1.1, "cm"))+
  theme_bw()

combined_plot <- ggarrange(gg1, gg2, ncol = 1, nrow = 2)  

ggsave(plot = combined_plot, filename = "99_Figures/Fig10_MapPerf.png",
       width = 10, height = 8, dpi = 300)

