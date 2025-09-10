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

setwd("/Users/cyrilthebault/Postdoc_Ucal/02_DATA/FUSE-Dynamic-Combination-Sim-Paper")

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

## rivers

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
# FUSE HQ
Eval_FUSE_long_HQ = loadRData(file.path(dir_FUSE, "1", "Eval_FUSE.Rdata"))

decision_medians_HQ <- Eval_FUSE_long_HQ %>%
  group_by(ModelDecisions) %>%
  summarise(median_Cal_KGE = median(`Cal : KGE :  1`, na.rm = TRUE))

best_decision_HQ <- decision_medians_HQ %>%
  filter(median_Cal_KGE == max(median_Cal_KGE)) %>%
  pull(ModelDecisions)

Eval_FUSE_best_HQ <- Eval_FUSE_long_HQ %>%
  filter(ModelDecisions == best_decision_HQ) %>%
  dplyr::select(Codes,`Eval : KGE :  1`, `Eval : KGE : -1`) %>%
  rename(`KGE(Q)` = `Eval : KGE :  1`, `KGE(1/Q)` = `Eval : KGE : -1`) %>%
  mutate(Type = "High-flow benchmark: KGE(Q)")



# FUSE LQ
Eval_FUSE_long_LQ = loadRData(file.path(dir_FUSE, "-1", "Eval_FUSE.Rdata"))

decision_medians_LQ <- Eval_FUSE_long_LQ %>%
  group_by(ModelDecisions) %>%
  summarise(median_Cal_KGEinv = median(`Cal : KGE : -1`, na.rm = TRUE))

best_decision_LQ <- decision_medians_LQ %>%
  filter(median_Cal_KGEinv == max(median_Cal_KGEinv)) %>%
  pull(ModelDecisions)

Eval_FUSE_best_LQ <- Eval_FUSE_long_LQ %>%
  filter(ModelDecisions == best_decision_LQ) %>%
  dplyr::select(Codes,`Eval : KGE :  1`, `Eval : KGE : -1`) %>%
  rename(`KGE(Q)` = `Eval : KGE :  1`, `KGE(1/Q)` = `Eval : KGE : -1`) %>%
  mutate(Type = "Low-flow benchmark: KGE(1/Q)")



# FUSE WA
Eval_FUSE_long_WA = loadRData(file.path(dir_FUSE, "HQLQ_WA_q", "Eval_FUSE.Rdata"))

best_decision_WA <- "16_17_5"

Eval_FUSE_best_WA <- Eval_FUSE_long_WA %>%
  filter(ModelDecisions == best_decision_WA) %>%
  dplyr::select(Codes,`Eval : KGE :  1`, `Eval : KGE : -1`) %>%
  rename(`KGE(Q)` = `Eval : KGE :  1`, `KGE(1/Q)` = `Eval : KGE : -1`) %>%
  mutate(Type = "Dynamic combination",
         `Difference with high-flow benchmark: KGE(Q)` = `KGE(Q)` - Eval_FUSE_best_HQ$`KGE(Q)`,
         `Difference with low-flow benchmark: KGE(1/Q)` = `KGE(1/Q)` - Eval_FUSE_best_LQ$`KGE(1/Q)`
  )



shp_outlets <- shp_outlets %>%
  left_join(Eval_FUSE_best_WA, by = "Codes")

###### Plot

shp_outlets = shp_outlets[! (is.na(shp_outlets$`KGE(Q)`)|
                               is.na(shp_outlets$`KGE(1/Q)`)|
                               is.na(shp_outlets$`Difference with high-flow benchmark: KGE(Q)`)|
                               is.na(shp_outlets$`Difference with low-flow benchmark: KGE(1/Q)`)
),]

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
  labs(title = "a) High-flow evaluation")+
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
  labs(title = "b) Low-flow evaluation")+
  annotation_scale(location = "bl", style = "bar") +
  annotation_north_arrow(location = "br", style = north_arrow_fancy_orienteering, height = unit(1.1, "cm"), width = unit(1.1, "cm"))+
  theme_bw()

gg3 <- ggplot() +
  geom_sf(data = NorthAm, fill = "grey90", color = "grey40") +
  geom_sf(
    data = shp_outlets,
    aes(fill = `Difference with high-flow benchmark: KGE(Q)`),
    size = 3, shape = 21, color = "black"
  ) +
  scale_fill_gradient2(
    low = "red",      # For negative values
    mid = "white",    # For zero
    high = "blue",    # For positive values
    midpoint = 0,
    limits = c(-0.5, 0.5),
    oob = scales::squish,  # Values beyond limits will be squished to end colors
    name = "Δ KGE(Q)\nwith high-flow\nbenchmark",     # Legend title (you can customize it)
    breaks = c(-0.5, 0, 0.5),
    labels = c("< -0.5", "0", "> 0.5")
  ) +
  labs(title = " ") +
  annotation_scale(location = "bl", style = "bar") +
  annotation_north_arrow(
    location = "br",
    style = north_arrow_fancy_orienteering,
    height = unit(1.1, "cm"),
    width = unit(1.1, "cm")
  ) +
  theme_bw()

gg4 <- ggplot() +
  geom_sf(data = NorthAm, fill = "grey90", color = "grey40") +
  geom_sf(
    data = shp_outlets,
    aes(fill = `Difference with low-flow benchmark: KGE(1/Q)`),
    size = 3, shape = 21, color = "black"
  ) +
  scale_fill_gradient2(
    low = "red",      # For negative values
    mid = "white",    # For zero
    high = "blue",    # For positive values
    midpoint = 0,
    limits = c(-0.5, 0.5),
    oob = scales::squish,  # Values beyond limits will be squished to end colors
    name = "Δ KGE(1/Q)\nwith low-flow\nbenchmark",     # Legend title (you can customize it)
    breaks = c(-0.5, 0, 0.5),
    labels = c("< -0.5", "0", "> 0.5")
  ) +
  labs(title = " ") +
  annotation_scale(location = "bl", style = "bar") +
  annotation_north_arrow(
    location = "br",
    style = north_arrow_fancy_orienteering,
    height = unit(1.1, "cm"),
    width = unit(1.1, "cm")
  ) +
  theme_bw()

# Define the layout
layout <- "
AB
CD
"

# Combine the plots
combined_plot <- gg1 + gg3 + gg2 + gg4 + 
  plot_layout(design = layout)


ggsave(plot = combined_plot, filename = "99_Figures/Fig8_MapPerf.png",
       width = 12, height = 6, dpi = 300)

