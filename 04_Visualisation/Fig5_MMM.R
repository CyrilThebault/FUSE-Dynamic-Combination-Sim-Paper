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

library(tidyverse)
library(sf)
library(ggplot2)
library(ggpubr)
library(viridis)

generate_random_colors <- function(num_colors, seed = 1997) {
  set.seed(seed)
  colors <- sample(colors(), num_colors)
  return(colors)
}

border_path <- "Shp/boundaries_p_2021_v3.shp"

# Load and process CONUS shapefile
gdf_borders <- st_read(border_path)
conus <- gdf_borders %>%
  filter(COUNTRY == "USA") %>%
  filter(!STATEABB %in% c("US-AK", "US-HI", "US-PR", "US-VI"))
conus <- st_transform(conus, 4326)

df_eval_gumboot = read.table("00_DATA/SamplingUncertainty/brief_analysis_modeling_results_with_gumboot_eval.csv",
               header = TRUE, sep = ",")

# Create spatial data
gdf_sub <- st_as_sf(df_eval_gumboot, coords = c("gauge_lon", "gauge_lat"), crs = 4326)
gdf_sub <- st_transform(gdf_sub, st_crs(conus))

# Remove basins without uncertainty info
gdf_gumboot <- gdf_sub %>% filter(score > -990)

list_names = gsub("_above_p5","",grep("_above_p5$", names(df_eval_gumboot), value = TRUE))

# Create a color dictionary
color_dict <- setNames(generate_random_colors(length(unique(list_names))),
                       unique(list_names))

gdf_gumboot <- gdf_gumboot %>%
  rowwise() %>%
  mutate(kgecomp_MMM_SU_eval = get(paste0(mod_need, "_val_kge"))) %>%
  ungroup()

p1 <-  ggplot() +
  geom_sf(data = conus, fill = "grey60", color = "white") +
  geom_sf(data = gdf_gumboot, aes(fill = kgecomp_MMM_SU_eval), size = 3, shape = 21, color = "black") +
  scale_fill_gradientn(
    limits = c(-0.42, 1),
    colors = c(rep("purple",2), viridis(80)),
    breaks = seq(-0.4,1,0.1),
    labels = c("< -0.4", seq(-3,10,1)/10),
    oob = scales::squish,
    guide = guide_colorbar(title = expression("KGE"[comp]), barwidth = 1.5, barheight = 15)
  ) +
  labs(x = "", y = "", title = "a)")+
  theme_bw()

outliers <- gdf_gumboot[!(gdf_gumboot$score >= -0.1), ]

p2 <- ggplot(gdf_gumboot, aes(x = kgecomp_MMM_SU_eval)) +
  geom_histogram(binwidth = 0.025, fill = "grey", color = "black") +
  xlim(-0.1, 1) +
  labs(x = expression("KGE"[comp]), y = "Basins", title = "")+
  geom_label(
    x = -0.1, y = Inf, 
    label = paste("Number of outliers: ", nrow(outliers)), 
    hjust = 0, vjust = 1.5, 
    size = 4, 
    color = "black", 
    fill = "white", 
    alpha = 0.7
  )+
  theme_bw()



models_needed <- names(sort(table(gdf_gumboot$mod_need), decreasing = TRUE))
models_needed_lab <-  gsub("_", "", models_needed)

basins_filled <- unname(sort(table(gdf_gumboot$mod_need), decreasing = TRUE))

# Loop over the models we identified and create incremental plots
ix =  length(models_needed) 
model <- models_needed[ix]

# Define bar plotting stats
bar_x <- models_needed # always the same but kept here for clarity
bar_y <- cumsum(basins_filled[1:(ix)])
bar_y <- c(bar_y, rep(-1, length(bar_x) - length(bar_y)))

# Initialize plot

gdf_plot = gdf_gumboot[gdf_gumboot$mod_need %in%  head(models_needed,ix),]

# Make the actual figure

p3 <- ggplot() +
  geom_sf(data = conus, fill = "grey60", color = "white", size = 0.5) +
  geom_sf(data = gdf_plot, aes(fill = factor(mod_need, levels = models_needed)), size = 3, shape = 21, color = "black") +
  scale_fill_manual(values = color_dict,
                    guide = 'none') +
  labs(x = "", y = "", title = "b)")+
  theme_bw()

p4 <- ggplot(data.frame(x = bar_x, y = bar_y), aes(x = factor(x, levels = models_needed), y = y, fill = factor(x, levels = models_needed))) +
  geom_bar(stat = "identity",  col = "black") +
  scale_fill_manual(values = color_dict,
                    guide = 'none') +
  scale_x_discrete(labels = models_needed_lab) +
  labs(title = "", x = "Models", y = "Catchments") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
  ylim(0, nrow(gdf_gumboot))


plot_kgecomp_eval <- ggarrange(
  p1, p2,p3,p4,
  ncol = 2, nrow = 2,
  widths = c(2, 1)  # Set width ratio
)

# Save the plot
ggsave(filename = "99_Figures/Fig5_MMM.png",
       plot = plot_kgecomp_eval, 
       width = 12, height = 8, dpi = 300, units = "in",
       bg = "white")
