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
library(abind)
library(gridExtra)

source(file = "Metrics.R")


##-----------------------------------------
##-------------- VARIABLES ----------------
##-----------------------------------------

dir_FUSE = '00_DATA'

su_array = loadRData(file.path(dir_FUSE,  "SamplingUncertainty.Rdata"))

su_numeric <- su_array[, , , setdiff(dimnames(su_array)[[4]], "GOF_stat")]
su_numeric <- apply(su_numeric, c(1,2,3,4), as.numeric)

# Force numeric conversion before subtraction
p95_numeric <- as.numeric(su_array[,,, "p95"])
p05_numeric <- as.numeric(su_array[,,, "p05"])

# Compute range
range_array <- su_numeric[,,, "p95"] - su_numeric[,,, "p05"]
dim(range_array) <- c(dim(range_array), 1)
dimnames(range_array)[[4]] <- "range"

# Attach using abind
su <- abind(su_numeric, range_array, along = 4)

# Extract relevant dimension names
catchments <- dimnames(su)[[1]]
decisions <- dimnames(su)[[2]]
metrics <- c("KGE", "KGEinv")
variables <- dimnames(su)[[4]]

# Identify index of WA and other decisions
wa_index <- which(decisions == "WA")
other_indices <- which(decisions == "HQ" | decisions == "LQ")

# Create an array to store the classification results
result_array <- array(NA, dim = c(length(catchments), length(metrics), length(other_indices)),
                      dimnames = list(catchments, metrics, paste("WA_vs", decisions[other_indices], sep = "_")))

# Classification function
classify_performance <- function(p05_M, p95_M, score_M, p05_W, p95_W, score_W) {
  if (is.na(p05_M) | is.na(p95_M) | is.na(p05_W) | is.na(p95_W)) {
    return(NA)
  }
  
  highest_score = max(c(score_M, score_W))
  lowest_score = min(c(score_M, score_W))
  
  p05 = ifelse(score_M == highest_score, p05_M, p05_W)
  p95 = ifelse(score_M == highest_score, p95_M, p95_W)
  
  if (score_W == highest_score & p05 > lowest_score) {
    return("Better")
  } else if (score_M == highest_score & p05 > lowest_score) {
    return("Worst")
  } else {
    return("Equivalent")
  }
}

# Loop through catchments, metrics, and other decisions
for (i in seq_along(catchments)) {
  for (j in seq_along(metrics)) {
    for (k in seq_along(other_indices)) {
      idx_other <- other_indices[k]
      
      # Extract required values
      p05_M <- su[i, idx_other, j, which(variables == "p05")]
      p95_M <- su[i, idx_other, j, which(variables == "p95")]
      score_M <- su[i, idx_other, j, which(variables == "p50")]
      
      p05_W <- su[i, wa_index, j, which(variables == "p05")]
      p95_W <- su[i, wa_index, j, which(variables == "p95")]
      score_W <- su[i, wa_index, j, which(variables == "p50")]
      
      # Apply classification
      result_array[i, j, k] <- classify_performance(p05_M, p95_M, score_M, p05_W, p95_W, score_W)
    }
  }
}

# Identify rows (catchments) that contain any NA across the full 3D array
na_catchments <- dimnames(result_array)[[1]][apply(result_array, 1, function(x) any(is.na(x)))]

# Drop those rows
result_array <- result_array[!(dimnames(result_array)[[1]] %in% na_catchments), , ]

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

shp_outlets = subset(shp_outlets_DB, Codes %in% dimnames(result_array)[[1]])

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

# Define fill colors
fill_colors <- c("Better" = "blue", "Equivalent" = "white", "Worst" = "red")

# Extract dimension names
comparisons <- dimnames(result_array)[[3]]
metrics <- dimnames(result_array)[[2]]

# Create a list to store plots
plot_list <- list()

# Loop over each metric and comparison
for (comparison in comparisons) {
  for (metric in metrics) {
    
    # Create a new column in shp_outlets for current comparison/metric
    new_col_name <- paste0(comparison, "_", metric)
    shp_outlets[[new_col_name]] <- result_array[, metric, comparison]
    
    # Generate map plot
    p <- ggplot() +
      geom_sf(data = NorthAm, fill = "grey90", color = "grey40") +
      geom_sf(data = shp_outlets, aes_string(fill = new_col_name), size = 3, shape = 21, color = "black") +
      scale_fill_manual(values = fill_colors, na.value = "grey80") +
      annotation_scale(location = "bl", style = "bar") +
      annotation_north_arrow(location = "br", style = north_arrow_fancy_orienteering,
                             height = unit(1.1, "cm"), width = unit(1.1, "cm")) +
      theme_bw() +
      theme(legend.position = "none") 
    
    # Save the plot to list
    plot_list[[paste(metric, comparison, sep = "_")]] <- p
  }
}

gg1 = grid.arrange(grobs = plot_list, ncol = 2)


ggsave(plot = gg1, filename = "99_Figures/Fig11_MapEqui.png",
       width = 12, height = 6, dpi = 300)
