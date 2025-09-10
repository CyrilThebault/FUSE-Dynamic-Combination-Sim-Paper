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

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(abind)

source(file = "Metrics.R")


##-----------------------------------------
##-------------- VARIABLES ----------------
##-----------------------------------------

dir_FUSE = file.path('00_DATA')

su_array = loadRData(file.path(dir_FUSE, "SamplingUncertainty.Rdata"))

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

##################################
##           Boxplot            ##
##################################

range_data <- su[,c("HQ","LQ", "WA"),c("KGE", "KGEinv"),"range"]

df_range <- melt(range_data, varnames = c("catchment", "model", "criterion"), value.name = "range")

df_range$model <- factor(df_range$model,
                         levels = rev(c("HQ", "LQ", "WA")),  
                         labels = rev(c("High-flow benchmark: KGE(Q)",
                                        "Low-flow benchmark: KGE(1/Q)",
                                        "Dynamic combination")))

df_range$criterion <- factor(df_range$criterion,
                             levels = c("KGE", "KGEinv"),  
                             labels = c("High-flow evaluation",
                                        "Low-flow evaluation"))


# Define color palette
color_dict <- c(
  "High-flow benchmark: KGE(Q)" = "dodgerblue4",
  "Low-flow benchmark: KGE(1/Q)" = "forestgreen",
  "Dynamic combination" = "firebrick4"
)


# Split the data
df_high <- subset(df_range, criterion == "High-flow evaluation")
df_low <- subset(df_range, criterion == "Low-flow evaluation")

# Create plot for High-flow
gg1 <- ggplot(df_high, aes(x = range, y = model, fill = model)) +
  stat_summary(
    geom = "boxplot",
    fun.data = function(x) {
      data.frame(
        y = median(x),
        ymin = quantile(x, 0.10),
        lower = quantile(x, 0.25),
        middle = median(x),
        upper = quantile(x, 0.75),
        ymax = quantile(x, 0.90)
      )
    },
    position = position_dodge2()
  ) +
  geom_vline(xintercept = 0, col = "black", linetype = "dashed") +
  scale_fill_manual(values = color_dict) +
  coord_cartesian(xlim = c(0, 0.6)) +
  labs(
    title = "a) High-flow evaluation",
    y = "",
    x = ""
  ) +
  guides(fill = "none") + 
  theme_bw()

# Create plot for Low-flow
gg2 <- ggplot(df_low, aes(x = range, y = model, fill = model)) +
  stat_summary(
    geom = "boxplot",
    fun.data = function(x) {
      data.frame(
        y = median(x),
        ymin = quantile(x, 0.10),
        lower = quantile(x, 0.25),
        middle = median(x),
        upper = quantile(x, 0.75),
        ymax = quantile(x, 0.90)
      )
    },
    position = position_dodge2()
  ) +
  geom_vline(xintercept = 0, col = "black", linetype = "dashed") +
  scale_fill_manual(values = color_dict) +
  coord_cartesian(xlim = c(0, 0.6)) +
  labs(
    title = "b) Low-flow evaluation",
    x = "Uncertainty range (p95 - p05 from gumboot)",
    y = ""
  ) +
  guides(fill = "none")+
  theme_bw()

# Combine plots vertically with patchwork
combined_plot <- gg1 / gg2

ggsave(plot = combined_plot, filename = "99_Figures/Fig10_BoxplotSU.png",
       width = 12, height = 6, dpi = 300)
