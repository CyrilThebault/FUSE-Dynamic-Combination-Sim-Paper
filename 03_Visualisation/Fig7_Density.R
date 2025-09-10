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

library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(purrr)
library(ggnewscale)
library(ggpubr)

source(file = "Metrics.R")

##-----------------------------------------
##---------------- MAIN ------------------
##-----------------------------------------

dir_FUSE = file.path("00_DATA")

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
  mutate(Type = "High-flow\nbenchmark:\nKGE(Q)")



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
  mutate(Type = "Low-flow\nbenchmark:\nKGE(1/Q)")



# FUSE WA
Eval_FUSE_long_WA = loadRData(file.path(dir_FUSE, "HQLQ_WA_q", "Eval_FUSE.Rdata"))

best_decision_WA <- "16_17_5"

Eval_FUSE_best_WA <- Eval_FUSE_long_WA %>%
  filter(ModelDecisions == best_decision_WA) %>%
  dplyr::select(Codes,`Eval : KGE :  1`, `Eval : KGE : -1`) %>%
  rename(`KGE(Q)` = `Eval : KGE :  1`, `KGE(1/Q)` = `Eval : KGE : -1`) %>%
  mutate(Type = "Dynamic\ncombination")



Eval_FUSE_best_HQ = Eval_FUSE_best_HQ[,c("Type", "KGE(Q)", "KGE(1/Q)")]
Eval_FUSE_best_LQ = Eval_FUSE_best_LQ[,c("Type", "KGE(Q)", "KGE(1/Q)")]
Eval_FUSE_best_WA = Eval_FUSE_best_WA[,c("Type", "KGE(Q)", "KGE(1/Q)")]

# Ensure the rownames are preserved as a column for merging
combined_df <-rbind(Eval_FUSE_best_HQ,Eval_FUSE_best_LQ,
                    Eval_FUSE_best_WA
)

##########################
#------Density plot------#
##########################

# Initialize the ggplot object
gg1 <- ggplot() + 
  theme_bw() + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(-2.02, 1), ylim = c(-2.02, 1))+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  theme(
    legend.position = 'none',
    # axis.text.x = element_blank(),  # Remove x-axis text
    # axis.text.y = element_blank(),  # Remove y-axis text
    # axis.ticks.x = element_blank(), # Remove x-axis ticks
    # axis.ticks.y = element_blank(), # Remove y-axis ticks
    axis.title.x = element_blank(), # Remove x-axis label
    axis.title.y = element_blank()  # Remove y-axis label
    # panel.spacing.x = unit(1.5, "lines"),
    # panel.spacing.y = unit(1.5, "lines"),
    # plot.margin = unit(c(1, 1.5, 1, 1), "lines")
  )

# Create a red color ramp with 10 levels of varying alpha
bins = 5
alpha_levels <- seq(0, 1, length.out = bins)
mycolors = c("dodgerblue4", "forestgreen", "firebrick4")


# Add data for each Type with distinct scales using ggnewscale
mytypes = unique(combined_df$Type)
# mytypes = c("FUSE - KGE(Q)", "FUSE - KGE(1/Q)")

mycount = 0
for (type in mytypes) {
  mycount = mycount+1
  
  mycolor_ramp <- sapply(alpha_levels, function(alpha) scales::alpha(mycolors[mycount], alpha))
  
  gg1 <- gg1 + 
    new_scale_fill() +  # Reset fill scale for each Type
    geom_density_2d_filled(
      data = combined_df %>% filter(Type == type),
      aes(x = `KGE(1/Q)`, y = `KGE(Q)`, fill = ..level..),  # Use `..level..` for fill
      bins = bins, 
      n = 3000,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = mycolor_ramp)
}




##########################
#------Boxplot HQ--------#
##########################


# Define color palette
color_dict <- c(
  "High-flow\nbenchmark:\nKGE(Q)" = "dodgerblue4",
  "Low-flow\nbenchmark:\nKGE(1/Q)" = "forestgreen",
  "Dynamic\ncombination" = "firebrick4"
)


# Convert Type to factor with ordered levels
combined_df$Type <- factor(combined_df$Type, 
                           levels = rev(c("Dynamic\ncombination", 
                                          "Low-flow\nbenchmark:\nKGE(1/Q)", 
                                          "High-flow\nbenchmark:\nKGE(Q)" )))

# Create gg1
gg2 <- ggplot(combined_df, aes(x = Type, y = `KGE(Q)`, fill = Type)) +
  geom_violin(alpha = 0.6) +
  stat_summary(geom = "boxplot", 
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
               width = 0.2) +
  scale_fill_manual(values = color_dict) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-2.02, 1))+
  geom_hline(yintercept = 1, col = "black", linetype = "dashed") +
  labs(x = NULL, y = "KGE(Q)") +
  guides(fill = guide_legend(title = "Modelling approach"))+
  theme_bw() +
  theme(legend.position = "none")


##########################
#------Boxplot LQ--------#
##########################

# Convert Type to factor with ordered levels
combined_df$Type <- factor(combined_df$Type, 
                           levels = c("Dynamic\ncombination", 
                                      "Low-flow\nbenchmark:\nKGE(1/Q)", 
                                      "High-flow\nbenchmark:\nKGE(Q)" ))

gg3 <- ggplot(combined_df, aes(x = `KGE(1/Q)`, y = Type, fill = Type)) +
  geom_violin(alpha = 0.6) +
  stat_summary(geom = "boxplot",
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
               width = 0.2) +
  scale_fill_manual(values = color_dict) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(-2.02, 1))+
  geom_vline(xintercept = 1, col = "black", linetype = "dashed") +
  labs(y = NULL, x = "KGE(1/Q)") +
  guides(fill = guide_legend(title = "Modelling approach"))+
  theme_bw() +
  theme(legend.position = "none")

##########################
#-------Final plot-------#
##########################

# Create an empty plot for D
emptyplot <- ggplot() + theme_void()

# Define the layout
layout <- "
BAA
BAA
DCC
"

# Combine the plots
combined_plot <- gg1 + gg2 + gg3 + emptyplot + 
  plot_layout(design = layout)


ggsave(plot = combined_plot, filename = "99_Figures/Fig7_Density.png",
       width = 8.2, height = 8.2, dpi = 300)
