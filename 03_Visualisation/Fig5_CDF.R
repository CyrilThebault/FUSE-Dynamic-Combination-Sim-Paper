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
library(ggpubr)

source(file = "Metrics.R")

##-----------------------------------------
##---------------- BEST MOD ---------------
##-----------------------------------------

dir_FUSE = file.path('00_DATA')

# FUSE HQ
Eval_FUSE_long_HQ = loadRData(file.path(dir_FUSE, "1", "Eval_FUSE.Rdata"))
Eval_FUSE_HQ <- data.frame(Eval_FUSE_long_HQ %>%
                             dplyr::select(Codes, ModelDecisions, `Cal : KGE :  1`) %>%  # dplyr::select only the relevant columns
                             pivot_wider(names_from = ModelDecisions, values_from = `Cal : KGE :  1`)) %>% 
  column_to_rownames(var = "Codes")

best_HQ = names(which.max(apply(Eval_FUSE_HQ, 2, median, na.rm = TRUE)))

# FUSE LQ
Eval_FUSE_long_LQ = loadRData(file.path(dir_FUSE, "-1", "Eval_FUSE.Rdata"))

Eval_FUSE_LQ <- data.frame(Eval_FUSE_long_LQ %>%
                             dplyr::select(Codes, ModelDecisions, `Cal : KGE : -1`) %>%  # dplyr::select only the relevant columns
                             pivot_wider(names_from = ModelDecisions, values_from = `Cal : KGE : -1`)) %>%
  column_to_rownames(var = "Codes")

best_LQ = names(which.max(apply(Eval_FUSE_LQ, 2, median, na.rm = TRUE)))


best_models <- c(
  "High-flow ensemble: KGE(Q)"    = paste0(best_HQ,".x"),
  "Low-flow ensemble: KGE(1/Q)"  = paste0(best_LQ,".y")
)


##-----------------------------------------
##---------------- Part1 ------------------
##-----------------------------------------

# FUSE HQ
Eval_FUSE_long_HQ = loadRData(file.path(dir_FUSE, "1", "Eval_FUSE.Rdata"))
Eval_FUSE_HQ <- data.frame(Eval_FUSE_long_HQ %>%
                             dplyr::select(Codes, ModelDecisions, `Eval : KGE :  1`) %>%  # dplyr::select only the relevant columns
                             pivot_wider(names_from = ModelDecisions, values_from = `Eval : KGE :  1`)) %>% 
  column_to_rownames(var = "Codes")

# FUSE LQ
Eval_FUSE_long_LQ = loadRData(file.path(dir_FUSE, "-1", "Eval_FUSE.Rdata"))

Eval_FUSE_LQ <- data.frame(Eval_FUSE_long_LQ %>%
                             dplyr::select(Codes, ModelDecisions, `Eval : KGE :  1`) %>%  # dplyr::select only the relevant columns
                             pivot_wider(names_from = ModelDecisions, values_from = `Eval : KGE :  1`)) %>%
  column_to_rownames(var = "Codes")

# Ensure the rownames are preserved as a column for merging
Eval_FUSE_HQ$rownames <- rownames(Eval_FUSE_HQ)
Eval_FUSE_LQ$rownames <- rownames(Eval_FUSE_LQ)


data_frames <- list(Eval_FUSE_HQ, Eval_FUSE_LQ)

# Merge the data frames by rownames
combined_df <- reduce(data_frames, full_join, by = "rownames")

colnames(combined_df)[colnames(combined_df) == "rownames"] = "Catchment"

# Pivot the data frame into a long format
long_df <- combined_df %>%
  pivot_longer(
    cols = -Catchment,            # Pivot all columns except 'Catchment'
    names_to = "Model",           # New column for the model names
    values_to = "KGE_Eval"        # New column for the model values
  )

long_df$KGE_Eval[is.na(long_df$KGE_Eval)] = -9999

# Create a CDF dataframe for each ModelDecision
cdf_combined <- long_df %>%
  group_by(Model) %>%
  arrange(as.numeric(KGE_Eval)) %>%
  mutate(CDF = seq(0, 1, length.out = n()),
         Model = as.character(Model)) %>%
  ungroup()



# Modify the color_dict and add a new column for legend grouping
cdf_combined <- cdf_combined %>%
  mutate(LegendGroup = case_when(
    grepl("^X[0-9]+.x$", Model) ~ "High-flow ensemble: KGE(Q)",
    grepl("^X[0-9]+.y$", Model) ~ "Low-flow ensemble: KGE(1/Q)"
  ))

cdf_combined$LegendGroup <- factor(cdf_combined$LegendGroup, 
                                   levels = c("High-flow ensemble: KGE(Q)","Low-flow ensemble: KGE(1/Q)"))

cdf_combined <- cdf_combined %>%
  mutate(IsBest = ifelse(Model %in% best_models, TRUE, FALSE))

# cdf_combined$LegendGroup[cdf_combined$Model %in% c("X138","X96","X170","X34")] = "FUSE Original Models"


color_dict <- c(
  "High-flow ensemble: KGE(Q)" = "dodgerblue4",
  "Low-flow ensemble: KGE(1/Q)" = "forestgreen"
)

targets = c(0.125, 0.25, 0.375, 0.5, 0.625 ,0.75, 0.875)

find_closest_cdf <- function(cdf_values, targets) {
  sapply(targets, function(target) {
    which.min(abs(cdf_values - target))  # Index of the closest CDF value
  })
}

CDF_values = unique(cdf_combined$CDF)

index = find_closest_cdf(CDF_values, targets)

approx_target = CDF_values[index]

cdf_boxplot_data <- cdf_combined %>%
  filter(CDF %in% approx_target)

# Create a small offset for the boxplots to avoid overlap (shift by LegendGroup)
cdf_boxplot_data <- cdf_boxplot_data %>%
  mutate(CDF_offset = case_when(
    LegendGroup == "High-flow ensemble: KGE(Q)" ~ CDF,
    LegendGroup == "Low-flow ensemble: KGE(1/Q)" ~ CDF
  ))

# Create the plot
gg1 <- ggplot(cdf_combined, aes(x = KGE_Eval, y = CDF, group = Model)) +
  # FUSE ensemble CDF
  geom_line(data = subset(cdf_combined, !IsBest), 
            aes(color = LegendGroup), linewidth = 0.4, alpha = 0.2) +
  
  # Best models CDF
  geom_line(data = subset(cdf_combined, IsBest), color = "black", linewidth = 1.4, alpha = 1) +
  geom_line(data = subset(cdf_combined, IsBest), 
            aes(color = LegendGroup), linewidth = 1.2, alpha = 1) +
  
  # Add vertical reference line
  geom_vline(xintercept = 1, col = "black", linetype = "dashed") +
  
  # Add boxplots at CDF = 0.25, 0.5, and 0.75 for each LegendGroup
  geom_boxplot(data = cdf_boxplot_data, aes(x = KGE_Eval, y = CDF_offset, group = interaction(CDF, LegendGroup),
                                            fill = LegendGroup),
               width = 0.05, alpha = 0.5, outlier.shape = NA) +  # `width` controls horizontal width
  
  # Labels and theme
  labs(x = "KGE(Q) - Evaluation Period", y = "Cumulative distribution function (CDF)", color = "Model", title = "a) High-flow evaluation") +
  theme_bw() +
  coord_cartesian(xlim = c(-2.02, 1.02), ylim = c(-0.02, 1.02), expand = FALSE) +
  
  # Customize colors
  scale_color_manual(values = color_dict, 
                     breaks = c("Composite ensemble: KGEcomp", "High-flow ensemble: KGE(Q)", "Low-flow ensemble: KGE(1/Q)"),
                     labels = c("Composite ensemble: KGEcomp", "High-flow ensemble: KGE(Q)", "Low-flow ensemble: KGE(1/Q)")) +
  scale_fill_manual(values = color_dict) +  # Match boxplot fill colors
  
  guides(
    color = guide_legend(
      title = NULL,  # remove title
      override.aes = list(alpha = 1, linewidth = 1)  # solid colors in legend
    ),
    fill = "none"  # remove the fill legend ("LegendGroup" one)
  ) +
  
  # Remove legend for cleaner visualization
  theme(  legend.position = c(0.02, 0.98),  # x and y coordinates (in [0,1] range)
          legend.justification = c(0, 1),  # aligns legend box to top-left corner
          legend.background = element_rect(fill = "white", color = "black"))   


##-----------------------------------------
##---------------- Part2 ------------------
##-----------------------------------------


# FUSE HQ
Eval_FUSE_long_HQ = loadRData(file.path(dir_FUSE, "1", "Eval_FUSE.Rdata"))
Eval_FUSE_HQ <- data.frame(Eval_FUSE_long_HQ %>%
                             dplyr::select(Codes, ModelDecisions, `Eval : KGE : -1`) %>%  # dplyr::select only the relevant columns
                             pivot_wider(names_from = ModelDecisions, values_from = `Eval : KGE : -1`)) %>% 
  column_to_rownames(var = "Codes")

# FUSE LQ
Eval_FUSE_long_LQ = loadRData(file.path(dir_FUSE, "-1", "Eval_FUSE.Rdata"))

Eval_FUSE_LQ <- data.frame(Eval_FUSE_long_LQ %>%
                             dplyr::select(Codes, ModelDecisions, `Eval : KGE : -1`) %>%  # dplyr::select only the relevant columns
                             pivot_wider(names_from = ModelDecisions, values_from = `Eval : KGE : -1`)) %>%
  column_to_rownames(var = "Codes")


# Ensure the rownames are preserved as a column for merging
Eval_FUSE_HQ$rownames <- rownames(Eval_FUSE_HQ)
Eval_FUSE_LQ$rownames <- rownames(Eval_FUSE_LQ)


data_frames <- list(Eval_FUSE_HQ, Eval_FUSE_LQ)


# Merge the data frames by rownames
combined_df <- reduce(data_frames, full_join, by = "rownames")

colnames(combined_df)[colnames(combined_df) == "rownames"] = "Catchment"

# Pivot the data frame into a long format
long_df <- combined_df %>%
  pivot_longer(
    cols = -Catchment,            # Pivot all columns except 'Catchment'
    names_to = "Model",           # New column for the model names
    values_to = "KGE_Eval"        # New column for the model values
  )

long_df$KGE_Eval[is.na(long_df$KGE_Eval)] = -9999

# Create a CDF dataframe for each ModelDecision
cdf_combined <- long_df %>%
  group_by(Model) %>%
  arrange(as.numeric(KGE_Eval)) %>%
  mutate(CDF = seq(0, 1, length.out = n()),
         Model = as.character(Model)) %>%
  ungroup()



# Modify the color_dict and add a new column for legend grouping
cdf_combined <- cdf_combined %>%
  mutate(LegendGroup = case_when(
    grepl("^X[0-9]+.x$", Model) ~ "High-flow ensemble: KGE(Q)",
    grepl("^X[0-9]+.y$", Model) ~ "Low-flow ensemble: KGE(1/Q)"
  ))

cdf_combined$LegendGroup <- factor(cdf_combined$LegendGroup, 
                                   levels = c("High-flow ensemble: KGE(Q)","Low-flow ensemble: KGE(1/Q)"))

cdf_combined <- cdf_combined %>%
  mutate(IsBest = ifelse(Model %in% best_models, TRUE, FALSE))


# cdf_combined$LegendGroup[cdf_combined$Model %in% c("X138","X96","X170","X34")] = "FUSE Original Models"

color_dict <- c(
  "High-flow ensemble: KGE(Q)" = "dodgerblue4",
  "Low-flow ensemble: KGE(1/Q)" = "forestgreen"
)

targets = c(0.125, 0.25, 0.375, 0.5, 0.625 ,0.75, 0.875)

find_closest_cdf <- function(cdf_values, targets) {
  sapply(targets, function(target) {
    which.min(abs(cdf_values - target))  # Index of the closest CDF value
  })
}

CDF_values = unique(cdf_combined$CDF)

index = find_closest_cdf(CDF_values, targets)

approx_target = CDF_values[index]

cdf_boxplot_data <- cdf_combined %>%
  filter(CDF %in% approx_target)

# Create a small offset for the boxplots to avoid overlap (shift by LegendGroup)
cdf_boxplot_data <- cdf_boxplot_data %>%
  mutate(CDF_offset = case_when(
    LegendGroup == "High-flow ensemble: KGE(Q)" ~ CDF,
    LegendGroup == "Low-flow ensemble: KGE(1/Q)" ~ CDF,
  ))

# Create the plot
gg2 <- ggplot(cdf_combined, aes(x = KGE_Eval, y = CDF, group = Model)) +
  # FUSE ensemble CDF
  geom_line(data = subset(cdf_combined, !IsBest), 
            aes(color = LegendGroup), linewidth = 0.4, alpha = 0.2) +
  
  # Best models CDF
  geom_line(data = subset(cdf_combined, IsBest), color = "black", linewidth = 1.4, alpha = 1) +
  geom_line(data = subset(cdf_combined, IsBest), 
            aes(color = LegendGroup), linewidth = 1.2, alpha = 1) +
  
  # Add vertical reference line
  geom_vline(xintercept = 1, col = "black", linetype = "dashed") +
  
  # Add boxplots at CDF = 0.25, 0.5, and 0.75 for each LegendGroup
  geom_boxplot(data = cdf_boxplot_data, aes(x = KGE_Eval, y = CDF_offset, group = interaction(CDF, LegendGroup),
                                            fill = LegendGroup),
               width = 0.05, alpha = 0.5, outlier.shape = NA) +  # `width` controls horizontal width
  
  # Labels and theme
  labs(x = "KGE(1/Q) - Evaluation Period", y = "Cumulative distribution function (CDF)", color = "Model", title = "b) Low-flow evaluation") +
  theme_bw() +
  coord_cartesian(xlim = c(-2.02, 1.02), ylim = c(-0.02, 1.02), expand = FALSE) +
  
  # Customize colors
  scale_color_manual(values = color_dict, 
                     breaks = c("High-flow ensemble: KGE(Q)", "Low-flow ensemble: KGE(1/Q)"),
                     labels = c("High-flow ensemble: KGE(Q)", "Low-flow ensemble: KGE(1/Q)")) +
  scale_fill_manual(values = color_dict) +  # Match boxplot fill colors
  
  # Remove legend for cleaner visualization
  theme(legend.position = "none")   


##-----------------------------------------
##---------------- Plot ------------------
##-----------------------------------------

combined_plot <- ggarrange(gg1, gg2, ncol = 2, nrow = 1)    

ggsave(plot = combined_plot, filename = "99_Figures/Fig5_CDF.png",
       width = 12, height = 6, dpi = 600)
