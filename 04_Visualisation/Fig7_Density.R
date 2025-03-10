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

# FUSE HQ
Eval_FUSE_long_HQ = loadRData(file.path(dir_FUSE, "1", "Eval_FUSE.Rdata")) %>%
  mutate(`Cal : KGEcomp` = (`Cal : KGE :  1` + `Cal : KGE : -1`) / 2,
         `Eval : KGEcomp` = (`Eval : KGE :  1` + `Eval : KGE : -1`) / 2)

decision_medians_HQ <- Eval_FUSE_long_HQ %>%
  group_by(ModelDecisions) %>%
  summarise(median_Cal_KGEcomp = median(`Cal : KGEcomp`, na.rm = TRUE))

best_decision_HQ <- decision_medians_HQ %>%
  filter(median_Cal_KGEcomp == max(median_Cal_KGEcomp)) %>%
  pull(ModelDecisions)

Eval_FUSE_best_HQ <- Eval_FUSE_long_HQ %>%
  filter(ModelDecisions == best_decision_HQ) %>%
  dplyr::select(Codes,`Eval : KGE :  1`, `Eval : KGE : -1`) %>%
  rename(`KGE(Q)` = `Eval : KGE :  1`, `KGE(1/Q)` = `Eval : KGE : -1`) %>%
  mutate(Type = "FUSE - KGE(Q)")



# FUSE LQ
Eval_FUSE_long_LQ = loadRData(file.path(dir_FUSE, "-1", "Eval_FUSE.Rdata"))%>%
  mutate(`Cal : KGEcomp` = (`Cal : KGE :  1` + `Cal : KGE : -1`) / 2,
         `Eval : KGEcomp` = (`Eval : KGE :  1` + `Eval : KGE : -1`) / 2)

decision_medians_LQ <- Eval_FUSE_long_LQ %>%
  group_by(ModelDecisions) %>%
  summarise(median_Cal_KGEcomp = median(`Cal : KGEcomp`, na.rm = TRUE))

best_decision_LQ <- decision_medians_LQ %>%
  filter(median_Cal_KGEcomp == max(median_Cal_KGEcomp)) %>%
  pull(ModelDecisions)

Eval_FUSE_best_LQ <- Eval_FUSE_long_LQ %>%
  filter(ModelDecisions == best_decision_LQ) %>%
  dplyr::select(Codes,`Eval : KGE :  1`, `Eval : KGE : -1`) %>%
  rename(`KGE(Q)` = `Eval : KGE :  1`, `KGE(1/Q)` = `Eval : KGE : -1`) %>%
  mutate(Type = "FUSE - KGE(1/Q)")



# FUSE comp
Eval_FUSE_long_Comp = loadRData(file.path(dir_FUSE, "Comp", "Eval_FUSE.Rdata"))%>%
  mutate(`Cal : KGEcomp` = (`Cal : KGE :  1` + `Cal : KGE : -1`) / 2,
         `Eval : KGEcomp` = (`Eval : KGE :  1` + `Eval : KGE : -1`) / 2)

decision_medians_Comp <- Eval_FUSE_long_Comp %>%
  group_by(ModelDecisions) %>%
  summarise(median_Cal_KGEcomp = median(`Cal : KGEcomp`, na.rm = TRUE))

best_decision_Comp <- decision_medians_Comp %>%
  filter(median_Cal_KGEcomp == max(median_Cal_KGEcomp)) %>%
  pull(ModelDecisions)

Eval_FUSE_best_Comp <- Eval_FUSE_long_Comp %>%
  filter(ModelDecisions == best_decision_Comp) %>%
  dplyr::select(Codes,`Eval : KGE :  1`, `Eval : KGE : -1`) %>%
  rename(`KGE(Q)` = `Eval : KGE :  1`, `KGE(1/Q)` = `Eval : KGE : -1`) %>%
  mutate(Type = "FUSE - KGEcomp")


# FUSE Mosa
df_eval_gumboot = read.table("00_DATA/SamplingUncertainty/brief_analysis_modeling_results_with_gumboot_eval.csv",
                             header = TRUE, sep = ",")
df_eval_gumboot$Codes = ifelse(nchar(df_eval_gumboot$gauge_id) == 7, paste0("USA_0", df_eval_gumboot$gauge_id), paste0("USA_", df_eval_gumboot$gauge_id))

Eval_FUSE_long_HQ_tmp <- Eval_FUSE_long_HQ %>%
  mutate(ModelDecisions = paste0("HF_", ModelDecisions))

Eval_FUSE_long_LQ_tmp <- Eval_FUSE_long_LQ %>%
  mutate(ModelDecisions = paste0("LF_", ModelDecisions))

Eval_FUSE_long_Mosa <- bind_rows(Eval_FUSE_long_HQ_tmp, Eval_FUSE_long_LQ_tmp)

df_lookup <- data.frame(Codes = df_eval_gumboot$Codes, mod_need = df_eval_gumboot$mod_need)

Eval_FUSE_best_Mosa <- Eval_FUSE_long_Mosa %>%
  left_join(df_lookup, by = "Codes") %>%
  group_by(Codes) %>%
  filter(if_else(is.na(mod_need), row_number() == 1, ModelDecisions == mod_need)) %>% 
  ungroup()

Eval_FUSE_best_Mosa[is.na(Eval_FUSE_best_Mosa$mod_need),3:12] = NA

Eval_FUSE_best_Mosa <- Eval_FUSE_best_Mosa %>%
  dplyr::select(Codes,`Eval : KGE :  1`, `Eval : KGE : -1`) %>%
  rename(`KGE(Q)` = `Eval : KGE :  1`, `KGE(1/Q)` = `Eval : KGE : -1`) %>%
  mutate(Type = "Multi Model Mosaic")


# FUSE WA
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




Eval_FUSE_best_HQ = Eval_FUSE_best_HQ[,c("Type", "KGE(Q)", "KGE(1/Q)")]
Eval_FUSE_best_LQ = Eval_FUSE_best_LQ[,c("Type", "KGE(Q)", "KGE(1/Q)")]
Eval_FUSE_best_Comp = Eval_FUSE_best_Comp[,c("Type", "KGE(Q)", "KGE(1/Q)")]
Eval_FUSE_best_Mosa = Eval_FUSE_best_Mosa[,c("Type", "KGE(Q)", "KGE(1/Q)")]
Eval_FUSE_best_WA = Eval_FUSE_best_WA[,c("Type", "KGE(Q)", "KGE(1/Q)")]


##########################
#------Density plot------#
##########################

# Ensure the rownames are preserved as a column for merging
combined_df <-rbind(Eval_FUSE_best_HQ,Eval_FUSE_best_LQ,
                    Eval_FUSE_best_Comp
                    # Eval_LowBenchmark_best_HQ,Eval_LowBenchmark_best_LQ,
                    # Eval_MidBenchmark_best_HQ,Eval_MidBenchmark_best_LQ,
                    # Benchmark,
                    # Eval_FUSE_best_WA
                    )

# Initialize the ggplot object
gg1 <- ggplot() + 
  theme_bw() + 
  labs(
    x = "KGE(1/Q)",
    y = "KGE(Q)",
    title = "a)"
  ) + 
  scale_x_continuous(limits = c(-0.41, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.41, 1), expand = c(0, 0)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  theme(
    legend.position = 'none'#,
    # panel.spacing.x = unit(1.5, "lines"),
    # panel.spacing.y = unit(1.5, "lines"),
    # plot.margin = unit(c(1, 1.5, 1, 1), "lines")
  )

# Create a red color ramp with 10 levels of varying alpha
bins = 5
alpha_levels <- seq(0, 1, length.out = bins)
mycolors = c("dodgerblue4", "forestgreen", "darkorange")


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
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = mycolor_ramp)
}

# Display the plot
gg1

##########################
#------Density plot------#
##########################

# Ensure the rownames are preserved as a column for merging
combined_df <-rbind(Eval_FUSE_best_HQ,Eval_FUSE_best_LQ,
                    # Eval_FUSE_best_Comp,
                    # Eval_LowBenchmark_best_HQ,Eval_LowBenchmark_best_LQ,
                    # Eval_MidBenchmark_best_HQ,Eval_MidBenchmark_best_LQ,
                    Eval_FUSE_best_Mosa
                    # Eval_FUSE_best_WA
)

# Initialize the ggplot object
gg2 <- ggplot() + 
  theme_bw() + 
  labs(
    x = "KGE(1/Q)",
    y = "KGE(Q)",
    title = "b)"
  ) + 
  scale_x_continuous(limits = c(-0.41, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.41, 1), expand = c(0, 0)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  theme(
    legend.position = 'none'#,
    # panel.spacing.x = unit(1.5, "lines"),
    # panel.spacing.y = unit(1.5, "lines"),
    # plot.margin = unit(c(1, 1.5, 1, 1), "lines")
  )

# Create a red color ramp with 10 levels of varying alpha
bins = 5
alpha_levels <- seq(0, 1, length.out = bins)
mycolors = c("dodgerblue4", "forestgreen", "lightpink")


# Add data for each Type with distinct scales using ggnewscale
mytypes = unique(combined_df$Type)
# mytypes = c("FUSE - KGE(Q)", "FUSE - KGE(1/Q)")

mycount = 0
for (type in mytypes) {
  mycount = mycount+1
  
  mycolor_ramp <- sapply(alpha_levels, function(alpha) scales::alpha(mycolors[mycount], alpha))
  
  gg2 <- gg2 + 
    new_scale_fill() +  # Reset fill scale for each Type
    geom_density_2d_filled(
      data = combined_df %>% filter(Type == type),
      aes(x = `KGE(1/Q)`, y = `KGE(Q)`, fill = ..level..),  # Use `..level..` for fill
      bins = bins, 
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = mycolor_ramp)
}

# Display the plot
gg2


##########################
#------Density plot------#
##########################

# Ensure the rownames are preserved as a column for merging
combined_df <-rbind(Eval_FUSE_best_HQ,Eval_FUSE_best_LQ,
                    # Eval_FUSE_best_Comp,
                    # Eval_LowBenchmark_best_HQ,Eval_LowBenchmark_best_LQ,
                    # Eval_MidBenchmark_best_HQ,Eval_MidBenchmark_best_LQ,
                    # Benchmark,
                    Eval_FUSE_best_WA
                    )

# Initialize the ggplot object
gg3 <- ggplot() + 
  theme_bw() + 
  labs(
    x = "KGE(1/Q)",
    y = "KGE(Q)",
    title = "c)"
  ) + 
  scale_x_continuous(limits = c(-0.41, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.41, 1), expand = c(0, 0)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  theme(
    legend.position = 'none'#,
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
  
  gg3 <- gg3 + 
    new_scale_fill() +  # Reset fill scale for each Type
    geom_density_2d_filled(
      data = combined_df %>% filter(Type == type),
      aes(x = `KGE(1/Q)`, y = `KGE(Q)`, fill = ..level..),  # Use `..level..` for fill
      bins = bins, 
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = mycolor_ramp)
}

# Display the plot
gg3


##########################
#-------Final plot-------#
##########################

combined_plot <- ggarrange(gg1, gg2, gg3, ncol = 3, nrow = 1)    

combined_plot
ggsave(plot = combined_plot, filename = "99_Figures/Fig7_Density.png",
       width = 12, height = 4, dpi = 300)
