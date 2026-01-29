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
library(RColorBrewer)
library(patchwork)

source(file = "Metrics.R")

##-----------------------------------------
##---------------- MAIN ------------------
##-----------------------------------------

dir_FUSE = file.path("00_DATA")

# Single model HQ
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
  mutate(Type = "SM-HF\n(benchmark)")



# Single model LQ
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
  mutate(Type = "SM-LF\n(benchmark)")



# Combi static time + static space HQ (top-3 single mod)
Eval_FUSE_long_SAT_HQ = loadRData(file.path(dir_FUSE, "SA_HQ", "Eval_FUSE.Rdata"))

top3_decisions_HQ <- decision_medians_HQ %>%
  arrange(desc(median_Cal_KGE)) %>%
  slice_head(n = 3) %>%
  pull(ModelDecisions)

mat <- Eval_FUSE_long_SAT_HQ["KGE", , "Eval", paste(sort(top3_decisions_HQ), collapse = "_"), ]
Eval_FUSE_best_SAT_HQ <- data.frame(
  Codes    = factor(colnames(mat)),
  `KGE(Q)`   = as.numeric(mat["1",  ]),
  `KGE(1/Q)` = as.numeric(mat["-1", ]),
  Type     = "MM1-HF"
)

names(Eval_FUSE_best_SAT_HQ)[names(Eval_FUSE_best_SAT_HQ) == "KGE.Q."]   <- "KGE(Q)"
names(Eval_FUSE_best_SAT_HQ)[names(Eval_FUSE_best_SAT_HQ) == "KGE.1.Q."] <- "KGE(1/Q)"


# Combi static time + static space LQ (top-3 single mod)
Eval_FUSE_long_SAT_LQ = loadRData(file.path(dir_FUSE, "SA_LQ", "Eval_FUSE.Rdata"))

top3_decisions_LQ <- decision_medians_LQ %>%
  arrange(desc(median_Cal_KGEinv)) %>%
  slice_head(n = 3) %>%
  pull(ModelDecisions)

mat <- Eval_FUSE_long_SAT_LQ["KGE", , "Eval", paste(sort(top3_decisions_LQ), collapse = "_"), ]
Eval_FUSE_best_SAT_LQ <- data.frame(
  Codes    = factor(colnames(mat)),
  `KGE(Q)`   = as.numeric(mat["1",  ]),
  `KGE(1/Q)` = as.numeric(mat["-1", ]),
  Type     = "MM1-LF")

names(Eval_FUSE_best_SAT_LQ)[names(Eval_FUSE_best_SAT_LQ) == "KGE.Q."]   <- "KGE(Q)"
names(Eval_FUSE_best_SAT_LQ)[names(Eval_FUSE_best_SAT_LQ) == "KGE.1.Q."] <- "KGE(1/Q)"


# Combi static time + variable space HQ (top-3 single mod)
# Eval_FUSE_long_SA_HQ = loadRData(file.path(dir_FUSE, "SA_HQ", "Eval_FUSE.Rdata"))
Eval_FUSE_long_SA_HQ = Eval_FUSE_long_SAT_HQ

top3_per_catchment_HQ <- Eval_FUSE_long_HQ %>%
  group_by(Codes, ModelDecisions) %>%
  summarise(median_Cal_KGE = median(`Cal : KGE :  1`, na.rm = TRUE), .groups = "drop") %>%
  group_by(Codes) %>%
  slice_max(median_Cal_KGE, n = 3, with_ties = FALSE) %>%
  arrange(Codes, desc(median_Cal_KGE)) %>%
  summarise(model_key = paste(sort(ModelDecisions), collapse = "_"), .groups = "drop")

mat <- mapply(
  FUN = function(code, key) {
    if (!key %in% dimnames(Eval_FUSE_long_SA_HQ)$model) return(c(`1` = NA_real_, `-1` = NA_real_))
    Eval_FUSE_long_SA_HQ["KGE", c("1", "-1"), "Eval", key, code]
  },
  code = as.character(top3_per_catchment_HQ$Codes),
  key  = top3_per_catchment_HQ$model_key,
  SIMPLIFY = TRUE
)

Eval_FUSE_best_SA_HQ <- data.frame(
  Codes    = factor(colnames(mat)),
  `KGE(Q)`   = as.numeric(mat["1",  ]),
  `KGE(1/Q)` = as.numeric(mat["-1", ]),
  Type     = "MM2-HF")

names(Eval_FUSE_best_SA_HQ)[names(Eval_FUSE_best_SA_HQ) == "KGE.Q."]   <- "KGE(Q)"
names(Eval_FUSE_best_SA_HQ)[names(Eval_FUSE_best_SA_HQ) == "KGE.1.Q."] <- "KGE(1/Q)"


# Combi static time + variable space LQ (top-3 single mod)
#Eval_FUSE_long_SA_LQ = loadRData(file.path(dir_FUSE, "SA_LQ", "Eval_FUSE.Rdata"))
Eval_FUSE_long_SA_LQ = Eval_FUSE_long_SAT_LQ

top3_per_catchment_LQ <- Eval_FUSE_long_LQ %>%
  group_by(Codes, ModelDecisions) %>%
  summarise(median_Cal_KGEinv = median(`Cal : KGE : -1`, na.rm = TRUE), .groups = "drop") %>%
  group_by(Codes) %>%
  slice_max(median_Cal_KGEinv, n = 3, with_ties = FALSE) %>%
  arrange(Codes, desc(median_Cal_KGEinv)) %>%
  summarise(model_key = paste(sort(ModelDecisions), collapse = "_"), .groups = "drop")

mat <- mapply(
  FUN = function(code, key) {
    if (!key %in% dimnames(Eval_FUSE_long_SA_LQ)$model) return(c(`1` = NA_real_, `-1` = NA_real_))
    Eval_FUSE_long_SA_LQ["KGE", c("1", "-1"), "Eval", key, code]
  },
  code = as.character(top3_per_catchment_LQ$Codes),
  key  = top3_per_catchment_LQ$model_key,
  SIMPLIFY = TRUE
)

Eval_FUSE_best_SA_LQ <- data.frame(
  Codes    = factor(colnames(mat)),
  `KGE(Q)`   = as.numeric(mat["1",  ]),
  `KGE(1/Q)` = as.numeric(mat["-1", ]),
  Type     = "MM2-LF")

names(Eval_FUSE_best_SA_LQ)[names(Eval_FUSE_best_SA_LQ) == "KGE.Q."]   <- "KGE(Q)"
names(Eval_FUSE_best_SA_LQ)[names(Eval_FUSE_best_SA_LQ) == "KGE.1.Q."] <- "KGE(1/Q)"



# Combi static time + variable space HQ (top combi of 3 models)
#Eval_FUSE_long_SAopti_HQ = loadRData(file.path(dir_FUSE, "SA_HQ", "Eval_FUSE.Rdata"))
Eval_FUSE_long_SAopti_HQ = Eval_FUSE_long_SAT_HQ

all_models <- dimnames(Eval_FUSE_long_SAopti_HQ)$model
models3 <- all_models[grepl("^[0-9]+(_[0-9]+){2}$", all_models)]  # e.g. "198_30_194"
codes <- dimnames(Eval_FUSE_long_SAopti_HQ)$code

best_combi_per_catchment <- sapply(codes, function(cd) {
  v <- Eval_FUSE_long_SAopti_HQ["KGE", "1", "Cal", models3, cd]  # restrict to models3
  models3[which.max(replace(v, is.na(v), -Inf))]             # ignore NA in max
})

mat <- mapply(
  FUN = function(cd, mdl) Eval_FUSE_long_SAopti_HQ["KGE", c("1","-1"), "Eval", mdl, cd],
  cd  = codes,
  mdl = best_combi_per_catchment,
  SIMPLIFY = TRUE
)

Eval_FUSE_best_SAopti_HQ <- data.frame(
  Codes      = factor(codes, levels = codes),
  `KGE(Q)`   = as.numeric(mat["1",  ]),
  `KGE(1/Q)` = as.numeric(mat["-1", ]),
  Type     = "MM3-HF",
  check.names = FALSE
)

# Combi static time + variable space LQ (top combi of 3 models)
#Eval_FUSE_long_SAopti_LQ = loadRData(file.path(dir_FUSE, "SA_LQ", "Eval_FUSE.Rdata"))
Eval_FUSE_long_SAopti_LQ = Eval_FUSE_long_SAT_LQ

all_models <- dimnames(Eval_FUSE_long_SAopti_LQ)$model
models3 <- all_models[grepl("^[0-9]+(_[0-9]+){2}$", all_models)]  # e.g. "198_30_194"
codes <- dimnames(Eval_FUSE_long_SAopti_LQ)$code

best_combi_per_catchment <- sapply(codes, function(cd) {
  v <- Eval_FUSE_long_SAopti_LQ["KGE", "-1", "Cal", models3, cd]  # restrict to models3
  models3[which.max(replace(v, is.na(v), -Inf))]             # ignore NA in max
})

mat <- mapply(
  FUN = function(cd, mdl) Eval_FUSE_long_SAopti_LQ["KGE", c("1","-1"), "Eval", mdl, cd],
  cd  = codes,
  mdl = best_combi_per_catchment,
  SIMPLIFY = TRUE
)

Eval_FUSE_best_SAopti_LQ <- data.frame(
  Codes      = factor(codes, levels = codes),
  `KGE(Q)`   = as.numeric(mat["1",  ]),
  `KGE(1/Q)` = as.numeric(mat["-1", ]),
  Type     = "MM3-LF",
  check.names = FALSE
)


#####

Eval_FUSE_best_HQ = Eval_FUSE_best_HQ[,c("Type", "KGE(Q)", "KGE(1/Q)")]
Eval_FUSE_best_LQ = Eval_FUSE_best_LQ[,c("Type", "KGE(Q)", "KGE(1/Q)")]
Eval_FUSE_best_SAT_HQ = Eval_FUSE_best_SAT_HQ[,c("Type", "KGE(Q)", "KGE(1/Q)")]
Eval_FUSE_best_SAT_LQ = Eval_FUSE_best_SAT_LQ[,c("Type", "KGE(Q)", "KGE(1/Q)")]
Eval_FUSE_best_SA_HQ = Eval_FUSE_best_SA_HQ[,c("Type", "KGE(Q)", "KGE(1/Q)")]
Eval_FUSE_best_SA_LQ = Eval_FUSE_best_SA_LQ[,c("Type", "KGE(Q)", "KGE(1/Q)")]
Eval_FUSE_best_SAopti_HQ = Eval_FUSE_best_SAopti_HQ[,c("Type", "KGE(Q)", "KGE(1/Q)")]
Eval_FUSE_best_SAopti_LQ = Eval_FUSE_best_SAopti_LQ[,c("Type", "KGE(Q)", "KGE(1/Q)")]

# Ensure the rownames are preserved as a column for merging
combined_df <-rbind(Eval_FUSE_best_HQ,Eval_FUSE_best_LQ,
                    Eval_FUSE_best_SAT_HQ,Eval_FUSE_best_SAT_LQ,
                    Eval_FUSE_best_SA_HQ,Eval_FUSE_best_SA_LQ,
                    Eval_FUSE_best_SAopti_HQ,Eval_FUSE_best_SAopti_LQ
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


# Add data for each Type with distinct scales using ggnewscale
mytypes = unique(combined_df$Type)

# Create a color ramp with 5 levels of varying alpha
bins = 5
alpha_levels <- seq(0, 1, length.out = bins)
mycolors <- brewer.pal(n = length(mytypes), name = "Paired")

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
color_dict <- setNames(mycolors, mytypes)


# Convert Type to factor with ordered levels
combined_df$Type <- factor(combined_df$Type, 
                           levels = rev(mytypes))

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
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
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
                           levels = rev(mytypes))

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
BBAA
BBAA
DDCC
DDCC
"

# Combine the plots
combined_plot <- gg1 + gg2 + gg3 + emptyplot + 
  plot_layout(design = layout)


ggsave(plot = combined_plot, filename = "99_Figures/FigD1_DensityBenchmarks.png",
       width = 10.2, height = 10.2, dpi = 300)
