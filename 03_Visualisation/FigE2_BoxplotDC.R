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
library(dplyr)
library(tidyr)
library(patchwork)

source(file = "Metrics.R")


##-----------------------------------------
##-------------- VARIABLES ----------------
##-----------------------------------------

dir_FUSE = file.path("00_DATA")


# Selected decision overall

Eval_FUSE_long_WA = loadRData(file.path(dir_FUSE, "HQLQ_WA_q", "Eval_FUSE.Rdata"))%>%
  mutate(`Cal : KGEcomp` = (`Cal : KGE :  1` + `Cal : KGE : -1`) / 2,
         `Eval : KGEcomp` = (`Eval : KGE :  1` + `Eval : KGE : -1`) / 2)

selected_decision_allcatch <- "16_17_5"

Eval_FUSE_best_WA_Exp1 <- Eval_FUSE_long_WA %>%
  filter(ModelDecisions == selected_decision_allcatch) %>%
  dplyr::select(Codes,`Eval : KGE :  1`, `Eval : KGE : -1`) %>%
  rename(`KGE(Q)` = `Eval : KGE :  1`, `KGE(1/Q)` = `Eval : KGE : -1`) %>%
  mutate(Type = "Manual selection across all catchments")


# Selected decision overall

decision_medians_WA <- Eval_FUSE_long_WA %>%
  filter(!startsWith(as.character(ModelDecisions), "1_")) %>%
  group_by(ModelDecisions) %>%
  summarise(median_Cal_KGEcomp = median(`Cal : KGEcomp`, na.rm = TRUE))

best_decision_allcatch <- decision_medians_WA %>%
  filter(median_Cal_KGEcomp == max(median_Cal_KGEcomp)) %>%
  pull(ModelDecisions)

Eval_FUSE_best_WA_Exp2 <- Eval_FUSE_long_WA %>%
  filter(ModelDecisions == best_decision_allcatch) %>%
  dplyr::select(Codes,`Eval : KGE :  1`, `Eval : KGE : -1`) %>%
  rename(`KGE(Q)` = `Eval : KGE :  1`, `KGE(1/Q)` = `Eval : KGE : -1`) %>%
  mutate(Type = "Automatic selection across all catchments")



# Best decision per-catchment

best_decision_percatch <- Eval_FUSE_long_WA %>%
  filter(!startsWith(as.character(ModelDecisions), "1_")) %>%
  group_by(Codes) %>%
  summarise(
    best_model = {
      vals <- `Cal : KGEcomp`
      if (all(is.na(vals))) NA_character_
      else as.character(ModelDecisions[which.max(vals)])
    },
    .groups = "drop"
  )%>%
  pull(best_model)


Eval_FUSE_best_WA_Exp3 <- tibble(
  Codes = unique(as.character(Eval_FUSE_long_WA$Codes)),
  ModelDecisions = best_decision_percatch
) %>%
  left_join(
    Eval_FUSE_long_WA %>%
      select(Codes, ModelDecisions, `Eval : KGE :  1`, `Eval : KGE : -1`),
    by = c("Codes", "ModelDecisions")
  ) %>%
  transmute(
    Codes = factor(Codes, levels = unique(as.character(Eval_FUSE_long_WA$Codes))),
    `KGE(Q)` = `Eval : KGE :  1`,
    `KGE(1/Q)` = `Eval : KGE : -1`,
    Type = "Automatic selection per catchments"
  )


# Selected decision overall

Eval_FUSE_long_WA_tauvar = loadRData(file.path(dir_FUSE, "HQLQ_WA_q_tauvar", "Eval_FUSE.Rdata"))%>%
  mutate(`Cal : KGEcomp` = (`Cal : KGE :  1` + `Cal : KGE : -1`) / 2,
         `Eval : KGEcomp` = (`Eval : KGE :  1` + `Eval : KGE : -1`) / 2)

best_decision_percatch_tauvar <- Eval_FUSE_long_WA_tauvar %>%
  group_by(Codes) %>%
  summarise(
    best_model = {
      vals <- `Cal : KGEcomp`
      if (all(is.na(vals))) NA_character_
      else as.character(ModelDecisions[which.max(vals)])
    },
    .groups = "drop"
  )%>%
  pull(best_model)


Eval_FUSE_best_WA_Exp4 <- tibble(
  Codes = unique(as.character(Eval_FUSE_long_WA_tauvar$Codes)),
  ModelDecisions = best_decision_percatch_tauvar
) %>%
  left_join(
    Eval_FUSE_long_WA_tauvar %>%
      select(Codes, ModelDecisions, `Eval : KGE :  1`, `Eval : KGE : -1`),
    by = c("Codes", "ModelDecisions")
  ) %>%
  transmute(
    Codes = factor(Codes, levels = unique(as.character(Eval_FUSE_long_WA_tauvar$Codes))),
    `KGE(Q)` = `Eval : KGE :  1`,
    `KGE(1/Q)` = `Eval : KGE : -1`,
    Type = "Automatic selection per catchments (with dynamic τ)"
  )


##################################
##           Boxplot            ##
##################################

median_ini_HQ = median(Eval_FUSE_best_WA_Exp1$`KGE(Q)`, na.rm = TRUE)
median_ini_LQ = median(Eval_FUSE_best_WA_Exp1$`KGE(1/Q)`, na.rm = TRUE)

Eval_FUSE_best_All <- bind_rows(
  Eval_FUSE_best_WA_Exp1,
  Eval_FUSE_best_WA_Exp2,
  Eval_FUSE_best_WA_Exp3,
  Eval_FUSE_best_WA_Exp4
)

df_long <- Eval_FUSE_best_All %>%
  mutate(Type = factor(Type)) %>%
  pivot_longer(
    cols = c(`KGE(Q)`, `KGE(1/Q)`),
    names_to  = "Metric",
    values_to = "KGE"
  )




# Convert Type to factor with ordered levels
df_long$Type <- factor(df_long$Type, 
                           levels = c("Automatic selection per catchments (with dynamic τ)",
                                      "Automatic selection per catchments", 
                                      "Automatic selection across all catchments", 
                                      "Manual selection across all catchments"))

# Create gg
gg_HQ <- ggplot(filter(df_long, Metric == "KGE(Q)"), aes(x = KGE, y = Type, fill = Type)) +
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
               position = position_dodge2()) +
  coord_cartesian(xlim = c(-0.41,1))+
  scale_fill_brewer(palette = "Dark2") +
  geom_vline(xintercept = 1, col = "deepskyblue4", size = 0.8, linetype = "dashed") +
  geom_vline(xintercept = median_ini_HQ, col = "black", size = 0.5, linetype = "dashed") +
  labs(y = NULL, x = "KGE(Q)", title = "a) High-flow evaluation") +
  guides(fill = guide_legend(title = "Modelling approach"))+
  theme_bw() +
  theme(legend.position = "none")
# theme(axis.text.x = element_blank())

gg_HQ


gg_LQ <- ggplot(filter(df_long, Metric == "KGE(1/Q)"), aes(x = KGE, y = Type, fill = Type)) +
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
               position = position_dodge2()) +
  coord_cartesian(xlim = c(-0.41,1))+
  scale_fill_brewer(palette = "Dark2") +
  geom_vline(xintercept = 1, col = "deepskyblue4", size = 0.8, linetype = "dashed") +
  geom_vline(xintercept = median_ini_LQ, col = "black", size = 0.5, linetype = "dashed") +
  labs(y = NULL, x = "KGE(1/Q)", title = "b) Low-flow evaluation") +
  guides(fill = guide_legend(title = "Modelling approach"))+
  theme_bw() +
  theme(legend.position = "none")
# theme(axis.text.x = element_blank())

gg_LQ



# Combine plots vertically with patchwork
combined_plot <- gg_HQ / gg_LQ

ggsave(plot = combined_plot, filename = "99_Figures/FigE2_BoxplotDC.png",
       width = 12, height = 6, dpi = 300)

