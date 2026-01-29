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


# Using MAE

Eval_FUSE_long_WA_MAE = loadRData(file.path(dir_FUSE, "HQLQ_WA_q", "Eval_FUSE.Rdata"))%>%
  mutate(`Cal : KGEcomp` = (`Cal : KGE :  1` + `Cal : KGE : -1`) / 2,
         `Eval : KGEcomp` = (`Eval : KGE :  1` + `Eval : KGE : -1`) / 2)

selected_decision_allcatch <- "16_17_5"

Eval_FUSE_best_WA_MAE <- Eval_FUSE_long_WA_MAE %>%
  filter(ModelDecisions == selected_decision_allcatch) %>%
  dplyr::select(Codes,`Eval : KGE :  1`, `Eval : KGE : -1`) %>%
  rename(`KGE(Q)` = `Eval : KGE :  1`, `KGE(1/Q)` = `Eval : KGE : -1`) %>%
  mutate(Type = "MAE")



# Using KGE(Q) and KGE(1/Q) depending of quantile

Eval_FUSE_long_WA_HQLQ = loadRData(file.path(dir_FUSE, "HQLQ_WA_q_HQLQ", "Eval_FUSE.Rdata"))%>%
  mutate(`Cal : KGEcomp` = (`Cal : KGE :  1` + `Cal : KGE : -1`) / 2,
         `Eval : KGEcomp` = (`Eval : KGE :  1` + `Eval : KGE : -1`) / 2)

selected_decision_allcatch <- "16_17_5"

Eval_FUSE_best_WA_HQLQ <- Eval_FUSE_long_WA_HQLQ %>%
  filter(ModelDecisions == selected_decision_allcatch) %>%
  dplyr::select(Codes,`Eval : KGE :  1`, `Eval : KGE : -1`) %>%
  rename(`KGE(Q)` = `Eval : KGE :  1`, `KGE(1/Q)` = `Eval : KGE : -1`) %>%
  mutate(Type =  "KGE(Q) and KGE(1/Q)")



# Using KGEcomp depending of quantile

Eval_FUSE_long_WA_KGEcomp = loadRData(file.path(dir_FUSE, "HQLQ_WA_q_KGEcomp", "Eval_FUSE.Rdata"))%>%
  mutate(`Cal : KGEcomp` = (`Cal : KGE :  1` + `Cal : KGE : -1`) / 2,
         `Eval : KGEcomp` = (`Eval : KGE :  1` + `Eval : KGE : -1`) / 2)

selected_decision_allcatch <- "16_17_5"

Eval_FUSE_best_WA_KGEcomp <- Eval_FUSE_long_WA_KGEcomp %>%
  filter(ModelDecisions == selected_decision_allcatch) %>%
  dplyr::select(Codes,`Eval : KGE :  1`, `Eval : KGE : -1`) %>%
  rename(`KGE(Q)` = `Eval : KGE :  1`, `KGE(1/Q)` = `Eval : KGE : -1`) %>%
  mutate(Type =  "KGEcomp")




##################################
##           Boxplot            ##
##################################

Eval_FUSE_best_All <- bind_rows(
  Eval_FUSE_best_WA_MAE,
  Eval_FUSE_best_WA_HQLQ,
  Eval_FUSE_best_WA_KGEcomp
)

# Split MAE and other types
df_mae <- Eval_FUSE_best_All %>%
  filter(Type == "MAE") %>%
  select(Codes,
         KGE_Q_MAE  = `KGE(Q)`,
         KGE_LQ_MAE = `KGE(1/Q)`)

df_other <- Eval_FUSE_best_All %>%
  filter(Type != "MAE")

# Join so MAE is always on x-axis
df_plot <- df_other %>%
  left_join(df_mae, by = "Codes")

p_high <- ggplot(df_plot,
                 aes(x = KGE_Q_MAE,
                     y = `KGE(Q)`,
                     color = Type)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_cartesian(xlim=c(-2,1), ylim=c(-2,1), expand = FALSE)+
  labs(
    title = "(a) High-flow evaluation",
    x = "KGE(Q)\nreference metric used in the dynamic combination (MAE)",
    y = "KGE(Q)\nalternative metric tested in the dynamic combination",
    color = "Type"
  ) +
  theme_bw()+
  theme(legend.title = element_blank())


p_low <- ggplot(df_plot,
                aes(x = KGE_LQ_MAE,
                    y = `KGE(1/Q)`,
                    color = Type)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_cartesian(xlim=c(-2,1), ylim=c(-2,1), expand = FALSE)+
  labs(
    title = "(b) Low-flow evaluation",
    x = "KGE(1/Q)\nreference metric used in the dynamic combination (MAE)",
    y = "KGE(1/Q)\nalternative metric tested in the dynamic combination",
    color = "Type"
  ) +
  theme_bw()+
  theme(legend.title = element_blank())


combined_plot <- (p_high + p_low + guide_area()) +
  plot_layout(
    ncol = 3,
    widths = c(1, 1, 0.5),
    guides = "collect"
  ) &
  theme(legend.position = "right")

ggsave(plot = combined_plot, filename = "99_Figures/FigF1_ScatterComp.png",
       width = 10, height = 4.8, dpi = 300)

