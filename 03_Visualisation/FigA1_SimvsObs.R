#! ---------------------------------------------------------------------------------------
#!
#! Description       :
#!
#! Authors           : Cyril Thebault <cyril.thebault@inrae.fr>
#!
#! Creation date     : 2022-08-08 16:54:28
#! Modification date :
#!
#! Comments          :
#!
#! ---------------------------------------------------------------------------------------


#! ----------------------------- path definition

#! -------------- sources

#! ----------------------------- workbook directory

setwd("/Users/cyrilthebault/Postdoc_Ucal/02_DATA/FUSE-Dynamic-Combination-Sim-Paper")

#! ----------------------------- package loading

library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(purrr)
library(ggnewscale)
library(ggpubr)
library(patchwork)

##-----------------------------------------
##---------------- MAIN ------------------
##-----------------------------------------

dir_FUSE = file.path('00_DATA')

# FUSE WA SIM
Eval_FUSE_long_WA_Sim = loadRData(file.path(dir_FUSE, "HQLQ_WA_q", "Eval_FUSE.Rdata"))

# FUSE WA SIM
Eval_FUSE_long_WA_Obs = loadRData(file.path(dir_FUSE, "HQLQ_Obs_q", "Eval_FUSE.Rdata"))


df_sim <- Eval_FUSE_long_WA_Sim %>%
  select(Codes, ModelDecisions, KGE_eval = `Eval : KGE :  1`, KGEinv_eval = `Eval : KGE : -1`) %>%
  mutate(Framework = "Sim")

df_obs <- Eval_FUSE_long_WA_Obs %>%
  select(Codes, ModelDecisions, KGE_eval = `Eval : KGE :  1`, KGEinv_eval = `Eval : KGE : -1`) %>%
  mutate(Framework = "Obs")

df_all <- bind_rows(df_sim, df_obs)

##-----------------------------------------
##---------------- PLOT -------------------
##-----------------------------------------


df_HQ <- df_all %>%
  select(Codes, ModelDecisions, KGE_eval, Framework) %>%
  tidyr::pivot_wider(names_from = Framework, values_from = KGE_eval)%>%
  mutate(Diff = Sim - Obs)


gg1 <- ggplot(df_HQ, aes(x = ModelDecisions, y = Diff)) +
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
    position = position_dodge2(),
    fill=NA,
    col = "grey70"
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 16,     
    size = 1,
    color = "red"
  ) +

  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Parameter sets of the dynamic combination",
    y = expression(Delta * "KGE(Q)"),
    title = "(a) Difference in high-flow performance using simulations or observations to define current window"
  ) +
  coord_cartesian(ylim = c(   -0.5, 0.25)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),  
        panel.grid.minor.x = element_blank()  
        )



df_LQ <- df_all %>%
  select(Codes, ModelDecisions, KGEinv_eval, Framework) %>%
  tidyr::pivot_wider(names_from = Framework, values_from = KGEinv_eval)%>%
  mutate(Diff = Sim - Obs)


gg2 <- ggplot(df_LQ, aes(x = ModelDecisions, y = Diff)) +
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
    position = position_dodge2(),
    fill=NA,
    col = "grey70"
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 16,     
    size = 1,
    color = "red"
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Parameter sets of the dynamic combination",
    y = expression(Delta * "KGE(1/Q)"),
    title = "(b) Difference in low-flow performance using simulations or observations to define current window"
  ) +
  coord_cartesian(ylim = c(   -0.5, 0.25)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),  
        panel.grid.minor.x = element_blank()  
  )


combined_plot <- gg1 / gg2

ggsave(plot = combined_plot, filename = "99_Figures/FigA1_SimvsObs.png",
       width = 12, height = 8, dpi = 600)