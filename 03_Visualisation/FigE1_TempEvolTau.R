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
library(reshape2)
library(ggpubr)
library(data.table)

source(file = "Metrics.R")

##-----------------------------------------
##---------------- PART 1 -----------------
##-----------------------------------------

dir_FUSE = '00_DATA'

decisions = read.table(file.path(dir_FUSE, "list_decision_78.txt"), sep = ";", header = TRUE)$ID

catchment = "USA_01022500"

# FUSE ensemble
QsimHQ = loadRData(file.path(dir_FUSE, "1", "SimArray.Rdata"))
QsimLQ = loadRData(file.path(dir_FUSE, "-1", "SimArray.Rdata"))

df_QsimHQ = QsimHQ[,,catchment]
df_QsimLQ = QsimLQ[,,catchment]

colnames(df_QsimHQ) = paste0("HQ",colnames(df_QsimHQ))
colnames(df_QsimLQ) = paste0("LQ",colnames(df_QsimLQ))

rm(QsimHQ)
rm(QsimLQ)


# Dynamic combination outputs
WA_RES = loadRData(file.path(dir_FUSE, "HQLQ_WA_q_tauvar", paste0(catchment,"_tauvar_9_7.Rdata"))) 

DatesR = WA_RES$dates

simulated = cbind(data.frame(Dates = DatesR, WA = WA_RES$weighted_avg), df_QsimHQ, df_QsimLQ)

# Obs
Obs = loadRData(file.path(dir_FUSE,"ObservedStreamflow.Rdata"))$Qobs_mm
observed = data.frame(Dates = DatesR, Q = Obs[,catchment])


# Tau 
tauvar = WA_RES$tauvar

tau <- data.frame(
  Dates = DatesR,
  tau  = tauvar
)

#------------- Plot

Ind_Plot = which(format(DatesR, "%Y") == "2005")

simulated_plot = simulated[Ind_Plot,]
simulated_plot_long = simulated_plot %>%
  pivot_longer(
    cols = -c(Dates, WA),
    names_to = "Simulation",
    values_to = "Q",
    names_transform = list(Simulation = as.character)  # Ensure names are character
  )

observed_plot = observed[Ind_Plot,]


# Adjust the hydrograph plot
hydrograph_plot <- ggplot() +
  geom_line(data = simulated_plot_long, aes(x = Dates, y = Q, group = Simulation), color = "grey75") +  #####ICI
  geom_line(data = simulated_plot, aes(x = Dates, y = WA), color = "firebrick4") + 
  geom_line(data = observed_plot, aes(x = Dates, y = Q), color = "black", linetype = "dashed") +  # Add observed data
  coord_cartesian(xlim = x_limits, ylim = c(0,max(c(simulated_plot_long$Q, observed_plot$Q))*1.1), expand = FALSE) +  
  theme_bw() +
  labs(
    x = "Date",
    y = "Streamflow [mm]",
    color = "Legend"
  ) +
  theme(
    plot.margin = margin(5, 5, 5, 5)
  ) +
  scale_x_date(
    breaks = "1 month",  # Show every month
    labels = scales::date_format("%b")  # Show only month names
  ) +
  guides(color = guide_legend(override.aes = list(linetype = c(1,1))))


tau_plot <- tau[Ind_Plot,]

gg2 <- ggplot(tau_plot, aes(x = Dates, y = tau)) +
  geom_line(col="royalblue4")+
  coord_cartesian(xlim = x_limits, ylim = c(0,30), expand = FALSE) +  
  theme_bw()+
  labs(x = "Date", y = expression(tau~"["*days*"]")) +
  theme(
    plot.margin = margin(5, 5, 5, 5)
  ) +
  scale_x_date(
    breaks = "1 month",  # Show every month
    labels = scales::date_format("%b")  # Show only month names
  ) 

gg2

# Combined plot
combined_plot <-hydrograph_plot / gg2

combined_plot

ggsave(plot = combined_plot, filename = "99_Figures/FigE1_TempEvolTau.png",
       width = 8, height = 5, dpi = 300)
