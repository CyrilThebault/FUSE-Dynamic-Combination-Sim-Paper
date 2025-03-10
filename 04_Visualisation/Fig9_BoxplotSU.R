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

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

source(file = "Metrics.R")


##-----------------------------------------
##-------------- VARIABLES ----------------
##-----------------------------------------

dir_FUSE = file.path('00_DATA')

su_HQ = loadRData(file.path(dir_FUSE, "SamplingUncertainty", "SU_HQ.Rdata"))
su_HQ <- su_HQ[order(rownames(su_HQ)), ]
su_HQ$range = su_HQ$p95 - su_HQ$p05

su_LQ = loadRData(file.path(dir_FUSE, "SamplingUncertainty", "SU_LQ.Rdata"))
su_LQ <- su_LQ[order(rownames(su_LQ)), ]
su_LQ$range = su_LQ$p95 - su_LQ$p05

su_WA = loadRData(file.path(dir_FUSE, "SamplingUncertainty", "SU_WA.Rdata"))
su_WA <- su_WA[order(rownames(su_WA)), ]
su_WA$range = su_WA$p95 - su_WA$p05

su_Mosa = loadRData(file.path(dir_FUSE, "SamplingUncertainty", "SU_Mosa.Rdata"))
su_Mosa <- su_Mosa[order(rownames(su_Mosa)), ]
su_Mosa$range = su_Mosa$p95 - su_Mosa$p05

su_Comp = loadRData(file.path(dir_FUSE, "SamplingUncertainty", "SU_Comp.Rdata"))
su_Comp <- su_Comp[order(rownames(su_Comp)), ]
su_Comp$range = su_Comp$p95 - su_Comp$p05


##################################
##           Boxplot            ##
##################################

# Create a dataframe for plotting
data_uncertainty <- data.frame(
  Approach = rep(c("FUSE - KGE(Q)","FUSE - KGE(1/Q)","FUSE - KGEcomp", "Multi-Model Mosaic", "Dynamic Combination"), 
                 times = c(length(na.omit(su_HQ$range)),
                           length(na.omit(su_LQ$range)),
                           length(na.omit(su_Comp$range)), 
                           length(na.omit(su_Mosa$range)), 
                           length(na.omit(su_WA$range)))),
  Range = c(na.omit(su_HQ$range), 
            na.omit(su_LQ$range), 
            na.omit(su_Comp$range), 
            na.omit(su_Mosa$range), 
            na.omit(su_WA$range))
)


data_uncertainty$Approach <- factor(data_uncertainty$Approach, 
                                    levels = c("Dynamic Combination",
                                               "Multi-Model Mosaic",
                                               "FUSE - KGEcomp",
                                               "FUSE - KGE(1/Q)", 
                                               "FUSE - KGE(Q)"))

# Create a boxplot without outliers
gg <- ggplot(data_uncertainty, aes(x = Approach, y = Range, fill = Approach)) +
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
  geom_hline(yintercept = 0, col = "black", linetype = "dashed") +
  scale_fill_manual(values = c("FUSE - KGE(Q)" = "dodgerblue4",
                               "FUSE - KGE(1/Q)" = "forestgreen",
                               "FUSE - KGEcomp" = "darkorange", 
                               "Multi-Model Mosaic" = "lightpink", 
                               "Dynamic Combination" = "firebrick4")) +
  coord_flip(ylim = c(0,0.75))+
  labs(title = "",
       y = "KGEcomp uncertainty range (p95 - p05 from gumboot)", 
       x = "") +
  theme_bw() +  
  theme(
    legend.position = "none"
  )


ggsave(plot = gg, filename = "99_Figures/Fig9_BoxplotSU.png",
       width = 10, height = 3, dpi = 300)
