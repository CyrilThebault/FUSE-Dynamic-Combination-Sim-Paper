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

library(ggpubr)

source(file = "Metrics.R")

##-----------------------------------------
##---------------- MAIN ------------------
##-----------------------------------------

# WA
dir_FUSE = file.path('00_DATA')
Eval_FUSE_long_WA = loadRData(file.path(dir_FUSE, "WA", "Eval_FUSE.Rdata"))
row_names_FUSE_WA <- as.character(unique(Eval_FUSE_long_WA$Codes))

Eval_FUSE_WA_LQ <- data.frame(Eval_FUSE_long_WA %>%
                             dplyr::select(Codes, ModelDecisions, `Eval : KGE : -1`) %>%  # Select only the relevant columns
                             pivot_wider(names_from = ModelDecisions, values_from = `Eval : KGE : -1`)) %>% 
  column_to_rownames(var = "Codes")

Eval_FUSE_WA_HQ <- data.frame(Eval_FUSE_long_WA %>%
                                dplyr::select(Codes, ModelDecisions, `Eval : KGE :  1`) %>%  # Select only the relevant columns
                                pivot_wider(names_from = ModelDecisions, values_from = `Eval : KGE :  1`)) %>% 
  column_to_rownames(var = "Codes")



# Prepare the data
df_median_LQ <- Eval_FUSE_WA_LQ %>%
  pivot_longer(cols = everything(), names_to = "param_combo", values_to = "value") %>%
  separate(param_combo, into = c("Param1", "Param2", "Param3"), sep = "_") %>%
  mutate(across(starts_with("Param"), ~ as.numeric(gsub("X", "", .)))) %>%
  group_by(Param1, Param2, Param3) %>%
  summarize(median_value = median(value, na.rm = TRUE), .groups = "drop")%>%
  mutate(Param1 = factor(Param1, labels = paste("\u03C4 =", unique(Param1), " days")))


df_median_HQ <- Eval_FUSE_WA_HQ %>%
  pivot_longer(cols = everything(), names_to = "param_combo", values_to = "value") %>%
  separate(param_combo, into = c("Param1", "Param2", "Param3"), sep = "_") %>%
  mutate(across(starts_with("Param"), ~ as.numeric(gsub("X", "", .)))) %>%
  group_by(Param1, Param2, Param3) %>%
  summarize(median_value = median(value, na.rm = TRUE), .groups = "drop")%>%
  mutate(Param1 = factor(Param1, labels = paste("\u03C4 =", unique(Param1), " days")))



# Create gg1 and gg2 as before
gg1 = ggplot(df_median_HQ, aes(x = Param2, y = Param3, fill = median_value)) +   
  geom_tile() +  
  facet_wrap(~ Param1, ncol = 5) +  
  scale_fill_gradientn(name = "Median\nKGE(Q)", 
                       colors = viridis::viridis(100), 
                       limits = c(0.65, 0.75),
                       breaks = c(0.65,0.70,0.75),
                       labels = c("< 0.65", "0.70", "0.75"),
                       oob = scales::squish,
                       guide = guide_colorbar(barwidth = 1, barheight = 6)) +  # Custom scale with limits
  labs(x = "m", y = "k", title = "a)") +   
  scale_x_continuous(breaks = seq(1, 19, by = 2), expand = c(0, 0)) +   
  scale_y_continuous(breaks = seq(1, 19, by = 2), expand = c(0, 0)) +   
  theme_bw() +   
  theme(strip.text = element_text(size = 8), axis.text.x = element_text(angle = 45, hjust = 1))

gg2 = ggplot(df_median_LQ, aes(x = Param2, y = Param3, fill = median_value)) +   
  geom_tile() +  
  facet_wrap(~ Param1, ncol = 5) +  
  scale_fill_gradientn(name = "Median\nKGE(1/Q)", 
                       colors = viridis::viridis(100), 
                       limits = c(0.57, 0.67),
                       breaks = c(0.57,0.62,0.67),
                       labels = c("< 0.57", "0.62", "0.67"),
                       oob = scales::squish,
                       guide = guide_colorbar(barwidth = 1, barheight = 6)) +  # Custom scale with limits
  labs(x = "m", y = "k", title = "b)") +   
  scale_x_continuous(breaks = seq(1, 19, by = 2), expand = c(0, 0)) +   
  scale_y_continuous(breaks = seq(1, 19, by = 2), expand = c(0, 0)) +   
  theme_bw() +   
  theme(strip.text = element_text(size = 8), axis.text.x = element_text(angle = 45, hjust = 1))



# Align the legends

gg = ggarrange(gg1, gg2, ncol = 1, nrow = 2)

ggsave(plot = gg, filename = "99_Figures/Fig6_ParamDC.png",
       width = 9, height = 8, dpi = 300)
