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

EvalStart = "1999-01-01"
EvalEnd = "2009-12-31"

DatesR <- seq(from = as.Date(EvalStart),
              to = as.Date(EvalEnd), 
              by = 'day')


# Dynamic combination outputs
WA_RES = loadRData(file.path(dir_FUSE, "WA", paste0(catchment,"_16_15_19.Rdata"))) 

df_weights = WA_RES$weights
colnames(df_weights) = 1:156
rownames(df_weights) = as.character(DatesR)

simulated = data.frame(Dates = DatesR, Q = WA_RES$weighted_avg)

# Obs
Obs = loadRData(file.path(dir_FUSE,"ObservedStreamflow.Rdata"))$Qobs_mm
Obs = Obs[as.Date(rownames(Obs)) %in% DatesR,]
observed = data.frame(Dates = DatesR, Q = Obs[,catchment])

#------------- Plot

Ind_Plot = which(format(DatesR, "%Y") == "2003")

# Extract column names to determine flow type
flow_type <- ifelse(as.numeric(colnames(df_weights))<=78, "highflow", "lowflow")

# Melt the matrix to long format
weights_melt <- melt(df_weights[Ind_Plot,])

# Map colors based on flow type
weights_melt$color <- ifelse(
  flow_type[weights_melt$Var2] == "highflow",
  ifelse(weights_melt$value > 0, "dodgerblue4", NA),
  ifelse(weights_melt$value > 0, "forestgreen", NA)
)

simulated_plot = simulated[Ind_Plot,]
observed_plot = observed[Ind_Plot,]



x_limits = range(simulated_plot$Dates)
y_limits <- c(1, ncol(df_weights))

# Adjust the heatmap plot
heatmap_plot <- ggplot(weights_melt, aes(x = as.Date(Var1), y = Var2, fill = color)) +
  geom_tile() +
  scale_fill_identity() +
  coord_cartesian(xlim = x_limits, ylim = y_limits, expand = FALSE) +  
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )+
  labs(
    title = "a)",
    y = "Model Decision"
  ) +
  scale_x_date(
    breaks = "1 month",  # Show every month
    labels = scales::date_format("%b")  # Show only month names
  )

# Adjust the hydrograph plot
hydrograph_plot <- ggplot() +
  geom_line(data = simulated_plot, aes(x = Dates, y = Q), color = "firebrick4") + 
  geom_line(data = observed_plot, aes(x = Dates, y = Q), color = "black", linetype = "dashed") +  # Add observed data
  coord_cartesian(xlim = x_limits, ylim = c(0,15), expand = FALSE) +  
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

# Combine the plots using patchwork
gg1 <- ggarrange(heatmap_plot, hydrograph_plot,
                           nrow = 2, heights = c(2, 1)) # Adjust height ratio if needed


##-----------------------------------------
##---------------- PART 2 -----------------
##-----------------------------------------

# Weights
myWeights = loadRData(file.path(dir_FUSE, "WA","Weights_16_15_19.Rdata")) 

Obs$Date <- rownames(Obs)
#------------- Plot

# Step 1: Compute Percentile Bins for Observed Streamflow
percentile_bins <- seq(0, 100, by = 1)  # 100 bins (0-1%, 1-2%, ... 99-100%)

# Convert myWeights to a data frame
weights_df <- as.data.frame.table(myWeights, responseName = "Weight")
colnames(weights_df) <- c("Catchment", "Date", "Decision", "Weight")
weights_df <- weights_df[weights_df$Weight != 0,] #Filter out zero weights (avoid counting irrelevant entries)

# Convert Obs to long format
obs_long <- melt(Obs, id.vars = "Date", variable.name = "Catchment", value.name = "Streamflow")

# Step 2: Assign each streamflow value to a percentile bin
obs_long <- obs_long %>%
  group_by(Catchment) %>%
  mutate(
    Percentile = {
      # Compute quantiles
      breaks <- unique(quantile(Streamflow, probs = percentile_bins / 100, na.rm = TRUE))
      
      # Ensure there are at least 2 unique breaks
      if (length(breaks) < 2) {
        breaks <- c(min(Streamflow, na.rm = TRUE), max(Streamflow, na.rm = TRUE) + 1)
      }
      
      # Assign percentiles
      cut(Streamflow, breaks = breaks, include.lowest = TRUE, labels = percentile_bins[2:length(breaks)])
    }
  )



# Step 3: Merge Observations with Weights
# Convert to data.table for faster operations
obs_long <- as.data.table(obs_long)
weights_df <- as.data.table(weights_df)

# Ensure Date is already in Date format
obs_long[, Date := as.character(Date)]
weights_df[, Date := as.character(Date)]

# Set keys for faster merging
setkey(obs_long, Catchment, Date)
setkey(weights_df, Catchment, Date)

# Perform optimized join
merged_data <- weights_df[obs_long, nomatch = 0]

# Step 4: Summarize HF & LF Proportions by Percentile
merged_data$Category <- ifelse(grepl("HF", merged_data$Decision), "HF", "LF")

merged_data <- as.data.table(merged_data)


percentile_summary <- merged_data[, .(
  `High-Flow` = sum(Category == "HF") / .N,
  `Low-Flow` = sum(Category == "LF") / .N
), by = Percentile]

percentile_summary_long <- melt(percentile_summary, id.vars = "Percentile", 
                                variable.name = "Flow_Category", value.name = "Proportion")


percentile_summary_long$Flow_Category <- factor(percentile_summary_long$Flow_Category, 
                                                levels = rev(levels(factor(percentile_summary_long$Flow_Category))))

percentile_summary_long = percentile_summary_long[!is.na(percentile_summary_long$Percentile),]

percentile_summary_long$Proportion = as.numeric(percentile_summary_long$Proportion)
percentile_summary_long$Percentile = as.numeric(percentile_summary_long$Percentile)

percentile_summary_long$Flow_Category <- factor(percentile_summary_long$Flow_Category, 
                                                levels = c("High-Flow", "Low-Flow"))

gg2 <- ggplot(percentile_summary_long, aes(x = Percentile, y = Proportion, fill = Flow_Category)) +
  geom_bar(stat = "identity", position = "stack", width = 1.1) +
  scale_fill_manual(values = c("High-Flow" = "#7998BB", "Low-Flow" = "#81B984"), 
                    limits = c("High-Flow", "Low-Flow")) +  # Set legend order
  scale_x_continuous(limits = c(1, 100), breaks = seq(0, 100, by = 5), expand = c(0, 0)) +  
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0),  
                     sec.axis = sec_axis(~ 1 - ., name = "Proportion of models calibrated\non high-flow - KGE(Q)")) +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.title.y = element_text(color = "forestgreen", size = 12),  
    axis.text.y = element_text(color = "forestgreen"),  
    axis.title.y.right = element_text(color = "dodgerblue4", size = 12),  
    axis.text.y.right = element_text(color = "dodgerblue4")  
  ) +
  labs(title = "b)",
       x = "Streamflow Percentile",
       y = "Proportion of models calibrated\non low-flow - KGE(1/Q)",
       fill = "")


##-----------------------------------------
##-------------  Combined -----------------
##-----------------------------------------

combined_plot <- ggarrange(gg1,gg2,nrow = 2)


ggsave(plot = combined_plot, filename = "99_Figures/Fig12_TempEvol.png",
       width = 8, height = 8, dpi = 300)
