#! ---------------------------------------------------------------------------------------
#!
#! Description       :
#!
#! Authors           : Cyril Thebault <cyril.thebault@ucalgary.ca>
#!
#! Creation date     : 2024-11-02 18:14:16
#! Modification date :
#!
#! Comments          :
#!
#! ---------------------------------------------------------------------------------------

#! ----------------------------- package loading


loadRData <- function(file_name) {
  load(file_name)
  get(ls()[ls() != "file_name"])
}

library(abind)

FUSE_path = "/work/comphyd_lab/users/cyril.thebault/Postdoc_Ucal/02_DATA/FUSE"

dataset = "CAMELS"
spatialisation = "Lumped"
inputdata = "daymet"
CritCal = "KGE"

# Sampling Uncertainty
df_gumboot = read.table(file = file.path(FUSE_path,dataset,spatialisation,inputdata,CritCal,"SamplingUncertainty", 
                                         "brief_analysis_modeling_results_with_gumboot_eval.csv"),
                        header = TRUE,
                        sep = ",")
df_gumboot$Code = ifelse(nchar(df_gumboot$gauge_id)==7, paste0("USA_0", df_gumboot$gauge_id), paste0("USA_", df_gumboot$gauge_id))

# Simulations
QsimHF = loadRData(file.path(FUSE_path,dataset,spatialisation,inputdata,CritCal,"1", "SimArray.Rdata"))
QsimLF = loadRData(file.path(FUSE_path,dataset,spatialisation,inputdata,CritCal,"-1", "SimArray.Rdata"))

Qsim <- abind(QsimHF, QsimLF, along = 4) # Combine along the 4th dimension

dimnames(Qsim)[[4]] <- c("1", "-1")

rm(QsimHF,QsimLF)

# Obervations
Qobs = loadRData(file = file.path(FUSE_path,dataset,spatialisation,inputdata,CritCal,"1","ObservedStreamflow.Rdata"))$Qobs

# Evaluations
EvalHF = loadRData(file.path(FUSE_path,dataset,spatialisation,inputdata,CritCal,"1", "Eval_FUSE.Rdata"))
EvalLF = loadRData(file.path(FUSE_path,dataset,spatialisation,inputdata,CritCal,"-1", "Eval_FUSE.Rdata"))


CalStart = "1989-01-01"
CalEnd = "1998-12-31"

EvalStart = "1999-01-01"
EvalEnd = "2009-12-31"

time_windows = seq(1,30,3)
k_models = seq(1,20,2)
k_neighbours = seq(1,20,2)


DatesR = as.Date(dimnames(Qsim)[[1]])
models = dimnames(Qsim)[[2]]
catchments = dimnames(Qsim)[[3]]
transfo = dimnames(Qsim)[[4]]

mydf = expand.grid(catchments, time_windows, k_models, k_neighbours)
colnames(mydf) = c("Code", "Tau", "Km", "Kn")
mydf = mydf[order(mydf$Code),]

############

args <- commandArgs(trailingOnly = TRUE)
id <- as.numeric(args[1])

catchment = mydf$Code[id]
tau = mydf$Tau[id]
km = mydf$Km[id]
knn = mydf$Kn[id]


Ind_Cal = which(DatesR >= CalStart & DatesR <= CalEnd)

Ind_Eval = which(DatesR >= EvalStart & DatesR <= EvalEnd)

# Subset simulated
simulatedHF = Qsim[,models,catchment, "1"]
colnames(simulatedHF) = paste0("HF_",colnames(simulatedHF))

simulatedLF = Qsim[,models,catchment, "-1"]
colnames(simulatedLF) = paste0("LF_",colnames(simulatedLF))

simulated = cbind(simulatedHF, simulatedLF)

# Subset observed
qmm = Qobs[,catchment]

# Subset evaluation
evaluationHF = EvalHF[EvalHF$ModelDecisions %in% models & EvalHF$Codes %in% catchment,]
evaluationHF$ModelDecisions <- ifelse(grepl("^HF_", evaluationHF$ModelDecisions), 
                                      evaluationHF$ModelDecisions, 
                                      paste0("HF_", evaluationHF$ModelDecisions))

evaluationLF = EvalLF[EvalLF$ModelDecisions %in% models & EvalLF$Codes %in% catchment,]
evaluationLF$ModelDecisions <- ifelse(grepl("^HF_", evaluationLF$ModelDecisions), 
                                      evaluationLF$ModelDecisions, 
                                      paste0("HF_", evaluationLF$ModelDecisions))

evaluation = rbind(evaluationHF, evaluationLF)

# Select best model over calibration period

best_model <- df_gumboot$mod_need[df_gumboot$Code == catchment]

# Set up initial weights (equal weights)
init_weights <- rep(1 / ncol(simulated), ncol(simulated))

# Create a matrix to store the dynamic weights for each time step and each simulation
dynamic_weights <- matrix(NA, nrow = nrow(simulated), ncol = ncol(simulated))

# Create a vector to store the weighted average for each time step
dynamic_weighted_avg <- numeric(nrow(simulated))

# Initialize the progress bar
# pb <- txtProgressBar(min = 0, max = length(Ind_Eval), style = 3)
for (t in Ind_Eval) {
  
  
  window_indices <- (t - tau):(t-1)
  
  
  simulated_window <- simulated[window_indices, , drop = FALSE]

  best_window <- simulated[window_indices, best_model]
  d_best_window <- diff(best_window)
  fv_best_window <- c(best_window, d_best_window)
  
  
  # Initialize a vector to store MAE windows and models results
  mae_windows <- rep(NA, length(Ind_Cal))
  mae_models <- matrix(NA, nrow = length(Ind_Cal), ncol = ncol(simulated))
  
  # Loop through the calibration period in windows
  inc = 0
  for (i in Ind_Cal) {
    inc = inc +1
    # Extract the current window of observed data
    simulated_window_cal <- simulated[(i-tau):(i-1), , drop = FALSE]
    observed_window_cal <- qmm[(i-tau):(i-1)]
    
    mae_models[inc,] <- colMeans(abs(simulated_window_cal - observed_window_cal))
    
    d_observed_window_cal <- diff(observed_window_cal)
    fv_observed_window_cal <- c(observed_window_cal, d_observed_window_cal)
    
    # Calculate MAE between the combined observed data and the current simulated window
    mae_windows[inc] <- mean(abs(fv_best_window - fv_observed_window_cal))
    
  }
  
  
  
  best_indices_windows <- order(mae_windows)[1:knn]
  
  window_weights <- numeric(length(mae_windows))
  window_weights[best_indices_windows] <- 1 / knn
  
  sub_mae_models = mae_models[best_indices_windows,]
  
  if(length(best_indices_windows) > 1){
    sub_mae_models_mean = colMeans(sub_mae_models)
  } else{
    sub_mae_models_mean = sub_mae_models
  }
  
  best_indices_mod <- order(sub_mae_models_mean)[1:km]
  
  mod_weight <- numeric(ncol(mae_models))
  mod_weight[best_indices_mod] <- 1 / km
  
  
  # Store the dynamic weights for the current time step
  dynamic_weights[t, ] <- mod_weight
  
  # Calculate the weighted average for the current time step using the dynamic weights
  dynamic_weighted_avg[t] <- sum(simulated[t, ] * mod_weight, na.rm = TRUE)
  
}

simulated_WA <- list(
  weighted_avg = dynamic_weighted_avg[Ind_Eval],
  weights = dynamic_weights[Ind_Eval,]
)


#####################
# Save results
#####################

dir_out = file.path(FUSE_path,dataset,spatialisation,inputdata,CritCal, "WA", catchment, "output")
if(! dir.exists(dir_out)){dir.create(dir_out, recursive = TRUE)}

save(simulated_WA, file = paste0(dir_out,"/", catchment,"_",tau,"_",km,"_",knn,".Rdata" ))

print(paste0("Saved file: ",dir_out,"/", catchment,"_",tau,"_",km,"_",knn,".Rdata" ))
