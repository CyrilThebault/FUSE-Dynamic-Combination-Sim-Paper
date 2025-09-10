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

outputfolder = "HQLQ_WA_q"

# Simulations
QsimHF = loadRData(file.path(FUSE_path,dataset,spatialisation,inputdata,CritCal,"1", "SimArray.Rdata"))
QsimLF = loadRData(file.path(FUSE_path,dataset,spatialisation,inputdata,CritCal,"-1", "SimArray.Rdata"))

Qsim <- abind(QsimHF, QsimLF, along = 4) # Combine along the 4th dimension

dimnames(Qsim)[[4]] <- c("1", "-1")

rm(QsimHF,QsimLF)

# Obervations
Qobs = loadRData(file = file.path(FUSE_path,dataset,spatialisation,inputdata,CritCal,"1","ObservedStreamflow.Rdata"))$Qobs


SearchStart = "1989-01-01"
SearchEnd = "1993-12-31"

SelecStart = "1994-01-01"
SelecEnd = "1998-12-31"

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


Ind_Search = which(DatesR >= SearchStart & DatesR <= SearchEnd)
Ind_Selec = which(DatesR >= SelecStart & DatesR <= SelecEnd)
Ind_Eval = which(DatesR >= EvalStart & DatesR <= EvalEnd)

# Subset simulated
simulatedHF = Qsim[,models,catchment, "1"]
colnames(simulatedHF) = paste0("HF_",colnames(simulatedHF))

simulatedLF = Qsim[,models,catchment, "-1"]
colnames(simulatedLF) = paste0("LF_",colnames(simulatedLF))

simulated = cbind(simulatedHF, simulatedLF)

# Subset observed
qmm = Qobs[,catchment]


# Create a matrix to store the dynamic weights for each time step and each simulation
dynamic_weights <- matrix(NA, nrow = nrow(simulated), ncol = ncol(simulated))

# Create a vector to store the weighted average for each time step
dynamic_weighted_avg <- rep(NA,nrow(simulated))

# Loop to create the dynamic combination
for (t in c(Ind_Selec, Ind_Eval)) {
  
  #current window
  window_indices <- (t - tau):(t-1)
  
  simulated_window <- simulated[window_indices, , drop = FALSE]
  observed_window <- qmm[window_indices]
  dynamic_window <- dynamic_weighted_avg[window_indices]
  
  best_window <- dynamic_window
  ind_na = is.na(best_window)
  
  # Warm-up of the dynamic combination:
  # Use observation (or median of simulation if not available).
  # Only on the first Tau time steps can be NA (otherwise dynamic combination is continuous)
  if(any(ind_na)){
    
    best_window[ind_na] = observed_window[ind_na] 
    
    ind_na = is.na(best_window)
    if(any(ind_na)){
      if(sum(ind_na) == 1){
        best_window[ind_na] = median(simulated[window_indices[ind_na], ], na.rm = TRUE)
      } else {
        best_window[ind_na] = apply(simulated[window_indices[ind_na], ], 1, median, na.rm = TRUE)
      }
      
    }
  }
  
  #Create current condition feature vector
  fv_best_window <- c(best_window)
  
  # Initialize a vector to store MAE windows and models results
  mae_windows <- rep(NA, length(Ind_Search))
  mae_models <- matrix(NA, nrow = length(Ind_Search), ncol = ncol(simulated))
  
  # Loop through the search period
  inc = 0
  for (i in Ind_Search) {
    inc = inc +1
    # Extract the window of observed data in the search period (for time step i)
    simulated_window_cal <- simulated[(i-tau):(i-1), , drop = FALSE]
    observed_window_cal <- qmm[(i-tau):(i-1)]
    
    #Create observed feature vector in the search period (for time step i)
    fv_observed_window_cal <- c(observed_window_cal)
    
    # Calculate MAE between the observed window (for time step i) and the current simulated window
    mae_windows[inc] <- mean(abs(fv_best_window - fv_observed_window_cal))
    
    # Compute MAE for every models in the search period (for time step i)
    mae_models[inc,] <- colMeans(abs(simulated_window_cal - observed_window_cal))
    
  }
  
  # select the nearest neighbors
  best_indices_windows <- order(mae_windows)[1:knn]
  
  # attribute equal weights to the selected neighbors
  window_weights <- numeric(length(mae_windows))
  window_weights[best_indices_windows] <- 1 / knn
  
  # subset the MAE models matrix to check performance within the selected neighbors
  sub_mae_models = mae_models[best_indices_windows,]
  
  if(length(best_indices_windows) > 1){
    sub_mae_models_mean = colMeans(sub_mae_models)
  } else{
    sub_mae_models_mean = sub_mae_models
  }
  
  # select the best models
  best_indices_mod <- order(sub_mae_models_mean)[1:km]
  
  # attribute equal weights to the selected models
  mod_weight <- numeric(ncol(mae_models))
  mod_weight[best_indices_mod] <- 1 / km
  
  
  # Store the dynamic weights for the current time step
  dynamic_weights[t, ] <- mod_weight
  
  # Calculate the weighted average for the current time step using the dynamic weights
  dynamic_weighted_avg[t] <- sum(simulated[t, ] * mod_weight, na.rm = TRUE)

}

# Store results in a list
simulated_WA <- list(
  dates = DatesR,
  weighted_avg = dynamic_weighted_avg,
  weights = dynamic_weights,
  searchPeriod = c(SearchStart, SearchEnd),
  runPeriod = c(SelecStart, EvalEnd)
)


#####################
# Save results
#####################

dir_out = file.path(FUSE_path,dataset,spatialisation,inputdata,CritCal, outputfolder, catchment, "output")
if(! dir.exists(dir_out)){dir.create(dir_out, recursive = TRUE)}

save(simulated_WA, file = paste0(dir_out,"/", catchment,"_",tau,"_",km,"_",knn,".Rdata" ))

print(paste0("Saved file: ",dir_out,"/", catchment,"_",tau,"_",km,"_",knn,".Rdata" ))
