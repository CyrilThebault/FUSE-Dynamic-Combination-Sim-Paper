# Define path to base directory
base_path <- "/work/comphyd_lab/users/cyril.thebault/Postdoc_Ucal/02_DATA/FUSE/CAMELS/Lumped/daymet/KGE/WA"

paramset_name <- "16_15_19"
  
# Get list of catchments
catchments <- list.dirs(base_path, recursive = FALSE, full.names = FALSE)

# Initialize an empty list to store results
catchment_data <- list()

# Iterate through catchments
for (catchment in catchments) {
  print(catchment)
  output_path <- file.path(base_path, catchment, "output")
  
  # Check if output directory exists
  if (!dir.exists(output_path)) next
  
  # Get list of paramset files
  paramset_file <- list.files(output_path, pattern = paste0(paramset_name,"\\.Rdata$"), full.names = TRUE)
  

  load(paramset_file)
    
  catchment_data[[catchment]] <- simulated_WA$weights
  
}


final_array <- array(NA, dim = c(length(catchments), 4018, 156))

# Fill the array
for (i in seq_along(catchments)) {
  catchment <- catchments[i]
  if (catchment %in% names(catchment_data)) {
    final_array[i, , ] <- catchment_data[[catchment]]
  }
}

EvalStart = "1999-01-01"
EvalEnd = "2009-12-31"

DatesR <- seq(from = as.Date(EvalStart),
              to = as.Date(EvalEnd),
              by = 'day')

models = read.table(paste0(base_path, "/../../../../../list_decision_78.txt"), sep = ";", header = TRUE)[,1]
models_tot = c(paste0("HF_", models), paste0("LF_", models))


catchment_names = names(catchment_data)

dimnames(final_array) <- list("Codes" = catchment_names, "Dates" = as.character(DatesR), "Decisions" = models_tot)


# Save the array and catchment names for reference
save(final_array, file = file.path(base_path,paste0("Weights_",paramset_name,".Rdata")))

print("Processing complete! Data saved as catchment_weights.Rdata")
