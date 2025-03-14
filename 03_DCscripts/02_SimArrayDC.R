#! ---------------------------------------------------------------------------------------
#!
#! Description       :
#!
#! Authors           : Cyril Thebault <cyril.thebault@ucalgary.ca>
#!
#! Creation date     : 2024-11-12 08:45:48
#! Modification date :
#!
#! Comments          :
#!
#! ---------------------------------------------------------------------------------------

library(ncdf4)

#! ----------------------------- path

FUSE_path <-  "/work/comphyd_lab/users/cyril.thebault/Postdoc_Ucal/02_DATA/FUSE"
BDD_path <- "/home/cyril.thebault/Postdoc_Ucal/02_DATA/BDD"

#! ----------------------------  function

loadRData <- function(file_name) {
  load(file_name)
  get(ls()[ls() != "file_name"])
}

#! ----------------------------  main

catchments <-  unlist(read.table(file.path(BDD_path, "liste_BV_CAMELS_559.txt")))

decisions <- apply(read.table(file.path(FUSE_path,"list_WA_1000.txt"), header = TRUE, sep = ";"), 1, function(row) paste(row, collapse = "_"))

decisionsID = as.character(decisions)

dataset = "CAMELS"
spatialisation = "Lumped"
inputdata = "daymet"
CritCal = "KGE"
TransfoCal="WA"

for(catchment in catchments){
  
  # print("############")
  print(catchment)
  # print("############")
  
  
  for(decision in decisionsID){
    
    # print(decision)
    
    SIM_path = file.path(FUSE_path, dataset, spatialisation, inputdata, CritCal, TransfoCal, catchment, "output")
    myfile = paste0(catchment, "_",decision,".Rdata")
    
    if(!file.exists(file.path(SIM_path, myfile))){next}
    
    Qsim = loadRData(file.path(SIM_path, myfile))$weighted_avg
    Qsim = as.data.frame(Qsim)
    Qsim$Date = seq(as.Date("1999-01-01"),as.Date("2009-12-31"), by = "day")
    if(!exists("Qsim_array")){
      
      DatesR = seq(as.Date("1987-01-01"),as.Date("2009-12-31"), by = "day")
      
      # Create an empty array with dimensions (dates, decisions, catchments)
      Qsim_array <- array(NA, dim = c(length(DatesR), length(decisionsID), length(catchments)),
                          dimnames = list(Dates = as.character(DatesR),
                                          Decisions = decisionsID,
                                          Catchments = catchments))
      
    }
    
    Qsim_array[as.character(Qsim$Date),decision,catchment] = Qsim$Qsim
    
  }
  
}

save(Qsim_array, file = file.path(FUSE_path, dataset, spatialisation, inputdata, CritCal, TransfoCal, "SimArray.Rdata"))