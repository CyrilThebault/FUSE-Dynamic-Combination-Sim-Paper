
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

#! ----------------------------- package loading

library(abind)
library(ncdf4)

source(file = "Metrics.R")

#! ---------------------------- bassins

FUSE_path = "/Users/cyrilthebault/Postdoc_Ucal/02_DATA/FUSE/01_Paper"
models = read.table(file.path(FUSE_path, "list_decision_78.txt"), sep = ";", header = TRUE)[,1]

QsimHQ = loadRData(file.path(FUSE_path,"1", "SimArray.Rdata"))
QsimLQ = loadRData(file.path(FUSE_path,"-1", "SimArray.Rdata"))

# Combine along the 4th dimension
Qsim <- abind(QsimHQ, QsimLQ, along = 4)
dimnames(Qsim)[[4]] <- c("1", "-1")

rm(QsimHQ,QsimLQ)

DatesR = as.Date(dimnames(Qsim)[[1]])
models = dimnames(Qsim)[[2]]
catchments = dimnames(Qsim)[[3]]
transfo = dimnames(Qsim)[[4]]

SimStart = "1989-01-01"
SimEnd = "2009-12-31"

DatesSim <- seq(from = as.Date(SimStart),
              to = as.Date(SimEnd), 
              by = 'day')

IDNumbers = as.numeric(substr(catchments, 5, nchar(catchments)))

##-----------------------------------------
##---------------- FUSE -------------------
##-----------------------------------------

for(trans in transfo){
    
  for(mod in models){
    
    myQsim = array(NA, dim = c(1, length(catchments), length(DatesSim)))
    myQsim[1,,] = t(Qsim[DatesR %in% DatesSim, mod, ,trans])
    
    
    dir_nc = file.path(FUSE_path, "SamplingUncertainty")
    prefix = ifelse(trans == "1", "HF", "LF")
    
    ncfname <- paste0(dir_nc, "/",prefix,"_",mod,"_fuse",prefix,mod,"_modelResults.nc")
    
    # define dimensions
    hrudim <- ncdim_def(name = "catchment_index", units = "", vals = 1:length(catchments), create_dimvar=FALSE)
    pardim <- ncdim_def(name = "par_index", units = "", vals = 1:6, create_dimvar=FALSE)
    objdim <- ncdim_def(name = "obj_fn_index", units = "", vals = 1:1, create_dimvar=FALSE)
    timedim <-  ncdim_def(name = "time", units = "", vals = 1:length(DatesSim), create_dimvar=FALSE)
    
    # define variables
    CalParam_def <- ncvar_def(name = "Calibrated_parameters", 
                              units = "See FUSE documentation",
                              dim = list(objdim, pardim, hrudim),
                              longname = "Parameter values after calibration",
                              prec = "double"
    )
    
    GaugeID_def <- ncvar_def(name = "Gauge_ID", 
                           units = "n/a",
                           dim = hrudim,
                           longname = "Gauge ID as used by USGS and in the CAMELS-US dataset",
                           prec = "double"
    )
    
    ObjFunCal_def <- ncvar_def(name = "Objective_function_cal", 
                               units = "[-]",
                               dim = list(objdim, hrudim),
                               longname = "Objective function obtained during model calibration",
                               prec = "double"
    )
    
    ObjFunEval_def <- ncvar_def(name = "Objective_function_eval", 
                                units = "[-]",
                                dim = list(objdim, hrudim),
                                longname = "Objective function obtained during model evaluation",
                                prec = "double"
    )
    
    SimEa_def <- ncvar_def(name = "Sim_ea", 
                            units = "[mm/d]",
                            dim = list(objdim, timedim, hrudim),
                            longname = "Time series of simulated total evaporation",
                            prec = "double"
    )
    
    SimFluxDinf_def <- ncvar_def(name = "Sim_flux_dinf", 
                           units = "[mm/d]",
                           dim = list(objdim, timedim, hrudim),
                           longname = "Time series of simulated flux dinf",
                           prec = "double"
    )
    
    SimFluxExc_def <- ncvar_def(name = "Sim_flux_exc", 
                           units = "[mm/d]",
                           dim = list(objdim, timedim, hrudim),
                           longname = "Time series of simulated flux exc",
                           prec = "double"
    )
    
    SimFluxFlow_def <- ncvar_def(name = "Sim_flux_flow", 
                           units = "[mm/d]",
                           dim = list(objdim, timedim, hrudim),
                           longname = "Time series of simulated flux flow",
                           prec = "double"
    )
    
    SimFluxGwf_def <- ncvar_def(name = "Sim_flux_gwf", 
                           units = "[mm/d]",
                           dim = list(objdim, timedim, hrudim),
                           longname = "Time series of simulated flux gwf",
                           prec = "double"
    )
    
    SimFluxInf_def <- ncvar_def(name = "Sim_flux_inf", 
                           units = "[mm/d]",
                           dim = list(objdim, timedim, hrudim),
                           longname = "Time series of simulated flux inf",
                           prec = "double"
    )
    
    SimFluxInt_def <- ncvar_def(name = "Sim_flux_int", 
                           units = "[mm/d]",
                           dim = list(objdim, timedim, hrudim),
                           longname = "Time series of simulated flux int",
                           prec = "double"
    )
    
    SimFluxRec_def <- ncvar_def(name = "Sim_flux_rec", 
                           units = "[mm/d]",
                           dim = list(objdim, timedim, hrudim),
                           longname = "Time series of simulated flux rec",
                           prec = "double"
    )
    
    SimFluxRun_def <- ncvar_def(name = "Sim_flux_run", 
                           units = "[mm/d]",
                           dim = list(objdim, timedim, hrudim),
                           longname = "Time series of simulated flux run",
                           prec = "double"
    )
    
    SimFluxSeep_def <- ncvar_def(name = "Sim_flux_seep", 
                           units = "[mm/d]",
                           dim = list(objdim, timedim, hrudim),
                           longname = "Time series of simulated flux seep",
                           prec = "double"
    )
    
    SimFluxSmf_def <- ncvar_def(name = "Sim_flux_smf", 
                           units = "[mm/d]",
                           dim = list(objdim, timedim, hrudim),
                           longname = "Time series of simulated flux smf",
                           prec = "double"
    )
    
    SimFluxSrun_def <- ncvar_def(name = "Sim_flux_srun", 
                           units = "[mm/d]",
                           dim = list(objdim, timedim, hrudim),
                           longname = "Time series of simulated flux srun",
                           prec = "double"
    )
    
    SimFluxTrap_def <- ncvar_def(name = "Sim_flux_trap", 
                           units = "[mm/d]",
                           dim = list(objdim, timedim, hrudim),
                           longname = "Time series of simulated flux trap",
                           prec = "double"
    )
    
    SimQ_def <- ncvar_def(name = "Sim_q", 
                           units = "[mm/d]",
                           dim = list(objdim, hrudim, timedim),
                           longname = "Time series of simulated flows",
                           prec = "double"
    )
    
    SimStoreS1_def <- ncvar_def(name = "Sim_store_S1", 
                           units = "[mm]",
                           dim = list(objdim, timedim, hrudim),
                           longname = "Time series of simulated store S1",
                           prec = "double"
    )
    
    SimStoreS2_def <- ncvar_def(name = "Sim_store_S2", 
                           units = "[mm]",
                           dim = list(objdim, timedim, hrudim),
                           longname = "Time series of simulated store S2",
                           prec = "double"
    )
    
    SimStoreS3_def <- ncvar_def(name = "Sim_store_S3", 
                           units = "[mm]",
                           dim = list(objdim, timedim, hrudim),
                           longname = "Time series of simulated store S3",
                           prec = "double"
    )
    
    SimStoreS4_def <- ncvar_def(name = "Sim_store_S4", 
                           units = "[mm]",
                           dim = list(objdim, timedim, hrudim),
                           longname = "Time series of simulated store S4",
                           prec = "double"
    )
    
    SimStoreS5_def <- ncvar_def(name = "Sim_store_S5", 
                           units = "[mm]",
                           dim = list(objdim, timedim, hrudim),
                           longname = "Time series of simulated store S5",
                           prec = "double"
    )
    
    SimStoreWarmup_def <- ncvar_def(name = "Store_warmup", 
                           units = "[-]",
                           dim = list(objdim, hrudim),
                           longname = "Number of times year 1 needed to be iterated before storage values stabilized (<1% difference) with the calibrated parameter set",
                           prec = "double"
    )
    
    
    ncout <- nc_create(ncfname,list(CalParam_def, GaugeID_def, ObjFunCal_def, ObjFunEval_def, SimEa_def, SimFluxDinf_def, SimFluxExc_def, SimFluxFlow_def, SimFluxGwf_def, SimFluxInf_def, SimFluxInt_def, SimFluxRec_def, SimFluxRun_def, SimFluxSeep_def, SimFluxSmf_def, SimFluxSrun_def, SimFluxTrap_def, SimQ_def, SimStoreS1_def, SimStoreS2_def, SimStoreS3_def, SimStoreS4_def, SimStoreS5_def, SimStoreWarmup_def), force_v4 = TRUE)
    
    # put variables
    # ncvar_put(ncout,CalParam_def, Param)
    ncvar_put(ncout,GaugeID_def, IDNumbers)
    # ncvar_put(ncout,ObjFunCal_def,EvalCal)
    # ncvar_put(ncout,ObjFunEval_def,EvalVal)
    ncvar_put(ncout,SimQ_def,myQsim)
    
    
    # add global attributes
    ncatt_put(ncout,0,"author", "Cyril Thébault")
    ncatt_put(ncout,0,"date",paste0("Created ", format(Sys.time(), format = "%Y/%m/%d %H:%M:%S")))
    ncatt_put(ncout,0,"institution",'University of Calgary')
    ncatt_put(ncout,0,"model_name", paste0('f_',mod,'_fuse',mod))
    ncatt_put(ncout,0,"model_number_stores",'see model decision')
    ncatt_put(ncout,0,"model_number_parameters",'see model decision')
    ncatt_put(ncout,0,"model_parameter_ranges_min",'na')
    ncatt_put(ncout,0,"model_parameter_ranges_max",'na')
    ncatt_put(ncout,0,"solver_name",'SCE')
    ncatt_put(ncout,0,"solver_resnorm_tolerance",'nq')
    ncatt_put(ncout,0,"solver_resnorm_maxiter",'nq')
    ncatt_put(ncout,0,"solver_delta_t_[d]",'nq')
    ncatt_put(ncout,0,"objective_function",ifelse(trans == "1", 'KGE(Q)', 'KGE(1/Q'))
    ncatt_put(ncout,0,"time_calibration_start",'1989-01-01')
    ncatt_put(ncout,0,"time_calibration_end",'1998-12-31')
    ncatt_put(ncout,0,"time_evaluation_start",'1999-01-01')
    ncatt_put(ncout,0,"time_evaluation_end",'2009-12-31')
    ncatt_put(ncout,0,"time_step_size",'1 [d]')
    ncatt_put(ncout,0,"algorithm",'na')
    
    
    
    nc_close(ncout)
    
    print(paste0("FUSE --- ",trans," --- ", mod))
  }
}
