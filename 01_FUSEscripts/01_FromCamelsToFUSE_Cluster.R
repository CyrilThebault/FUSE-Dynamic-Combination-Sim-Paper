#! ---------------------------------------------------------------------------------------
#!
#! Description       :
#!
#! Authors           : Cyril Thebault <cyril.thebault@ucalgary.ca>
#!
#! Creation date     : 2024-10-05 17:24:18
#! Modification date :
#!
#! Comments          :
#!
#! ---------------------------------------------------------------------------------------

#! ----------------------------- package loading

library(dplyr)
library(tidyr)
library(ncdf4)
library(sf)
library(airGR)


#! ----------------------------  File manager

FUSE_path = "/work/comphyd_lab/users/cyril.thebault/Postdoc_Ucal/02_DATA/FUSE"   # Path of the workflow
CAMELS_path = "/home/cyril.thebault/Postdoc_Ucal/02_DATA/CAMELS/basin_timeseries_v1p2_metForcing_obsFlow/basin_dataset_public_v1p2" # Path of the dataset       # (02) Path of the CAMELS-spat database\
dir_BV = "/home/cyril.thebault/Postdoc_Ucal/02_DATA/BDD" # Path of the list of catchments

models = unlist(read.table(file.path(FUSE_path, "list_decision_78.txt"), sep = ";", header = TRUE)$ID)

catchments = unlist(read.table(file.path(dir_BV, "liste_BV_CAMELS_559.txt")))

SimStart = "1989-01-01"
SimEnd = "2009-12-31"

CalStart = "1989-01-01"
CalEnd = "1998-12-31"

WU = 2

CritCal = "KGE"
TransfoCal = "1"
dataset = "CAMELS"
spatialisation = "Lumped"
inputdata = "daymet"

tsCal = 'daily'

#! ----------------------------  catchments

args <- commandArgs(trailingOnly = TRUE)
id <- as.numeric(args[1])

catchment = catchments[id]

#! ----------------------------  create folder
mydemo = "fuse_template"
DemoFolder = file.path(FUSE_path, mydemo)
CatchmentFolder = file.path(FUSE_path, dataset, spatialisation, inputdata, CritCal, TransfoCal, catchment)

#create folder tree
if(!dir.exists(CatchmentFolder)){
  dir.create(path = file.path(CatchmentFolder, "input"), recursive = TRUE)
  dir.create(path = file.path(CatchmentFolder, "output"), recursive = TRUE)
  dir.create(path = file.path(CatchmentFolder, "settings/fuse_zDecisions"), recursive = TRUE)
}

if(length(list.files(file.path(CatchmentFolder, "settings")))==1){
  #copy settings files
  settings_files = list.files(file.path(DemoFolder,"settings"))
  
  invisible(file.copy(from = paste0(DemoFolder,"/settings/", settings_files), 
                      to = file.path(CatchmentFolder,"settings"),
                      recursive = FALSE)
  )
  
}

DecisionFiles = paste0("fuse_zDecisions_",models,".txt")
if(!any(DecisionFiles %in% list.files(file.path(CatchmentFolder, "settings/fuse_zDecisions")))){
  
  invisible(file.copy(from = paste0(DemoFolder,"/settings/fuse_zDecisions/", DecisionFiles), 
                      to = paste0(CatchmentFolder,"/settings/fuse_zDecisions"), recursive = TRUE)
  )
  
}

#! ----------------------------  get input data

if((!paste0(catchment,"_elev_bands.nc") %in% list.files(file.path(CatchmentFolder, "input")))|
   (!paste0(catchment,"_input.nc") %in% list.files(file.path(CatchmentFolder, "input")))){
  
  ######## Metadata
  df_metadata <- read.csv2(file.path(CAMELS_path, "camels-spat-metadata.csv"), sep = ",")
  
  Country = strsplit(catchment, "_")[[1]][1]
  IDnumber = strsplit(catchment, "_")[[1]][2]
  
  category_value <- df_metadata$subset_category[df_metadata$Country == Country & df_metadata$Station_id == IDnumber]
  catchment_area <- as.numeric(df_metadata$Basin_area_km2[df_metadata$Country == Country & df_metadata$Station_id == IDnumber])
  
  ######## Shapefiles
  shp_catchment = read_sf(paste0(CAMELS_path,"/shapefiles/",category_value,"/shapes-", tolower(spatialisation),"/", catchment, "/",catchment,"_",tolower(spatialisation),".shp"))
  
  if(!paste0(catchment,"_input.nc") %in% list.files(file.path(CatchmentFolder, "input"))){
    
    ####### Dates
    DatesR <- seq(from = as.Date(gsub(format(as.Date(min(SimStart, CalStart)), "%Y"),
                                      as.numeric(format(as.Date(min(SimStart, CalStart)), "%Y"))-WU,
                                      min(SimStart, CalStart))),
                  to = as.Date(max(SimEnd, CalEnd)), 
                  by = 'day')
    
    ###### Precipitation, temperature & potential evapotranspiration
    
    # Basin mean
    
    forcing_path = paste0(CAMELS_path,"/forcing/",category_value,"/",inputdata,"/",inputdata,"-", tolower(spatialisation))
    forcing_file = paste0(catchment,"_",inputdata,"_",tolower(spatialisation),".nc")
    if(inputdata == "em-earth"){
      forcing_file = gsub("em-earth", "em_earth", forcing_file)
    }
    
    nc_forcing = nc_open(file.path(forcing_path,forcing_file))
    
    cf <- CFtime(ncatt_get(nc_forcing,"time","units")$value, calendar = ncatt_get(nc_forcing,"time","calendar")$value, ncvar_get(nc_forcing,"time")) # convert time to CFtime class
    timestamps <- as_timestamp(cf) # get character-string times
    
    
    if(inputdata == "daymet"){
      
      time_LT <- as.POSIXct(timestamps, format="%Y-%m-%d %H:%M", tz="UTC")
      
      forcing = data.frame(time_LT = time_LT) 
      forcing$Date = as.Date(time_LT)
      
      lat = as.numeric(unique(ncvar_get(nc_forcing, "latitude")))
      lon = as.numeric(unique(ncvar_get(nc_forcing, "longitude")))
      
      forcing[,'Tmean[C]'] = rowMeans(data.frame(Tmin = ncvar_get(nc_forcing, "tmin"), Tmax = ncvar_get(nc_forcing, "tmax")))
      forcing[,'P[mm/d]'] = ncvar_get(nc_forcing, "prcp")
      forcing[,'P[mm/d]'][forcing[,'P[mm/d]'] < 0] = 0
      forcing[,'PET_Oudin[mm/d]'] <- PE_Oudin(JD = as.POSIXlt(forcing[,"Date"])$yday + 1,
                                              Temp = forcing[,'Tmean[C]'],
                                              Lat = lat, LatUnit = "deg")
      
      subforcing = forcing[forcing$Date %in% DatesR,]
    } else if(inputdata == "em-earth"){
      
      time_LT <- as.POSIXct(timestamps, format="%Y-%m-%d %H:%M:%S", tz="UTC")
      
      forcing_h = data.frame(time_LT = time_LT) 
      forcing_h$Date = as.Date(time_LT)
      
      lat = as.numeric(unique(ncvar_get(nc_forcing, "latitude")))
      lon = as.numeric(unique(ncvar_get(nc_forcing, "longitude")))
      
      forcing_h[,'Tmean[C]'] = ncvar_get(nc_forcing, "tmean") - 273.15
      forcing_h[,'P[mm/h]'] = ncvar_get(nc_forcing, "prcp")*3600
      forcing_h[,'P[mm/h]'][forcing_h[,'P[mm/h]'] < 0] = 0
      
      forcing <- forcing_h %>%
        group_by(Date) %>%
        summarise(
          `Tmean[C]` = mean(`Tmean[C]`, na.rm = TRUE),
          `P[mm/d]` = sum(`P[mm/h]`, na.rm = TRUE)
        )
      
      
      forcing[,'PET_Oudin[mm/d]'] <- PE_Oudin(JD = as.POSIXlt(forcing$Date)$yday + 1,
                                              Temp = as.numeric(forcing[['Tmean[C]']]),
                                              Lat = lat, LatUnit = "deg")
      
      subforcing = forcing[match(DatesR, forcing$Date),]
    }
    
    
    nc_close(nc_forcing)
    
    ###### Observed streamflow
    
    observed_path = paste0(CAMELS_path, "/observations/",category_value,"/obs-",tsCal)
    observed_file = paste0(catchment, "_", tsCal, "_flow_observations.nc")
    nc_observed = nc_open(file.path(observed_path, observed_file))
    
    cf <- CFtime(ncatt_get(nc_observed,"time","units")$value, calendar = ncatt_get(nc_observed,"time","calendar")$value, ncvar_get(nc_observed,"time")) # convert time to CFtime class
    timestamps <- as_timestamp(cf) # get character-string times
    time_LT <- as.POSIXct(timestamps, format="%Y-%m-%d", tz="UTC")
    
    observed = data.frame(time_LT = time_LT) 
    observed$Date = as.Date(time_LT)
    
    observed$qm3s =  ncvar_get(nc_observed, "q_obs")
    
    observed$qmm = observed$qm3s*86400/(catchment_area*10^3) 
    
    qmm = observed$qmm[match(DatesR, observed$Date)]
    
    qmm[is.nan(qmm)] = NA
    qmm[qmm< 0] = NA
    
    subforcing[,'Qobs[mm/d]'] = qmm
    
    nc_close(nc_observed)
    #! ----------------------------  create input file
    
    inputname <- paste0(CatchmentFolder, "/input/",catchment,"_input.nc")
    
    
    # define dimensions
    latdim <- ncdim_def(name = "latitude", units = "degreesN", vals = lat, longname = "latitude", create_dimvar=TRUE)
    londim <- ncdim_def(name = "longitude", units = "degreesE", vals = lon, longname = "longitude", create_dimvar=TRUE)
    timedim <-  ncdim_def(name = "time", units = paste0("days since ", min(DatesR)), vals = (1:nrow(subforcing))-1, longname = "time",create_dimvar=TRUE, unlim = TRUE)
    
    # define variables
    PET_def <- ncvar_def(name = "pet", 
                         units = "mm/day",
                         dim = list(londim, latdim, timedim),
                         longname = "Potential evaportanspiration estimated using Oudin et al., 2005, JoH",
                         prec = "double"
    )
    
    P_def <- ncvar_def(name = "pr", 
                       units = "mm/day",
                       dim = list(londim, latdim, timedim),
                       longname = "Mean daily precipitation",
                       prec = "double"
    )
    
    Qobs_def <- ncvar_def(name = "q_obs", 
                          units = "mm/day",
                          dim = list(londim, latdim, timedim),
                          longname = "Mean observed daily discharge",
                          prec = "double"
    )
    
    T_def <- ncvar_def(name = "temp", 
                       units = "degC",
                       dim = list(londim, latdim, timedim),
                       longname = "Mean daily temperature",
                       prec = "double"
    )
    
    
    ncinput <- nc_create(inputname,list(PET_def, P_def, Qobs_def, T_def), force_v4 = TRUE)
    
    # put variables
    
    pet_array <- array(subforcing[["PET_Oudin[mm/d]"]], dim = c(1, 1, nrow(subforcing)))
    p_array    <- array(subforcing[["P[mm/d]"]], dim = c(1, 1, nrow(subforcing)))
    qobs_array <- array(subforcing[["Qobs[mm/d]"]], dim = c(1, 1, nrow(subforcing)))
    t_array    <- array(subforcing[["Tmean[C]"]], dim = c(1, 1, nrow(subforcing)))
    
    ncvar_put(ncinput,PET_def, pet_array)
    ncvar_put(ncinput,P_def, p_array)
    ncvar_put(ncinput,Qobs_def,qobs_array)
    ncvar_put(ncinput,T_def,t_array)
    
    
    # add global attributes
    ncatt_put(ncinput,0,"author", "Cyril Thébault")
    ncatt_put(ncinput,0,"date",paste0("Created ", format(Sys.time(), format = "%Y/%m/%d %H:%M:%S")))
    ncatt_put(ncinput,0,"institution",'University of Calgary')
    
    nc_close(ncinput)
  }
  
  #! ----------------------------  get elevation band data
  
  if(!paste0(catchment,"_elev_bands.nc") %in% list.files(file.path(CatchmentFolder, "input"))){
    
    # Elevation 
    elev_path = paste0(CAMELS_path, "/geospatial/",category_value,"/merit/", catchment)
    elev_file = paste0(catchment, "_merit_hydro_elv.tif")
    elev = raster(file.path(elev_path, elev_file))
    
    elev_vals <- getValues(elev)
    elev_vals <- elev_vals[!is.na(elev_vals) & is.finite(elev_vals)]
    
    band_breaks <- seq(floor(min(elev_vals, na.rm = TRUE) / 100) * 100, ceiling(max(elev_vals, na.rm = TRUE) / 100) * 100, by = 100)
    elev_band <- cut(elev_vals, breaks = band_breaks, right = FALSE)
    
    band_fraction <- as.numeric(prop.table(table(elev_band)))
    band_midpoints <- (band_breaks[-1] + band_breaks[-length(band_breaks)]) / 2
    
    #! ----------------------------  create elevation band file
    
    elevname <- paste0(CatchmentFolder, "/input/",catchment,"_elev_bands.nc")
    
    # define dimensions
    latdim <- ncdim_def(name = "latitude", units = "degreesN", vals = lat, longname = "latitude", create_dimvar=TRUE)
    londim <- ncdim_def(name = "longitude", units = "degreesE", vals = lon, longname = "longitude", create_dimvar=TRUE)
    elevdim <-  ncdim_def(name = "elevation_band", units = "-", vals = 1:length(band_fraction), longname = "elevation_band",create_dimvar=TRUE)
    
    
    # define variables
    AreaFrac_def <- ncvar_def(name = "area_frac", 
                              units = "-",
                              dim = list(londim, latdim, elevdim),
                              longname = "Fraction of the catchment covered by each elevation band",
                              prec = "double"
    )
    
    MeanElev_def <- ncvar_def(name = "mean_elev", 
                              units = "m asl",
                              dim = list(londim, latdim, elevdim),
                              longname = "Mid-point elevation of each elevation band",
                              prec = "double"
    )
    
    PrecFrac_def <- ncvar_def(name = "prec_frac", 
                              units = "-",
                              dim = list(londim, latdim, elevdim),
                              longname = "Fraction of catchment precipitation that falls on each elevation band - same as area_frac",
                              prec = "double"
    )
    
    
    ncelev <- nc_create(elevname,list(AreaFrac_def, MeanElev_def, PrecFrac_def), force_v4 = TRUE)
    
    # put variables
    areafrac_array    <- array(as.numeric(band_fraction), dim = c(1, 1, length(band_fraction)))
    meanelev_array <- array(as.numeric(band_midpoints), dim = c(1, 1, length(band_fraction)))
    precfrac_array    <- array(as.numeric(band_fraction), dim = c(1, 1, length(band_fraction)))
    
    ncvar_put(ncelev,AreaFrac_def, areafrac_array)
    ncvar_put(ncelev,MeanElev_def, meanelev_array)
    ncvar_put(ncelev,PrecFrac_def, precfrac_array)
    
    # add global attributes
    ncatt_put(ncelev,0,"author", "Cyril Thébault")
    ncatt_put(ncelev,0,"date",paste0("Created ", format(Sys.time(), format = "%Y/%m/%d %H:%M:%S")))
    ncatt_put(ncelev,0,"institution",'University of Calgary')
    
    nc_close(ncelev)
  }
}

#! ----------------------------  fm_catch file

for(mod in models){
  
  FM = readLines( file.path(DemoFolder, "fm_catch.txt"))
  
  FM[which(grepl("! SETNGS_PATH", FM))] = gsub( "/my/path/to/fuse/catchment/settings/" , file.path("/dev/shm",catchment,"settings/") , FM[which(grepl("! SETNGS_PATH", FM))] )
  FM[which(grepl("! INPUT_PATH", FM))] = gsub( "/my/path/to/fuse/catchment/input/" , file.path("/dev/shm",catchment,"input/") , FM[which(grepl("! INPUT_PATH", FM))] )
  FM[which(grepl("! OUTPUT_PATH", FM))] = gsub( "/my/path/to/fuse/catchment/output/" , file.path("/dev/shm",catchment,"output/") , FM[which(grepl("! OUTPUT_PATH", FM))] )
  
  FM[which(grepl("! Q_ONLY", FM))] = gsub( "FALSE" , "TRUE" , FM[which(grepl("! Q_ONLY", FM))] )
  
  FM[which(grepl("! M_DECISIONS", FM))] = gsub( "902" , mod , FM[which(grepl("! M_DECISIONS", FM))] )
  FM[which(grepl("! FMODEL_ID", FM))] = gsub( "902" , mod , FM[which(grepl("! FMODEL_ID", FM))] )
  
  FM[which(grepl("! date_start_sim", FM))] = gsub( "2000-10-01" , 
                                                   as.Date(gsub(format(as.Date(min(SimStart, CalStart)), "%Y"),
                                                                as.numeric(format(as.Date(min(SimStart, CalStart)), "%Y"))-WU,
                                                                min(SimStart, CalStart))), 
                                                   FM[which(grepl("! date_start_sim", FM))] )
  FM[which(grepl("! date_end_sim", FM))] = gsub( "2005-09-30" , SimEnd , FM[which(grepl("! date_end_sim", FM))] )
  FM[which(grepl("! date_start_eval", FM))] = gsub( "2001-10-01" , CalStart , FM[which(grepl("! date_start_eval", FM))] )
  FM[which(grepl("! date_end_eval", FM))] = gsub( "2005-09-30" , CalEnd , FM[which(grepl("! date_end_eval", FM))] )
  
  FM[which(grepl("! METRIC", FM))] = gsub( "RMSE" , CritCal , FM[which(grepl("! METRIC", FM))] )
  FM[which(grepl("! TRANSFO", FM))] = gsub( "1" , TransfoCal , FM[which(grepl("! TRANSFO", FM))] )
  
  FM[which(grepl("! MAXN", FM))] = gsub( "20" , 10000 , FM[which(grepl("! MAXN", FM))] )
  
  writeLines( FM , paste0(CatchmentFolder,"/",catchment,"_",mod,".txt") )
  
}

print(paste0(catchment, " DONE"))
