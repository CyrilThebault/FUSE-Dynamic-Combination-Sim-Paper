# This R script was created to reproduce the results from Knoben et al. (2025)
#
# Knoben, W. J. M., Raman, A., Gründemann, G. J., Kumar, M., Pietroniro, A., Shen, C.,
# Song, Y., Thébault, C., van Werkhoven, K., Wood, A. W., & Clark, M. P. (2024). 
# Technical note: How many models do we need to simulate hydrologic processes across 
# large geographical domains? Hydrology and Earth System Sciences Discussions, 1–21.
# https://doi.org/10.5194/hess-2024-279



#! ---------------------------------------------------------------------------------------
#!
#! Description       :
#!
#! Authors           : Cyril Thebault <cyril.thebault@inrae.fr>
#!
#! Creation date     : 2024-09-14 13:37:55
#! Modification date :
#!
#! Comments          :
#!
#! ---------------------------------------------------------------------------------------


# Load required libraries
library(tidyverse)
library(sf)
library(gumboot)
library(lubridate)
library(ncdf4)
library(lpSolve)

##########################
#   FUNCTIONS
##########################

# Helper functions
kge <- function(obs, sim) {
  
  suis_NA <- is.na(obs) | is.na(sim)
  
  qsim = sim[!suis_NA]
  qobs = obs[!suis_NA]
  
  sigma_Obs <- sd(qobs)
  sigma_Sim <- sd(qsim)
  alpha <- sigma_Sim/sigma_Obs
  ere <- cor(qobs, 
             qsim, method = "pearson")
  beta <- sum(qsim)/sum(qobs)
  return(list(KGE = 1-sqrt((ere-1)^2+(beta-1)^2+(alpha-1)^2),
              alpha = alpha,
              beta = beta,
              ere = ere))
}

extract_letter_number <- function(label) {
  parts <- strsplit(label, "_")[[1]]
  list(letter = tolower(parts[1]), number = as.integer(parts[2]))
}

generate_random_colors <- function(num_colors, seed = 1997) {
  set.seed(seed)
  colors <- sample(colors(), num_colors)
  return(colors)
}

# Function to get CAMELS observations
get_camels_obs_in_mm_d <- function(qobs_files, basin, area) {
  qobs_file <- qobs_files[grep(sprintf("/%08d_", as.integer(basin)), qobs_files)]
  
  if (length(qobs_file) != 1) {
    stop("Case not implemented")
  }
  
  tmp <- read.table(qobs_file, header = FALSE, col.names = c("gauge_id", "year", "month", "day", "qobs", "qc"))
  tmp <- tmp %>% select(-gauge_id, -qc)
  
  tmp$date <- as.Date(paste(tmp$year, tmp$month, tmp$day, sep = "-"))
  tmp <- tmp %>% select(date, qobs) %>% arrange(date)
  
  tmp$qobs[tmp$qobs < 0] <- NA
  
  m3_per_ft3 <- 0.028317
  s_per_d <- 86400
  m2_per_km2 <- 1000000
  mm_per_m <- 1000
  area_m2 <- area * m2_per_km2
  tmp$qobs_mm_d <- tmp$qobs * m3_per_ft3 * s_per_d / area_m2 * mm_per_m
  
  return(tmp)
}

# Function to get MARRMoT simulations
get_marrmot_sim_in_mm_d <- function(model_netcdf_names, basin=NULL, model) {
  
  if(! any(grepl(paste0(model,"_"),model_netcdf_names))){
    print("Model not found, return NA")
    return(NA)
  } 
  
  model_netcdf <- nc_open(model_netcdf_names[grepl(paste0(model,"_"),model_netcdf_names)])
  
  # Extract calibration and validation dates from netcdf attributes
  cal_s <- ncatt_get(model_netcdf,0, "time_calibration_start")$value
  cal_e <- ncatt_get(model_netcdf,0, "time_calibration_end")$value
  val_s <- ncatt_get(model_netcdf,0, "time_evaluation_start")$value
  val_e <- ncatt_get(model_netcdf,0, "time_evaluation_end")$value
  
  # Get simulations for all basins
  if(strsplit(model, "_")[[1]][1] == "m"){
    
    qsim <- data.frame(t(ncvar_get(model_netcdf, "Sim_q")[,,1]))
  } else {
    
    qsim <- data.frame(t(ncvar_get(model_netcdf, "Sim_q")))
    
  }
  
  colnames(qsim) = ncvar_get(model_netcdf, "Gauge_ID")
  time <- seq(from=as.Date(cal_s), to=as.Date(val_e), by="day") 
  
  nc_close(model_netcdf)
  
  if(! is.null(basin)){
    qsim = data.frame(qsim[,as.character(basin)])
    colnames(qsim) = basin
  }
  
  result <- list(
    simulation_data = cbind(date = time, qsim),
    calibration_start = cal_s,
    calibration_end = cal_e,
    validation_start = val_s,
    validation_end = val_e
  )
  
  return(result)
  
  
}


# Calculate KGE scores for each model
calculate_kge_scores <- function(model_netcdf_names, qobs_all) {
  results <- list()
  
  for (model_netcdf_name in model_netcdf_names) {
    model_name <- tools::file_path_sans_ext(basename(model_netcdf_name))
    model_id <- paste(strsplit(model_name, "_")[[1]][1:2], collapse = "_")
    
    print(paste("Running model", model_id, ":", model_name))
    
    sim_info <- get_marrmot_sim_in_mm_d(model_netcdf_name, NULL, model_id)
    
    qsim_all <- sim_info[["simulation_data"]]
    
    kges_c <- c()
    kges_v <- c()
    basins <- c()
    
    for (basin in colnames(qsim_all)[-1]) {  # Skip the date column
      qobs <- qobs_all %>% select(date, !!basin)
      qsim <- qsim_all %>% select(date, !!basin)
      
      # Subset to cal and val periods
      qobs_cal <- qobs %>% filter(date >= sim_info[["calibration_start"]], date <= sim_info[["calibration_end"]])
      qsim_cal <- qsim %>% filter(date >= sim_info[["calibration_start"]], date <= sim_info[["calibration_end"]])
      qobs_val <- qobs %>% filter(date >= sim_info[["validation_start"]], date <= sim_info[["validation_end"]])
      qsim_val <- qsim %>% filter(date >= sim_info[["validation_start"]], date <= sim_info[["validation_end"]])
      
      qobs_full <- qobs %>% filter(date >= "1987-01-01", date <= "2009-12-31")
      qsim_full <- qsim %>% filter(date >= "1987-01-01", date <= "2009-12-31")
      
      # Remove NAs
      cal_data <- na.omit(inner_join(qobs_cal, qsim_cal, by = "date"))
      val_data <- na.omit(inner_join(qobs_val, qsim_val, by = "date"))
      
      # Calculate KGE
      kge_cal <- kge(cal_data[[paste0(basin, ".x")]], cal_data[[paste0(basin, ".y")]])$KGE
      kge_val <- kge(val_data[[paste0(basin, ".x")]], val_data[[paste0(basin, ".y")]])$KGE
      
      # Calculate KGEinv
      epsilon <- mean(cal_data[[paste0(basin, ".x")]], na.rm = TRUE)/100
      kge_cal_inv <- kge(obs = (epsilon+cal_data[[paste0(basin, ".x")]])^-1, sim =(epsilon+cal_data[[paste0(basin, ".y")]])^-1)$KGE
      kge_val_inv <- kge(obs = (epsilon+val_data[[paste0(basin, ".x")]])^-1, sim =(epsilon+val_data[[paste0(basin, ".y")]])^-1)$KGE
      
      # Calculate KGEcomp
      kge_cal_comp <- (kge_cal+kge_cal_inv)/2
      kge_val_comp <- (kge_val+kge_val_inv)/2
      
      basins <- c(basins, basin)
      kges_c <- c(kges_c, kge_cal_comp)
      kges_v <- c(kges_v, kge_val_comp)
      
    }
    
    results[[model_id]] <- data.frame(
      gauge_id = basins,
      cal_kge = kges_c,
      val_kge = kges_v
    )
    
  }
  
  return(results)
}

# Modified bootjack function from gumboot package (Clark et al., 2021) to use KGEcomp score
#
# Clark, M. P., Vogel, R. M., Lamontagne, J. R., Mizukami, N., Knoben, W. J. M., Tang, G., 
# Gharari, S., Freer, J. E., Whitfield, P. H., Shook, K. R., & Papalexiou, S. M. (2021). 
# The Abuse of Popular Performance Metrics in Hydrologic Modeling. Water Resources Research, 57(9), 
# e2020WR029001. https://doi.org/10.1029/2020WR029001

bootjack_v2 <- function (flows, GOF_stat = c("NSE", "KGE", "KGEcomp"), nSample = 1000, 
                         waterYearMonth = 10, startYear = NULL, endYear = NULL, minDays = 100, 
                         minYears = 10, returnSamples = FALSE, seed = NULL, bootYearFile = NULL) {
  if ("KGE" %in% GOF_stat) {
    KGE_is_present <- TRUE
  } else {
    KGE_is_present <- FALSE
  }
  
  if ("NSE" %in% GOF_stat) {
    NSE_is_present <- TRUE
  } else {
    NSE_is_present <- FALSE
  }
  
  if ("KGEcomp" %in% GOF_stat) {
    KGEcomp_is_present <- TRUE
  } else {
    KGEcomp_is_present <- FALSE
  }
  
  options(dplyr.summarise.inform = F)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  write_bootyears <- FALSE
  read_bootyears <- FALSE
  if (!is.null(bootYearFile)) {
    if (!file.exists(bootYearFile)){ 
      write_bootyears <- TRUE
    } else {
      read_bootyears <- TRUE
    }
  }
  flows$year <- as.numeric(format(flows$date, format = "%Y"))
  flows$month <- as.numeric(format(flows$date, format = "%m"))
  flows$day <- as.numeric(format(flows$date, format = "%d"))
  iyUnique <- unique(flows$year)
  nTrials <- length(endYear)
  flows$iyWater <- ifelse(flows$month >= waterYearMonth, flows$year + 
                            1, flows$year)
  iyWater <- ifelse(flows$month >= waterYearMonth, flows$year + 
                      1, flows$year)
  nYears <- length(unique(iyWater))
  
  
  yearsJack <- array(-9999, c(nYears))
  yearsBoot <- array(-9999, c(nYears, nSample))
  zeroVal <- -1e-10
  ixValid <- which((flows$obs > zeroVal) & (flows$sim > zeroVal))
  good_flows <- flows[ixValid, ]
  valid_days <- good_flows %>% group_by(iyWater) %>% summarise(good_days = n_distinct(date))
  valid_years <- valid_days[valid_days$good_days > minDays, ]
  if (!is.null(startYear)) {
    valid_years <- valid_years[valid_years$iyWater >= startYear, ]
  }
  if (!is.null(endYear)) {
    valid_years <- valid_years[valid_years$iyWater <= endYear, ]
  }
  nyValid <- nrow(valid_years)
  if (nyValid < minYears) {
    errorStats <- data.frame(GOF_stat = "", seJack = NA_real_, 
                             seBoot = NA_real_, p05 = NA_real_, p50 = NA_real_, 
                             p95 = NA_real_, score = NA_real_, biasJack = NA_real_, 
                             biasBoot = NA_real_, seJab = NA_real_)
    return(errorStats)
  }
  izUnique <- valid_years$iyWater
  samplingStrategy <- c("jack", "boot")
  n_sampling_strategies <- length(samplingStrategy)
  mSample <- c(nyValid + 1, nSample + 1)
  
  statsJack <- lapply(1:3, function(x) data.frame(matrix(NA, nrow = mSample[1], ncol = 6)))
  names(statsJack) <- GOF_stat
  colnames_list <- c("meanSim", "meanObs", "varSim", "varObs", "rProd", "score")
  statsJack <- lapply(statsJack, function(df) { colnames(df) <- colnames_list; df })
  
  statsBoot <- lapply(1:3, function(x) data.frame(matrix(NA, nrow = mSample[2], ncol = 6)))
  names(statsBoot) <- GOF_stat
  colnames_list <- c("meanSim", "meanObs", "varSim", "varObs", "rProd", "score")
  statsBoot <- lapply(statsBoot, function(df) { colnames(df) <- colnames_list; df })
  
  if("KGEcomp" %in% GOF_stat){
    
    statsJack$KGEcomp$kge <- NA
    statsJack$KGEcomp$meanSiminv <- NA
    statsJack$KGEcomp$meanObsinv <- NA
    statsJack$KGEcomp$varSiminv <- NA
    statsJack$KGEcomp$varObsinv <- NA
    statsJack$KGEcomp$rProdinv <- NA
    statsJack$KGEcomp$kgeinv <- NA
    
    statsJack$KGEcomp <- statsJack$KGEcomp[, c(
      "meanSim", "meanObs", "varSim", "varObs", "rProd",
      "kge", "meanSiminv", "meanObsinv", "varSiminv", "varObsinv", "rProdinv", "kgeinv",
      "score"
    )]
    
    
    statsBoot$KGEcomp$kge <- NA
    statsBoot$KGEcomp$meanSiminv <- NA
    statsBoot$KGEcomp$meanObsinv <- NA
    statsBoot$KGEcomp$varSiminv <- NA
    statsBoot$KGEcomp$varObsinv <- NA
    statsBoot$KGEcomp$rProdinv <- NA
    statsBoot$KGEcomp$kgeinv <- NA
    
    statsBoot$KGEcomp <- statsBoot$KGEcomp[, c(
      "meanSim", "meanObs", "varSim", "varObs", "rProd",
      "kge", "meanSiminv", "meanObsinv", "varSiminv", "varObsinv", "rProdinv", "kgeinv",
      "score"
    )]
  }
  
  
  
  if (read_bootyears) {
    bootYears <- read.csv(bootYearFile, header = FALSE)
  }
  for (iStrategy in 1:n_sampling_strategies) {
    for (iSample in 1:mSample[iStrategy]) {
      if (iSample == 1) {
        jxValid <- ixValid
      } else {
        if (samplingStrategy[iStrategy] == "jack") {
          yearsJack[iSample - 1] = izUnique[iSample - 
                                              1]
          iyIndex <- which(iyWater[ixValid] != izUnique[iSample - 
                                                          1])
          jxValid <- ixValid[iyIndex]
        }
        else {
          if (!read_bootyears) {
            uRand <- runif(nyValid)
            ixYear <- floor(uRand * nyValid) + 1
            iyYear <- izUnique[ixYear]
          }
          else {
            iyYear <- bootYears[, iSample - 1]
          }
          yearsBoot[1:nyValid, iSample - 1] <- iyYear
          iyIndex <- which(iyWater[ixValid] == iyYear[1])
          for (iYear in 2:nyValid) iyIndex <- c(iyIndex, 
                                                which(iyWater[ixValid] == iyYear[iYear]))
          jxValid <- ixValid[iyIndex]
        }
      }
      
      
      if (NSE_is_present) {
        
        qSimValid <- flows$sim[jxValid]
        qObsValid <- flows$obs[jxValid]
        wyValid = iyWater[jxValid]
        meanSim <- mean(qSimValid, na.rm = TRUE)
        meanObs <- mean(qObsValid, na.rm = TRUE)
        varSim <- var(qSimValid, na.rm = TRUE)
        varObs <- var(qObsValid, na.rm = TRUE)
        rProd <- cor(qSimValid, qObsValid)
        
        yBeta = (meanObs - meanSim)/sqrt(varObs)
        alpha = sqrt(varSim)/sqrt(varObs)
        nse = 2 * alpha * rProd - yBeta^2 - alpha^2
        if (samplingStrategy[iStrategy] == "jack") {
          statsJack$NSE[iSample,] = c(meanSim, meanObs,
                                      varSim, varObs,
                                      rProd,
                                      nse)
        }
        if (samplingStrategy[iStrategy] == "boot") {
          statsBoot$NSE[iSample,] = c(meanSim, meanObs,
                                      varSim, varObs,
                                      rProd,
                                      nse)
        }
      }
      
      if (KGE_is_present) {
        
        qSimValid <- flows$sim[jxValid]
        qObsValid <- flows$obs[jxValid]
        wyValid = iyWater[jxValid]
        meanSim <- mean(qSimValid, na.rm = TRUE)
        meanObs <- mean(qObsValid, na.rm = TRUE)
        varSim <- var(qSimValid, na.rm = TRUE)
        varObs <- var(qObsValid, na.rm = TRUE)
        rProd <- cor(qSimValid, qObsValid)
        
        xBeta <- meanSim/meanObs
        alpha <- sqrt(varSim)/sqrt(varObs)
        kge <- 1 - sqrt((xBeta - 1)^2 + (alpha - 1)^2 + 
                          (rProd - 1)^2)
        
        if (samplingStrategy[iStrategy] == "jack") {
          statsJack$KGE[iSample,] = c(meanSim, meanObs,
                                      varSim, varObs,
                                      rProd,
                                      kge)
        }
        if (samplingStrategy[iStrategy] == "boot") {
          statsBoot$KGE[iSample,] = c(meanSim, meanObs,
                                      varSim, varObs,
                                      rProd,
                                      kge)
        }
      }
      
      if (KGEcomp_is_present) {
        
        epsilon <- mean(flows$obs, na.rm = TRUE )/100
        transfo <- -1
        qSimValidinv <- (epsilon+flows$sim[jxValid])^transfo
        qObsValidinv <- (epsilon+flows$obs[jxValid])^transfo
        wyValid = iyWater[jxValid]
        meanSiminv <- mean(qSimValidinv, na.rm = TRUE)
        meanObsinv <- mean(qObsValidinv, na.rm = TRUE)
        varSiminv <- var(qSimValidinv, na.rm = TRUE)
        varObsinv <- var(qObsValidinv, na.rm = TRUE)
        rProdinv <- cor(qSimValidinv, qObsValidinv)
        
        xBetainv <- meanSiminv/meanObsinv
        alphainv <- sqrt(varSiminv)/sqrt(varObsinv)
        kgeinv <- 1 - sqrt((xBetainv - 1)^2 + (alphainv - 1)^2 + 
                             (rProdinv - 1)^2)
        
        
        qSimValid <- flows$sim[jxValid]
        qObsValid <- flows$obs[jxValid]
        wyValid = iyWater[jxValid]
        meanSim <- mean(qSimValid, na.rm = TRUE)
        meanObs <- mean(qObsValid, na.rm = TRUE)
        varSim <- var(qSimValid, na.rm = TRUE)
        varObs <- var(qObsValid, na.rm = TRUE)
        rProd <- cor(qSimValid, qObsValid)
        
        xBeta <- meanSim/meanObs
        alpha <- sqrt(varSim)/sqrt(varObs)
        kge <- 1 - sqrt((xBeta - 1)^2 + (alpha - 1)^2 + 
                          (rProd - 1)^2)
        
        
        kgecomp = (kge+kgeinv)/2
        
        
        
        if (samplingStrategy[iStrategy] == "jack") {
          statsJack$KGEcomp[iSample,] = c(meanSim, meanObs,
                                          varSim, varObs,
                                          rProd,
                                          kge,
                                          meanSiminv, meanObsinv,
                                          varSiminv, varObsinv,
                                          rProdinv,
                                          kgeinv,
                                          kgecomp)
        }
        if (samplingStrategy[iStrategy] == "boot") {
          statsBoot$KGEcomp[iSample,] = c(meanSim, meanObs,
                                          varSim, varObs,
                                          rProd,
                                          kge,
                                          meanSiminv, meanObsinv,
                                          varSiminv, varObsinv,
                                          rProdinv,
                                          kgeinv,
                                          kgecomp)
        }
        
      }
      
      
    }
  }
  if (write_bootyears & (samplingStrategy[iStrategy] == "boot")) {
    write.table(yearsBoot, file = bootYearFile, row.names = FALSE, 
                col.names = FALSE, sep = ",")
  }
  if (returnSamples) {
    return_vals <- list(statsBoot = statsBoot, statsJack = statsJack)
    return(return_vals)
  }
  errorStats <- data.frame(GOF_stat = GOF_stat, seJack = NA_real_, 
                           seBoot = NA_real_, p05 = NA_real_, p50 = NA_real_, p95 = NA_real_, 
                           score = NA_real_, biasJack = NA_real_, biasBoot = NA_real_, 
                           seJab = NA_real_)
  # if (any(sapply(statsJack, function(df) any(df <= -9998, na.rm = TRUE)))) {
  #   return(errorStats)
  # }
  numstats <- length(GOF_stat)
  colnames <- names(statsJack)
  if (KGE_is_present) {
    kge_col <- "KGE"
  }
  if (NSE_is_present) {
    nse_col <- "NSE"
  }
  if (KGEcomp_is_present) {
    kgecomp_col <- "KGEcomp"
  }
  
  for (iPlot in 1:numstats) {
    if (GOF_stat[iPlot] == "NSE") {
      ixPos <- nse_col
    }
    if (GOF_stat[iPlot] == "KGE") {
      ixPos <- kge_col
    }
    if (GOF_stat[iPlot] == "KGEcomp") {
      ixPos <- kgecomp_col
    }
    xJack <- statsJack[[ixPos]][, "score"]
    xBoot <- statsBoot[[ixPos]][, "score"]
    iSort <- order(xJack)
    score <- xJack[1]
    zJack <- xJack[2:(nYears + 1)]
    zBoot <- xBoot[2:(nSample + 1)]
    ixJack <- which(zJack > -9998 & (!is.na(zJack)))
    nJack <- length(ixJack)
    jackMean <- mean(zJack, na.rm = TRUE)
    jackScore <- (nJack * score) - (nJack - 1) * jackMean
    sumSqErr <- (nJack - 1) * sum((jackMean - zJack[ixJack])^2)
    seJack <- sqrt(sumSqErr/nJack)
    ySample <- zBoot[order(zBoot)]
    seBoot <- sd(zBoot)
    p05 <- ySample[floor(0.05 * nSample) + 1]
    p50 <- ySample[floor(0.5 * nSample) + 1]
    p95 <- ySample[floor(0.95 * nSample) + 1]
    biasJack <- (nJack - 1) * (jackMean - score)
    biasBoot <- mean(zBoot) - score
    jabData <- vector("numeric", nYears)
    for (iYear in 2:(nYears + 1)) {
      matchYear <- vector("integer", nSample)
      for (iSample in 1:nSample) {
        ixMatch <- which(yearsBoot[, iSample] == iyUnique[iYear])
        nMatch <- length(ixMatch)
        matchYear[iSample] <- nMatch
      }
      ixMissing <- which(matchYear == 0)
      nMissing <- length(ixMissing)
      xSample <- zBoot[ixMissing]
      ySample <- xSample[order(xSample)]
      p05jack_R <- quantile(xSample, 0.05, type = 3)
      p95jack_R <- quantile(xSample, 0.95, type = 3)
      p05jack <- ySample[floor(0.05 * nMissing) + 1]
      p95jack <- ySample[floor(0.95 * nMissing) + 1]
      jabData[iYear - 1] <- p95jack - p05jack
    }
    jabMean <- mean(jabData)
    sumSqErr <- (nYears - 1) * sum((jabMean - jabData)^2)
    seJab <- sqrt(sumSqErr/nYears)
    errorStats[iPlot, ] <- c(GOF_stat[iPlot], seJack, seBoot, 
                             p05, p50, p95, score, biasJack, biasBoot, seJab)
  }
  errorStats$seJack <- as.numeric(errorStats$seJack)
  errorStats$seBoot <- as.numeric(errorStats$seBoot)
  errorStats$p05 <- as.numeric(errorStats$p05)
  errorStats$p50 <- as.numeric(errorStats$p50)
  errorStats$p95 <- as.numeric(errorStats$p95)
  errorStats$score <- as.numeric(errorStats$score)
  errorStats$biasJack <- as.numeric(errorStats$biasJack)
  errorStats$biasBoot <- as.numeric(errorStats$biasBoot)
  errorStats$seJab <- as.numeric(errorStats$seJab)
  return(errorStats)
}


min_set_cover <- function(constraint_matrix) {
  
  obj <- rep(1, ncol(constraint_matrix))
  
  constraint_matrix[is.na(constraint_matrix)] <- FALSE
  
  constraint_matrix <- as.data.frame(lapply(constraint_matrix, as.numeric))
  
  # Check if any elements are not covered by any subset
  uncovered_elements <- as.character(df_eval_gumboot$gauge_id[rowSums(constraint_matrix) == 0])
  
  # Remove rows corresponding to uncovered elements
  if (length(uncovered_elements) > 0) {
    constraint_matrix <- constraint_matrix[rowSums(constraint_matrix) > 0, , drop = FALSE]
  }
  
  # Set up the linear programming problem
  direction <- rep(">=", nrow(constraint_matrix))
  rhs <- rep(1, nrow(constraint_matrix))
  
  # Solve the problem
  result <- lp("min", obj, constraint_matrix, direction, rhs, all.bin = TRUE)
  
  # Check if a solution was found
  if (result$status != 0) {
    stop("No feasible solution found")
  }
  
  covered_elements = as.character(df_eval_gumboot$gauge_id[!as.character(df_eval_gumboot$gauge_id) %in% uncovered_elements])
  # Return the results
  list(
    num_selected = sum(result$solution),
    selected_indices = which(result$solution == 1),
    covered_elements = covered_elements,
    uncovered_elements = uncovered_elements,
    coverage_percentage = round(length(covered_elements) / nrow(df_eval_gumboot) * 100, 2)
  )
}


sort_selected_subsets_by_coverage <- function(all_basins, selected_subsets, selected_models, match_output_to_existing = FALSE) {
  
  # Prep
  remaining_basins <- all_basins # Start with a full list of basins
  remaining_subsets <- selected_subsets # Start with all selected subsets
  remaining_models <- selected_models # Start with all selected models
  sorted_models <- c()
  basins_filled <- c()
  
  # Sorted by coverage
  while (length(remaining_basins) > 0) {
    
    # Find how many of the remaining basins each remaining subset covers
    basins_covered <- sapply(remaining_subsets, function(subset) {
      length(intersect(remaining_basins, subset))
    })
    
    # Find the model that covers the most of these basins
    num_basins <- max(basins_covered)
    idx_to_use <- which.max(basins_covered)
    model <- remaining_models[idx_to_use]
    
    # Update the "remaining_" variables
    remaining_basins <- setdiff(remaining_basins, remaining_subsets[[idx_to_use]])
    remaining_models <- remaining_models[-idx_to_use]
    remaining_subsets <- remaining_subsets[-idx_to_use]
    
    # Update the outputs
    sorted_models <- c(sorted_models, model)
    basins_filled <- c(basins_filled, num_basins)
  }
  
  # Fix up output if needed
  if (match_output_to_existing) {
    sorted_models <- paste0(sorted_models, "_above_p5")
    basins_filled <- as.integer(basins_filled)
  }
  
  return(list(sorted_models = sorted_models, basins_filled = basins_filled))
}

##########################
#   MAIN
##########################

# File paths
camels_topo <- "/Users/cyrilthebault/Postdoc_Ucal/02_DATA/CAMELS/camels_topo.txt"
camels_flow <- "/Users/cyrilthebault/Postdoc_Ucal/02_DATA/CAMELS/basin_timeseries_v1p2_metForcing_obsFlow/basin_dataset_public_v1p2/usgs_streamflow/"
brief_folder <- "/Users/cyrilthebault/Postdoc_Ucal/02_DATA/FUSE/01_Paper/SamplingUncertainty"


# Load gauge metadata
df <- read.csv(camels_topo, sep = ";") %>%
  select(gauge_id, gauge_lat, gauge_lon, area_gages2, area_geospa_fabric)



# Find streamflow files
print("---- Get observed streamflow time series")
qobs_files <- list.files(path = camels_flow, pattern = "*.txt", recursive = TRUE, full.names = TRUE)

# Initialize an empty dataframe for merging
qobs_all <- NULL

# Track last progress update
last_reported_progress <- 0

# Loop through each row of the dataframe (assuming 'df' is your dataframe)
total_iterations = nrow(df)
df$area_contrib = NA
for (i in 1:total_iterations) {
  basin <- df[i, 'gauge_id']
  
  IDnumber = ifelse(nchar(basin)==7,paste0("0",basin), as.character(basin))
  BM_forcing_path = file.path(camels_flow,"..","basin_mean_forcing", "daymet")
  BM_forcing_file = NULL
  for(folder in list.files(BM_forcing_path)){
    BM_forcing_file = list.files(file.path(BM_forcing_path, folder), pattern = IDnumber)
    if(length(BM_forcing_file) != 0){break}
  }
  
  area = as.numeric(readLines(paste0(BM_forcing_path,"/",folder,"/", BM_forcing_file))[3])*10^-6
  
  df[i,'area_contrib'] = area
  
  # Get the observed data for the basin
  qobs <- get_camels_obs_in_mm_d(qobs_files, basin, area)
  
  # Drop the 'qobs' column and rename 'qobs_mm_d' to the basin name
  qobs <- qobs[, !colnames(qobs) %in% 'qobs']
  colnames(qobs)[colnames(qobs) == 'qobs_mm_d'] <- basin
  
  # Merge dataframes by a common column (e.g., time), use 'all = TRUE' to handle differing row numbers
  if (is.null(qobs_all)) {
    qobs_all <- qobs
  } else {
    qobs_all <- merge(qobs_all, qobs, by = "date", all = TRUE)
  }
  
  # Print progress in increments of 10%
  progress <- floor((i / total_iterations) * 100)
  if (progress >= last_reported_progress + 10) {
    print(paste("Progress:", progress, "%"))
    last_reported_progress <- progress
  }
  
}

print("---- Get simulated streamflow time series and calculate KGE")
# Find the model result files
model_netcdf_names <- list.files(path = brief_folder, pattern = "*.nc", full.names = TRUE)

kge_results <- calculate_kge_scores(model_netcdf_names, qobs_all)

# Get the names of the elements in the list
list_names <- names(kge_results)

# Iterate over the rest of the list
for (i in 1:length(kge_results)) {
  
  df_name <- list_names[i]
  
  if(i ==1){
    df_tmp <- kge_results[[i]]
    colnames(df_tmp)[2:3] = paste0(df_name,"_",colnames(df_tmp)[2:3])
    df_results <- df_tmp
  } else {
    df_tmp <- kge_results[[i]]
    colnames(df_tmp)[2:3] = paste0(df_name,"_",colnames(df_tmp)[2:3])
    df_results <- cbind(df_results, df_tmp[,2:3])
  }
  
}

# Merge with gauge metadata
df_ini <- merge(df, df_results, by = "gauge_id")

# Save results
write.csv(df_ini, file.path(brief_folder, "img/brief_analysis_updated_modeling_results.csv"), row.names = FALSE)

print("---- Find the best model")
# Create df_eval
df_eval <- df_ini

# Find best model for each basin
df_eval$best_cal_column <- apply(df_eval %>% select(ends_with("cal_kge")), 1, function(x) {
  # Check if x is not empty and contains valid values
  if (length(x) == 0 || all(is.na(x))) {
    return(NA)
  }
  
  # Find the name of the maximum value
  max_index <- which.max(x)
  if (length(max_index) == 0 || is.na(max_index)) {
    return(NA)
  }
  
  max_name <- names(x)[max_index]
  if (is.na(max_name)) {
    return(NA)
  } else {
    return(max_name)
  }
})

df_eval$best_cal_model <- sapply(strsplit(df_eval$best_cal_column, "_"), function(x) {
  if(NA %in% x){
    return(NA)
  } else{
    return(paste(x[1:2], collapse = "_"))
  }
  
})

df_eval$best_cal_score <- apply(df_eval %>% select(ends_with("cal_kge")), 1, function(x){
  if(all(is.na(x))){
    return(NA)
  } else{
    return(max(x, na.rm = TRUE))
  }
})

# Save results
write.csv(df_eval, file.path(brief_folder, "img/brief_analysis_updated_modeling_results_eval.csv"), row.names = FALSE)


print("---- Sampling uncertainty (gumboot)")
# Create df_gumboot
df_gumboot = df_eval

# Prepare for gumboot analysis
gumboot_cols <- c("GOF_stat", "seJack", "seBoot", "p05", "p50", "p95", "score", "biasJack", "biasBoot", "seJab")
for (col in gumboot_cols) {
  df_gumboot[,col] <- NA
}

# Gumboot analysis
for (i in 1:nrow(df_gumboot)) {
  basin <- df_gumboot[i,"gauge_id"]
  print(paste("Running", basin))
  
  qobs <- get_camels_obs_in_mm_d(qobs_files, basin, df_gumboot$area_contrib[i])
  
  sim_info <- get_marrmot_sim_in_mm_d(model_netcdf_names, basin, df_gumboot$best_cal_model[i])
  if(length(sim_info) != 5){next}
  qsim <- sim_info[["simulation_data"]]
  
  qobs_sub <- qobs %>% filter(date >= sim_info[["calibration_start"]], date <= sim_info[["calibration_end"]])
  qsim_sub <- qsim %>% filter(date >= sim_info[["calibration_start"]], date <= sim_info[["calibration_end"]])
  
  foot <- inner_join(qobs_sub, qsim_sub, by = "date") %>%
    select(date, obs = qobs_mm_d, sim = as.character(basin)) %>%
    na.omit()
  
  clean_foot <- na.omit(foot)
  
  # Run gumboot, catching errors if needed
  result <- try({
    result <- bootjack_v2(clean_foot, GOF_stat = "KGEcomp", seed = 1)
  })
  
  # Check if we encountered an error, and return all NaNs if so
  if (inherits(result, "try-error")) {
    print("An error occurred, returning all NA")
    next
  } else if (all(is.na(result))) {
    print("Not enough data, minimum 10 years required")
    next
  }
  
  df_gumboot[i,colnames(result)] = result
}

# Save results with gumboot analysis
write.csv(df_gumboot, file.path(brief_folder, "img/brief_analysis_modeling_results_with_gumboot.csv"), row.names = FALSE)

print("---- Analysis of models within uncertainty bounds")

df_eval_gumboot = df_gumboot

# Calculate range
df_eval_gumboot$range_5_95 <- df_eval_gumboot$p95 - df_eval_gumboot$p05

# Analysis of models within uncertainty bounds
columns_cal <- grep("_cal_kge$", names(df_eval_gumboot), value = TRUE)

for (column_cal in columns_cal) {
  model <- paste(strsplit(column_cal, "_")[[1]][1:2], collapse = "_")
  flag_column <- paste0(model, "_above_p5")
  df_eval_gumboot[[flag_column]] <- df_eval_gumboot[[column_cal]] >= df_eval_gumboot$p05
}

columns_flags = grep("_above_p5$", names(df_eval_gumboot), value = TRUE)
df_eval_gumboot$similar_model_count =  rowSums(df_eval_gumboot[,columns_flags], na.rm = TRUE)

constraint_matrix <- df_eval_gumboot[,columns_flags]

# Function to solve the set cover problem
result <- min_set_cover(constraint_matrix)
num_selected <- result$num_selected
selected_indices <- result$selected_indices
selected_models <- list_names[selected_indices]
covered_basins <- colSums(df_eval_gumboot[,columns_flags[selected_indices]], na.rm = TRUE)

# Print results
cat("Number of subsets selected:", num_selected, "\n")
cat("Indices of selected subsets:", selected_indices, "\n")
cat("Selected models:", selected_models, "\n")
cat("Basins covered:", covered_basins, "\n")


df_eval_gumboot$mod_need = NA

tmp_df <- df_eval_gumboot[! df_eval_gumboot$similar_model_count == 0, paste0(selected_models, "_above_p5")]
basins_left <- nrow(df)

while (basins_left > 0) {
  # Perform a greedy selection for the basins we have left
  model_needed <- names(which.max(colSums(tmp_df, na.rm = TRUE)))
  
  # Subset the dataframe to keep only those basins for which we have no model yet
  
  ind_tmp <- tmp_df[[model_needed]]
  ind_tmp[is.na(ind_tmp)] <- FALSE
  tmp_df <- tmp_df[!ind_tmp, ]
  
  ind_gumboot = is.na(df_eval_gumboot$mod_need) & df_eval_gumboot[[model_needed]]
  ind_gumboot[is.na(ind_gumboot)] <- FALSE
  df_eval_gumboot$mod_need[ind_gumboot] <- sapply(strsplit(model_needed, "_"), function(x) paste(x[1:2], collapse = "_"))
  
  # Update the count of basins we still need - this will terminate the loop
  basins_left <- nrow(tmp_df)
}

print(sort(table(df_eval_gumboot$mod_need), decreasing = TRUE))

# Save results with gumboot analysis
write.csv(df_eval_gumboot, file.path(brief_folder, "img/brief_analysis_modeling_results_with_gumboot_eval.csv"), row.names = FALSE)

