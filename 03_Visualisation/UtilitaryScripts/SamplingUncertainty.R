
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


#! ----------------------------- path definition

#! -------------- sources
setwd("/Users/cyrilthebault/Postdoc_Ucal/02_DATA/FUSE-Dynamic-Combination-Sim-Paper")


#! ----------------------------- package loading

library(sf)
library(dplyr)
library(dplyr)
library(tidyr)
library(abind)

source("Metrics.R")

# /!\ Adapted from gumboot package to add KGE(1/Q)
# Uses jackknife and bootstrap methods to quantify the sampling uncertainty in 
# goodness-of-fit statistics. Full details are in Clark et al. (2021), 
# "The abuse of popular performance metrics in hydrologic modeling", 
# Water Resources Research, <doi:10.1029/2020WR029001>. 
bootjack_v2 <- function (flows, GOF_stat = c("NSE", "KGE", "KGEinv"), nSample = 1000, 
                         waterYearMonth = 10, startYear = NULL, endYear = NULL, minDays = 100, 
                         minYears = 10, returnSamples = FALSE, seed = NULL, bootYearFile = NULL) 
{
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
  
  if ("KGEinv" %in% GOF_stat) {
    KGEinv_is_present <- TRUE
  } else {
    KGEinv_is_present <- FALSE
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
      
      
      if (KGEinv_is_present) {
        
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

        
        if (samplingStrategy[iStrategy] == "jack") {
          statsJack$KGEinv[iSample,] = c(meanSiminv, meanObsinv,
                                          varSiminv, varObsinv,
                                          rProdinv,
                                          kgeinv)
        }
        if (samplingStrategy[iStrategy] == "boot") {
          statsBoot$KGEinv[iSample,] = c(meanSiminv, meanObsinv,
                                          varSiminv, varObsinv,
                                          rProdinv,
                                          kgeinv)
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
  if (KGEinv_is_present) {
    kgeinv_col <- "KGEinv"
  }
  
  for (iPlot in 1:numstats) {
    if (GOF_stat[iPlot] == "NSE") {
      ixPos <- nse_col
    }
    if (GOF_stat[iPlot] == "KGE") {
      ixPos <- kge_col
    }
    if (GOF_stat[iPlot] == "KGEinv") {
      ixPos <- kgeinv_col
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


dir_FUSE = file.path('00_DATA')

##-----------------------------------------
##---------------- Obs --------------------
##-----------------------------------------

Obs = loadRData(file.path(dir_FUSE, "ObservedStreamflow.Rdata"))$Qobs_mm

##-----------------------------------------
##---------------- KGE(Q) -----------------
##-----------------------------------------

Eval_FUSE_long_HQ = loadRData(file.path(dir_FUSE, "1", "Eval_FUSE.Rdata"))

decision_medians_HQ <- Eval_FUSE_long_HQ %>%
  group_by(ModelDecisions) %>%
  summarise(median_Cal_KGE = median(`Cal : KGE :  1`, na.rm = TRUE))

best_decision_HQ <- decision_medians_HQ %>%
  filter(median_Cal_KGE == max(median_Cal_KGE)) %>%
  pull(ModelDecisions)

Sim_FUSE_HQ = loadRData(file.path(dir_FUSE,"1", "SimArray.Rdata"))

Sim_FUSE_HQ_sub = Sim_FUSE_HQ[,as.character(best_decision_HQ),]

rm(Sim_FUSE_HQ)

##-----------------------------------------
##---------------- KGE(1/Q) ---------------
##-----------------------------------------

Eval_FUSE_long_LQ = loadRData(file.path(dir_FUSE,"-1", "Eval_FUSE.Rdata"))

decision_medians_LQ <- Eval_FUSE_long_LQ %>%
  group_by(ModelDecisions) %>%
  summarise(median_Cal_KGEinv = median(`Cal : KGE : -1`, na.rm = TRUE))

best_decision_LQ <- decision_medians_LQ %>%
  filter(median_Cal_KGEinv == max(median_Cal_KGEinv)) %>%
  pull(ModelDecisions)

Sim_FUSE_LQ = loadRData(file.path(dir_FUSE,"-1", "SimArray.Rdata"))

Sim_FUSE_LQ_sub = Sim_FUSE_LQ[,as.character(best_decision_LQ),]

rm(Sim_FUSE_LQ)


##-----------------------------------------
##------------------ WA -------------------
##-----------------------------------------

Eval_FUSE_long_WA = loadRData(file.path(dir_FUSE, "HQLQ_WA_q", "Eval_FUSE.Rdata"))

best_decision_WA <- "16_17_5"

# FUSE WA
Sim_FUSE_WA = loadRData(file.path(dir_FUSE, "HQLQ_WA_q", "SimArray.Rdata"))

Sim_FUSE_WA_sub = Sim_FUSE_WA[,as.character(best_decision_WA),]

rm(Sim_FUSE_WA)

gc()

##-----------------------------------------
##----------------- SU --------------------
##-----------------------------------------
Ind_Eval = which(as.Date(rownames(Sim_FUSE_WA_sub)) >= as.Date("1999-01-01"))
catchments = colnames(Sim_FUSE_WA_sub)

# Define model variants you want to process
criteria <- c("KGE", "KGEinv")
model_names <- c("HQ", "LQ", "WA")
gumboot_columns <- c("GOF_stat", "seJack", "seBoot", "p05", "p50", "p95", "score", "biasJack", "biasBoot", "seJab")

# Build named list of simulation sub-data
sim_data_list <- list(
  HQ = Sim_FUSE_HQ_sub,
  LQ = Sim_FUSE_LQ_sub,
  WA = Sim_FUSE_WA_sub
)


# Initialize array: [catchment, model, criterion, metric]
su_array <- array(NA,
                  dim = c(length(catchments), length(model_names), length(criteria), length(gumboot_columns)),
                  dimnames = list(catchments, model_names, criteria, gumboot_columns))

# -----------------------------------------
# Main loop
# -----------------------------------------

for (catchment in catchments) {
  print(catchment)
  
  for (model in model_names) {
    sim_mat <- sim_data_list[[model]]
    
    flow_df <- data.frame(
      date = as.Date(rownames(sim_mat)[Ind_Eval]),
      obs = Obs[Ind_Eval, catchment],
      sim = sim_mat[Ind_Eval, catchment]
    )
    
    if(all(flow_df$sim == 0, na.rm = TRUE)){
      print(paste0("No data in ", catchment, " - ", model))
      next
    }
    
    result_df <- bootjack_v2(flow_df, GOF_stat = criteria, seed = 1)
    
    su_array[catchment, model, , ] <- as.matrix(result_df)
  }
}

# -----------------------------------------
# Save all result data.frames
# -----------------------------------------

save(su_array, file = file.path(dir_FUSE, "SamplingUncertainty.Rdata"))
