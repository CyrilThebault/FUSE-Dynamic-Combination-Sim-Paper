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

path_res = "../02_DATA"

#! ----------------------------- package loading

library(dplyr)
library(tidyr)
library(abind)

source(file = "Metrics.R")

bootjack_v2 <- function (flows, GOF_stat = c("NSE", "KGE", "KGEcomp"), nSample = 1000, 
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



dir_FUSE = file.path(path_res, 'FUSE')

##-----------------------------------------
##---------------- Obs --------------------
##-----------------------------------------

Obs = loadRData(file.path(dir_FUSE,"01_Paper", "ObservedStreamflow.Rdata"))$Qobs_mm

##-----------------------------------------
##---------------- KGE(Q) -----------------
##-----------------------------------------

Eval_FUSE_long_HQ = loadRData(file.path(dir_FUSE,"01_Paper", "1", "Eval_FUSE.Rdata")) %>%
  mutate(`Cal : KGEcomp` = (`Cal : KGE :  1` + `Cal : KGE : -1`) / 2,
         `Eval : KGEcomp` = (`Eval : KGE :  1` + `Eval : KGE : -1`) / 2)

decision_medians_HQ <- Eval_FUSE_long_HQ %>%
  group_by(ModelDecisions) %>%
  summarise(median_Cal_KGEcomp = median(`Cal : KGEcomp`, na.rm = TRUE))

best_decision_HQ <- decision_medians_HQ %>%
  filter(median_Cal_KGEcomp == max(median_Cal_KGEcomp)) %>%
  pull(ModelDecisions)

Sim_FUSE_HQ = loadRData(file.path(dir_FUSE,"01_Paper","1", "SimArray.Rdata"))

Sim_FUSE_HQ_sub = Sim_FUSE_HQ[,as.character(best_decision_HQ),]

rm(Sim_FUSE_HQ)

##-----------------------------------------
##---------------- KGE(1/Q) ---------------
##-----------------------------------------

Eval_FUSE_long_LQ = loadRData(file.path(dir_FUSE,"01_Paper", "1", "Eval_FUSE.Rdata")) %>%
  mutate(`Cal : KGEcomp` = (`Cal : KGE :  1` + `Cal : KGE : -1`) / 2,
         `Eval : KGEcomp` = (`Eval : KGE :  1` + `Eval : KGE : -1`) / 2)

decision_medians_LQ <- Eval_FUSE_long_LQ %>%
  group_by(ModelDecisions) %>%
  summarise(median_Cal_KGEcomp = median(`Cal : KGEcomp`, na.rm = TRUE))

best_decision_LQ <- decision_medians_LQ %>%
  filter(median_Cal_KGEcomp == max(median_Cal_KGEcomp)) %>%
  pull(ModelDecisions)

Sim_FUSE_LQ = loadRData(file.path(dir_FUSE,"01_Paper","-1", "SimArray.Rdata"))

Sim_FUSE_LQ_sub = Sim_FUSE_LQ[,as.character(best_decision_LQ),]

rm(Sim_FUSE_LQ)

##-----------------------------------------
##---------------- Comp -------------------
##-----------------------------------------

Eval_FUSE_long_Comp = loadRData(file.path(dir_FUSE,"01_Paper", "Comp", "Eval_FUSE.Rdata"))%>%
  mutate(`Cal : KGEcomp` = (`Cal : KGE :  1` + `Cal : KGE : -1`) / 2,
         `Eval : KGEcomp` = (`Eval : KGE :  1` + `Eval : KGE : -1`) / 2)

decision_medians_Comp <- Eval_FUSE_long_Comp %>%
  group_by(ModelDecisions) %>%
  summarise(median_Cal_KGEcomp = median(`Cal : KGEcomp`, na.rm = TRUE))

best_decision_Comp <- decision_medians_Comp %>%
  filter(median_Cal_KGEcomp == max(median_Cal_KGEcomp)) %>%
  pull(ModelDecisions)

Sim_FUSE_Comp = loadRData(file.path(dir_FUSE,"01_Paper","Comp", "SimArray.Rdata"))

Sim_FUSE_Comp_sub = Sim_FUSE_Comp[,as.character(best_decision_Comp),]

rm(Sim_FUSE_Comp)

##-----------------------------------------
##------------------ WA -------------------
##-----------------------------------------

Eval_FUSE_long_WA = loadRData(file.path(dir_FUSE,"01_Paper", "WA", "Eval_FUSE.Rdata"))%>%
  mutate(`Cal : KGEcomp` = (`Cal : KGE :  1` + `Cal : KGE : -1`) / 2,
         `Eval : KGEcomp` = (`Eval : KGE :  1` + `Eval : KGE : -1`) / 2)

decision_medians_WA <- Eval_FUSE_long_WA %>%
  group_by(ModelDecisions) %>%
  summarise(median_Eval_KGEcomp = median(`Eval : KGEcomp`, na.rm = TRUE))

best_decision_WA <- decision_medians_WA %>%
  filter(median_Eval_KGEcomp == max(median_Eval_KGEcomp)) %>%
  pull(ModelDecisions)

# FUSE WA
Sim_FUSE_WA = loadRData(file.path(dir_FUSE,"01_Paper", "WA", "SimArray.Rdata"))

Sim_FUSE_WA_sub = Sim_FUSE_WA[,as.character(best_decision_WA),]

rm(Sim_FUSE_WA)


##-----------------------------------------
##---------------- Mosa -------------------
##-----------------------------------------
df_eval_gumboot = read.table(file.path(dir_FUSE,"01_Paper","SamplingUncertainty","img","brief_analysis_modeling_results_with_gumboot_eval.csv"),
                             header = TRUE, sep = ",")
df_eval_gumboot$Codes = ifelse(nchar(df_eval_gumboot$gauge_id) == 7, paste0("USA_0", df_eval_gumboot$gauge_id), paste0("USA_", df_eval_gumboot$gauge_id))


dir_FUSE = file.path(path_res, 'FUSE')

Sim_FUSE_HF = loadRData(file.path(dir_FUSE,"01_Paper","1", "SimArray.Rdata"))
dimnames(Sim_FUSE_HF)[["Decisions"]] = paste0("HF_", dimnames(Sim_FUSE_HF)[["Decisions"]])

Sim_FUSE_LF = loadRData(file.path(dir_FUSE,"01_Paper","-1", "SimArray.Rdata"))
dimnames(Sim_FUSE_LF)[["Decisions"]] = paste0("LF_", dimnames(Sim_FUSE_LF)[["Decisions"]])

Sim_FUSE <- abind(Sim_FUSE_HF, Sim_FUSE_LF, along = 2)

rm(Sim_FUSE_HF,Sim_FUSE_LF)

df_eval_gumboot <- df_eval_gumboot[match(dimnames(Sim_FUSE)[[3]], df_eval_gumboot$Codes), ]

mod_need_sorted <- df_eval_gumboot$mod_need
mod_need_indices <- match(mod_need_sorted, dimnames(Sim_FUSE)[[2]])

Sim_Mosa <- matrix(NA, nrow = length(dimnames(Sim_FUSE)[[1]]), ncol = length(dimnames(Sim_FUSE)[[3]]))
rownames(Sim_Mosa) = dimnames(Sim_FUSE)[[1]]
colnames(Sim_Mosa) = dimnames(Sim_FUSE)[[3]]

valid_indices <- which(!is.na(mod_need_indices))

for (i in valid_indices) {
  Sim_Mosa[, i] <- Sim_FUSE[, mod_need_indices[i], i]
}


##-----------------------------------------
##----------------- SU --------------------
##-----------------------------------------
Ind_Eval = which(as.Date(rownames(Sim_Mosa)) >= as.Date("1999-01-01"))
catchments = colnames(Sim_Mosa)

su_HQ = data.frame(matrix(nrow = length(catchments), ncol = 10))
rownames(su_HQ) = catchments
colnames(su_HQ) = c("GOF_stat", "seJack", "seBoot", "p05", "p50", "p95", "score", "biasJack", "biasBoot", "seJab")

su_LQ = data.frame(matrix(nrow = length(catchments), ncol = 10))
rownames(su_LQ) = catchments
colnames(su_LQ) = c("GOF_stat", "seJack", "seBoot", "p05", "p50", "p95", "score", "biasJack", "biasBoot", "seJab")

su_WA = data.frame(matrix(nrow = length(catchments), ncol = 10))
rownames(su_WA) = catchments
colnames(su_WA) = c("GOF_stat", "seJack", "seBoot", "p05", "p50", "p95", "score", "biasJack", "biasBoot", "seJab")

su_Mosa = data.frame(matrix(nrow = length(catchments), ncol = 10))
rownames(su_Mosa) = catchments
colnames(su_Mosa) = c("GOF_stat", "seJack", "seBoot", "p05", "p50", "p95", "score", "biasJack", "biasBoot", "seJab")

su_Comp = data.frame(matrix(nrow = length(catchments), ncol = 10))
rownames(su_Comp) = catchments
colnames(su_Comp) = c("GOF_stat", "seJack", "seBoot", "p05", "p50", "p95", "score", "biasJack", "biasBoot", "seJab")

for(catchment in catchments){
  
  print(catchment)
  
  flow_HQ <- data.frame(date = as.Date(rownames(Sim_FUSE_HQ_sub)[Ind_Eval]), obs = Obs[Ind_Eval,catchment], sim = Sim_FUSE_HQ_sub[Ind_Eval, catchment])
  flow_LQ <- data.frame(date = as.Date(rownames(Sim_FUSE_LQ_sub)[Ind_Eval]), obs = Obs[Ind_Eval,catchment], sim = Sim_FUSE_LQ_sub[Ind_Eval, catchment])
  flow_WA <- data.frame(date = as.Date(rownames(Sim_FUSE_WA_sub)[Ind_Eval]), obs = Obs[Ind_Eval,catchment], sim = Sim_FUSE_WA_sub[Ind_Eval, catchment])
  flow_Mosa <- data.frame(date = as.Date(rownames(Sim_Mosa)[Ind_Eval]), obs = Obs[Ind_Eval,catchment], sim = Sim_Mosa[Ind_Eval, catchment])
  flow_Comp <- data.frame(date = as.Date(rownames(Sim_FUSE_Comp_sub)[Ind_Eval]), obs = Obs[Ind_Eval,catchment], sim = Sim_FUSE_Comp_sub[Ind_Eval, catchment])
  
  su_HQ[catchment,]  <- bootjack_v2(flow_HQ, GOF_stat = "KGEcomp", seed = 1)
  su_LQ[catchment,]  <- bootjack_v2(flow_LQ, GOF_stat = "KGEcomp", seed = 1)
  su_WA[catchment,] <- bootjack_v2(flow_WA, GOF_stat = "KGEcomp", seed = 1)
  su_Mosa[catchment,]  <- bootjack_v2(flow_Mosa, GOF_stat = "KGEcomp", seed = 1)
  su_Comp[catchment,]  <- bootjack_v2(flow_Comp, GOF_stat = "KGEcomp", seed = 1)
  
}

save(su_HQ, file = file.path(dir_FUSE, "01_Paper", "SU_HQ.Rdata"))
save(su_LQ, file = file.path(dir_FUSE, "01_Paper", "SU_LQ.Rdata"))
save(su_WA, file = file.path(dir_FUSE, "01_Paper", "SU_WA.Rdata"))
save(su_Mosa, file = file.path(dir_FUSE, "01_Paper", "SU_Mosa.Rdata"))
save(su_Comp, file = file.path(dir_FUSE, "01_Paper", "SU_Comp.Rdata"))

