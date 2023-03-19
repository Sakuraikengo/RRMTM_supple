machine <- c("M100", "M100", "P4M")
year <- c("2019", "2020", "2021")
source("scripts/R/1.functionCode.R")
source("scripts/R/1.MTM_2.R")

#' # 0. Read packages
options(stringsAsFactors = FALSE)
library(stringr)
library(psych)
library(ggplot2)
library(corrplot)
library(tagcloud)
library(doParallel)
library(asreml)

# make the folder
baseFolder <- "midstream/5.2.RRMMTMyear/"
# baseFolder <- "C:/Users/biometrics/Desktop/5.2.RRMMTMyear/"
if (!dir.exists(baseFolder)) {
  dir.create(baseFolder)
}

# read the genome data
amat0 <- as.matrix(read.csv("data/genome/amat173583SNP.csv", row.names = 1, header = T))
colnames(amat0)[colnames(amat0) == "X5002T"] <- "5002T"
colnames(amat0)[colnames(amat0) == "HOUJAKU_KUWAZU"] <- "Houjaku Kuwazu"
rownames(amat0)[rownames(amat0) == "HOUJAKU_KUWAZU"] <- "Houjaku Kuwazu"

# read the CV of FreshWeight
cv <- read.csv("midstream/2.0.freshWeightCV/freshWeightCV.csv", 
               header = T, row.names = 1)
colnames(cv) <- "CV"
cv <- scale(cv)

# read the varieties which will be used
selVariety0 <- read.csv("midstream/2.0.freshWeightCV/selVariety.csv", 
                        row.names = 1)
selVariety <- selVariety0[, 1]

# set the root folder
rootFolder <- "data/MSdata/0.3.dataBind/"

# save or read the seed
seedIndCsv <- paste0("midstream/seedInd.csv")
if (file.exists(seedIndCsv)) {
  seedInd <- read.csv(paste0(seedIndCsv), row.names = 1, header = T)
  seedInd <- c(as.matrix(seedInd))
} else {
  seedInd <- sample(1:500, 10, replace = F)
  write.csv(x = seedInd, file = paste0(seedIndCsv))
}

#### read the legendre parameters all days ####
pathList <- list.files(path = rootFolder, pattern = ".csv", recursive = T)
treatment <- c("W1", "W2", "W3", "W4")

# choose the VIs name
VIsName <- c("NDVI", "NDRE")

matList <- lapply(pathList, function(eachPath) {
  # eachPath <- pathList[[1]]
  # get the file data
  yearEach <- str_split(eachPath, pattern = "/")[[1]][1]
  dayEach <- str_sub(eachPath, start = -20, end = -17)
  eachDf <- read.csv(paste0(rootFolder, eachPath), header = T, row.names = 1)
  treatmentRep <- rep(treatment, length(VIsName))
  VIsRep <- rep(VIsName, each = length(treatment))
  
  colName <- paste0(yearEach, "_", 
                    dayEach, "_", 
                    treatmentRep, "_", 
                    VIsRep)
  
  # make the matrix to input the data
  eachMat <- matrix(NA, nrow = length(selVariety), 
                    ncol = length(colName))
  rownames(eachMat) <- selVariety
  colnames(eachMat) <- colName
  
  # fill the values
  for (eachTraitInd in 1:length(colName)) {
    # eachTraitInd <- 1
    treatmentNow <- treatmentRep[eachTraitInd]
    VIsNow <- VIsRep[eachTraitInd]
    
    eachDfSel <- eachDf[eachDf$treatment == treatmentNow, c("variety", VIsNow)]
    eachDfSel <- na.omit(eachDfSel)
    rownames(eachDfSel) <- eachDfSel$variety
    eachMat[, eachTraitInd] <- eachDfSel[selVariety, VIsNow]
  }
  
  if (yearEach == 2019) {
    colnames(eachMat) <- gsub(pattern = "W2", replacement = "W5", x = colnames(eachMat))
    colnames(eachMat) <- gsub(pattern = "W3", replacement = "W2", x = colnames(eachMat))
  }
  
  return(eachMat)
})

valueMat <- do.call(what = cbind, args = matList)
colnames(valueMat) <- gsub(pattern = "W1", replacement = "C", x = colnames(valueMat))
colnames(valueMat) <- gsub(pattern = "W2", replacement = "D10", x = colnames(valueMat))
colnames(valueMat) <- gsub(pattern = "W3", replacement = "W10", x = colnames(valueMat))
colnames(valueMat) <- gsub(pattern = "W4", replacement = "D", x = colnames(valueMat))
# valueMat <- scale(valueMat)

resultIndex <- c("Correlation", "R2", "RMSE")

# choose the VIs name
VIsName <- c("NDVI")

# set the treatment
treatment <- c("W10", "D10", "D")

sowingDayVec <- c(as.Date("2019-07-10"), 
                  as.Date("2020-07-08"), 
                  as.Date("2021-07-06"))

dataSelVec <- c(as.Date("2019-08-25"), as.Date("2019-08-31"), 
                as.Date("2019-09-03"), as.Date("2019-09-04"), 
                as.Date("2020-08-22"), as.Date("2020-08-24"), 
                as.Date("2020-08-27"), as.Date("2020-08-29"), 
                as.Date("2020-09-01"), as.Date("2020-09-04"), 
                as.Date("2021-08-22"), as.Date("2021-08-25"), 
                as.Date("2021-08-27"), as.Date("2021-08-29"), 
                as.Date("2021-08-30"), as.Date("2021-08-31"), 
                as.Date("2021-09-05"))

# select the use day
dataSel0 <- c("2019_0825", "2019_0831", "2019_0903", "2019_0904", 
              "2020_0822", "2020_0824", "2020_0827", "2020_0829", "2020_0901", "2020_0904", 
              "2021_0822", "2021_0825", "2021_0827", "2021_0829", 
              "2021_0830", "2021_0831", "2021_0905")
dataSel <- paste(dataSel0, collapse = "|")


#### calculate the MTM ######
if (length(list.files(baseFolder)) < 4) {
  for (treatmentInd in 1:length(treatment)) {
    # treatmentInd <- 1
    # make the save folder
    eachTreatment <- treatment[treatmentInd]
    matNow <- valueMat[, grepl(pattern = paste0("_", eachTreatment, "_"), x = colnames(valueMat))]
    # colnames(matNow)
    
    saveFolder <- paste0(baseFolder, eachTreatment, "/")
    if (!dir.exists(saveFolder)) {
      dir.create(saveFolder, recursive = T)
    }
    
    # select the VI
    for (ViInd in 1:length(VIsName)) {
      # ViInd <- 1
      eachVi <- VIsName[ViInd]
      selMat <- matNow[, grepl(pattern = eachVi, x = colnames(matNow))]
      # colnames(selMat)
      selMat <- selMat[, grepl(pattern = dataSel, x = colnames(selMat))]
      
      # make the matrix to input the result
      # row means the test year, col means the training year
      resMat <- matrix(NA, nrow = length(year), ncol = length(year))
      rownames(resMat) <- colnames(resMat) <- year
      
      yearUse <- unique(str_sub(colnames(selMat), 1, 4))
      
      for (trainYear in yearUse) {
        # trainYear <- yearUse[3]
        trainMat <- selMat[, grepl(pattern = trainYear, x = colnames(selMat)), drop = F]
        sowingDayTrain <- sowingDayVec[grep(pattern = trainYear, x = sowingDayVec)]
        useDayTrain <- dataSelVec[grep(pattern = trainYear, x = dataSelVec)]
        dayTrain <- as.numeric(useDayTrain - sowingDayTrain)
        for (testYear in yearUse) {
          # testYear <- yearUse[2]
          if (trainYear == testYear) {
            next
          }
          testMat <- selMat[, grepl(pattern = testYear, x = colnames(selMat)), drop = F]
          sowingDayTest <- sowingDayVec[grep(pattern = testYear, x = sowingDayVec)]
          useDayTest <- dataSelVec[grep(pattern = testYear, x = dataSelVec)]
          dayTest <- as.numeric(useDayTest - sowingDayTest)
          
          # prepare the data for MTM
          # set the GRM
          amat <- amat0[rownames(selMat), rownames(selMat)]
          
          # structure of the variance-covariance matrix of genotypic values
          n_trait <- 3
          MTM_K <- list(list(
            K = amat,
            COV = list(type = 'UN',
                       df0 = n_trait,
                       S0 = diag(n_trait))
          ))
          
          # structure of the variance-covariance matrix of residual effects
          MTM_resCov <- list(type = 'UN',
                             df0 = n_trait,
                             S0 = diag(n_trait))
          
          # set the phenotype data (only CV)
          MTM_Y <- matrix(NA, nrow = nrow(selMat), ncol = 3)
          rownames(MTM_Y) <- rownames(selMat)
          colnames(MTM_Y) <- c("CV", paste0("L_", 0:1))
          MTM_Y[, "CV"] <- cv[rownames(MTM_Y), "CV"]
          target <- MTM_Y[, "CV"]
          
          
          # make the case for the result
          resultEachSeed <- matrix(NA, ncol = length(resultIndex), nrow = length(seedInd))
          colnames(resultEachSeed) <- resultIndex
          rownames(resultEachSeed) <- seedInd
          
          rep10 <- rep(1:10, length(target))[1:length(target)]
          
          # make all trainning and test data
          modelDataList <- NULL
          for (seedIndEach in 1:length(seedInd)) {
            # seedIndEach <- 1
            set.seed(seedInd[seedIndEach])
            crossCvInd <- sample(rep10, length(target), replace = F)
            
            for (times in 1:10) {
              # times <- 1
              
              # make the RR model (training data)
              trainMatNow <- trainMat[crossCvInd != times, ]
              trainValue <- MakeRRMLinear(matWide = trainMatNow, day = dayTrain, k = 1)
              
              testMatNow <- testMat[crossCvInd == times, ]
              testValue <- MakeRRMLinear(matWide = testMatNow, day = dayTest, k = 1)
              modelDataList[[((seedIndEach - 1) * 10) + times]] <- list(trainValue, testValue)
            }
          }
          length(modelDataList)
          cl <- makeCluster(10)
          registerDoParallel(cl)
          resList <- foreach(seedIndEach = 1:length(seedInd)) %dopar% {
            # for (seedIndEach in 1:length(seedInd)) {
            # seedIndEach <- 1
            set.seed(seedInd[seedIndEach])
            crossCvInd <- sample(rep10, length(target), replace = F)
            predictionData <- rep(NA, length(target))
            names(predictionData) <- rownames(trainMat)
            
            for (times in 1:10) {
              # times <- 1
              MTM_XF <- NULL
              dataInd <- ((seedIndEach - 1) * 10) + times
              trainValue <- modelDataList[[dataInd]][[1]]
              testValue <- modelDataList[[dataInd]][[2]]
              
              # enter the trainnig data
              MTM_Y_NA <- MTM_Y
              MTM_Y_NA[rownames(trainValue), 2:ncol(MTM_Y_NA)] <- trainValue
              
              # enter the test data
              MTM_Y_NA[rownames(testValue), 2:ncol(MTM_Y_NA)] <- testValue
              MTM_Y_NA[crossCvInd == times, "CV"] <- NA
              
              # # delete the .dat files
              # datFiles <- list.files(saveFolder, pattern = ".dat", full.names = T)
              # if (length(datFiles) != 0) {
              #   file.remove(datFiles)
              # }
              
              # estimation (for RMarkdown)
              MTM_output <- capture.output(
                res_MTM <- MTM_2(XF = MTM_XF, Y = MTM_Y_NA, K = MTM_K, resCov = MTM_resCov,
                                 nIter = 12000, burnIn = 2000, thin = 20,
                                 saveAt = paste0(saveFolder, 
                                                 "train", trainYear, "_", 
                                                 "test", testYear, "_", 
                                                 eachVi, "_", 
                                                 seedIndEach, "_", 
                                                 times, "_"))
              )
              
              # put into the predicted value
              # predictionData[crossCvInd == times] <- res_MTM$YHat[, "target"][crossCvInd == times]
              if (is.null(MTM_XF)) {
                predictEach <- (res_MTM$mu + res_MTM$K[[1]]$U)[, 1]
              } else {
                predictEach <- (res_MTM$mu + MTM_XF %*% res_MTM$B.f + res_MTM$K[[1]]$U)[, 1]
              }
              predictionData[crossCvInd == times] <- predictEach[crossCvInd == times]
            }
            
            predictData <- predictionData
            obsData <- target
            
            # calculate the R2 and RMSE
            correlation <- cor(obsData, predictData)
            R2 <- 1 - sum((obsData - predictData) ^ 2) / sum((obsData - mean(obsData)) ^ 2)
            RMSE <- sqrt(sum((obsData - predictData) ^ 2) / length(obsData))
            
            # # input the result
            # resultEachSeed[seedIndEach, ] <- c(correlation, R2, RMSE)
            # 
            xlim <- ylim <- range(predictData, obsData)
            png(paste0(saveFolder, "train", trainYear, "_", 
                       "test", testYear, "_", 
                       eachVi, 
                       "_", seedIndEach, ".png"))
            plot(obsData, predictData, 
                 main = paste0("train", trainYear, "_", 
                               "test", testYear, "_", eachVi, " r = ", 
                               correlation), 
                 xlim = xlim, ylim = ylim)
            abline(0, 1, col = 2, lty = 2)
            dev.off()
            
            return(c(correlation, R2, RMSE))
          }
          stopCluster(cl)
          
          resultEachSeed <- do.call(what = rbind, args = resList)
          write.csv(resultEachSeed, file = paste0(saveFolder, "all_res_train", 
                                                  trainYear, "test", testYear, "_", eachVi, ".csv"))
          resMat[trainYear, testYear] <- apply(resultEachSeed, 2, mean)[1]
        }
        write.csv(resMat, file = paste0(saveFolder, "result_", eachVi, ".csv"))
      }
    }
  }
}


######## visualize the each year results ############
# read the result of G model
resultPathG <- paste0("midstream/2.1.genomicPrediction/result.csv")
resultG <- read.csv(resultPathG, header = T, row.names = 1)[1, 1]

# read the result of All day model
resultAllDayPath <- paste0("midstream/4.2.MTMyear/dfAll.csv")
dfAllDay <- read.csv(resultAllDayPath, header = T, row.names = 1)
dfAllDay <- na.omit(dfAllDay)

# read the result of RRM within each year
resultWithinPath <- paste0("midstream/5.1.RRMMTM/dfAll.csv")
resultWithin <- read.csv(resultWithinPath, header = T, row.names = 1)
resultWithin <- na.omit(resultWithin)
resultWithin <- resultWithin[!(resultWithin$Treatment == "W1"),  ]
resultWithin <- resultWithin[resultWithin$VI == "NDVI",  ]
resultWithin[(resultWithin$Treatment == "W2" & resultWithin$Year == 2019),  "Treatment"] <- "W5"
resultWithin[(resultWithin$Treatment == "W3" & resultWithin$Year == 2019),  "Treatment"] <- "W2"

resultWithin$Treatment[resultWithin$Treatment == "W2"] <- "D10"
resultWithin$Treatment[resultWithin$Treatment == "W3"] <- "W10"
resultWithin$Treatment[resultWithin$Treatment == "W4"] <- "D"

dfLegWithin <- data.frame(resultWithin, resultWithin$Year)
colnames(dfLegWithin) <- c("Accuracy", "SE", "Treatment", "VI", "Train", "Test")

dfLegCross <- NULL
for (eachTreatmentInd in 1:length(treatment)) {
  # eachTreatmentInd <- 1
  # make the save folder
  eachTreatment <- treatment[eachTreatmentInd]
  treatmentFolder <- paste0(baseFolder, eachTreatment, "/")
  
  
  # extract the Name of VIs
  resList0 <- list.files(treatmentFolder, pattern = "all_res")
  trainYear <- str_sub(resList0, start = 14, end = 17)
  testYear <- str_sub(resList0, start = 22, end = 25)
  viNameEach <- str_sub(resList0, start = 27, end = -5)
  
  resList <- paste0(treatmentFolder, resList0)
  
  corList <- lapply(resList, function(eachResList) {
    # eachResList <- resList[[1]]
    eachResMat <- read.csv(eachResList, header = T, row.names = 1)
    correlation <- mean(eachResMat[, 1])
    standardDeviation <- sd(eachResMat[, 1])
    standardError <- standardDeviation / sqrt(nrow(eachResMat))
    return(c(correlation, standardError))
  })
  resMat <- do.call(what = rbind, args = corList)
  dfLong <- data.frame(resMat, 
                       Treatment = eachTreatment, 
                       VI = viNameEach, 
                       Train = trainYear, 
                       Test = testYear)
  colnames(dfLong)[1:2] <- c("Accuracy", "SE")
  write.csv(dfLong, file = paste0(treatmentFolder, "corDf.csv"))
  
  # bind the all data
  dfNew <- data.frame(dfLong)
  dfLegCross <- rbind(dfLegCross, dfNew)
}

dfLeg <- rbind(dfLegCross, dfLegWithin)

write.csv(dfLeg, paste0(baseFolder, "dfAll.csv"))
dfLeg$Treatment <- factor(dfLeg$Treatment, levels = c("W5", "W10", "D10", "D"))

# calculate the ratio, model using all day data, model using legendre parameters
if (all(dfLeg[, 3:6] == dfAllDay[, 3:6])) {
  ratio <- dfLeg[, 1] / dfAllDay[, 1]
  percentage <- (ratio - 1) * 100
  dfRatio <- cbind(Ratio = percentage, 
                   dfLeg[, 3:6])
}
dataList <- list(list(dfEach = dfLeg, key = "Leg"), 
                 list(dfEach = dfRatio, key = "Ratio"))

lapply(dataList, function(dfList) {
  # dfList <- dataList[[1]]
  dfEach <- dfList$dfEach
  key <- dfList$key
  
  # extract the year cross validation
  dfCross0 <- dfEach[dfEach$Train != dfEach$Test, ]
  trainTestX <- paste0(dfCross0$Train, "/", dfCross0$Test)
  trainTest <- 1:nrow(dfCross0)
  dfCross <- data.frame(dfCross0, trainTest, trainTestX)
  
  if (key == "Leg") {
    g <- ggplot(dfEach, aes(x = Train, y = Accuracy, fill = Treatment)) + 
      geom_bar(stat = "identity", position = "dodge") + 
      scale_fill_manual(values = c("olivedrab3", "lightskyblue1", "darkgoldenrod1", "chocolate2")) + 
      facet_wrap(~Test, nrow = 3) + 
      geom_errorbar(aes(ymin = Accuracy - SE, ymax = Accuracy + SE), 
                    position = position_dodge(0.9), width = .3) + 
      geom_hline(aes(yintercept = resultG), linetype = "dashed") + 
      theme(text = element_text(size = 15))
    png(paste0(baseFolder, key, "_NDVI.png"), 
        width = 1440, height = 1440, res = 214)
    print(g)
    dev.off()
    
    gCross <- ggplot(dfCross, aes(x = trainTest, y = Accuracy, fill = Treatment)) + 
      geom_bar(stat = "identity", position = "dodge") + 
      scale_fill_manual(values = c("lightskyblue1", "darkgoldenrod1", "chocolate2")) + 
      geom_errorbar(aes(ymin = Accuracy - SE, ymax = Accuracy + SE), 
                    position = position_dodge(0.9), width = .3) + 
      scale_x_continuous(labels = NULL, breaks = 1:nrow(dfCross)) + 
      xlab(NULL) + 
      geom_hline(aes(yintercept = resultG), linetype = "dashed") + 
      theme(text = element_text(size = 15))
    png(paste0(baseFolder, key, "_NDVI_cross.png"), 
        width = 1440, height = 1440, res = 214)
    print(gCross)
    dev.off()
    
    gCrossX <- ggplot(dfCross, aes(x = trainTestX, y = Accuracy, fill = Treatment)) + 
      geom_bar(stat = "identity", position = "dodge") + 
      scale_fill_manual(values = c("lightskyblue1", "darkgoldenrod1", "chocolate2")) + 
      geom_errorbar(aes(ymin = Accuracy - SE, ymax = Accuracy + SE), 
                    position = position_dodge(0.9), width = .3) + 
      xlab(NULL) + 
      geom_hline(aes(yintercept = resultG), linetype = "dashed") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1), 
            text = element_text(size = 15))  
    png(paste0(baseFolder, key, "_NDVI_crossX.png"), 
        width = 1440, height = 1440, res = 214)
    print(gCrossX)
    dev.off()
    
  } else if (key == "Ratio") {
    g <- ggplot(dfEach, aes(x = Train, y = Ratio, fill = Treatment)) + 
      geom_bar(stat = "identity", position = "dodge") + 
      scale_fill_manual(values = c("olivedrab3", "lightskyblue1", "darkgoldenrod1", "chocolate2")) + 
      facet_wrap(~Test, nrow = 3) + 
      labs(y = "Proportion of improvement (%)") + 
      geom_hline(aes(yintercept = 1), linetype = "dashed") + 
      theme(text = element_text(size = 15))
    png(paste0(baseFolder, key, "_NDVI.png"), 
        width = 1440, height = 1440, res = 214)
    print(g)
    dev.off()
    
    gCross <- ggplot(dfCross, aes(x = trainTest, y = Ratio, fill = Treatment)) + 
      geom_bar(stat = "identity", position = "dodge") + 
      scale_fill_manual(values = c("lightskyblue1", "darkgoldenrod1", "chocolate2")) + 
      scale_x_continuous(labels = NULL, breaks = 1:nrow(dfCross)) + 
      xlab(NULL) + 
      labs(y = "Proportion of improvement (%)") + 
      geom_hline(aes(yintercept = 0), linetype = "dashed") + 
      theme(text = element_text(size = 15))
    png(paste0(baseFolder, key, "_NDVI_cross.png"), 
        width = 1440, height = 1440, res = 214)
    print(gCross)
    dev.off()
    
    gCrossX <- ggplot(dfCross, aes(x = trainTestX, y = Ratio, fill = Treatment)) + 
      geom_bar(stat = "identity", position = "dodge") + 
      scale_fill_manual(values = c("lightskyblue1", "darkgoldenrod1", "chocolate2")) + 
      xlab(NULL) + 
      labs(y = "Proportion of improvement (%)") + 
      geom_hline(aes(yintercept = 0), linetype = "dashed") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1), 
            text = element_text(size = 15))  
    png(paste0(baseFolder, key, "_NDVI_crossX.png"), 
        width = 1440, height = 1440, res = 214)
    print(gCrossX)
    dev.off()
    
  }
})

# W102020_2021 <- dfLeg[dfLeg$Treatment == "W10" & dfLeg$Train == 2020 & dfLeg$Test == 2021, "Accuracy"]
# W102021_2020 <- dfLeg[dfLeg$Treatment == "W10" & dfLeg$Train == 2021 & dfLeg$Test == 2020, "Accuracy"]
# D102020_2021 <- dfLeg[dfLeg$Treatment == "D10" & dfLeg$Train == 2020 & dfLeg$Test == 2021, "Accuracy"]
# D102021_2020 <- dfLeg[dfLeg$Treatment == "D10" & dfLeg$Train == 2021 & dfLeg$Test == 2020, "Accuracy"]
# D2020_2021 <- dfLeg[dfLeg$Treatment == "D" & dfLeg$Train == 2020 & dfLeg$Test == 2021, "Accuracy"]
# D2021_2020 <- dfLeg[dfLeg$Treatment == "D" & dfLeg$Train == 2021 & dfLeg$Test == 2020, "Accuracy"]
# 
# (W102020_2021 + W102021_2020) / (D2020_2021 + D2021_2020)
# (D102020_2021 + D102021_2020) / (D2020_2021 + D2021_2020)
# 
# mean(c(W102020_2021, W102021_2020))
# mean(c(D102020_2021, D102021_2020))
# 
# dfRatioCross <- dfRatio[dfRatio$Train != dfRatio$Test, ]
# tapply(dfRatioCross$Ratio, dfRatioCross$Treatment, mean)
