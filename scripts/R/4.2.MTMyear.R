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

# make the folder
baseFolder <- "midstream/4.2.MTMyear/"
# baseFolder <- "C:/Users/biometrics/Desktop/4.2.MTMyear/"
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
VIsName <- c("NDVI", "NDRE")
VIsName <- c("NDVI")

# set the treatment
treatment <- c("C", "W10", "D10", "D")
treatment <- c("W10", "D10", "D")

# select the use day
useDayList <- list(c("2019_0831|2020_0827|2021_0827"), 
                   c("2019_0903|2020_0829|2021_0829"), 
                   c("2019_0904|2020_0901|2021_0830"))

dataSel0 <- c("2019_0831", "2019_0903", "2019_0904", 
              "2020_0827", "2020_0829", "2020_0901", 
              "2021_0827", "2021_0829", "2021_0830")
dataSel <- paste(dataSel0, collapse = "|")

#### calculate the MTM ######
if (length(list.files(baseFolder)) < 9) {
  resultAll <- foreach(treatmentInd = 1:length(treatment), .packages = c("doParallel")) %do% {
    # treatmentInd <- 1
    # make the save folder
    eachTreatment <- treatment[treatmentInd]
    matNow <- valueMat[, grepl(pattern = paste0("_", eachTreatment, "_"), x = colnames(valueMat))]
    # colnames(matNow)
    
    # select the VI
    for (ViInd in 1:length(VIsName)) {
      # ViInd <- 1
      eachVi <- VIsName[ViInd]
      selMat <- matNow[, grepl(pattern = eachVi, x = colnames(matNow))]
      # colnames(selMat)
      selMat <- selMat[, grepl(pattern = dataSel, x = colnames(selMat))]
      
      # save the correlation between each day"s VI value
      corMat <- cor(selMat)
      write.csv(x = corMat, file = paste0(baseFolder, "corMat", eachTreatment, "_", eachVi, ".csv"))
      colnames(corMat) <- rownames(corMat) <- NULL
      
      png(paste0(baseFolder, "corPlot", eachTreatment, "_", eachVi, ".png"), 
          width = 1440, height = 1440, res = 214)
      corrplot(corMat)
      dev.off()
      
      selMatNow <- selMat
      
      # make the matrix to input the result
      # row means the test year, col means the training year
      resMat <- matrix(NA, nrow = length(year), ncol = length(year))
      rownames(resMat) <- colnames(resMat) <- year
      
      yearUse <- unique(str_sub(colnames(selMatNow), 1, 4))
      
      saveFolder <- paste0(baseFolder, eachTreatment, "/All/")
      if (!dir.exists(saveFolder)) {
        dir.create(saveFolder, recursive = T)
      }
      
      
      for (trainYear in yearUse) {
        # trainYear <- yearUse[3]
        trainMat <- selMatNow[, grepl(pattern = trainYear, x = colnames(selMatNow)), drop = F]
        # trainMat <- scale(trainMat)
        # colnames(trainMat)
        for (testYear in yearUse) {
          # testYear <- yearUse[2]
          if (trainYear == testYear) {
            next
          }
          testMat <- selMatNow[, grepl(pattern = testYear, x = colnames(selMatNow)), drop = F]
          # testMat <- scale(testMat)
          # colnames(testMat)
          
          # prepare the data for MTM
          # set the GRM
          amat <- amat0[rownames(selMatNow), rownames(selMatNow)]
          
          # structure of the variance-covariance matrix of genotypic values
          n_trait <- ncol(trainMat) + 1
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
          MTM_Y <- matrix(NA, nrow = nrow(selMatNow), ncol = (ncol(trainMat) + 1))
          rownames(MTM_Y) <- rownames(selMatNow)
          colnames(MTM_Y) <- c("CV", paste0("Day", 1:3))
          MTM_Y[, "CV"] <- cv[rownames(MTM_Y), "CV"]
          target <- MTM_Y[, "CV"]
          # pairs(cbind(testMat, trainMat))
          
          # make the case for the result
          resultEachSeed <- matrix(NA, ncol = length(resultIndex), nrow = length(seedInd))
          colnames(resultEachSeed) <- resultIndex
          rownames(resultEachSeed) <- seedInd
          
          rep10 <- rep(1:10, length(target))[1:length(target)]
          # rep5 <- rep(1:5, length(target))[1:length(target)]
          
          cl <- makeCluster(10)
          registerDoParallel(cl)
          resList <- foreach(seedIndEach = 1:length(seedInd)) %dopar% {
            # for (seedIndEach in 1:length(seedInd)) {
            # seedIndEach <- 1
            set.seed(seedInd[seedIndEach])
            crossCVInd <- sample(rep10, length(target), replace = F)
            predictionData <- rep(NA, length(target))
            names(predictionData) <- rownames(trainMat)
            
            for (times in 1:10) {
              # times <- 1
              MTM_XF <- NULL
              
              # enter the trainnig data
              MTM_Y_NA <- MTM_Y
              MTM_Y_NA[crossCVInd != times, 2:ncol(MTM_Y_NA)] <- scale(trainMat[crossCVInd != times, ])
              
              # enter the test data
              MTM_Y_NA[crossCVInd == times, 2:ncol(MTM_Y_NA)] <- scale(testMat[crossCVInd == times, ])
              MTM_Y_NA[crossCVInd == times, "CV"] <- NA
              
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
              # predictionData[crossCVInd == times] <- res_MTM$YHat[, "target"][crossCVInd == times]
              if (is.null(MTM_XF)) {
                predictEach <- (res_MTM$mu + res_MTM$K[[1]]$U)[, 1]
              } else {
                predictEach <- (res_MTM$mu + MTM_XF %*% res_MTM$B.f + res_MTM$K[[1]]$U)[, 1]
              }
              predictionData[crossCVInd == times] <- predictEach[crossCVInd == times]
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

# read the result of All day model within each year
resultWithinPath <- paste0("midstream/4.1.MTM/dfAll.csv")
resultWithin <- read.csv(resultWithinPath, header = T, row.names = 1)
resultWithin <- na.omit(resultWithin)
resultWithin <- resultWithin[!(resultWithin$Treatment == "W1"),  ]
resultWithin <- resultWithin[resultWithin$VI == "NDVI",  ]
resultWithin[(resultWithin$Treatment == "W2" & resultWithin$Year == 2019),  "Treatment"] <- "W5"
resultWithin[(resultWithin$Treatment == "W3" & resultWithin$Year == 2019),  "Treatment"] <- "W2"

resultWithin$Treatment[resultWithin$Treatment == "W2"] <- "D10"
resultWithin$Treatment[resultWithin$Treatment == "W3"] <- "W10"
resultWithin$Treatment[resultWithin$Treatment == "W4"] <- "D"

resultOld <- data.frame(resultWithin, resultWithin$Year)
colnames(resultOld) <- c("Accuracy", "SE", "Treatment", "VI", "Train", "Test")

df <- NULL
for (eachTreatmentInd in 1:length(treatment)) {
  # eachTreatmentInd <- 1
  # make the save folder
  eachTreatment <- treatment[eachTreatmentInd]
  treatmentFolder <- paste0(baseFolder, eachTreatment, "/All/")
  
  
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
  df <- rbind(df, dfNew)
}

dfRR <- rbind(df, resultOld)
dfRR <- na.omit(dfRR)
write.csv(dfRR, paste0(baseFolder, "dfAll.csv"))
dfRR$Treatment <- factor(dfRR$Treatment, levels = c("W5", "W10", "D10", "D"))

g <- ggplot(dfRR, aes(x = Train, y = Accuracy, fill = Treatment)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c("olivedrab3", "lightskyblue1", "darkgoldenrod1", "chocolate2")) + 
  facet_wrap(~Test, nrow = 3) + 
  geom_hline(aes(yintercept = resultG), linetype = "dashed")
png(paste0(baseFolder, "resNDVI.png"), 
    width = 1440, height = 1440, res = 214)
print(g)
dev.off()
