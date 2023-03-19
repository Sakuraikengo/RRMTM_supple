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
baseFolder <- "midstream/4.3.MTMsmallData/"
# baseFolder <- "C:/Users/biometrics/Desktop/4.3.MTMsmallData/"
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

# set the seed for this code
seedIndCsv2 <- paste0("midstream/seedInd2.csv")
if (file.exists(seedIndCsv2)) {
  seedInd2 <- read.csv(paste0(seedIndCsv2), row.names = 1, header = T)
  seedInd2 <- c(as.matrix(seedInd2))
} else {
  seedInd2 <- sample(1:5000, 5, replace = F)
  write.csv(x = seedInd2, file = paste0(seedIndCsv2))
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
  return(eachMat)
})

valueMat <- do.call(what = cbind, args = matList)
# valueMat <- scale(valueMat)

resultIndex <- c("Correlation", "R2", "RMSE")

# choose the VIs name
VIsName <- c("NDVI")

# set the treatment
treatment <- c("W2", "W3", "W4")

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
valueMat <- valueMat[, grepl(pattern = dataSel, x = colnames(valueMat))]

# the percentage of lines which have train data
dataRatioVec <- c(10, 25, 50)


#### calculate the MTM ######
if (length(list.files(baseFolder)) < 3) {
  for (eachDataInd in 1:length(year)) {
    # eachDataInd <- 3
    eachYear <- year[eachDataInd]
    eachMachine <- machine[eachDataInd]
    
    matEachYear <- valueMat[, grepl(pattern = eachYear, x = colnames(valueMat))]
    # colnames(matEachYear)
    
    for (dataRatioInd in 1:length(dataRatioVec)) {
      # dataRatioInd <- 1
      dataRatio <- dataRatioVec[dataRatioInd]
      # make the save folder
      saveFolder <- paste0(baseFolder, eachYear, "/", eachMachine, "/train", dataRatio, "/")
      if (!dir.exists(saveFolder)) {
        dir.create(saveFolder, recursive = T)
      }
      
      # select the VI
      for (ViInd in 1:length(VIsName)) {
        # ViInd <- 1
        eachVi <- VIsName[ViInd]
        matEachVi <- matEachYear[, grepl(pattern = eachVi, x = colnames(matEachYear))]
        # colnames(matEachVi)
        
        # make the matrix to input the result
        resMat <- matrix(NA, nrow = length(treatment), ncol = length(resultIndex))
        rownames(resMat) <- treatment
        colnames(resMat) <- resultIndex
        
        for (treatmentInd in 1:length(treatment)) {
          # treatmentInd <- 2
          eachTreatment <- treatment[treatmentInd]
          resFile <- paste0(saveFolder, "all_res_", eachTreatment, "_", eachVi, ".csv")
          if (file.exists(resFile)) {
            resMat[eachTreatment, ] <- apply(read.csv(resFile, row.names = 1), 2, mean)
            next
          }
          
          matNow <-  matEachVi[, grepl(pattern = paste0("_", eachTreatment, "_"), x = colnames(matEachVi))]
          
          # colnames(matNow)
          
          # prepare the data for MTM
          # set the GRM
          amat <- amat0[rownames(matNow), rownames(matNow)]
          
          # structure of the variance-covariance matrix of genotypic values
          n_trait <- ncol(matNow) + 1
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
          MTM_Y <- matrix(NA, nrow = nrow(matNow), ncol = n_trait)
          rownames(MTM_Y) <- rownames(matNow)
          colnames(MTM_Y) <- c("CV", str_sub(colnames(matNow), 6, 9))
          MTM_Y[, "CV"] <- cv[rownames(MTM_Y), "CV"]
          target <- MTM_Y[, "CV"]
          
          
          rep10 <- rep(1:10, length(target))[1:length(target)]
          
          cl <- makeCluster(10)
          registerDoParallel(cl)
          resList <- foreach(seedIndEach = 1:length(seedInd)) %dopar% {
            # for (seedIndEach in 1:length(seedInd)) {
            # seedIndEach <- 1
            set.seed(seedInd[seedIndEach])
            crossCvInd <- sample(rep10, length(target), replace = F)
            predictionData <- rep(NA, length(target))
            names(predictionData) <- rownames(matNow)
            
            # make the case for the result
            resultEachSeed <- matrix(NA, nrow = length(resultIndex), 
                                     ncol = length(seedInd2))
            rownames(resultEachSeed) <- resultIndex
            colnames(resultEachSeed) <- paste0(seedInd[seedIndEach], "_", 
                                               seedInd2)
            
            for (seedIndEach2 in 1:length(seedInd2)) {
              # seedIndEach2 <- 1
              
              for (times in 1:10) {
                # times <- 1
                MTM_XF <- NULL
                
                # reduce the number of secondary trait data 
                trainMatNow <- matNow[crossCvInd != times, ]
                set.seed(seedInd2[seedIndEach2])
                useDataInd <- sample(1:(100 / dataRatio), nrow(trainMatNow), replace = T)
                trainMatUse <- trainMatNow[useDataInd == 1, ]
                
                # enter the trainnig data
                MTM_Y_NA <- MTM_Y
                MTM_Y_NA[rownames(trainMatUse), 2:ncol(MTM_Y_NA)] <- scale(trainMatUse)
                
                # enter the test data
                MTM_Y_NA[crossCvInd == times, 2:ncol(MTM_Y_NA)] <- scale(matNow)[crossCvInd == times, ]
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
                                                   eachVi, "_", 
                                                   eachTreatment, "_", 
                                                   seedIndEach, "_", 
                                                   seedIndEach2, "_", 
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
              
              # input the result
              resultEachSeed[, seedIndEach2] <- c(correlation, R2, RMSE)
              
              if (seedIndEach2 == 1) {
                xlim <- ylim <- range(predictData, obsData)
                png(paste0(saveFolder, eachVi, 
                           "_", eachTreatment, "_", seedIndEach, ".png"))
                plot(obsData, predictData, 
                     main = paste0(eachVi, " ", eachTreatment, " r = ", 
                                   correlation), 
                     xlim = xlim, ylim = ylim)
                abline(0, 1, col = 2, lty = 2)
                dev.off()
              }
              write.csv(resultEachSeed, file = paste0(saveFolder, eachVi, 
                                                      "_", eachTreatment, "_", seedIndEach, ".csv"))
            }
            return(c(correlation, R2, RMSE))
          }
          stopCluster(cl)
          
          resultAll <- do.call(what = rbind, args = resList)
          write.csv(resultAll, file = paste0(saveFolder, "all_res_", 
                                                  eachTreatment, "_", eachVi, ".csv"))
          resMat[eachTreatment, ] <- apply(resultAll, 2, mean)
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

# read the result of model using all data
resultPathAll <- paste0("midstream/4.1.MTM/dfAll.csv")
dfAllData0 <- read.csv(resultPathAll, header = T, row.names = 1)
dfAllData0 <- na.omit(dfAllData0)
dfAllData <- data.frame(dfAllData0, Ratio = 100)

df <- NULL
for (eachDataInd in 1:length(machine)) {
  # eachDataInd <- 1
  # make the save folder
  eachYear <- year[eachDataInd]
  eachMachine <- machine[eachDataInd]
  yearFolder <- paste0(baseFolder, eachYear, "/", eachMachine, "/")
  
  
  # extract the Name of VIs
  resList0 <- list.files(yearFolder, pattern = "all_res", recursive = T)
  treatmentEach <- str_sub(resList0, start = 17, end = -10)
  viNameEach <- str_sub(resList0, start = 20, end = -5)
  dataRatio <- str_sub(resList0, start = 6, end = 7)
  
  resList <- paste0(yearFolder, resList0)
  
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
                       Treatment = treatmentEach, 
                       VI = viNameEach, 
                       Year = eachYear, 
                       Ratio = dataRatio)
  colnames(dfLong)[1:2] <- c("Accuracy", "SE")
  write.csv(dfLong, file = paste0(yearFolder, "corDf.csv"))
  
  # bind the all data
  df <- rbind(df, dfLong)
}

write.csv(df, paste0(baseFolder, "dfAll.csv"))
df <- na.omit(df)

# calculate the ratio, model using all day data, model using legendre parameters
dfAll <- rbind(df, dfAllData)

dfAll <- dfAll[dfAll$VI == "NDVI", ]
dfAll <- dfAll[!(dfAll$Treatment == "W1"),  ]
dfAll[(dfAll$Treatment == "W2" & dfAll$Year == 2019),  "Treatment"] <- "W5"
dfAll[(dfAll$Treatment == "W3" & dfAll$Year == 2019),  "Treatment"] <- "W2"

dfAll$Treatment[dfAll$Treatment == "W2"] <- "D10"
dfAll$Treatment[dfAll$Treatment == "W3"] <- "W10"
dfAll$Treatment[dfAll$Treatment == "W4"] <- "D"
dfAll$Treatment <- factor(dfAll$Treatment, levels = c("W5", "W10", "D10", "D"))
dfAll$Ratio <- as.numeric(dfAll$Ratio)

g <- ggplot(dfAll, aes(x = Ratio, y = Accuracy, color = Treatment)) + 
  geom_line(size = 1) + 
  geom_errorbar(aes(ymin = Accuracy - SE, ymax = Accuracy + SE), 
                width = .5) + 
  scale_color_manual(values = c("olivedrab3", "lightskyblue1", "darkgoldenrod1", "chocolate2")) + 
  facet_wrap(~Year, nrow = 3) + 
  scale_y_continuous(limits = c(0.27, 0.53)) + 
  labs(x = "Ratio (%)") + 
  geom_hline(aes(yintercept = resultG), linetype = "dashed")
png(paste0(baseFolder, "diffRatio_NDVI.png"), 
    width = 1440, height = 1440, res = 214)
print(g)
dev.off()

# dfAll$Accuracy[dfAll$Treatment == "W10" & dfAll$Ratio == 100] / dfAll$Accuracy[dfAll$Treatment == "W10" & dfAll$Ratio == 10]
# dfAll$Accuracy[dfAll$Treatment == "D10" & dfAll$Ratio == 100] / dfAll$Accuracy[dfAll$Treatment == "D10" & dfAll$Ratio == 10]
# 
# dfAll[dfAll$Treatment == "D10" & dfAll$Ratio == 100, ]
