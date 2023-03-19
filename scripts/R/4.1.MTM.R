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
baseFolder <- "midstream/4.1.MTM/"
# baseFolder <- "C:/Users/biometrics/Desktop/4.1.MTM/"
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

# set the root folder
rootFolder <- "data/MSdata/0.3.dataBind/"

# select the day based on the weather (only sunny & 11:00-13:00)
selDayList <- list(c("0825", "0831", "0903", "0904"), 
                   c("0822", "0824", "0827", "0829", "0901", "0904"), 
                   c("0822", "0825", "0827", "0829", "0830", "0831", "0905"))


# save or read the seed
seedIndCsv <- paste0("midstream/seedInd.csv")
if (file.exists(seedIndCsv)) {
  seedInd <- read.csv(paste0(seedIndCsv), row.names = 1, header = T)
  seedInd <- c(as.matrix(seedInd))
} else {
  seedInd <- sample(1:500, 10, replace = F)
  write.csv(x = seedInd, file = paste0(seedIndCsv))
}

treatment <- c("W1", "W2", "W3", "W4")
resultIndex <- c("Correlation", "R2", "RMSE")

# choose the VIs name
VIsName <- c("NDVI", "NDRE")

# set the plot number
condition <- c("C", "W5", "W10", "D")

#### calculate the MTM ######
if (length(list.files(baseFolder)) < 3) {
  resultAll <- foreach(eachDataInd = 1:length(machine), .packages = c("doParallel")) %do% {
    # eachDataInd <- 2
    # make the save folder
    eachYear <- year[eachDataInd]
    eachMachine <- machine[eachDataInd]
    saveFolder <- paste0(baseFolder, eachYear, "/", eachMachine, "/")
    if (!dir.exists(saveFolder)) {
      dir.create(saveFolder, recursive = T)
    }
    
    # read the MS data
    dataPath <- paste0(rootFolder, eachYear, "/", eachMachine, "/")
    dataFiles0 <- list.files(dataPath, full.names = T, pattern = ".csv")
    selDays <- selDayList[[eachDataInd]]
    selDaysVec <- paste(selDays, collapse = "|")
    dataFiles <- dataFiles0[grep(pattern = selDaysVec, x = dataFiles0)]
    
    # read the all day files
    dataList <- lapply(dataFiles, function(eachDataFile) {
      # eachDataFile <- dataFiles[1]
      # extract the day 
      eachDayInd <- regexpr("0[0-9]{3}", eachDataFile)[1]
      eachDay <- str_sub(eachDataFile, start = eachDayInd, end = (eachDayInd + 3))
      
      # read the data
      eachData <- read.csv(eachDataFile, header = T, row.names = 1)
      colnames(eachData)[4:ncol(eachData)] <- paste0(colnames(eachData)[4:ncol(eachData)], "_", eachDay)
      return(eachData)
    })
    dataAll <- do.call(cbind, args = dataList)
    
    # select the VI
    for (ViInd in 1:length(VIsName)) {
      # ViInd <- 1
      eachVi <- VIsName[ViInd]
      
      # make the matrix to input the result
      resMat <- matrix(NA, nrow = length(treatment), ncol = length(resultIndex))
      rownames(resMat) <- treatment
      colnames(resMat) <- resultIndex
      
      for (treatmentInd in 1:length(treatment)) {
        # treatmentInd <- 1
        eachTreatment <- treatment[treatmentInd]
        eachCondition <- condition[treatmentInd]
        selDataInd <- paste(eachVi, selDays, sep = "_")
        selDataInd <- selDataInd[selDataInd %in% colnames(dataAll)]
        selData <- dataAll[dataAll$treatment == eachTreatment, 
                           c("variety", selDataInd)]
        selData <- na.omit(selData)
        rownames(selData) <- selData$variety
        if (dim(selData)[1] < 10) {
          next
        }
        
        df <- data.frame(CV = cv, 
                         selData[rownames(cv), 2:ncol(selData)])
        n_trait <- ncol(df)
        # set the GRM
        amat <- amat0[rownames(df), rownames(df)]
        
        # Y
        MTM_Y <- as.matrix(df)
        MTM_Y <- scale(MTM_Y)
        
        # structure of the variance-covariance matrix of genotypic values
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
        
        target <- MTM_Y[, "CV"]
        
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
          
          for (times in 1:10) {
            # times <- 1
            MTM_XF <- NULL
            # enter NA
            MTM_Y_Na <- MTM_Y
            MTM_Y_Na[crossCVInd == times, "CV"] <- NA
            
            # # delete the .dat files
            # datFiles <- list.files(saveFolder, pattern = ".dat", full.names = T)
            # if (length(datFiles) != 0) {
            #   file.remove(datFiles)
            # }
            
            # estimation (for RMarkdown)
            MTM_output <- capture.output(
              res_MTM <- MTM_2(XF = MTM_XF, Y = MTM_Y_Na, K = MTM_K, resCov = MTM_resCov,
                               nIter = 12000, burnIn = 2000, thin = 20,
                               saveAt = paste0(saveFolder, 
                                               eachTreatment, "_", 
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
          png(paste0(saveFolder, "/predictionPlot_", eachTreatment, "_", eachVi, 
                     "_", seedIndEach, ".png"))
          plot(obsData, predictData, 
               main = paste0(eachTreatment, "_", eachVi, " r = ", 
                             correlation), 
               xlim = xlim, ylim = ylim)
          abline(0, 1, col = 2, lty = 2)
          dev.off()
          
          return(c(correlation, R2, RMSE))
        }
        stopCluster(cl)
        
        resultEachSeed <- do.call(what = rbind, args = resList)
        write.csv(resultEachSeed, file = paste0(saveFolder, "all_res_", 
                                                eachTreatment, "_", eachVi, ".csv"))
        resMat[eachTreatment, ] <- apply(resultEachSeed, 2, mean)
      }
      write.csv(resMat, file = paste0(saveFolder, "result_", eachVi, ".csv"))
    }
  }
}

######## visualize the each year results ############
# read the result of G model
resultPathG <- paste0("midstream/2.1.genomicPrediction/result.csv")
resultG <- read.csv(resultPathG, header = T, row.names = 1)[1, 1]

df <- NULL
for (eachDataInd in 1:length(machine)) {
  # eachDataInd <- 1
  # make the save folder
  eachYear <- year[eachDataInd]
  eachMachine <- machine[eachDataInd]
  yearFolder <- paste0(baseFolder, eachYear, "/", eachMachine, "/")
  
  
  # extract the Name of VIs
  resList0 <- list.files(yearFolder, pattern = "all_res")
  treatmentEach <- str_sub(resList0, start = 9, end = -10)
  viNameEach <- str_sub(resList0, start = 12, end = -5)
  
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
                       VI = viNameEach)
  colnames(dfLong)[1:2] <- c("Accuracy", "SE")
  write.csv(dfLong, file = paste0(yearFolder, "corDf.csv"))
  
  pdf(paste0(yearFolder, "barPlot.pdf"))
  g <- ggplot(dfLong, aes(x = VI, y = Accuracy, fill = Treatment)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    scale_fill_manual(values = c("blue", "lightskyblue1", "darkgoldenrod1", "chocolate2"))
  print(g)
  
  dfLongDiv <- dfLong
  dfLongDiv$Accuracy <- dfLongDiv$Accuracy / resultG
  colnames(dfLongDiv)[1] <- "Ratio"
  
  g <- ggplot(dfLongDiv, aes(x = VI, y = Ratio, fill = Treatment)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    scale_fill_manual(values = c("blue", "lightskyblue1", "darkgoldenrod1", "chocolate2")) + 
    geom_hline(aes(yintercept = 1), linetype = "dashed")
  print(g)
  
  dev.off()
  
  # bind the all data
  dfNew <- data.frame(dfLong, Year = eachYear)
  df <- rbind(df, dfNew)
}

write.csv(df, paste0(baseFolder, "dfAll.csv"))
dfAllDay <- na.omit(df)
dfAllDay <- dfAllDay[!(dfAllDay$Treatment == "W1"),  ]
dfAllDay[(dfAllDay$Treatment == "W2" & dfAllDay$Year == 2019),  "Treatment"] <- "W5"
dfAllDay[(dfAllDay$Treatment == "W3" & dfAllDay$Year == 2019),  "Treatment"] <- "W2"

dfAllDay$Treatment[dfAllDay$Treatment == "W2"] <- "D10"
dfAllDay$Treatment[dfAllDay$Treatment == "W3"] <- "W10"
dfAllDay$Treatment[dfAllDay$Treatment == "W4"] <- "D"
dfAllDay$Treatment <- factor(dfAllDay$Treatment, levels = c("W5", "W10", "D10", "D"))

dfAllDay_NDVI <- dfAllDay[dfAllDay$VI == "NDVI", ]
dfAllDay_NDRE <- dfAllDay[dfAllDay$VI == "NDRE", ]


g <- ggplot(dfAllDay_NDVI, aes(x = Year, y = Accuracy, fill = Treatment)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = Accuracy - SE, ymax = Accuracy + SE), 
                position = position_dodge(0.9), width = .3) + 
  scale_fill_manual(values = c("olivedrab3", "lightskyblue1", "darkgoldenrod1", "chocolate2")) + 
  geom_hline(aes(yintercept = resultG), linetype = "dashed")
png(paste0(baseFolder, "result_NDVI.png"), 
    width = 1440, height = 1440, res = 214)
print(g)
dev.off()

g <- ggplot(dfAllDay_NDRE, aes(x = Year, y = Accuracy, fill = Treatment)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = Accuracy - SE, ymax = Accuracy + SE), 
                position = position_dodge(0.9), width = .3) + 
  scale_fill_manual(values = c("olivedrab3", "lightskyblue1", "darkgoldenrod1", "chocolate2")) + 
  geom_hline(aes(yintercept = resultG), linetype = "dashed")
png(paste0(baseFolder, "result_NDRE.png"), 
    width = 1440, height = 1440, res = 214)
print(g)
dev.off()
