machine <- c("M100", "M100", "P4M")
year <- c("2019", "2020", "2021")
source("scripts/R/1.functionCode.R")

#' # 0. Read packages
options(stringsAsFactors = FALSE)
library(stringr)
library(psych)
library(ggplot2)
library(corrplot)
library(asreml)
library(tagcloud)
library(RAINBOWR)

# make the folder
baseFolder <- "midstream/3.1.RRMlinear/"
if (!dir.exists(baseFolder)) {
  dir.create(baseFolder)
}

# read the genome data
amat0 <- as.matrix(read.csv("data/genome/amat173583SNP.csv", row.names = 1, header = T))
colnames(amat0)[colnames(amat0) == "X5002T"] <- "5002T"
colnames(amat0)[colnames(amat0) == "HOUJAKU_KUWAZU"] <- "Houjaku Kuwazu"
rownames(amat0)[rownames(amat0) == "HOUJAKU_KUWAZU"] <- "Houjaku Kuwazu"

# read the Mean of FreshWeight
cv <- read.csv("midstream/2.0.freshWeightCV/freshWeightCV.csv", 
               header = T, row.names = 1)
colnames(cv) <- "CV"

# set the root folder
rootFolder <- "data/MSdata/0.4.timeSeriesMS/"

treatment <- c("W1", "W2", "W3", "W4")
# choose the VIs name
VIsName <- c("NDVI", "NDRE")

# read the selected varieties
useVariety <- read.csv("midstream/2.0.freshWeightCV/selVariety.csv", 
                       row.names = 1, check.names = F)[, 1]

sowingDayList <- list(as.Date("2019-07-10"), 
                      as.Date("2020-07-08"), 
                      as.Date("2021-07-06"))

# clarify the day when change the irrigation management in W2, W3
changeDayList <- list(as.Date("2019-09-03"), 
                      as.Date("2020-08-30"), 
                      as.Date("2021-08-29"))

dataPeriodList <- list(c(as.Date("2019-08-25"), as.Date("2019-09-04")), 
                       c(as.Date("2020-08-22"), as.Date("2020-09-04")), 
                       c(as.Date("2021-08-22"), as.Date("2021-09-05")))

dataSelList <- list(c(as.Date("2019-08-25"), as.Date("2019-08-31"), 
                      as.Date("2019-09-03"), as.Date("2019-09-04")), 
                    c(as.Date("2020-08-22"), as.Date("2020-08-24"), 
                      as.Date("2020-08-27"), as.Date("2020-08-29"), 
                      as.Date("2020-09-01"), as.Date("2020-09-04")), 
                    c(as.Date("2021-08-22"), as.Date("2021-08-25"), 
                      as.Date("2021-08-27"), as.Date("2021-08-29"), 
                      as.Date("2021-08-30"), as.Date("2021-08-31"), 
                      as.Date("2021-09-05")))

# set the plot number
condition <- c("C", "W5", "W10", "D")


for (eachDataInd in 1:length(machine)) {
  # eachDataInd <- 1
  # make the save folder
  eachYear <- year[eachDataInd]
  eachMachine <- machine[eachDataInd]
  saveFolder <- paste0(baseFolder, eachYear, "/", eachMachine, "/")
  if (!dir.exists(saveFolder)) {
    dir.create(saveFolder, recursive = T)
  }
  
  # make the figure folder
  figFolder <- paste0(saveFolder, "figure/")
  if (!dir.exists(figFolder)) {
    dir.create(figFolder, recursive = T)
  }
  
  
  # data folder
  dataFolder <- paste0(rootFolder, eachYear, "/", eachMachine, "/")
  msNameList <- list.files(dataFolder, pattern = "timeSeries")
  msDataList <- list.files(dataFolder, pattern = "timeSeries", full.names = T)
  
  
  # choose the date period which will be used the below analysis
  dataRange <- dataPeriodList[[eachDataInd]]
  dayDiff <- as.numeric(dataRange[2] - dataRange[1])
  changeInd <- as.numeric(changeDayList[[eachDataInd]] - dataRange[1]) + 1
  selDayInd <- dataSelList[[eachDataInd]] - sowingDayList[[eachDataInd]]
  selDayVec <- paste(selDayInd, collapse = "|")
  
  for (viInd in 1:length(msDataList)) {
    # viInd <- 1
    viName0 <- msNameList[[viInd]]
    viName <- str_sub(viName0, start = 12, end = (nchar(viName0) - 4))
    
    # read the ms data
    msDataAll <- read.csv(msDataList[[viInd]], 
                          header = T, row.names = 1, check.names = F)
    
    # select the date which we use
    colInd <- grep(pattern = selDayVec, 
                   x = colnames(msDataAll))
    
    msData <- msDataAll[, c(1:2, colInd)]
    msData <- msData[msData$variety %in% rownames(amat0), ]
    msData <- msData[msData$variety %in% useVariety, ]
    
    for (treatmentInd in 1:length(treatment)) {
      # treatmentInd <- 3
      treatmentEach <- treatment[treatmentInd]
      conditionEach <- condition[treatmentInd]
      
      msEach <- msData[msData$treatment == treatmentEach, ]
      msWide <- msEach[, 3:ncol(msEach)]
      rownames(msWide) <- msEach$variety
      
      # remove the day of NA
      usedDay <- apply(is.na(msWide), 2, sum) < 100
      msWide <- msWide[, usedDay, drop = F]
      if (ncol(msWide) < 2) {
        next
      }
      
      # scale the ms data 
      # only standardize the mean based on the first day value
      meanFirst <- mean(msWide[, 1])
      msWide <- msWide - meanFirst
      
      # change df wide type to df long type
      msLong <- WideToLong(msWide, 
                           dfRowInd = rownames(msWide), 
                           dfColInd = colnames(msWide))
      
      colnames(msLong) <- c("Value", "Day", "Variety")
      
      # scale the value
      msLong$Value <- scale(msLong$Value)
      
      msLong$Day <- as.numeric(msLong$Day)
      dayMean <- mean(msLong$Day)
      msLong$Day <- msLong$Day - dayMean
      dataNum <- length(unique(msLong$Day))
      msLong$Variety <- as.factor(msLong$Variety)
      # dim(msLong)
      amat <- amat0[msEach$variety, msEach$variety]
      nVariety <- nrow(amat)
      
      # add the daySq (day^2)
      msLongSq <- data.frame(msLong[, 1:2], 
                             DaySq = (msLong[, 2] ^ 2), 
                             Variety = msLong[, 3])
      
      # set the parameter for legendre
      k <- 1
      res <- asreml(fixed = Value ~ as.factor(Day),
                    random = ~leg(Day, k):id(Variety),
                    data = msLongSq, maxit = 200)
      
      # extract the random effect (intercept, legendre term)
      randomValue <- summary(res, coef = T)$coef.random
      u <- randomValue[grep(pattern = "leg", x = rownames(randomValue)), 1]
      # u <- randomValue[, 1]
      # length(u)
      
      # make the rowname
      varietyList0 <- rownames(randomValue)[1:nVariety]
      varietyList <- str_split(varietyList0, pattern = "_")
      varietyName <- sapply(varietyList, function(eachVariety) {
        # eachVariety <- varietyList[[1]]
        return(eachVariety[3])
      })
      
      # invert vector to matrix
      uMat <- matrix(u, byrow = F, ncol = (k + 1))
      rownames(uMat) <- varietyName
      colnames(uMat) <- paste("L", seq(0, k), sep = "_")
      write.csv(uMat, file = paste0(saveFolder, 
                                    "legendre_", treatmentEach, "_", viName, ".csv"))
      
      
      # bind the "Rank" of fresh weight and the parameters
      df <- na.omit(cbind(uMat, CV = cv[rownames(uMat), ]))
      
      pdf(file = paste0(figFolder, 
                        "legendre_", treatmentEach, "_", viName, ".pdf"))
      plot(res)
      pairs.panels(df)
      
      # set the coefficient matrix of legendre
      d <- 0:dayDiff
      phai <- PhaiLegendre(d = d, k = (k + 1))
      
      
      # visualize the change of VI
      matplot(t(df[, 1:(k + 1)] %*% t(phai)), 
              type = "l", lty = "solid")
      abline(v = changeInd, col="gray40", lty = 2)
      dev.off()
      
      # coloring based on the cv(fresh weight)
      colorEachInd <- smoothPalette(-df[, "CV"], pal = "RdBu")
      
      png(paste0(figFolder,  "matplotAll_", treatmentEach, "_", viName, ".png"), 
          width = 1440, height = 1440, res = 214)
      matplot(t(df[, 1:(k + 1)] %*% t(phai)), 
              type = "l", lty = "solid", col = colorEachInd)
      abline(v = changeInd, col="gray40", lty = 2)
      dev.off()
      # select the top & worst 10 data
      selVariety <- c(rownames(df)[order(df[, "CV"], decreasing = T)][1:10], 
                      rownames(df)[order(df[, "CV"], decreasing = F)][1:10])
      uMatSel <- uMat[selVariety, ]
      
      png(paste0(figFolder,  "matplotSel_", treatmentEach, "_", viName, ".png"), 
          width = 1440, height = 1440, res = 214)
      matplot(t(uMatSel[, 1:(k + 1)] %*% t(phai)), 
              type = "l", lty = "solid", lwd = 2, 
              col = rep(c("Tomato", "SkyBlue"), each = 10), 
              main = "Red:high CV, Blue:low CV")
      abline(v = changeInd, col="gray40", lty = 2)
      dev.off()
      
      # visualize the raw data
      ms0 <- apply(msWide[rownames(df), ], 2, function(eachDayData) {
        # eachDayData <- msWide[rownames(df), 1]
        msDayMean <- mean(eachDayData)
        eachDayData0 <- eachDayData - msDayMean
        return(eachDayData0)
      })
      png(paste0(figFolder,  "matplotRaw_", treatmentEach, "_", viName, ".png"), 
          width = 1440, height = 1440, res = 214)
      matplot(t(ms0), 
              type = "l", lty = "solid", col = colorEachInd)
      dev.off()
      
      # calculate the genomic heritability
      h2Mat <- matrix(NA, nrow = 1, ncol = ncol(uMat))
      colnames(h2Mat) <- colnames(uMat)
      rownames(h2Mat) <- "h2"
      
      # set the GRM
      amatEach <- amat[rownames(uMat), rownames(uMat)]
      
      # make amat list for each trait data
      amatZ <- design.Z(pheno.labels = rownames(uMat),
                        geno.names = rownames(amatEach))
      amatList <- list(Z = amatZ,
                       K = amatEach)
      
      # set each ZETA 
      ZETA <- list(amatList = amatList)
      
      for (eachTraitInd in 1:ncol(uMat)) {
        # eachTraitInd <- 1
        eachTrait <- uMat[, eachTraitInd]
        mmfit <- EMM.cpp(y = scale(eachTrait), ZETA = ZETA)
        h2 <- mmfit$Vu / (mmfit$Vu + mmfit$Ve)
        h2Mat[1, eachTraitInd] <- h2
        # print(mmfit$Vu / (mmfit$Vu + mmfit$Ve))
        # plot(fwMean, mmfit$u)
      }
      write.csv(h2Mat, paste0(saveFolder, 
                              "h2_", treatmentEach, "_", viName, ".csv"))
      
    }
  }
}

# read the genotypic correlation
band_name <- c("NDVI", "NDRE")
index <- c("L0", "L1")
treatment <- c("W5", "W10", "D10", "D")
yearIndex <- paste0(rep(year, each = 2), "_", rep(index, 2))
# make array to input the gCor(resCor) between dryWeight and other phenotypes
h2ArrayAll<- array(NA, dim = c(length(yearIndex),
                               length(treatment),
                               length(band_name)))
dimnames(h2ArrayAll) <- list(yearIndex,
                             treatment,
                             band_name)

fileList <- list.files(path = baseFolder, pattern = "h2", recursive = T)
for (eachFile in fileList) {
  # eachFile <- fileList[[1]]
  
  if (!grepl(pattern = "/", eachFile)) {
    break
  }
  
  # get the file data
  yearEach <- str_split(eachFile, pattern = "/")[[1]][1]
  treatmentEach <- str_split(eachFile, pattern = "_")[[1]][2]
  VIeach0 <- str_split(eachFile, pattern = "_")[[1]][3]
  VIeach <- str_sub(VIeach0, start = 1, end = -5)
  
  # read the data
  h2 <- read.csv(paste0(baseFolder, eachFile), header = T, row.names = 1)
  h2 <- round(as.matrix(h2), 2)
  
  # input the genetic correlation value to the array
  rowInd <- grep(yearEach, yearIndex)
  if (yearEach == 2019 & treatmentEach == "W2") {
    h2ArrayAll[rowInd, 1, VIeach] <- h2
  } else if (yearEach == 2019 & treatmentEach == "W3") {
    h2ArrayAll[rowInd, 3, VIeach] <- h2
  } else if (treatmentEach == "W2") {
    h2ArrayAll[rowInd, 3, VIeach] <- h2
  } else if (treatmentEach == "W3") {
    h2ArrayAll[rowInd, 2, VIeach] <- h2
  } else if (treatmentEach == "W4") {
    h2ArrayAll[rowInd, 4, VIeach] <- h2
  }
}

for (eachData in 1:length(band_name)) {
  write.csv(h2ArrayAll[, , eachData], 
            paste0(baseFolder, "h2", band_name[eachData], ".csv"))
}


#### calculate the correlation among each legendre parameter ####

pathList0 <- list.files(path = baseFolder, pattern = "legendre_.._NDVI", recursive = T)
pathList <- pathList0[!grepl(pattern = "figure", x = pathList0)]

dfList <- lapply(pathList, function(eachPath) {
  # eachPath <- pathList[[1]]
  
  # get the file data
  yearEach <- str_split(eachPath, pattern = "/")[[1]][1]
  treatmentEach <- str_split(eachPath, pattern = "_")[[1]][2]
  if (yearEach == 2019 & treatmentEach == "W2") {
    treatmentEach <- "W5"
  } else if (yearEach == 2019 & treatmentEach == "W3") {
    treatmentEach <- "W2"
  }
  VIeach0 <- str_split(eachPath, pattern = "_")[[1]][3]
  VIeach <- str_sub(VIeach0, start = 1, end = -5)
  eachDf <- read.csv(paste0(baseFolder, eachPath), header = T, row.names = 1)
  colnames(eachDf) <- paste0(yearEach, "_", 
                             treatmentEach, "_", 
                             VIeach, "_", 
                             colnames(eachDf))
  return(eachDf)
})

valueDf <- do.call(what = cbind, args = dfList)
dim(valueDf)
treatment <- c("W1", "W2", "W3", "W4")

for (eachTreatment in treatment) {
  # eachTreatment <- treatment[2]
  valueDfEach <- valueDf[, grepl(pattern = eachTreatment, colnames(valueDf))]
  corrplot(cor(valueDfEach))
}
