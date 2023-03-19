machine <- c("M100", "M100", "P4M")
year <- c("2019", "2020", "2021")
source("scripts/R/1.functionCode.R")

#' # 0. Read packages
options(stringsAsFactors = FALSE)
library(stringr)
library(psych)
library(ggplot2)
library(corrplot)
library(tagcloud)
library(RAINBOWR)
library(asreml)

# make the folder
baseFolder <- "midstream/3.0.RRM/"
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

# make the data frame to bind the result of RR model for each parameter
dfRR <- NULL

# make the array to input the genomic heritability of each day data
band_name <- c("NDVI", "NDRE")
condition <- c("W5", "W10", "D10", "D")
selDayList <- list(c("0825", "0831", "0903", "0904"), 
                   c("0822", "0824", "0827", "0829", "0901", "0904"), 
                   c("0822", "0825", "0827", "0829", "0830", "0831", "0905"))
day <- unlist(selDayList)
yearDay <- paste0(rep(year, sapply(selDayList, length)), "_", day)
# make array to input the gCor(resCor) between dryWeight and other phenotypes
h2ArrayAll <- array(NA, dim = c(length(yearDay),
                                length(condition),
                                length(band_name)))
dimnames(h2ArrayAll) <- list(yearDay,
                             condition,
                             band_name)


for (eachDataInd in 1:length(machine)) {
  # eachDataInd <- 1
  # make the save folder
  eachYear <- year[eachDataInd]
  eachMachine <- machine[eachDataInd]
  
  # data folder
  dataFolder <- paste0(rootFolder, eachYear, "/", eachMachine, "/")
  msNameList <- list.files(dataFolder, pattern = "timeSeries")
  msDataList <- list.files(dataFolder, pattern = "timeSeries", full.names = T)
  
  
  # choose the date period which will be used the below analysis
  dataRange <- dataPeriodList[[eachDataInd]]
  dayDiff <- as.numeric(dataRange[2] - dataRange[1])
  changeInd <- as.numeric(changeDayList[[eachDataInd]] - dataRange[1])
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
      
      msEach <- msData[msData$treatment == treatmentEach, ]
      msWide <- msEach[, 3:ncol(msEach)]
      rownames(msWide) <- msEach$variety
      
      # remove the day of NA
      usedDay <- apply(is.na(msWide), 2, sum) < 100
      msWide <- msWide[, usedDay, drop = F]
      if (ncol(msWide) < 2) {
        next
      }
      
      # calculating the genomic heritability
      amat <- amat0[rownames(msWide), rownames(msWide)]
      
      # make amat list for each trait data
      amatZ <- design.Z(pheno.labels = rownames(msWide),
                        geno.names = rownames(amat))
      amatList <- list(Z = amatZ,
                       K = amat)
      
      # set each ZETA 
      ZETA <- list(amatList = amatList)
      
      for (eachDayInd in 1:ncol(msWide)) {
        # eachDayInd <- 1
        eachTrait <- msWide[, eachDayInd]
        mmfit <- EMM.cpp(y = scale(eachTrait), ZETA = ZETA)
        h2 <- mmfit$Vu / (mmfit$Vu + mmfit$Ve)
        h2 <- round(h2, 2)
        
        # input the genetic correlation value to the array
        rowInd <- grep(eachYear, yearDay)[1] + eachDayInd - 1
        if (eachYear == 2019 & treatmentEach == "W2") {
          h2ArrayAll[rowInd, 1, viName] <- h2
        } else if (eachYear == 2019 & treatmentEach == "W3") {
          h2ArrayAll[rowInd, 3, viName] <- h2
        } else if (treatmentEach == "W2") {
          h2ArrayAll[rowInd, 3, viName] <- h2
        } else if (treatmentEach == "W3") {
          h2ArrayAll[rowInd, 2, viName] <- h2
        } else if (treatmentEach == "W4") {
          h2ArrayAll[rowInd, 4, viName] <- h2
        }
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
      # calculating the loglik and AIC for each parameter
      for (k in 0:2) {
        # k <- 0
        res <- asreml(fixed = Value ~ as.factor(Day),
                      random = ~leg(Day, k):id(Variety),
                      data = msLongSq, maxit = 200)
        loglik <- res$loglik
        # parameters of "fixed effects", "random effects"
        nPara <- length(selDayInd) + ((k + 1) ^ 2 - (k + 1)) / 2 + (k + 1)
        AIC <- -2 * loglik + 2 * nPara
        dfRReach <- data.frame(Year = eachYear, Treatment = treatmentEach, 
                               VI = viName, k = k + 1, Loglik = loglik, AIC = AIC, 
                               p = nPara)
        dfRR <- rbind(dfRR, dfRReach)
      }
    }
  }
}

# save the result of RRM usingeach parameter 
dfRR$Treatment[dfRR$Treatment == "W2" & dfRR$Year == "2019"] <- "W5"
dfRR$Treatment[dfRR$Treatment == "W3" & dfRR$Year == "2019"] <- "W2"
dfRR$Treatment[dfRR$Treatment == "W1"] <- "C"
dfRR$Treatment[dfRR$Treatment == "W2"] <- "D10"
dfRR$Treatment[dfRR$Treatment == "W3"] <- "W10"
dfRR$Treatment[dfRR$Treatment == "W4"] <- "D"

write.csv(dfRR, file = paste0(baseFolder, "resultRRM.csv"))

# save the result of genomic heritability
# save the result
for (i in 1:dim(h2ArrayAll)[3]) {
  VI <- dimnames(h2ArrayAll)[[3]][i]
  write.csv(h2ArrayAll[, , i], file = paste0(baseFolder, "h2", VI, ".csv"))
}
