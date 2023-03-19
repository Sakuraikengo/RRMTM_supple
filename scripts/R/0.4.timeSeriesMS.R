machine <- c("M100", "M100", "P4M")
year <- c("2019", "2020", "2021")
source("scripts/R/1.functionCode.R")
options(stringsAsFactors = FALSE)

#' # 0. Read packages
library(BGLR)
library(doParallel)
library(tagcloud)
library(lme4)
library(RAINBOWR)
library(stringr)
library(psych)
library(ggplot2)
library(corrplot)

# make the folder
baseFolder <- "data/MSdata/0.4.timeSeriesMS/"
if (!dir.exists(baseFolder)) {
  dir.create(baseFolder)
}

# read the genome data
amat0 <- as.matrix(read.csv("data/genome/amat173583SNP.csv", row.names = 1, header = T))
colnames(amat0)[colnames(amat0) == "X5002T"] <- "5002T"
colnames(amat0)[colnames(amat0) == "HOUJAKU_KUWAZU"] <- "Houjaku Kuwazu"
rownames(amat0)[rownames(amat0) == "HOUJAKU_KUWAZU"] <- "Houjaku Kuwazu"


# set the root folder
rootFolder <- "data/MSdata/0.3.dataBind/"

treatment <- c("W1", "W2", "W3", "W4")
# choose the VIs name
VIsName <- c("NDVI", "NDRE")

sowingDayList <- list(as.Date("2019-07-10"), 
                      as.Date("2020-07-08"), 
                      as.Date("2021-07-06"))

# set the plot number
condition <- c("W1", "W2", "W3", "W4")
treatment <- rep(condition, each = 200)
plotNumberAll <- paste0(treatment, "-", 
                        rep(formatC(1:200, width = 3, flag = "0"), length(condition)))

for (eachDataInd in 1:length(machine)) {
  # make the save folder
  eachYear <- year[eachDataInd]
  eachMachine <- machine[eachDataInd]
  saveFolder <- paste0(baseFolder, eachYear, "/", eachMachine, "/")
  if (!dir.exists(saveFolder)) {
    dir.create(saveFolder, recursive = T)
  }
  
  # data folder
  dataFolder <- paste0(rootFolder, eachYear, "/", eachMachine, "/")
  dataDay0 <- list.files(dataFolder, pattern = ".csv")
  dataDay <- unique(str_sub(dataDay0, start = 1, end = 4))
  
  # change the data day to date
  dataDate <- as.Date(paste(eachYear, 
                            str_sub(dataDay, 1, 2), 
                            str_sub(dataDay, 3, 4), sep = "-"))
  
  # calculate the days after sowing (DTF)
  DTF <- dataDate - sowingDayList[[eachDataInd]]
  dateSeq <- seq(dataDate[1], dataDate[length(dataDate)], by = "days")
  dateSeq <- as.character(dateSeq)
  
  # make the array to input the MS data
  msArray <- array(data = NA, dim = c(length(plotNumberAll),
                                      length(DTF), 
                                      length(VIsName)))
  dimnames(msArray) <- list(plotNumberAll, DTF, VIsName)
  
  for (dayInd in 1:length(dataDay)) {
    # dayInd <- 9
    theDay <- dataDay[dayInd]
    dataAll <- read.csv(paste0(dataFolder, theDay, "_medianPheno.csv"), 
                        header = T, row.names = 1)
    dataMs <- dataAll[, VIsName]
    msArray[, dayInd, ] <- as.matrix(dataMs)
  }
  
  
  # write the csv
  for (msInd in 1:length(VIsName)) {
    # msInd <- 1
    msMat0 <- msArray[, , msInd]
    msMat <- cbind(dataAll[, 1:2], msMat0)
    write.csv(msMat, paste0(saveFolder, "timeSeries_", VIsName[msInd], ".csv"))
  }
}
