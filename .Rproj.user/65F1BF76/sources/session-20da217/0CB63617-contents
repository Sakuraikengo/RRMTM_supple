machine <- c("M100", "M100", "P4M")
year <- c("2019", "2020", "2021")
source("scripts/R/1.functionCode.R")

options(stringsAsFactors = FALSE)
library(readr)
library(stringr)
library(doParallel)

# make the folder
baseFolder <- "data/MSdata/0.2.spectralValue/"
if (!dir.exists(baseFolder)) {
  dir.create(baseFolder)
}

# set the root folder
rootFolder <- "raw_data/"

# set the condition
condition <- c("W1", "W2", "W3", "W4")

# choose the VIs name
VIsName <- c("GRVI", 
             "NDVI", 
             "GNDVI", 
             "BNDVI", 
             "NDRE", 
             "CIgreen", 
             "CIRE", 
             "MSR", 
             "NDI", 
             "ARI", 
             "mARI", 
             "RGBVI", 
             "RTVI")

for (eachDataInd in 1:length(machine)) {
  # eachDataInd <- 2
  
  # make the save folder
  eachYear <- year[eachDataInd]
  eachMachine <- machine[eachDataInd]
  saveFolder <- paste0(baseFolder, eachYear, "/", eachMachine, "/")
  if (!dir.exists(saveFolder)) {
    dir.create(saveFolder, recursive = T)
  }
  
  # root folder
  dataFolder <- paste0(rootFolder, eachYear, "/", eachMachine, "/MSdata/")
  dataDay <- list.files(dataFolder)
  
  cl <- makeCluster(10)
  registerDoParallel(cl)
  allResult <- foreach(dayInd = 1:length(dataDay), .packages = c("doParallel", "stringr")) %dopar% {
    # dayInd <- 3
    the_day <- dataDay[dayInd]
    print(the_day)
    spectralEachDayFolder <- paste0(dataFolder, the_day)
    
    msEachCondition <- foreach(conditionEach = condition) %do% {
      # conditionEach <- condition[2]
      spectralCsvList0 <- list.files(spectralEachDayFolder, 
                                     full.names = T, pattern = conditionEach)
      
      # read the gray scale file
      refBoardList <- list.files(spectralEachDayFolder, 
                                 pattern = "ref", full.names = T)
      refBoardFile <- refBoardList[grep(conditionEach, refBoardList)]
      refBoardValue <- read.csv(refBoardFile, header = T, row.names = 1)
      if (eachMachine == "M100") {
        colnames(refBoardValue) <- c("475550850_475", "475550850_550", "475550850_850", 
                                     "550660850_550", "550660850_660", "550660850_850", 
                                     "725")
      } else if (eachMachine == "P4M") {
        colnames(refBoardValue) <- c("blue", "green", "nir", 
                                     "red", "red edge")
      }
      if (length(grep("ref", spectralCsvList0)) > 0) {
        spectralCsvList0 <- spectralCsvList0[-grep("ref", spectralCsvList0)]
      }
      # get VIs value
      allMsValue <- lapply(spectralCsvList0, function(eachPlot) {
        # eachPlot <- spectralCsvList0[[1]]
        eachPlotValue <- GetMsData(eachPlot, VIsName, 
                                   machine = eachMachine, 
                                   refBoardValue = refBoardValue)
        if ((ncol(eachPlotValue) < 10) || (is.null(dim(eachPlotValue)))) { 
          eachMsData <- rep(NA, length(VIsName))
          names(eachMsData) <- rownames(eachPlotValue)
        } else {
          eachMsData <- apply(eachPlotValue, 1, median)
        }
        eachPlotMedian <- eachMsData
        
        return(eachPlotMedian)
      })
      
      msMedian <- do.call(what = rbind, args = allMsValue)

      nameLen <- nchar(spectralCsvList0[1])
      plotNumeber <- str_sub(spectralCsvList0, start = (nameLen - 9), 
                             end = (nameLen - 4))
      rownames(msMedian) <- plotNumeber
      print(conditionEach)
      return(msMedian)
    }
    
    msMedianEachDay <- do.call(what = rbind, args = msEachCondition)
    
    eachDay <- str_sub(the_day, start = 5, end = 8)
    write.csv(msMedianEachDay, paste0(saveFolder, eachDay, "_median.csv"))
  }
  stopCluster(cl)
}

