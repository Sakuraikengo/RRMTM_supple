machine <- c("M100", "M100", "P4M")
year <- c("2019", "2020", "2021")
source("scripts/R/1.functionCode.R")
options(stringsAsFactors = FALSE)

library(readr)
library(stringr)
library(psych)
library(ggplot2)
library(cluster)
# library(gaston)

# make the folder
baseFolder <- "data/MSdata/0.3.dataBind/"
if (!dir.exists(baseFolder)) {
  dir.create(baseFolder)
}

# set the root folder
rootFolder <- "data/MSdata/0.2.spectralValue/"
rawFolder <- "raw_data/"

# fill NA irregular plot
irrPlotList <- list(c("W1-036", "W1-085", "W1-115",
                      "W2-005", "W2-038", "W2-099", "W2-149",
                      "W3-014", "W3-036"),
                    c("W1-154"),
                    c("W2-068", "W2-171", "W3-046", "W3-094", "W4-056"))

sowingDayList <- list(as.Date("2019-07-10"), 
                      as.Date("2020-07-08"), 
                      as.Date("2021-07-06"))

# select the day based on the weather (only sunny & 11:00-13:00)
selDayList <- list(c("0802", "0810", "0817", "0818", 
                     "0824", "0825", "0831", "0903", "0904"), 
                   c("0801", "0803", "0809", "0811", "0813", 
                     "0815", "0818", "0820", "0822", "0824", 
                     "0827", "0829", "0901", "0904"), 
                   c("0724", "0725", "0726", "0730", "0731", 
                     "0802", "0806", "0815", "0822", "0825", 
                     "0827", "0829", "0830", "0831", "0905"))

condition <- c("W1", "W2", "W3", "W4")
treatment <- rep(condition, each = 200)
plotNumberAll <- paste0(treatment, "-", 
                        rep(formatC(1:200, width = 3, flag = "0"), length(condition)))
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
  dataFolder <- paste0(rootFolder, eachYear, "/", eachMachine, "/")
  dataDay <- selDayList[[eachDataInd]]
  selDayInd <- paste(dataDay, collapse = "|")
  # dataDay0 <- list.files(dataFolder)
  # dataDay <- unique(str_sub(dataDay0, start = 1, end = 4))
  
  # read the dry weight
  phenoFolder <- paste0(rawFolder, eachYear, "/field_data/")
  phenoFile <- list.files(phenoFolder, pattern = "ShootPhenotype", full.names = T)
  plantPhenotype0 <- read.csv(phenoFile, header = TRUE)
  plantPhenotype <- plantPhenotype0[!(plantPhenotype0$block %in% c("RC", "RD")), ]
  
  # calculate the each plot mean
  plotNumber <- paste0(plantPhenotype$block, "-", formatC(plantPhenotype$plot, width = 3, flag = "0"))
  plantPhenotype <- cbind(plantPhenotype, plotNumber)
  
  # dryLeaves <- aggregate(freshWeight_Leaves_g~plotNumber, data = plantPhenotype, FUN = mean)
  freshShoot <- aggregate(FreshWeight_Shoot_g~plotNumber, data = plantPhenotype, FUN = mean)
  
  # dryLeavesInd <- match(x = rownames(plantAreaAll), table = dryLeaves$plotNumber)
  freshShootInd <- match(x = plotNumberAll, table = freshShoot$plotNumber)
  freshWeight0 <- freshShoot[freshShootInd, "FreshWeight_Shoot_g"]
  
  freshWeight <- matrix(NA, nrow = length(plotNumberAll), ncol = 1)
  rownames(freshWeight) <- plotNumberAll
  freshWeight[, 1] <- freshWeight0
  
  # read the plotID and variety name
  plotIDfile <- list.files(phenoFolder, pattern = "PlotIDtable", full.names = T)
  plotIDraw <- read.csv(plotIDfile, header = T)
  plotID <- c(plotIDraw$Name[order(plotIDraw$plotW1)], plotIDraw$Name[order(plotIDraw$plotW2)], 
              plotIDraw$Name[order(plotIDraw$plotW3)], plotIDraw$Name[order(plotIDraw$plotW4)])
  
  medianCsvList0 <- list.files(dataFolder,
                               pattern = "median", full.names = T)
  medianCsvList <- medianCsvList0[grep(selDayInd, medianCsvList0)]
  eachDayDfList <- lapply(medianCsvList, function(eachdayCsv) {
    # eachdayCsv <- medianCsvList[[2]]
    
    # read the plant area each day
    theDay0 <- str_split(eachdayCsv, pattern = "/")[[1]][6]
    theDay <- str_sub(theDay0, 1, 4)
    
    # read the each day MS data
    eachDayDf0 <- read.csv(eachdayCsv, header = T, row.names = 1)
    
    # remove the 300s plot
    plotInd <- as.numeric(str_sub(rownames(eachDayDf0), start = 4, end = 6))
    eachDayDf <- eachDayDf0[plotInd <= 200, ]
    eachDayDf[abs(eachDayDf) == 1] <- NA
    
    eachDayPhenoDf <- cbind(variety = plotID, 
                            treatment = treatment, 
                            freshWeight = freshWeight[, 1], 
                            eachDayDf)
    eachDayPhenoDf[rownames(eachDayPhenoDf) %in% irrPlotList[[eachDataInd]], ] <- NA
    
    pcaData <- na.omit(eachDayPhenoDf[, 4:ncol(eachDayPhenoDf)])
    # pcaData <- scale(pcaData)
    
    if (nrow(pcaData) > 100) {
      # to decide the number of clusters
      gap <- clusGap(x = t(pcaData), FUNcluster = kmeans, K.max = 10, B = 10, verbose = F)
      pdf(paste0(saveFolder, theDay, "_PCA.pdf"))
      plot(gap, main = "clusGap")
      
      hc <- hclust(dist(t(scale(pcaData))))
      plot(hc)
      
      # do the PCA
      res <- prcomp(pcaData, scale = T)
      # summary(res)
      plot(res)
      
      op <- par(mfrow = c(2,2))
      # make biplot
      biplot(res, choices = 1:2)
      biplot(res, choices = 3:4)
      
      factor.loadings <- cor(pcaData, res$x[,1:4])
      theta <- 2 * pi * (0:100 / 100)
      x <- cos(theta)
      y <- sin(theta)
      plot(factor.loadings[,1:2], xlim = c(-1,1), ylim = c(-1,1), pch = " ")
      text(factor.loadings[,1:2], rownames(factor.loadings), col = "red")
      lines(x, y, col = "gray")
      abline(v = 0, h = 0)
      plot(factor.loadings[,3:4], xlim = c(-1,1), ylim = c(-1,1), pch = " ")
      text(factor.loadings[,3:4], rownames(factor.loadings), col = "red")
      lines(x, y, col = "gray")
      abline(v = 0, h = 0)
      par(op)
      
      dev.off()
    }
    return(eachDayPhenoDf)
  })
  
  # read the flowering date
  flowerDay0 <- read.csv(list.files(phenoFolder, pattern = "Flowering", full.names = T), 
                         header = T)
  flowerDay0 <- flowerDay0[flowerDay0$plot <= 200, ]
  flowerDay <- flowerDay0[, c("variety", "block", "FloweringDate")]
  flowerDate <- FlowerDataAsDate(flowerDate = flowerDay)
  rownames(flowerDate) <- plotNumberAll
  dayInd <- paste(eachYear, str_sub(dataDay, 1, 2), str_sub(dataDay, 3, 4), sep = "-")
  dayInd <- as.Date(dayInd)
  
  flowerDateList <- lapply(dayInd, function(eachDayInd) {
    # eachDayInd <- dayInd[1]
    flowerOrNot <- flowerDate[, 3] - eachDayInd
    # 0 means flowering already, 1 means not flowering, NA is also not flowering
    flowerOrNot[flowerOrNot > 0] <- 1
    flowerOrNot[flowerOrNot <= 0] <- 0
    flowerOrNot[is.na(flowerOrNot)] <- 1
    return(flowerOrNot)
  })
  
  ########## make matplot #########
  # we want to set the "variety" as the first column
  refDf <- cbind(eachDayDfList[[1]][, 1:2], 
                 flower = flowerDateList[[1]], 
                 eachDayDfList[[1]][, 3:ncol(eachDayDfList[[1]])])
  allDayArray <- array(NA, dim = c(nrow(refDf), (ncol(refDf) - 2), length(dataDay)))
  dimnames(allDayArray) <- list(rownames(refDf), colnames(refDf)[3:ncol(refDf)], dataDay)
  # data day split to the each group
  for (i in 1:length(dataDay)) {
    # i <- 1
    dayDf <- cbind(eachDayDfList[[i]][, 1:2], 
                   flower = as.numeric(flowerDateList[[i]]), 
                   eachDayDfList[[i]][, 3:ncol(eachDayDfList[[i]])])
    allDayArray[, , i] <- as.matrix(dayDf[, 3:ncol(dayDf)])
    write.csv(dayDf, 
              file = paste0(saveFolder, dataDay[i], "_medianPheno.csv"))
  }
  
  
  ########## make matplot #########
  
  condition <- c("W1", "W2", "W3", "W4")
  pdf(paste0(saveFolder, "matplot.pdf"))
  for (eachBandInd in 1:dim(allDayArray)[2]) {
    # eachBandInd <- 1
    eachBandValue <- allDayArray[, eachBandInd, ]
    dayNA <- apply(is.na(eachBandValue), 2, sum) > 500
    eachBandValue <- eachBandValue[, !dayNA]
    # eachBandValue <- na.omit(eachBandValue)
    # if (eachBandInd > 1) {
    #   eachBandValue <- scale(eachBandValue)
    # }
    eachBandName <- colnames(allDayArray)[eachBandInd]
    ylim <- range(eachBandValue, na.rm = T)
    matplot(t(eachBandValue), type = "l", main = paste0(eachBandName), 
            ylab = "value", xlab = "day", xaxt = "n", ylim = ylim)
    axis(1, at = c(1:ncol(eachBandValue)), labels = colnames(eachBandValue))
    
    eachCondition <- str_sub(rownames(eachBandValue), 1, 2)
    for (eachConditionInd in condition) {
      # eachConditionInd <- condition[2]
      eachConditionValue <- eachBandValue[eachCondition %in% eachConditionInd, ]
      matplot(t(eachConditionValue), type = "l", main = paste0(eachConditionInd, " ", eachBandName), 
              ylab = "value", xlab = "day", xaxt = "n", ylim = ylim)
      axis(1, at = c(1:ncol(eachBandValue)), labels = colnames(eachBandValue))
    }
  }
  dev.off()
  
  ########### visualize flowering time###########
  sowingDay <- sowingDayList[[eachDataInd]]
  measuringDay <- dayInd
  
  DTF <- flowerDate$FloweringDate - sowingDay
  flowerDf <- cbind(flowerDate, DTF)
  flowerDf <- na.omit(flowerDf)
  g <- ggplot(flowerDf, aes(x = FloweringDate)) + 
    facet_wrap(~block) + 
    geom_histogram(binwidth = 1) + 
    geom_vline(xintercept = measuringDay, linetype = "dashed") +
    scale_x_date(date_breaks = "1 months") + 
    xlab("Flowering Day") + 
    ylab("Count") + 
    theme(plot.title = element_text(size = 16, face = "bold"), 
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12,face = "bold"), 
          legend.text = element_text(size = 16))
  
  png(paste0(saveFolder, "floweringDay.png"), 
      height = 1440, width = 1440, res = 216)
  print(g)
  dev.off()
}
