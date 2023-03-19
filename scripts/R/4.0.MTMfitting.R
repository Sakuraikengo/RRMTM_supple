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
baseFolder <- "midstream/4.0.MTMfitting/"
# baseFolder <- "C:/Users/biometrics/Desktop/4.0.MTMfitting/"
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
cl <- makeCluster(3)
registerDoParallel(cl)
resultAll <- foreach(eachDataInd = 1:length(year), .packages = c("stringr", "MTM")) %dopar% {
  # eachDataInd <- 3
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
  
  # prepare the list to input the correlation among parameters and CV
  corMatList <- NULL
  
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
      
      df <- data.frame(selData[rownames(cv), 2:ncol(selData)], 
                       CV = cv)
      rownames(df) <- rownames(cv)
      n_trait <- ncol(df)
      
      # append the correlation to the list
      corMatList[[paste(eachTreatment, 
                        eachVi, 
                        eachYear, 
                        sep = "_")]] <- round(cor(df), 2)
      
      # When there is the result already, skip the below calculation
      fileBaseName <- paste0(saveFolder, 
                             eachTreatment, "_", 
                             eachVi)
      
      if (file.exists(paste0(fileBaseName, "_MTM_G.txt"))) {
        next
      }
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
      
      # estimation (for RMarkdown)
      MTM_output <- capture.output(
        res_MTM <- MTM_2(Y = MTM_Y, K = MTM_K, resCov = MTM_resCov,
                         nIter = 12000, burnIn = 2000, thin = 20,
                         saveAt = paste0(fileBaseName, "_"))
      )
      
      write.table(res_MTM$YHat, paste0(fileBaseName, "_MTM_fitting.txt"), quote = F)
      write.table(res_MTM$resCov$R, paste0(fileBaseName, "_MTM_resCov.txt"), quote = F)
      write.table(res_MTM$K[[1]]$G, paste0(fileBaseName, "_MTM_G.txt"), quote = F)
      write.table(res_MTM$K[[1]]$U, paste0(fileBaseName, "_MTM_U.txt"), quote = F)
      
      # diagnosis (log-likelihood)
      MTM_logLik <- read.table(paste0(fileBaseName, "_logLik.dat"))
      png(paste0(fileBaseName, "_logLik.png"))
      plot(MTM_logLik[[1]], cex = 0.5,
           main = paste('Log-likelihood', eachYear, eachTreatment, eachVi, collapse = " "),
           ylab = 'Log-likelihood')
      dev.off()
      
      # diagnosis (mu)
      MTM_mu <- readLines(paste0(fileBaseName, "_mu.dat"))
      MTM_mu_1 <- sapply(MTM_mu, strsplit, split = ' ')
      names(MTM_mu_1) <- NULL
      MTM_mu_2 <- sapply(MTM_mu_1, function(x) {
        length_max <- max(sapply(MTM_mu_1, length))
        y <- rep(NA, length_max)
        y[1:length(x)] <- x
        as.numeric(y)
      })
      MTM_mu_3 <- t(as.data.frame(MTM_mu_2))
      png(paste0(fileBaseName, "_mu.png"))
      matplot(MTM_mu_3, cex = 0.5, pch = 1:ncol(MTM_mu_3),
              main = paste('mu', eachYear, eachTreatment, eachVi, collapse = " "),
              ylab = 'mu')
      dev.off()
    }
  }
  return(corMatList)
}
stopCluster(cl)


# read the genotypic correlation
band_name <- c("NDVI", "NDRE")
treatment <- c("W1", "W2", "W3", "W4")

for (eachDataInd in 1:length(year)) {
  # eachDataInd <- 2
  # set working directory
  eachYear <- year[eachDataInd]
  eachMachine <- machine[eachDataInd]
  saveFolder <- paste0(baseFolder, eachYear, "/", eachMachine, "/")
  
  # make the folder to save the figures
  figFolder <- paste0(saveFolder, "figures/")
  if (!dir.exists(figFolder)) {
    dir.create(figFolder)
  }
  
  # read the cov files
  gCovFiles <- list.files(saveFolder, pattern = "MTM_G")
  resCovFiles <- list.files(saveFolder, pattern = "MTM_resCov")
  
  for (fileIndex in 1:length(gCovFiles)) {
    # fileIndex <- 4
    # get file info
    gCovRead <- gCovFiles[fileIndex]
    treatmentEach <- str_split(gCovRead, pattern = "_")[[1]][1]
    VIeach <- str_split(gCovRead, pattern = "_")[[1]][2]
    
    # extract the phenotypic correaltion
    dataInd <- paste(treatmentEach, VIeach, eachYear, sep = "_")
    corValuePheno <- resultAll[[eachDataInd]][[dataInd]]
    if (is.null(corValuePheno)) {
      next
    }
    dataDay <- str_sub(colnames(corValuePheno)[1:(ncol(corValuePheno) - 1)], 
                       start = 6, end = 9)
    
    # make cor plot "g"
    gCovPath <- paste0(saveFolder, "/", gCovRead)
    corValueG <- MakeCorMat(gCovPath, c(dataDay, "CV"))
    # corrplot(corValueG, method = "shade", tl.srt = 45, addCoef.col = "black",
    #          tl.col = "black", number.cex = .7)
    # mtext(paste0(the_day, " G cor ")
    #       , at = 3, line = 2.5, cex = 1)
    
    # get file info
    resCovRead <- resCovFiles[fileIndex]
    # treatment_res <- GetFileInfo(resCovRead)[[1]]
    # data_type_res <- GetFileInfo(resCovRead)[[2]]
    
    # make cor plot "g"
    resCovPath <- paste0(saveFolder, "/", resCovRead)
    CorValueRes <- MakeCorMat(resCovPath, c(dataDay, "CV"))
    
    # make the matrix 
    # (upper.tri = genotypic correlation, lower.tri = phenotypic correlation)
    corMat <- matrix(NA, nrow = nrow(corValuePheno), ncol = ncol(corValuePheno))
    rownames(corMat) <- colnames(corMat) <- rownames(corValueG)
    diag(corMat) <- 1
    corMat[upper.tri(corMat)] <- corValueG[upper.tri(corMat)]
    corMat[lower.tri(corMat)] <- corValuePheno[lower.tri(corMat)]
    
    corMatSort0 <- corMat[c(2:ncol(corMat), 1), ]
    corMatSort <- corMatSort0[, c(2:ncol(corMatSort0), 1)]
    # save the figures
    png(paste0(figFolder, dataInd, ".png"), 
        width = 1440, height = 1440, res = 214)
    corrplot(corMat, method = "shade", addCoef.col = "black",
             tl.col = "black", tl.pos = "d")
    dev.off()
    
    png(paste0(figFolder, dataInd, "_residual.png"), 
        width = 1440, height = 1440, res = 214)
    corrplot(CorValueRes, method = "shade", addCoef.col = "black",
             tl.col = "black", tl.pos = "d", type = "upper")
    dev.off()
  } 
}

# read the genotypic correlation
band_name <- c("NDVI", "NDRE")
treatment <- c("W5", "W10", "D10", "D")
day <- unlist(selDayList)
yearDay <- paste0(rep(year, sapply(selDayList, length)), "_", day)
# make array to input the gCor(resCor) between dryWeight and other phenotypes
gCorArrayAll<- array(NA, dim = c(length(yearDay),
                                 length(treatment),
                                 length(band_name)))
dimnames(gCorArrayAll) <- list(yearDay,
                               treatment,
                               band_name)

for (eachDataInd in 1:length(year)) {
  # eachDataInd <- 1
  # set working directory
  eachYear <- year[eachDataInd]
  eachMachine <- machine[eachDataInd]
  eachDayVec <- selDayList[[eachDataInd]]
  eachDayLen <- c(0, sapply(selDayList, length))
  
  saveFolder <- paste0(baseFolder, eachYear, "/", eachMachine, "/")
  
  # read the cov files
  gCovFiles <- list.files(saveFolder, pattern = "MTM_G")
  
  for (fileIndex in 1:length(gCovFiles)) {
    # fileIndex <- 6
    
    # get file info
    gCovRead <- gCovFiles[fileIndex]
    treatmentEach <- str_split(gCovRead, pattern = "_")[[1]][1]
    VIeach <- str_split(gCovRead, pattern = "_")[[1]][2]
    
    # make cor plot "g"
    gCovPath <- paste0(saveFolder, "/", gCovRead)
    corValueG <- MakeCorMat(gCovPath, c(eachDayVec, "CV"))
    
    # input the genetic correlation value to the array
    startInd <- sum(eachDayLen[1:eachDataInd]) + 1
    rowInd <- startInd:(startInd + eachDayLen[eachDataInd + 1] - 1)
    if (eachDataInd == 1 & treatmentEach == "W2") {
      gCorArrayAll[rowInd, 1, VIeach] <- corValueG["CV", 1:(ncol(corValueG) - 1)]
    } else if (eachDataInd == 1 & treatmentEach == "W3") {
      gCorArrayAll[rowInd, 3, VIeach] <- corValueG["CV", 1:(ncol(corValueG) - 1)]
    } else if (treatmentEach == "W2") {
      gCorArrayAll[rowInd, 3, VIeach] <- corValueG["CV", 1:(ncol(corValueG) - 1)]
    } else if (treatmentEach == "W3") {
      gCorArrayAll[rowInd, 2, VIeach] <- corValueG["CV", 1:(ncol(corValueG) - 1)]
    } else if (treatmentEach == "W4") {
      gCorArrayAll[rowInd, 4, VIeach] <- corValueG["CV", 1:(ncol(corValueG) - 1)]
    }
  } 
}

# save the result
for(i in 1:dim(gCorArrayAll)[3]) {
  VI <- dimnames(gCorArrayAll)[[3]][i]
  write.csv(gCorArrayAll[, , i], file = paste0(baseFolder, "gCor", VI, ".csv"))
}
