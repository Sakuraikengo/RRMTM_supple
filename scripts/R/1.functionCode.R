##### 0.0.soilMoisture.R ########
WideToLong <- function(baseDf, dfRowInd, dfColInd) {
  corMat <- as.matrix(baseDf)
  colnames(corMat) <- dfColInd
  value <- c(corMat)
  colInd <- rep(dfColInd, each = nrow(corMat))
  rowInd <- rep(dfRowInd, ncol(corMat))
  longValue <- cbind(value, colInd, rowInd)
  longValue <- data.frame(longValue)
  longValue$value <- as.numeric(longValue$value)
  return(longValue)
}

########## 0.2.spectralValue.R ############
GetMsData <- function(plotCsv, VIsName, machine, refBoardValue = NULL) {
  plotMat0 <- as.matrix(read.csv(plotCsv, header = T, row.names = 1))
  if (machine == "M100" & is.null(refBoardValue)) {
    refBoardValue <- matrix(1, nrow = 1, ncol = 7)
    colnames(refBoardValue) <- c("475550850_475", "475550850_550", "475550850_850", 
                                 "550660850_550", "550660850_660", "550660850_850", 
                                 "725")
  } else if (machine == "P4M" & is.null(refBoardValue)) {
    refBoardValue <- matrix(1, nrow = 1, ncol = 6)
    colnames(refBoardValue) <- c("blue", "green", "nir", 
                                 "red", "red edge")
  }
  
  if (machine == "M100") {
    spectral_blue <- plotMat0["475550850_475", ] / refBoardValue[, "475550850_475"]
    spectral_green <- plotMat0["550660850_550", ] / refBoardValue[, "550660850_550"]
    spectral_red <- plotMat0["550660850_660", ] / refBoardValue[, "550660850_660"]
    spectral_redge <- plotMat0["725", ] / refBoardValue[, "725"]
    spectral_nir <- plotMat0["550660850_850", ] / refBoardValue[, "550660850_850"]
  } else if (machine == "P4M") {
    spectral_blue <- plotMat0["blue", ] / refBoardValue[, "blue"]
    spectral_green <- plotMat0["green", ] / refBoardValue[, "green"]
    spectral_red <- plotMat0["red", ] / refBoardValue[, "red"]
    spectral_redge <- plotMat0["red edge", ] / refBoardValue[, "red edge"]
    spectral_nir <- plotMat0["nir", ] / refBoardValue[, "nir"]
  }
  
  # add VIs
  GRVI <- (spectral_green - spectral_red) / (spectral_green + spectral_red)
  NDVI <- (spectral_nir - spectral_red) / (spectral_nir + spectral_red)
  GNDVI <- (spectral_nir - spectral_green) / (spectral_nir + spectral_green)
  BNDVI <- (spectral_nir - spectral_blue) / (spectral_nir + spectral_blue)
  NDRE <- (spectral_nir - spectral_redge) / (spectral_nir + spectral_redge)
  CIgreen <- (spectral_nir / spectral_green) - 1
  CIRE <- (spectral_nir / spectral_redge) - 1
  MSR <- ((spectral_nir / spectral_red) - 1) / sqrt((spectral_nir / spectral_red) + 1)
  MSRRE <- ((spectral_nir / spectral_redge) - 1) / sqrt((spectral_nir / spectral_redge) + 1)
  NDI <- (spectral_redge - spectral_red) / (spectral_redge + spectral_red)
  BGVI <- (spectral_green - spectral_blue) / (spectral_green + spectral_blue)
  EVI <- 2.5 * (spectral_nir - spectral_red) / (spectral_nir + 6 * spectral_red - 7.5 * spectral_blue + 1)
  VARI <- (spectral_green - spectral_red) / (spectral_green + spectral_red - spectral_blue)
  ARI <- (1 / spectral_green) - (1 / spectral_redge)
  mARI <- ((1 / spectral_green) - (1 / spectral_redge)) * spectral_nir
  
  BRVI <- (spectral_blue - spectral_red) / (spectral_blue + spectral_red)
  BRENDI <- (spectral_redge - spectral_blue) / (spectral_redge + spectral_blue)
  GRENDI <- (spectral_redge - spectral_green) / (spectral_redge + spectral_green)
  
  RGBVI <- (spectral_green^2 - (spectral_blue * spectral_red)) / (spectral_green^2 + (spectral_blue * spectral_red))
  
  RTVI <- 100 * (spectral_nir - spectral_redge) - 10 * (spectral_nir - spectral_green)
  
  matVIs <- rbind(GRVI = GRVI,
                  NDVI = NDVI, 
                  GNDVI = GNDVI, 
                  BNDVI = BNDVI, 
                  NDRE = NDRE, 
                  CIgreen = CIgreen, 
                  CIRE = CIRE, 
                  MSR = MSR, 
                  MSRRE = MSRRE, 
                  NDI = NDI, 
                  BGVI = BGVI, 
                  EVI = EVI, 
                  VARI = VARI, 
                  ARI = ARI, 
                  mARI = mARI, 
                  BRVI = BRVI, 
                  BRENDI = BRENDI, 
                  GRENDI = GRENDI, 
                  RGBVI = RGBVI, 
                  RTVI = RTVI)
  plotMat <- matVIs[VIsName, ]
  # plotMat <- rbind(plotMat0, matVIsSel)
  plotMat[is.infinite(plotMat)] <- 0
  return(plotMat)
}

##### 0.3.dataBind.R ########
FlowerDataAsDate <- function(flowerDate) { 
  yearMonthDay <- str_split(flowerDate$FloweringDate, pattern = "/")
  flowerDate$FloweringDate <- sapply(yearMonthDay, function(eachYearMonthDay) {
    # eachYearMonthDay <- yearMonthDay[[1]]
    eachDate <- paste0(eachYearMonthDay[1], "-", 
                       str_sub(paste0("0", eachYearMonthDay[2]), start = -2), "-", 
                       str_sub(paste0("0", eachYearMonthDay[3]), start = -2))
    return(eachDate)
  })
  flowerDate$FloweringDate <- as.Date(flowerDate$FloweringDate)
  return(flowerDate)
}


##### 3.1.RRMlinear.R #####
PhaiLegendre <- function(d, k) {
  a <- -1 + (2 * (d - min(d))) / (max(d) - min(d))
  M <- 1
  for (i in 1:(k - 1)) {
    M <- cbind(M, a ^ i)
  }
  lambda1 <- sqrt(1 / 2) * c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)
  lambda2 <- sqrt(3 / 2) * c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)
  lambda3 <- sqrt(5 / 2) * c(0, 0, 0, 0, 0, 0, 0, 0, 3, 0, -1) / 2
  lambda4 <- sqrt(7 / 2) * c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1) / 8
  lambda5 <- sqrt(9 / 2) * c(0, 0, 0, 0, 0, 0, 35, 0, -30, 0, 3) / 8
  lambda6 <- sqrt(11 / 2) * c(0, 0, 0, 0, 0, 63, 0, -70, 0, 15, 0) / 8
  lambda7 <- sqrt(13 / 2) * c(0, 0, 0, 0, 231, 0, -315, 0, 105, 0, -5) / 16
  lambda8 <- sqrt(15 / 2) * c(0, 0, 0, 429, 0, -693, 0, 315, 0, -35, 0) / 16
  lambda9 <- sqrt(17 / 2) * c(0, 0, 6435, 0, -12012, 0, 6930, 0, -1260, 0, 35) / 128
  lambda10 <- sqrt(19 / 2) * c(0, 12155, 0, -25740, 0, 18018, 0, -4620, 0, 315, 0) / 128
  lambda11 <- sqrt(21 / 2) * c(46189, 0, -109395, 0, 90090, 0, -30030, 0, 3465, 0, -63) / 256
  Lambda <- cbind(rev(lambda1), rev(lambda2), rev(lambda3), 
                  rev(lambda4), rev(lambda5), rev(lambda6), 
                  rev(lambda7), rev(lambda8), rev(lambda9), 
                  rev(lambda10), rev(lambda11))
  
  Phai <- M %*% Lambda[1:k, 1:k]
  # print(image(cor(t(Phai))))
  return(Phai)
}


##### 4.0.RRMMTMfitting.R ######
MakeCorMat <- function(cov_file, index_names) {
  #cov_file <- list.files(data_folder)[5]
  cov_df <- read.table(cov_file)
  cor_value <- cov2cor(as.matrix(cov_df))
  cor_value <- round(cor_value, digits = 2)
  colnames(cor_value) <- index_names
  rownames(cor_value) <- index_names
  return(cor_value)
}

##### 4.2.RRMMTMyear.R ######
MakeRRM <- function(matWide, day, k) {
  matLong <- WideToLong(matWide, 
                        dfRowInd = rownames(matWide), 
                        dfColInd = day)
  
  colnames(matLong) <- c("Value", "Day", "Variety")
  
  # scale the value
  matLong$Value <- scale(matLong$Value)
  matLong$Day <- as.numeric(matLong$Day)
  dayMean <- mean(matLong$Day)
  matLong$Day <- matLong$Day - dayMean
  matLong$Variety <- as.factor(matLong$Variety)
  nVariety <- nrow(matWide)
  
  # add the daySq (day^2)
  matLongSq <- data.frame(matLong[, 1:2], 
                          DaySq = (matLong[, 2] ^ 2), 
                          Variety = matLong[, 3])
  
  # make the model
  res <- asreml(fixed = Value ~ as.factor(Day),
                random = ~leg(Day, 2):id(Variety),
                data = matLongSq, maxit = 200, trace = F)
  
  # extract the random effect (intercept, legendre term)
  randomValue <- summary(res, coef = T)$coef.random
  u <- randomValue[grep(pattern = "leg", x = rownames(randomValue)), 1]
  
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
  return(uMat)
}

MakeRRMLinear <- function(matWide, day, k) {
  matLong <- WideToLong(matWide, 
                        dfRowInd = rownames(matWide), 
                        dfColInd = day)
  
  colnames(matLong) <- c("Value", "Day", "Variety")
  
  # scale the value
  matLong$Value <- scale(matLong$Value)
  matLong$Day <- as.numeric(matLong$Day)
  dayMean <- mean(matLong$Day)
  matLong$Day <- matLong$Day - dayMean
  matLong$Variety <- as.factor(matLong$Variety)
  nVariety <- nrow(matWide)
  
  # add the daySq (day^2)
  matLongSq <- data.frame(matLong[, 1:2], 
                          DaySq = (matLong[, 2] ^ 2), 
                          Variety = matLong[, 3])
  
  # make the model
  res <- asreml(fixed = Value ~ as.factor(Day),
                random = ~leg(Day, 1):id(Variety),
                data = matLongSq, maxit = 200, trace = F)
  
  # extract the random effect (intercept, legendre term)
  randomValue <- summary(res, coef = T)$coef.random
  u <- randomValue[grep(pattern = "leg", x = rownames(randomValue)), 1]
  
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
  return(uMat)
}


