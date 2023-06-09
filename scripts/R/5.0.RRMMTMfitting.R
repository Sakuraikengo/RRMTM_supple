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
library(doParallel)
library(MTM)

# make the folder
baseFolder <- "midstream/5.0.RRMMTMfitting/"
# baseFolder <- "C:/Users/biometrics/Desktop/5.0.RRMMTMfitting/"
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
rootFolder <- "midstream/3.1.RRMlinear/"

treatment <- c("W1", "W2", "W3", "W4")
resultIndex <- c("Correlation", "R2", "RMSE")

# choose the VIs name
VIsName <- c("NDVI", "NDRE")

# set the plot number
condition <- c("C", "W5", "W10", "D")

#### calculate the MTM ######
cl <- makeCluster(3)
registerDoParallel(cl)
resultAll <- foreach(eachDataInd = 1:length(year), .packages = c("doParallel", "MTM")) %dopar% {
  # eachDataInd <- 2
  # make the save folder
  eachYear <- year[eachDataInd]
  eachMachine <- machine[eachDataInd]
  saveFolder <- paste0(baseFolder, eachYear, "/", eachMachine, "/")
  if (!dir.exists(saveFolder)) {
    dir.create(saveFolder, recursive = T)
  }
  
  # k <- 1
  # data folder
  dataFolder <- paste0(rootFolder, eachYear, "/", eachMachine, "/")
  msFileList0 <- list.files(dataFolder, pattern = "^legendre", full.names = T)
  msFileList <- msFileList0[grep(pattern = ".csv", x = msFileList0)]
  
  # prepare the list to input the correlation among parameters and CV
  corMatList <- NULL
  
  # select the VI
  for (ViInd in 1:length(VIsName)) {
    # ViInd <- 2
    eachVi <- VIsName[ViInd]
    selDataList <- msFileList[grep(pattern = eachVi, x = msFileList)]
    
    # make the matrix to input the result
    resMat <- matrix(NA, nrow = length(treatment), ncol = length(resultIndex))
    rownames(resMat) <- treatment
    colnames(resMat) <- resultIndex
    
    for (treatmentInd in 1:length(treatment)) {
      # treatmentInd <- 2
      eachTreatment <- treatment[treatmentInd]
      eachCondition <- condition[treatmentInd]
      selData <- selDataList[grep(pattern = eachTreatment, x = selDataList)]
      if (length(selData) < 1) {
        next
      }
      dfMs <- read.csv(selData, header = T, row.names = 1)
      
      df <- data.frame(dfMs[rownames(cv), ], CV = cv)
      rownames(df) <- rownames(cv)
      n_trait <- ncol(df)
      
      # append the correlation to the list
      corMatList[[paste(eachTreatment, eachVi, eachYear, sep = "_")]] <- round(cor(df), 2)
      
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
        res_MTM <- MTM(Y = MTM_Y, K = MTM_K, resCov = MTM_resCov,
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

###### Visualize the genetic and phenotypic correlation ##########

# read the genotypic correlation
band_name <- c("NDVI", "NDRE")
index <- c("L0", "L1")
treatment <- c("W5", "W10", "D10", "D")
yearIndex <- paste0(rep(year, each = 2), "_", rep(index, 2))
# make array to input the gCor(resCor) between dryWeight and other phenotypes
gCorArrayAll<- array(NA, dim = c(length(yearIndex),
                                 length(treatment),
                                 length(band_name)))
dimnames(gCorArrayAll) <- list(yearIndex,
                               treatment,
                               band_name)

for (eachDataInd in 1:length(year)) {
  # eachDataInd <- 1
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
    # fileIndex <- 1
    
    # get file info
    gCovRead <- gCovFiles[fileIndex]
    treatmentEach <- str_split(gCovRead, pattern = "_")[[1]][1]
    VIeach <- str_split(gCovRead, pattern = "_")[[1]][2]
    
    # make cor plot "g"
    gCovPath <- paste0(saveFolder, "/", gCovRead)
    corValueG <- MakeCorMat(gCovPath, c(index, "CV"))
    
    # input the genetic correlation value to the array
    rowInd <- (2*(eachDataInd - 1) + 1):(2*eachDataInd)
    if (eachDataInd == 1 & treatmentEach == "W2") {
      gCorArrayAll[rowInd, 1, VIeach] <- corValueG["CV", 1:2]
    } else if (eachDataInd == 1 & treatmentEach == "W3") {
      gCorArrayAll[rowInd, 3, VIeach] <- corValueG["CV", 1:2]
    } else if (treatmentEach == "W2") {
      gCorArrayAll[rowInd, 3, VIeach] <- corValueG["CV", 1:2]
    } else if (treatmentEach == "W3") {
      gCorArrayAll[rowInd, 2, VIeach] <- corValueG["CV", 1:2]
    } else if (treatmentEach == "W4") {
      gCorArrayAll[rowInd, 4, VIeach] <- corValueG["CV", 1:2]
    }
    
    corValueG <- round(corValueG, digits = 2)
    
    # get file info
    resCovRead <- resCovFiles[fileIndex]
    
    # make cor plot "g"
    resCovPath <- paste0(saveFolder, "/", resCovRead)
    CorValueRes <- MakeCorMat(resCovPath, c(index, "CV"))
    CorValueRes <- round(CorValueRes, digits = 2)
    
    # extract the phenotypic correaltion
    dataInd <- paste(treatmentEach, VIeach, eachYear, sep = "_")
    corValuePheno <- resultAll[[eachDataInd]][[dataInd]]
    
    # make the matrix 
    corMat <- matrix(NA, nrow = nrow(corValuePheno), ncol = ncol(corValuePheno))
    rownames(corMat) <- colnames(corMat) <- rownames(corValueG)
    diag(corMat) <- 1
    corMat[upper.tri(corMat)] <- corValueG[upper.tri(corMat)]
    corMat[lower.tri(corMat)] <- corValuePheno[lower.tri(corMat)]
    
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

# save the result
for (i in 1:dim(gCorArrayAll)[3]) {
  VI <- dimnames(gCorArrayAll)[[3]][i]
  write.csv(round(gCorArrayAll[, , i], digits = 2), 
            file = paste0(baseFolder, "gCor", VI, ".csv"))
}

# correlation test under bonferroni method for genetic correlation
cor2.test <- function(n, # the number of samples
                      r, # the estimated correlation
                      number,  # the number of test
                      sig.level = 0.05) {
  t <- abs(r) * sqrt((n - 2)/(1 - r^2))
  df <- n - 2
  p <- pt(t, df, lower.tail = FALSE) * 2
  significance <- (p < sig.level / number)
  return(significance)
}

gCorMat <- gCorArrayAll[, , 1]
sigMat <- matrix(NA, nrow = length(yearIndex), ncol = length(treatment))
rownames(sigMat) <- yearIndex
colnames(sigMat) <- treatment

n <- nrow(cv)
# number <- sum(!is.na(c(gCorMat)))
number <- 1
for (eachSigLevel in c(0.05, 0.01, 0.001)) {
  for (i in 1:nrow(sigMat)) {
    # i <- 1
    for (l in 1:ncol(sigMat)) {
      # l <- 1
      r <- gCorMat[i, l]
      if (is.na(r)) {
        next
      }
      significance <- cor2.test(n = n, r = r, number = number, sig.level = eachSigLevel)
      sigMat[i, l] <- significance
    }
  }
  write.csv(sigMat, 
            file = paste0(baseFolder, "sigTest", eachSigLevel, ".csv"))
}
