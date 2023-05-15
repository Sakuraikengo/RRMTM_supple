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

# make the folder
baseFolder <- "midstream/2.1.genomicPrediction/"
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


# save or read the seed
seedIndCsv <- paste0("midstream/seedInd.csv")
if (file.exists(seedIndCsv)) {
  seedInd <- read.csv(paste0(seedIndCsv), row.names = 1, header = T)
  seedInd <- c(as.matrix(seedInd))
} else {
  seedInd <- sample(1:500, 10, replace = F)
  write.csv(x = seedInd, file = paste0(seedIndCsv))
}

resultIndex <- c("Correlation", "R2", "RMSE")


# set the target trait
target <- scale(cv)

# make amat
amat <- amat0[rownames(target), rownames(target)]

# make amat list
amatZ <- design.Z(pheno.labels = rownames(target),
                  geno.names = rownames(amat))

amatList <- list(Z = amatZ,
                 K = amat)

# set the parameter
ZETA <- list(amat = amatList)
X0 <- NULL

# calculate the genomic heritability
fitEM3 <- EM3.cpp(y = target,
                  X0 = X0, n.core = 1, 
                  ZETA = ZETA)
h2 <- fitEM3$Vu / (fitEM3$Vu + fitEM3$Ve)
h2All <- c(fitEM3$Vu, fitEM3$Ve, h2)
names(h2All) <- c("Vu", "Ve", "h2")
write.csv(h2All, file = paste0(baseFolder, "h2.csv"))

# set the seed and cross-validation index
# rep5 <- rep(1:5, 1000)
# rep5 <- rep5[1:nrow(cv)]
rep10 <- rep(1:10, 1000)
rep10 <- rep10[1:nrow(cv)]


# make the case for the result
resultEachSeed <- matrix(NA, ncol = length(resultIndex), nrow = length(seedInd))
colnames(resultEachSeed) <- resultIndex
rownames(resultEachSeed) <- seedInd


for (seedIndEach in 1:length(seedInd)) {
  # seedIndEach <- 1
  set.seed(seedInd[seedIndEach])
  crosscvInd <- sample(rep10, nrow(cv), replace = F)
  
  # prepare the prediction box
  predictionDataRAINBOW <- rep(NA, nrow(cv))
  for (times in 1:10) {
    # times <- 1
    testInd <- crosscvInd == times
    
    targetNa <- target
    targetNa[testInd] <- NA
    
    # use RAINBOWR
    resEM3 <- EM3.cpp(y = targetNa,
                      X0 = X0, n.core = 1, 
                      ZETA = ZETA)
    
    # input the predicted data
    predictionDataRAINBOW[testInd] <- resEM3$y.pred[testInd]
  }
  # resEM3$weights
  predictData <- predictionDataRAINBOW
  obsData <- target
  # calculate the R2 and RMSE
  correlation <- cor(obsData, predictData)
  R2 <- 1 - sum((obsData - predictData) ^ 2) / sum((obsData - mean(obsData)) ^ 2)
  RMSE <- sqrt(sum((obsData - predictData) ^ 2) / length(obsData))
  
  # input the result
  resultEachSeed[seedIndEach, ] <- c(correlation, R2, RMSE)
  # arrayResultEachDay[ZETA, ] <- c(correlation, R2, RMSE)
  
  # make the plot
  xlim <- ylim <- range(predictData, obsData)
  
  png(paste0(baseFolder, "predictionPlot",
             seedIndEach, ".png"))
  plot((obsData), (predictData), 
       xlim = xlim, ylim = ylim, 
       main = paste0("RAINBOW prediction", " r = ", round(correlation, 2)))
  abline(0, 1, col = 2, lty = 2)
  dev.off()
}

resMat <- apply(resultEachSeed, 2, mean)
write.csv(resultEachSeed, file = paste0(baseFolder, "res_all.csv"))
write.csv(resMat, file = paste0(baseFolder, "result.csv"))
