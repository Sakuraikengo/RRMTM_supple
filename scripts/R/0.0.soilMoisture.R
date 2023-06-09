machine <- c("M100", "M100", "P4M")
year <- c("2019", "2020", "2021")
source("scripts/R/1.functionCode.R")

options(stringsAsFactors = FALSE)

library(readr)
library(stringr)
library(ggplot2)
library(ggsci)
library(readxl)

# make the folder
baseFolder <- "midstream/0.0.soilMoisture/"
if (!dir.exists(baseFolder)) {
  dir.create(baseFolder, recursive = T)
}

# set the each year data
# the days when MS collected
sowingDayList <- list(as.Date("2019-07-10"), 
                      as.Date("2020-07-08"), 
                      as.Date("2021-07-06"))

DestructiveDayList <- list(as.Date("2019-09-10"), 
                           as.Date("2020-09-08"), 
                           as.Date("2021-09-06"))

thinningDayList <- list(as.Date("2019-07-24"), 
                        as.Date("2020-07-20"), 
                        as.Date("2021-07-21"))


MeasuringDayList <- list(c("08-25", "08-31", "09-03", "09-04"), 
                         c("08-22", "08-24", 
                           "08-27", "08-29", "09-01", "09-04"), 
                         c("08-22", "08-25", 
                           "08-27", "08-29", "08-30", "08-31", "09-05"))

soilMoistureList <- list.files(path = "reference", 
                               full.names = T)

for (eachYearInd in 1:length(year)) {
  # eachYearInd <- 1
  eachYear <- year[eachYearInd]
  
  ######load_soilMoisture#####
  soilMoisture <- data.frame(read_xlsx(path = soilMoistureList[[eachYearInd]], 
                                       skip = 1))
  
  # remove the data from October
  soilMoisture <- soilMoisture[1:100, ]
  day <- as.Date(soilMoisture[, 1])
  rownames(soilMoisture) <- day
  
  # # colnames(soilMoisture)
  # # split the data into the 4 ridge(block)
  # ridge1 <- soilMoisture[, 5:16]
  # ridge2 <- soilMoisture[, 17:28]
  # ridge3 <- soilMoisture[, 29:40]
  # ridge4 <- soilMoisture[, 41:52]
  
  # split the data into the 4 ridge(block) only block for 200 accessions
  ridge1 <- soilMoisture[, 7:16]
  ridge2 <- soilMoisture[, 17:26]
  ridge3 <- soilMoisture[, 31:40]
  ridge4 <- soilMoisture[, 41:50]
  
  ridgeList <- list(ridge1, ridge2, ridge3, ridge4)
  
  # calculate the mean of soil moisture in each treatment
  ridgeResList <- lapply(ridgeList, function(eachRidge) {
    # eachRidge <- ridgeList[[1]]
    moistureMean <- rowMeans(eachRidge, na.rm = T)
    moistureSd <- apply(eachRidge, 1, sd, na.rm = T)
    return(list(moistureMean, moistureSd))
  })
  
  if (eachYear == "2019") {
    moisture_condition_W <- data.frame(C = ridgeResList[[1]][[1]], 
                                       W5 = ridgeResList[[3]][[1]], 
                                       D10 = ridgeResList[[4]][[1]], 
                                       D = ridgeResList[[2]][[1]])
    moisture_condition_sd <- data.frame(C = ridgeResList[[1]][[2]], 
                                        W5 = ridgeResList[[3]][[2]], 
                                        D10 = ridgeResList[[4]][[2]], 
                                        D = ridgeResList[[2]][[2]])
  } else if (eachYear == "2020") {
    moisture_condition_W <- data.frame(C = ridgeResList[[4]][[1]], 
                                       D10 = ridgeResList[[2]][[1]], 
                                       W10 = ridgeResList[[3]][[1]], 
                                       D = ridgeResList[[1]][[1]])
    moisture_condition_sd <- data.frame(C = ridgeResList[[4]][[2]], 
                                        D10 = ridgeResList[[2]][[2]], 
                                        W10 = ridgeResList[[3]][[2]], 
                                        D = ridgeResList[[1]][[2]])
  } else if (eachYear == "2021") {
    moisture_condition_W <- data.frame(C = ridgeResList[[3]][[1]], 
                                       D10 = ridgeResList[[2]][[1]], 
                                       W10 = ridgeResList[[1]][[1]], 
                                       D = ridgeResList[[4]][[1]])
    moisture_condition_sd <- data.frame(C = ridgeResList[[3]][[2]], 
                                        D10 = ridgeResList[[2]][[2]], 
                                        W10 = ridgeResList[[1]][[2]], 
                                        D = ridgeResList[[4]][[2]])
  }
  moistureDfMean <- WideToLong(baseDf = moisture_condition_W, 
                               dfRowInd = day, 
                               dfColInd = colnames(moisture_condition_W))
  moistureDfSd <- WideToLong(baseDf = moisture_condition_sd, 
                             dfRowInd = day, 
                             dfColInd = colnames(moisture_condition_sd))
  
  
  moistureDf <- cbind(moistureDfSd$value, moistureDfMean)
  
  
  colnames(moistureDf) <- c("sd", "mean", "condition", "Day")
  moistureDf$Day <- as.Date(day)
  moistureDf <- na.omit(moistureDf)
  
  measuringDay <- as.Date(paste0(eachYear, "-", MeasuringDayList[[eachYearInd]]))
  
  thinningDay <- as.Date(thinningDayList[[eachYearInd]])
  moistureDf$condition <- factor(moistureDf$condition, 
                                 levels = c("C", "W5", "W10", "D10", "D"))
  
  # set the color of each treatment
  if (eachYear == "2019") {
    col <- c("blue3", "olivedrab3", "darkgoldenrod1", "chocolate2")
  } else {
    col <- c("blue3", "lightskyblue1", "darkgoldenrod1", "chocolate2")
  }
  g_W <- ggplot(moistureDf, aes(x = Day, y = mean, group = condition, color = condition)) + 
    geom_line(aes(color = condition), size = 1.2) + 
    geom_vline(xintercept = measuringDay, linetype = "dashed") +
    # geom_vline(xintercept = measuringDayPheno, linetype = "dashed") + 
    geom_vline(xintercept = thinningDay, linetype = "solid") + 
    scale_color_manual(values = col) + 
    ylim(c(2, 7.5)) + 
    # geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, width = 1)) + 
    labs(x = "Day", y = "Soil Moisture (%)", colour = "Treatment") + 
    theme(text = element_text(size = 16))
  png(paste0(baseFolder, "/soilMoisture_", eachYear, ".png"), 
      height = 1440, width = 2160, res = 216)
  print(g_W)
  dev.off()
  
  # figure which have only drought treatment
  moistureDfDrought <- moistureDf[moistureDf$condition != "C", ]
  g_D <- ggplot(moistureDfDrought, aes(x = Day, y = mean, group = condition, color = condition)) + 
    geom_line(aes(color = condition), size = 1.2) + 
    geom_vline(xintercept = measuringDay, linetype = "dashed") +
    # geom_vline(xintercept = measuringDayPheno, linetype = "dashed") + 
    geom_vline(xintercept = thinningDay, linetype = "solid") + 
    scale_color_manual(values = col[2:4]) + 
    ylim(c(2, 7.5)) + 
    # geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, width = 1)) + 
    labs(x = "Day", y = "Soil Moisture (%)", colour = "Treatment") + 
    theme(text = element_text(size = 16))
  png(paste0(baseFolder, "/soilMoistureDrought_", eachYear, ".png"), 
      height = 1440, width = 2160, res = 216)
  print(g_D)
  dev.off()
}