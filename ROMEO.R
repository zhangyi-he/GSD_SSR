#' @title ROMEO
#' @author Zhangyi He
#' @version 0.1

# set the directory
setwd("~/Dropbox/Jeffery He/iResearch/Trials/ROMEO")

# install and load required packages
if (!require("rpact")) {
  install.packages("rpact")
  library("rpact")
}

# call required functions
source("./Code/GSD_SSR/getConditionalPower.R")
source("./Code/GSD_SSR/getPredictivePower.R")
source("./Code/GSD_SSR/getTargetSampleSize.R")
source("./Code/GSD_SSR/getTestStatistic.R")
# source("./Code/GSD_SSR/getVirtualTrials.R")
# source("./Code/GSD_SSR/getOperatingCharacteristics.R")

################################################################################

#' calculate the hazard rates given cumulative incidence or cumulative survival in the presence of competing risks
calculateHazardRates <- function(cumulativeIncidence,
                                 cumulativeSurvival,
                                 fixedTimePoint) {
  hazardRates <- matrix(data = NA, nrow = 2, ncol = 2)
  if (is.null(cumulativeIncidence)) {
    hazardRates <- -log(cumulativeSurvival) / fixedTimePoint
  } else {
    hazardRates[1, ] <- cumulativeIncidence[1, ] *
      (-log(1 - cumulativeIncidence[1, ] - cumulativeIncidence[2, ]) /
         (fixedTimePoint * (cumulativeIncidence[1, ] + cumulativeIncidence[2, ])))
    hazardRates[2, ] <- cumulativeIncidence[2, ] *
      (-log(1 - cumulativeIncidence[1, ] - cumulativeIncidence[2, ]) /
         (fixedTimePoint * (cumulativeIncidence[1, ] + cumulativeIncidence[2, ])))
  }

  return(hazardRates)
}

# # validation with Example 1 from PASS Chapter 716
# cumulativeIncidence <- matrix(c(0.10, 0.05, 0.65, 0.65),
#                               nrow = 2, ncol = 2, byrow = TRUE)
# cumulativeSurvival <- NULL
# fixedTimePoint <- 3
#
# calculateHazardRates(
#   cumulativeIncidence = cumulativeIncidence,
#   cumulativeSurvival = cumulativeSurvival,
#   fixedTimePoint = fixedTimePoint)

# # validation with ROMEO from PASS
# cumulativeIncidence <- matrix(c(0.70, 0.7916, 0.26, 0.196),
#                               nrow = 2, ncol = 2, byrow = TRUE)
# cumulativeSurvival <- NULL
# fixedTimePoint <- 2 # months (60 days)
#
# calculateHazardRates(
#   cumulativeIncidence = cumulativeIncidence,
#   cumulativeSurvival = cumulativeSurvival,
#   fixedTimePoint = fixedTimePoint)
#
# cumulativeIncidence <- matrix(c(0.70, 0.7916, 0.26, 0.196),
#                               nrow = 2, ncol = 2, byrow = TRUE)
# cumulativeSurvival <- NULL
# fixedTimePoint <- 60 # days
#
# calculateHazardRates(
#   cumulativeIncidence = cumulativeIncidence,
#   cumulativeSurvival = cumulativeSurvival,
#   fixedTimePoint = fixedTimePoint)

########################################

#' calculate the number of events given the number of subjects in the presence of competing risks
calculateNumberOfEvents <- function(numberOfSubjects,
                                    hazardRates,
                                    maxStudyDuration,
                                    followUpTime,
                                    allocationRatioPlanned) {
  lambda1ev <- hazardRates[1, 1]
  lambda2ev <- hazardRates[1, 2]
  lambda1cr <- hazardRates[2, 1]
  lambda2cr <- hazardRates[2, 2]

  pi1ev <- (lambda1ev / (lambda1ev + lambda1cr)) *
    (1 - (exp(-followUpTime * (lambda1ev + lambda1cr)) -
            exp(-(maxStudyDuration * (lambda1ev + lambda1cr)))) /
       ((maxStudyDuration - followUpTime) * (lambda1ev + lambda1cr)))
  pi2ev <- (lambda2ev / (lambda2ev + lambda2cr)) *
    (1 - (exp(-followUpTime * (lambda2ev + lambda2cr)) -
            exp(-(maxStudyDuration * (lambda2ev + lambda2cr)))) /
       ((maxStudyDuration - followUpTime) * (lambda2ev + lambda2cr)))

  numberOfEvents <- numberOfSubjects *
    (pi1ev * allocationRatioPlanned / (allocationRatioPlanned + 1) +
       pi2ev * 1 / (allocationRatioPlanned + 1))

  # return(floor(numberOfEvents))
  return(numberOfEvents)
}

# # validation with Example 1 from PASS Chapter 716
# cumulativeIncidence <- matrix(c(0.10, 0.05, 0.65, 0.65),
#                               nrow = 2, ncol = 2, byrow = TRUE)
# cumulativeSurvival <- NULL
# fixedTimePoint <- 3
# hazardRates <- calculateHazardRates(cumulativeIncidence,
#                                     cumulativeSurvival,
#                                     fixedTimePoint)
# numberOfSubjects <- 100 * (1 - 0.1)
# maxStudyDuration <- 7
# followUpTime <- 3
# allocationRatioPlanned <- 1
#
# calculateNumberOfEvents(
#   numberOfSubjects = numberOfSubjects,
#   hazardRates = hazardRates,
#   maxStudyDuration = maxStudyDuration,
#   followUpTime = followUpTime,
#   allocationRatioPlanned = allocationRatioPlanned)

# # validation with ROMEO from PASS
# numberOfSubjects <- 550
# hazardRates <- matrix(c(1.1735, 1.7603, 0.4359, 0.4359),
#                       nrow = 2, ncol = 2, byrow = TRUE)
# maxStudyDuration <- 39 + 12
# followUpTime <- 12
# allocationRatioPlanned <- 1
# calculateNumberOfEvents(
#   numberOfSubjects = numberOfSubjects,
#   hazardRates = hazardRates,
#   maxStudyDuration = maxStudyDuration,
#   followUpTime = followUpTime,
#   allocationRatioPlanned = allocationRatioPlanned)

########################################

#' calculate the number of subjects given the number of events in the presence of competing risks
calculateNumberOfSubjects <- function(numberOfEvents,
                                      hazardRates,
                                      maxStudyDuration,
                                      followUpTime,
                                      allocationRatioPlanned) {
  lambda1ev <- hazardRates[1, 1]
  lambda2ev <- hazardRates[1, 2]
  lambda1cr <- hazardRates[2, 1]
  lambda2cr <- hazardRates[2, 2]

  pi1ev <- (lambda1ev / (lambda1ev + lambda1cr)) *
    (1 - (exp(-followUpTime * (lambda1ev + lambda1cr)) -
            exp(-(maxStudyDuration * (lambda1ev + lambda1cr)))) /
       ((maxStudyDuration - followUpTime) * (lambda1ev + lambda1cr)))
  pi2ev <- (lambda2ev / (lambda2ev + lambda2cr)) *
    (1 - (exp(-followUpTime * (lambda2ev + lambda2cr)) -
            exp(-(maxStudyDuration * (lambda2ev + lambda2cr)))) /
       ((maxStudyDuration - followUpTime) * (lambda2ev + lambda2cr)))

  numberOfSubjects <- numberOfEvents /
    (pi1ev * allocationRatioPlanned / (allocationRatioPlanned + 1) +
       pi2ev * 1 / (allocationRatioPlanned + 1))

  # return(ceiling(numberOfSubjects))
  return(numberOfSubjects)
}

# # validation with Example 1 from PASS Chapter 716
# cumulativeIncidence <- matrix(c(0.10, 0.05, 0.65, 0.65),
#                               nrow = 2, ncol = 2, byrow = TRUE)
# cumulativeSurvival <- NULL
# fixedTimePoint <- 3
# hazardRates <- calculateHazardRates(cumulativeIncidence,
#                                     cumulativeSurvival,
#                                     fixedTimePoint)
# numberOfEvents <- 8
# maxStudyDuration <- 7
# followUpTime <- 3
# allocationRatioPlanned <- 1
#
# calculateNumberOfSubjects(
#   numberOfEvents = numberOfEvents,
#   hazardRates = hazardRates,
#   maxStudyDuration = maxStudyDuration,
#   followUpTime = followUpTime,
#   allocationRatioPlanned = allocationRatioPlanned)

# # validation with ROMEO from PASS
# numberOfEvents <- 420
# hazardRates <- matrix(c(1.1735, 1.7603, 0.4359, 0.4359),
#                       nrow = 2, ncol = 2, byrow = TRUE)
# maxStudyDuration <- 39 + 12
# followUpTime <- 12
# allocationRatioPlanned <- 1
# calculateNumberOfSubjects(
#   numberOfEvents = numberOfEvents,
#   hazardRates = hazardRates,
#   maxStudyDuration = maxStudyDuration,
#   followUpTime = followUpTime,
#   allocationRatioPlanned = allocationRatioPlanned)

################################################################################

# get the group sequential design with O'Brien Flemming stopping boundaries (beta spending)
kMax <- 3
alpha <- 0.025
beta = 0.1
sided = 1
informationRates = c(0.4, 0.6, 1)
typeOfDesign = "noEarlyEfficacy"
typeBetaSpending = "bsOF"

design <- getDesignGroupSequential(
  kMax = kMax,
  alpha = alpha,
  beta = beta,
  sided = sided,
  informationRates = informationRates,
  typeOfDesign = typeOfDesign,
  typeBetaSpending = typeBetaSpending)
design

# get the group sequential design plan
design <- design
lambda1 <- 1.7603 # the hazard rate for the event of interest in the treatment group (from PASS)
lambda2 <- 1.1735 # the hazard rate for the event of interest in the control group (from PASS)
allocationRatioPlanned <- 1
accrualTime <- c(0, 39)
followUpTime <- 12

designPlan <- getSampleSizeSurvival(
  design = design,
  lambda1 = lambda1,
  lambda2 = lambda2,
  allocationRatioPlanned = allocationRatioPlanned,
  accrualTime = accrualTime,
  followUpTime = followUpTime)
designPlan

# get the upper bound for test statistic for sample size reestimation
targetPower <- 1 - beta
targetPowerType <- "conditionalPower"
alternative <- "greater"
alpha <- tail(x = design$stageLevels, n = 1)
targetNumberOfEvents <- ceiling(x = designPlan$maxNumberOfEvents)
actualNumberOfEvents <- ceiling(x = designPlan$cumulativeEventsPerStage[kMax - 1, 1])
allocationRatioPlanned <- 1
hazardRatio <- lambda1 / lambda2

testStatisticUpperBound <- getTestStatisticSurvival(
  targetPower = targetPower,
  targetPowerType = targetPowerType,
  alternative = alternative,
  alpha = alpha,
  targetNumberOfEvents = targetNumberOfEvents,
  actualNumberOfEvents = actualNumberOfEvents,
  allocationRatioPlanned = allocationRatioPlanned,
  hazardRatio = hazardRatio)
testStatisticUpperBound

# get the target number of events from the target number of subjects
targetNumberOfSubjects <- 550
hazardRates <- matrix(c(1.1735, 1.7603, 0.4359, 0.4359),
                      nrow = 2, ncol = 2, byrow = TRUE)
maxStudyDuration <- 39 + 12
followUpTime <- 12
allocationRatioPlanned <- 1

targetNumberOfEvents <- calculateNumberOfEvents(
  numberOfSubjects = targetNumberOfSubjects,
  hazardRates = hazardRates,
  maxStudyDuration = maxStudyDuration,
  followUpTime = followUpTime,
  allocationRatioPlanned = allocationRatioPlanned)
targetNumberOfEvents <- floor(x = targetNumberOfEvents)
targetNumberOfEvents

# get the lower bound for test statistic for sample size reestimation
targetPower <- 0.8
targetPowerType <- "conditionalPower"
alternative <- "greater"
alpha <- tail(x = design$stageLevels, n = 1)
targetNumberOfEvents <- targetNumberOfEvents
actualNumberOfEvents <- ceiling(x = designPlan$cumulativeEventsPerStage[kMax - 1, 1])
allocationRatioPlanned <- 1
hazardRatio <- 1.33 # minimum clinical important difference

testStatisticLowerBound <- getTestStatisticSurvival(
  targetPower = targetPower,
  targetPowerType = targetPowerType,
  alternative = alternative,
  alpha = alpha,
  targetNumberOfEvents = targetNumberOfEvents,
  actualNumberOfEvents = actualNumberOfEvents,
  allocationRatioPlanned = allocationRatioPlanned,
  hazardRatio = hazardRatio)
testStatisticLowerBound

jpeg(file = "./Output/Output v1.0/GSD_SSR_SampleSizeRule1_CriticalValue.jpeg", width = 8, height = 6, units = "in", res = 500)
maxNumberOfEvents <- c(ceiling(x = designPlan$cumulativeEventsPerStage[kMax - 1, 1]),
                       ceiling(x = designPlan$maxNumberOfEvents),
                       targetNumberOfEvents,
                       ceiling(x = designPlan$maxNumberOfEvents),
                       ceiling(x = designPlan$maxNumberOfEvents))

testStatistic <- c(qnorm(p = 0.60, mean = 0, sd = 1),
                   tail(x = design$futilityBounds, n = 1),
                   testStatisticLowerBound,
                   testStatisticUpperBound,
                   qnorm(p = 0.99, mean = 0, sd = 1))
plot(x = testStatistic, y = maxNumberOfEvents, type = "s", lwd = 1,
     main = "Sample size rule (stage 2)", xlab = "Critical value", ylab = "Number of events")
dev.off()

jpeg(file = "./Output/Output v1.0/GSD_SSR_SampleSizeRule1_ConditionalPower.jpeg", width = 8, height = 6, units = "in", res = 500)
conditionalPower <- getConditionalPowerSurvival(
  testStatistic = testStatistic,
  alternative = "greater",
  alpha = tail(x = design$stageLevels, n = 1),
  targetNumberOfEvents = ceiling(x = designPlan$maxNumberOfEvents),
  actualNumberOfEvents = ceiling(x = designPlan$cumulativeEventsPerStage[kMax - 1, 1]),
  allocationRatioPlanned = 1,
  hazardRatio = 1.5)
plot(x = conditionalPower, y = maxNumberOfEvents, type = "s", lwd = 1,
     main = "Sample size rule (stage 2)", xlab = "Conditional power", ylab = "Number of events")
dev.off()


maxNumberOfEvents <- targetNumberOfEvents:ceiling(x = designPlan$maxNumberOfEvents)
testStatistic <- NULL
for (i in 1:length(maxNumberOfEvents)) {
  testStatistic <-
    c(testStatistic,
      getTestStatisticSurvival(
        targetPower = targetPower,
        targetPowerType = targetPowerType,
        alternative = alternative,
        alpha = alpha,
        targetNumberOfEvents = maxNumberOfEvents[i],
        actualNumberOfEvents = actualNumberOfEvents,
        allocationRatioPlanned = allocationRatioPlanned,
        hazardRatio = 1.33))
}

jpeg(file = "./Output/Output v1.0/GSD_SSR_SampleSizeRule2_CriticalValue.jpeg", width = 8, height = 6, units = "in", res = 500)
maxNumberOfEvents <- c(actualNumberOfEvents,
                       ceiling(x = designPlan$maxNumberOfEvents),
                       maxNumberOfEvents,
                       ceiling(x = designPlan$maxNumberOfEvents))

testStatistic <- c(qnorm(p = 0.60, mean = 0, sd = 1),
                   tail(x = design$futilityBounds, n = 1),
                   testStatistic,
                   qnorm(p = 0.99, mean = 0, sd = 1))
plot(x = testStatistic, y = maxNumberOfEvents, type = "s", lwd = 1,
     main = "Sample size rule (stage 2)", xlab = "Critical value", ylab = "Number of events")
dev.off()

jpeg(file = "./Output/Output v1.0/GSD_SSR_SampleSizeRule2_ConditionalPower.jpeg", width = 8, height = 6, units = "in", res = 500)
conditionalPower <- getConditionalPowerSurvival(
  testStatistic = testStatistic,
  alternative = "greater",
  alpha = tail(x = design$stageLevels, n = 1),
  targetNumberOfEvents = ceiling(x = designPlan$maxNumberOfEvents),
  actualNumberOfEvents = ceiling(x = designPlan$cumulativeEventsPerStage[kMax - 1, 1]),
  allocationRatioPlanned = 1,
  hazardRatio = 1.5)
plot(x = conditionalPower, y = maxNumberOfEvents, type = "s", lwd = 1,
     main = "Sample size rule (stage 2)", xlab = "Critical value", ylab = "Number of events")
dev.off()

################################################################################

