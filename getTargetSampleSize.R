#' @title Get target sample size
#' @author Zhangyi He

if (!require("MASS")) {
  install.packages("MASS")
  library("MASS")
}

source("./getConditionalPower.R")
source("./getPreditivePower.R")

################################################################################

#' return the target sample size given the target power for testing means in two samples
getTargetSampleSizeMeans <- function(targetPower,
                                     targetPowerType,
                                     testStatistic,
                                     sided,
                                     alpha,
                                     actualNumberOfSubjects,
                                     allocationRatioPlanned,
                                     mu,
                                     sigma,
                                     ...,
                                     alternative = NULL) {
  if (targetPowerType == "conditionalPower") {
    targetNumberOfSubjects <- actualNumberOfSubjects
    repeat {
      targetNumberOfSubjects <- targetNumberOfSubjects + 
        (1:500) * ifelse(length(strsplit(attr(fractions(allocationRatioPlanned), "fracs"), "/")[[1]]) == 2, 
                         sum(as.integer(strsplit(attr(fractions(allocationRatioPlanned), "fracs"), "/")[[1]])),
                         allocationRatioPlanned + 1)
      conditionalPower <- getConditionalPowerMeans(
        testStatistic = testStatistic,
        sided = sided,
        alpha = alpha,
        actualNumberOfSubjects = actualNumberOfSubjects,
        targetNumberOfSubjects = targetNumberOfSubjects,
        allocationRatioPlanned = allocationRatioPlanned,
        mu = mu,
        sigma = sigma,
        alternative = alternative)
      
      if (max(conditionalPower) > targetPower) {
        break
      } else {
        next
      }
    }
    
    targetNumberOfSubjects <- 
      targetNumberOfSubjects[which.min(abs(conditionalPower - targetPower))]
  } else {
    targetNumberOfSubjects <- actualNumberOfSubjects
    repeat {
      targetNumberOfSubjects <- targetNumberOfSubjects + 
        (1:500) * ifelse(length(strsplit(attr(fractions(allocationRatioPlanned), "fracs"), "/")[[1]]) == 2, 
                         sum(as.integer(strsplit(attr(fractions(allocationRatioPlanned), "fracs"), "/")[[1]])),
                         allocationRatioPlanned + 1)
      predictivePower <- getPredictivePowerMeans(
        testStatistic = testStatistic,
        sided = sided,
        alpha = alpha,
        actualNumberOfSubjects = actualNumberOfSubjects,
        targetNumberOfSubjects = targetNumberOfSubjects,
        allocationRatioPlanned = allocationRatioPlanned,
        sigma = sigma,
        alternative = alternative)
      
      if (max(predictivePower) > targetPower) {
        break
      } else {
        next
      }
    }
    
    targetNumberOfSubjects <- 
      targetNumberOfSubjects[which.min(abs(predictivePower - targetPower))]
  }
  
  return(targetNumberOfSubjects)
}

# # validation with Example 3 from PASS Chapter 433
# targetPower <- 0.8
# targetPowerType <- "conditionalPower"
# testStatistic <- 2.12
# sided <- 2
# alpha <- 0.05
# actualNumberOfSubjects <- 60
# allocationRatioPlanned <- 1
# mu <- c(0, 1.5)
# sigma <- c(6.7, 6.7)
# 
# getTargetSampleSizeMeans(
#   targetPower = targetPower,
#   targetPowerType = targetPowerType,
#   testStatistic = testStatistic,
#   sided = sided,
#   alpha = alpha,
#   actualNumberOfSubjects = actualNumberOfSubjects,
#   allocationRatioPlanned = allocationRatioPlanned,
#   mu = mu,
#   sigma = sigma)

########################################

#' return the target sample size given the target power for testing rates in two samples
getTargetSampleSizeRates <- function(targetPower,
                                     targetPowerType,
                                     testStatistic,
                                     sided,
                                     alpha,
                                     actualNumberOfSubjects,
                                     allocationRatioPlanned,
                                     rate,
                                     ...,
                                     alternative = NULL) {
  if (targetPowerType == "conditionalPower") {
    targetNumberOfSubjects <- actualNumberOfSubjects
    repeat {
      targetNumberOfSubjects <- targetNumberOfSubjects + 
        (1:500) * ifelse(length(strsplit(attr(fractions(allocationRatioPlanned), "fracs"), "/")[[1]]) == 2, 
                         sum(as.integer(strsplit(attr(fractions(allocationRatioPlanned), "fracs"), "/")[[1]])),
                         allocationRatioPlanned + 1)
      conditionalPower <- getConditionalPowerRates(
        testStatistic = testStatistic,
        sided = sided,
        alpha = alpha,
        actualNumberOfSubjects = actualNumberOfSubjects,
        targetNumberOfSubjects = targetNumberOfSubjects,
        allocationRatioPlanned = allocationRatioPlanned,
        rate = rate,
        alternative = alternative)
      
      if (max(conditionalPower) > targetPower) {
        break
      } else {
        next
      }
    }
    
    targetNumberOfSubjects <- 
      targetNumberOfSubjects[which.min(abs(conditionalPower - targetPower))]
  } else {
    targetNumberOfSubjects <- actualNumberOfSubjects
    repeat {
      targetNumberOfSubjects <- targetNumberOfSubjects + 
        (1:500) * ifelse(length(strsplit(attr(fractions(allocationRatioPlanned), "fracs"), "/")[[1]]) == 2, 
                         sum(as.integer(strsplit(attr(fractions(allocationRatioPlanned), "fracs"), "/")[[1]])),
                         allocationRatioPlanned + 1)
      predictivePower <- getPredictivePowerRates(
        testStatistic = testStatistic,
        sided = sided,
        alpha = alpha,
        actualNumberOfSubjects = actualNumberOfSubjects,
        targetNumberOfSubjects = targetNumberOfSubjects,
        allocationRatioPlanned = allocationRatioPlanned,
        rate = rate,
        alternative = alternative)
      
      if (max(predictivePower) > targetPower) {
        break
      } else {
        next
      }
    }
    
    targetNumberOfSubjects <- 
      targetNumberOfSubjects[which.min(abs(predictivePower - targetPower))]
  }
  
  return(targetNumberOfSubjects)
}

# # validation with Example 3 from PASS Chapter 202
# targetPower <- 0.8
# targetPowerType <- "conditionalPower"
# testStatistic <- 2.12
# sided <- 1
# alpha <- 0.025
# actualNumberOfSubjects <- 60
# allocationRatioPlanned <- 1
# rate <- c(0.643, 0.743)
# alternative <- "greater"
# 
# getTargetSampleSizeRates(
#   targetPower = targetPower,
#   targetPowerType = targetPowerType,
#   testStatistic = testStatistic,
#   sided = sided,
#   alpha = alpha,
#   actualNumberOfSubjects = actualNumberOfSubjects,
#   allocationRatioPlanned = allocationRatioPlanned,
#   rate = rate,
#   alternative = alternative)

########################################

#' return the target sample size given the target power for testing the hazard ratio in two samples
getTargetSampleSizeSurvival <- function(targetPower,
                                        targetPowerType,
                                        testStatistic,
                                        sided,
                                        alpha,
                                        actualNumberOfEvents,
                                        allocationRatioPlanned,
                                        hazardRate,
                                        ...,
                                        hazardRatio = NULL,
                                        alternative = NULL) {
  if (targetPowerType == "conditionalPower") {
    targetNumberOfEvents <- actualNumberOfEvents
    repeat {
      targetNumberOfEvents <- targetNumberOfEvents + (1:1000)
      conditionalPower <- getConditionalPowerSurvival(
        testStatistic = testStatistic,
        sided = sided,
        alpha = alpha,
        actualNumberOfEvents = actualNumberOfEvents,
        targetNumberOfEvents = targetNumberOfEvents,
        allocationRatioPlanned = allocationRatioPlanned,
        hazardRate = hazardRate,
        hazardRatio = hazardRatio,
        alternative = alternative)
      
      if (max(conditionalPower) > targetPower) {
        break
      } else {
        next
      }
    }
    
    targetNumberOfEvents <- 
      targetNumberOfEvents[which.min(abs(conditionalPower - targetPower))]
  } else {
    targetNumberOfEvents <- actualNumberOfEvents
    repeat {
      targetNumberOfEvents <- targetNumberOfEvents + (1:1000)
      predictivePower <- getPredictivePowerSurvival(
        testStatistic = testStatistic,
        sided = sided,
        alpha = alpha,
        actualNumberOfEvents = actualNumberOfEvents,
        targetNumberOfEvents = targetNumberOfEvents,
        allocationRatioPlanned = allocationRatioPlanned,
        alternative = alternative)
      
      if (max(predictivePower) > targetPower) {
        break
      } else {
        next
      }
    }
    
    targetNumberOfEvents <- 
      targetNumberOfEvents[which.min(abs(predictivePower - targetPower))]
  }
  
  return(targetNumberOfEvents)
}

# # validation with Example 3 from PASS Chapter 701
# targetPower <- 0.8
# targetPowerType <- "conditionalPower"
# testStatistic <- -2.12
# sided <- 1
# alpha <- 0.025
# actualNumberOfEvents <- 100
# allocationRatioPlanned <- 1
# hazardRatio <- 0.8
# alternative <- "less"
# 
# getTargetSampleSizeSurvival(
#   targetPower = targetPower,
#   targetPowerType = targetPowerType,
#   testStatistic = testStatistic,
#   sided = sided,
#   alpha = alpha,
#   actualNumberOfEvents = actualNumberOfEvents,
#   allocationRatioPlanned = allocationRatioPlanned,
#   hazardRatio = hazardRatio,
#   alternative = alternative)

################################################################################
