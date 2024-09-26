#' @title Get the sample size
#' @author Zhangyi He

#' version 0.1

# # install and load required packages
# if (!require("MASS")) {
#   install.packages("MASS")
#   library("MASS")
# }

# # call the getConditionalPower functions
# source("./getConditionalPower.R")
# # call the getPreditivePower functions
# source("./getPreditivePower.R")

################################################################################

#' return the sample size given the target power for testing means in two samples
getSampleSizeMeans <- function(targetPower,
                               targetPowerType = c("conditional", "predictive"),
                               testStatistic,
                               alternative = c("two.sided", "less", "greater"),
                               alpha,
                               actualNumberOfSubjects,
                               allocationRatio,
                               mu,
                               sigma) {
  if (targetPowerType == "conditional") {
    requiredNumberOfSubjects <- actualNumberOfSubjects
    repeat {
      # requiredNumberOfSubjects <- requiredNumberOfSubjects +
      #   (1:500) * ifelse(length(strsplit(attr(fractions(allocationRatio), "fracs"), "/")[[1]]) == 2,
      #                    sum(as.integer(strsplit(attr(fractions(allocationRatio), "fracs"), "/")[[1]])),
      #                    allocationRatio + 1)
      requiredNumberOfSubjects <- requiredNumberOfSubjects + 1:1000
      conditionalPower <- getConditionalPowerMeans(
        testStatistic = testStatistic,
        alternative = alternative,
        alpha = alpha,
        actualNumberOfSubjects = actualNumberOfSubjects,
        maxNumberOfSubjects = requiredNumberOfSubjects,
        allocationRatio = allocationRatio,
        mu = mu,
        sigma = sigma)

      if (all(conditionalPower < targetPower)) {
        next
      } else if (all(conditionalPower > targetPower)) {
        requiredNumberOfSubjects <- NA
        break
      } else {
        requiredNumberOfSubjects <-
          requiredNumberOfSubjects[which.min(abs(conditionalPower - targetPower))]
        break
      }
    }
  } else {
    requiredNumberOfSubjects <- actualNumberOfSubjects
    repeat {
      # requiredNumberOfSubjects <- requiredNumberOfSubjects +
      #   (1:500) * ifelse(length(strsplit(attr(fractions(allocationRatio), "fracs"), "/")[[1]]) == 2,
      #                    sum(as.integer(strsplit(attr(fractions(allocationRatio), "fracs"), "/")[[1]])),
      #                    allocationRatio + 1)
      requiredNumberOfSubjects <- requiredNumberOfSubjects + 1:1000
      predictivePower <- getPredictivePowerMeans(
        testStatistic = testStatistic,
        alternative = alternative,
        alpha = alpha,
        actualNumberOfSubjects = actualNumberOfSubjects,
        maxNumberOfSubjects = requiredNumberOfSubjects,
        allocationRatio = allocationRatio,
        sigma = sigma)

      if (all(predictivePower < targetPower)) {
        next
      } else if (all(predictivePower > targetPower)) {
        requiredNumberOfSubjects <- NA
        break
      } else {
        requiredNumberOfSubjects <-
          requiredNumberOfSubjects[which.min(abs(predictivePower - targetPower))]
        break
      }
    }
  }

  return(requiredNumberOfSubjects)
}

# # validation with Example 3 from PASS Chapter 433
# targetPower <- 0.8
# targetPowerType <- "conditional"
# testStatistic <- 2.12
# alternative <- "two.sided"
# alpha <- 0.05
# actualNumberOfSubjects <- 60
# allocationRatio <- c(1, 1)
# mu <- c(0, 1.5)
# sigma <- c(6.7, 6.7)
#
# getSampleSizeMeans(
#   targetPower = targetPower,
#   targetPowerType = targetPowerType,
#   testStatistic = testStatistic,
#   alternative = alternative,
#   alpha = alpha,
#   actualNumberOfSubjects = actualNumberOfSubjects,
#   allocationRatio = allocationRatio,
#   mu = mu,
#   sigma = sigma)

########################################

#' return the sample size given the target power for testing rates in two samples
getSampleSizeRates <- function(targetPower,
                               targetPowerType = c("conditional", "predictive"),
                               testStatistic,
                               alternative = c("two.sided", "less", "greater"),
                               alpha,
                               actualNumberOfSubjects,
                               allocationRatio,
                               pi) {
  if (targetPowerType == "conditional") {
    requiredNumberOfSubjects <- actualNumberOfSubjects
    repeat {
      # requiredNumberOfSubjects <- requiredNumberOfSubjects +
      #   (1:500) * ifelse(length(strsplit(attr(fractions(allocationRatio), "fracs"), "/")[[1]]) == 2,
      #                    sum(as.integer(strsplit(attr(fractions(allocationRatio), "fracs"), "/")[[1]])),
      #                    allocationRatio + 1)
      requiredNumberOfSubjects <- requiredNumberOfSubjects + 1:1000
      conditionalPower <- getConditionalPowerRates(
        testStatistic = testStatistic,
        alternative = alternative,
        alpha = alpha,
        actualNumberOfSubjects = actualNumberOfSubjects,
        maxNumberOfSubjects = requiredNumberOfSubjects,
        allocationRatio = allocationRatio,
        pi = pi)

      if (all(conditionalPower < targetPower)) {
        next
      } else if (all(conditionalPower > targetPower)) {
        requiredNumberOfSubjects <- NA
        break
      } else {
        requiredNumberOfSubjects <-
          requiredNumberOfSubjects[which.min(abs(conditionalPower - targetPower))]
        break
      }
    }
  } else {
    requiredNumberOfSubjects <- actualNumberOfSubjects
    repeat {
      # requiredNumberOfSubjects <- requiredNumberOfSubjects +
      #   (1:500) * ifelse(length(strsplit(attr(fractions(allocationRatio), "fracs"), "/")[[1]]) == 2,
      #                    sum(as.integer(strsplit(attr(fractions(allocationRatio), "fracs"), "/")[[1]])),
      #                    allocationRatio + 1)
      requiredNumberOfSubjects <- requiredNumberOfSubjects + 1:1000
      predictivePower <- getPredictivePowerRates(
        testStatistic = testStatistic,
        alternative = alternative,
        alpha = alpha,
        actualNumberOfSubjects = actualNumberOfSubjects,
        maxNumberOfSubjects = requiredNumberOfSubjects,
        allocationRatio = allocationRatio,
        pi = pi)

      if (all(predictivePower < targetPower)) {
        next
      } else if (all(predictivePower > targetPower)) {
        requiredNumberOfSubjects <- NA
        break
      } else {
        requiredNumberOfSubjects <-
          requiredNumberOfSubjects[which.min(abs(predictivePower - targetPower))]
        break
      }
    }
  }

  return(requiredNumberOfSubjects)
}

# # validation with Example 3 from PASS Chapter 202
# targetPower <- 0.8
# targetPowerType <- "conditional"
# testStatistic <- 2.12
# alternative <- "greater"
# alpha <- 0.025
# actualNumberOfSubjects <- 60
# allocationRatio <- c(1, 1)
# pi <- c(0.643, 0.743)
#
# getSampleSizeRates(
#   targetPower = targetPower,
#   targetPowerType = targetPowerType,
#   testStatistic = testStatistic,
#   alternative = alternative,
#   alpha = alpha,
#   actualNumberOfSubjects = actualNumberOfSubjects,
#   allocationRatio = allocationRatio,
#   pi = pi)

########################################

#' return the sample size given the target power for testing the hazard ratio in two samples
getSampleSizeSurvival <- function(targetPower,
                                  targetPowerType = c("conditional", "predictive"),
                                  testStatistic,
                                  alternative = c("two.sided", "less", "greater"),
                                  alpha,
                                  actualNumberOfEvents,
                                  allocationRatio,
                                  lambda,
                                  ...,
                                  hazardRatio = NULL) {
  if (targetPowerType == "conditional") {
    requiredNumberOfEvents <- actualNumberOfEvents
    repeat {
      requiredNumberOfEvents <- requiredNumberOfEvents + (1:1000)
      conditionalPower <- getConditionalPowerSurvival(
        testStatistic = testStatistic,
        alternative = alternative,
        alpha = alpha,
        actualNumberOfEvents = actualNumberOfEvents,
        maxNumberOfEvents = requiredNumberOfEvents,
        allocationRatio = allocationRatio,
        lambda = lambda,
        hazardRatio = hazardRatio)

      if (all(conditionalPower < targetPower)) {
        next
      } else if (all(conditionalPower > targetPower)) {
        requiredNumberOfEvents <- NA
        break
      } else {
        requiredNumberOfEvents <-
          requiredNumberOfEvents[which.min(abs(conditionalPower - targetPower))]
        break
      }
    }
  } else {
    requiredNumberOfEvents <- actualNumberOfEvents
    repeat {
      requiredNumberOfEvents <- requiredNumberOfEvents + (1:1000)
      predictivePower <- getPredictivePowerSurvival(
        testStatistic = testStatistic,
        alternative = alternative,
        alpha = alpha,
        actualNumberOfEvents = actualNumberOfEvents,
        maxNumberOfEvents = requiredNumberOfEvents,
        allocationRatio = allocationRatio)

      if (all(predictivePower < targetPower)) {
        next
      } else if (all(predictivePower > targetPower)) {
        requiredNumberOfEvents <- NA
        break
      } else {
        requiredNumberOfEvents <-
          requiredNumberOfEvents[which.min(abs(predictivePower - targetPower))]
        break
      }
    }
  }

  return(requiredNumberOfEvents)
}

# # validation with Example 3 from PASS Chapter 701
# targetPower <- 0.8
# targetPowerType <- "conditional"
# testStatistic <- -2.12
# alternative <- "less"
# alpha <- 0.025
# actualNumberOfEvents <- 100
# allocationRatio <- c(1, 1)
# hazardRatio <- 0.8
#
# getSampleSizeSurvival(
#   targetPower = targetPower,
#   targetPowerType = targetPowerType,
#   testStatistic = testStatistic,
#   alternative = alternative,
#   alpha = alpha,
#   actualNumberOfEvents = actualNumberOfEvents,
#   allocationRatio = allocationRatio,
#   hazardRatio = hazardRatio)

################################################################################
