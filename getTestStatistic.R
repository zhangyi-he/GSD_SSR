#' @title Get test statistic
#' @author Zhangyi He

source("./getConditionalPower.R")
source("./getPreditivePower.R")

################################################################################

#' return the test statistic given the target power and sample size for testing means in two samples
getTestStatisticMeans <- function(targetPower,
                                  targetPowerType,
                                  sided,
                                  alpha,
                                  targetNumberOfSubjects,
                                  actualNumberOfSubjects,
                                  allocationRatioPlanned,
                                  mu,
                                  sigma,
                                  ...,
                                  alternative = NULL) {
  if (targetPowerType == "conditionalPower") {
    # calculate the parameter being tested by the hypothesis
    theta <- mu[2] - mu[1]

    # calculate the actual information level
    actualInformationLevel <- 1 /
      (sigma[1]^2 / (actualNumberOfSubjects * (allocationRatioPlanned / (allocationRatioPlanned + 1))) +
         sigma[2]^2 / (actualNumberOfSubjects * (1 / (allocationRatioPlanned + 1))))
    # calculate the target information level
    targetInformationLevel <- 1 /
      (sigma[1]^2 / (targetNumberOfSubjects * (allocationRatioPlanned / (allocationRatioPlanned + 1))) +
         sigma[2]^2 / (targetNumberOfSubjects * (1 / (allocationRatioPlanned + 1))))

    if (sided == 1) {
      if (is.null(alternative)) {
        testStatistic <- 
          (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(targetInformationLevel) +
             theta * (targetInformationLevel - actualInformationLevel) +
             sqrt(targetInformationLevel - actualInformationLevel) *
             qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(actualInformationLevel))
        testStatistic <- 
          cbind(testStatistic,
                (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(targetInformationLevel) -
                   theta * (targetInformationLevel - actualInformationLevel) +
                   sqrt(targetInformationLevel - actualInformationLevel) *
                   qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(actualInformationLevel),
                deparse.level = 0)
      } else {
        if (alternative == "less") {
          testStatistic <- 
            (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(targetInformationLevel) +
               theta * (targetInformationLevel - actualInformationLevel) +
               sqrt(targetInformationLevel - actualInformationLevel) *
               qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(actualInformationLevel))
        } else {
          testStatistic <- 
            (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(targetInformationLevel) -
               theta * (targetInformationLevel - actualInformationLevel) +
               sqrt(targetInformationLevel - actualInformationLevel) *
               qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(actualInformationLevel)
        }
      }
    } else {
      testStatistic <- 
        (qnorm(p = 1 - alpha / 2, mean = 0, sd = 1) * sqrt(targetInformationLevel) +
           theta * (targetInformationLevel - actualInformationLevel) +
           sqrt(targetInformationLevel - actualInformationLevel) *
           qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(actualInformationLevel))
      testStatistic <- 
        cbind(testStatistic,
              (qnorm(p = 1 - alpha / 2, mean = 0, sd = 1) * sqrt(targetInformationLevel) -
                 theta * (targetInformationLevel - actualInformationLevel) +
                 sqrt(targetInformationLevel - actualInformationLevel) *
                 qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(actualInformationLevel),
              deparse.level = 0)
      testStatistic <- seq(from = min(testStatistic), to = max(testStatistic), length.out = 1e+06)

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

      testStatistic <- testStatistic[which(conditionalPower >= targetPower)]
      testStatistic <- c(min(testStatistic), max(testStatistic))
    }
  } else {
    # calculate the actual information level
    actualInformationLevel <- 1 /
      (sigma[1]^2 / (actualNumberOfSubjects * (allocationRatioPlanned / (allocationRatioPlanned + 1))) +
         sigma[2]^2 / (actualNumberOfSubjects * (1 / (allocationRatioPlanned + 1))))
    # calculate the target information level
    targetInformationLevel <- 1 /
      (sigma[1]^2 / (targetNumberOfSubjects * (allocationRatioPlanned / (allocationRatioPlanned + 1))) +
         sigma[2]^2 / (targetNumberOfSubjects * (1 / (allocationRatioPlanned + 1))))

    if (sided == 1) {
      if (is.null(alternative)) {
        testStatistic <- 
          (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
             sqrt(targetInformationLevel - actualInformationLevel) *
             qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(targetInformationLevel))
        testStatistic <- 
          cbind(testStatistic,
                (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
                   sqrt(targetInformationLevel - actualInformationLevel) *
                   qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(targetInformationLevel),
                deparse.level = 0)

      } else {
        if (alternative == "less") {
          testStatistic <- 
            (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
               sqrt(targetInformationLevel - actualInformationLevel) *
               qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(targetInformationLevel))
        } else {
          testStatistic <- 
            (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
               sqrt(targetInformationLevel - actualInformationLevel) *
               qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(targetInformationLevel)
        }
      }
    } else {
      testStatistic <- 
        (qnorm(p = 1 - alpha / 2, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
           sqrt(targetInformationLevel - actualInformationLevel) *
           qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(targetInformationLevel))
      testStatistic <- 
        cbind(testStatistic,
              (qnorm(p = 1 - alpha / 2, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
                 sqrt(targetInformationLevel - actualInformationLevel) *
                 qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(targetInformationLevel),
              deparse.level = 0)
      testStatistic <- seq(from = min(testStatistic), to = max(testStatistic), length.out = 1e+06)

      predictivePower <- getPredictivePowerMeans(
        testStatistic = testStatistic,
        sided = sided,
        alpha = alpha,
        actualNumberOfSubjects = actualNumberOfSubjects,
        targetNumberOfSubjects = targetNumberOfSubjects,
        allocationRatioPlanned = allocationRatioPlanned,
        sigma = sigma,
        alternative = alternative)

      testStatistic <- testStatistic[which(predictivePower >= targetPower)]
      testStatistic <- c(min(testStatistic), max(testStatistic))
    }
  }

  return(testStatistic)
}

# # validation with Example 2 from PASS Chapter 433
# targetPower <- 0.43342
# targetPowerType <- "conditionalPower"
# sided <- 2
# alpha <- 0.05
# targetNumberOfSubjects <- 120
# actualNumberOfSubjects <- 60
# allocationRatioPlanned <- 1
# mu <- c(0, 0.5)
# sigma <- c(4, 4)
# 
# getTestStatisticMeans(
#   targetPower = targetPower,
#   targetPowerType = targetPowerType,
#   sided = sided,
#   alpha = alpha,
#   targetNumberOfSubjects = targetNumberOfSubjects,
#   actualNumberOfSubjects = actualNumberOfSubjects,
#   allocationRatioPlanned = allocationRatioPlanned,
#   mu = mu,
#   sigma = sigma)
# 
# targetPower <- 0.8504
# targetPowerType <- "predictivePower"
# sided <- 2
# alpha <- 0.05
# targetNumberOfSubjects <- 120
# actualNumberOfSubjects <- 60
# allocationRatioPlanned <- 1
# mu <- c(0, 0.5)
# sigma <- c(4, 4)
# 
# getTestStatisticMeans(
#   targetPower = targetPower,
#   targetPowerType = targetPowerType,
#   sided = sided,
#   alpha = alpha,
#   targetNumberOfSubjects = targetNumberOfSubjects,
#   actualNumberOfSubjects = actualNumberOfSubjects,
#   allocationRatioPlanned = allocationRatioPlanned,
#   mu = mu,
#   sigma = sigma)

########################################

#' return the test statistic given the target power and sample size for testing rates in two samples
getTestStatisticRates <- function(targetPower,
                                  targetPowerType,
                                  sided,
                                  alpha,
                                  targetNumberOfSubjects,
                                  actualNumberOfSubjects,
                                  allocationRatioPlanned,
                                  rate,
                                  ...,
                                  alternative = NULL) {
  if (targetPowerType == "conditionalPower") {
    # calculate the parameter being tested by the hypothesis
    theta <- rate[2] - rate[1]

    # calculate the actual information level
    actualInformationLevel <- 1 / (((rate[2] + rate[1]) / 2) * (1 - ((rate[2] + rate[1]) / 2))) * 1 /
      (1 / (actualNumberOfSubjects * (allocationRatioPlanned / (allocationRatioPlanned + 1))) +
         1 / (actualNumberOfSubjects * (1 / (allocationRatioPlanned + 1))))
    # calculate the target information level
    targetInformationLevel <- 1 / (((rate[2] + rate[1]) / 2) * (1 - ((rate[2] + rate[1]) / 2))) * 1 /
      (1 / (targetNumberOfSubjects * (allocationRatioPlanned / (allocationRatioPlanned + 1))) +
         1 / (targetNumberOfSubjects * (1 / (allocationRatioPlanned + 1))))

    if (sided == 1) {
      if (is.null(alternative)) {
        testStatistic <- 
          (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(targetInformationLevel) +
             theta * (targetInformationLevel - actualInformationLevel) +
             sqrt(targetInformationLevel - actualInformationLevel) *
             qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(actualInformationLevel))
        testStatistic <- 
          cbind(testStatistic,
                (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(targetInformationLevel) -
                   theta * (targetInformationLevel - actualInformationLevel) +
                   sqrt(targetInformationLevel - actualInformationLevel) *
                   qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(actualInformationLevel),
                deparse.level = 0)
      } else {
        if (alternative == "less") {
          testStatistic <- 
            (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(targetInformationLevel) +
               theta * (targetInformationLevel - actualInformationLevel) +
               sqrt(targetInformationLevel - actualInformationLevel) *
               qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(actualInformationLevel))
        } else {
          testStatistic <- 
            (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(targetInformationLevel) -
               theta * (targetInformationLevel - actualInformationLevel) +
               sqrt(targetInformationLevel - actualInformationLevel) *
               qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(actualInformationLevel)
        }
      }
    } else {
      testStatistic <- 
        (qnorm(p = 1 - alpha / 2, mean = 0, sd = 1) * sqrt(targetInformationLevel) +
           theta * (targetInformationLevel - actualInformationLevel) +
           sqrt(targetInformationLevel - actualInformationLevel) *
           qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(actualInformationLevel))
      testStatistic <- 
        cbind(testStatistic,
              (qnorm(p = 1 - alpha / 2, mean = 0, sd = 1) * sqrt(targetInformationLevel) -
                 theta * (targetInformationLevel - actualInformationLevel) +
                 sqrt(targetInformationLevel - actualInformationLevel) *
                 qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(actualInformationLevel),
              deparse.level = 0)
      testStatistic <- seq(from = min(testStatistic), to = max(testStatistic), length.out = 1e+06)

      conditionalPower <- getConditionalPowerRates(
        testStatistic = testStatistic,
        sided = sided,
        alpha = alpha,
        actualNumberOfSubjects = actualNumberOfSubjects,
        targetNumberOfSubjects = targetNumberOfSubjects,
        allocationRatioPlanned = allocationRatioPlanned,
        rate = rate,
        alternative = alternative)

      testStatistic <- testStatistic[which(conditionalPower >= targetPower)]
      testStatistic <- c(min(testStatistic), max(testStatistic))
    }
  } else {
    # calculate the actual information level
    actualInformationLevel <- 1 / (((rate[2] + rate[1]) / 2) * (1 - ((rate[2] + rate[1]) / 2))) * 1 /
      (1 / (actualNumberOfSubjects * (allocationRatioPlanned / (allocationRatioPlanned + 1))) +
         1 / (actualNumberOfSubjects * (1 / (allocationRatioPlanned + 1))))
    # calculate the target information level
    targetInformationLevel <- 1 / (((rate[2] + rate[1]) / 2) * (1 - ((rate[2] + rate[1]) / 2))) * 1 /
      (1 / (targetNumberOfSubjects * (allocationRatioPlanned / (allocationRatioPlanned + 1))) +
         1 / (targetNumberOfSubjects * (1 / (allocationRatioPlanned + 1))))

    if (sided == 1) {
      if (is.null(alternative)) {
        testStatistic <- 
          (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
             sqrt(targetInformationLevel - actualInformationLevel) *
             qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(targetInformationLevel))
        testStatistic <- 
          cbind(testStatistic,
                (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
                   sqrt(targetInformationLevel - actualInformationLevel) *
                   qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(targetInformationLevel),
                deparse.level = 0)
      } else {
        if (alternative == "less") {
          testStatistic <- 
            (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
               sqrt(targetInformationLevel - actualInformationLevel) *
               qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(targetInformationLevel))
        } else {
          testStatistic <- 
            (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
               sqrt(targetInformationLevel - actualInformationLevel) *
               qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(targetInformationLevel)
        }
      }
    } else {
      testStatistic <- 
        (qnorm(p = 1 - alpha / 2, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
           sqrt(targetInformationLevel - actualInformationLevel) *
           qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(targetInformationLevel))
      testStatistic <- 
        cbind(testStatistic,
              (qnorm(p = 1 - alpha / 2, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
                 sqrt(targetInformationLevel - actualInformationLevel) *
                 qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(targetInformationLevel),
              deparse.level = 0)
      testStatistic <- seq(from = min(testStatistic), to = max(testStatistic), length.out = 1e+06)

      predictivePower <- getPredictivePowerRates(
        testStatistic = testStatistic,
        sided = sided,
        alpha = alpha,
        actualNumberOfSubjects = actualNumberOfSubjects,
        targetNumberOfSubjects = targetNumberOfSubjects,
        allocationRatioPlanned = allocationRatioPlanned,
        rate = rate,
        alternative = alternative)

      testStatistic <- testStatistic[which(predictivePower >= targetPower)]
      testStatistic <- c(min(testStatistic), max(testStatistic))
    }
  }

  return(testStatistic)
}

# # validation with Example 2 from PASS Chapter 202
# targetPower <- 0.16858
# targetPowerType <- "conditionalPower"
# sided <- 1
# alpha <- 0.025
# targetNumberOfSubjects <- 120
# actualNumberOfSubjects <- 60
# allocationRatioPlanned <- 1
# rate <- c(0.6, 0.7)
# alternative <- "greater"
# 
# getTestStatisticRates(
#   targetPower = targetPower,
#   targetPowerType = targetPowerType,
#   sided = sided,
#   alpha = alpha,
#   targetNumberOfSubjects = targetNumberOfSubjects,
#   actualNumberOfSubjects = actualNumberOfSubjects,
#   allocationRatioPlanned = allocationRatioPlanned,
#   rate = rate,
#   alternative = alternative)
# 
# targetPower <- 0.29262
# targetPowerType <- "predictivePower"
# sided <- 1
# alpha <- 0.025
# targetNumberOfSubjects <- 120
# actualNumberOfSubjects <- 60
# allocationRatioPlanned <- 1
# rate <- c(0.6, 0.7)
# alternative <- "greater"
# 
# getTestStatisticRates(
#   targetPower = targetPower,
#   targetPowerType = targetPowerType,
#   sided = sided,
#   alpha = alpha,
#   targetNumberOfSubjects = targetNumberOfSubjects,
#   actualNumberOfSubjects = actualNumberOfSubjects,
#   allocationRatioPlanned = allocationRatioPlanned,
#   rate = rate,
#   alternative = alternative)

########################################

#' return the test statistic given the target power and sample size for testing the hazard ratio in two samples
getTestStatisticSurvival <- function(targetPower,
                                     targetPowerType,
                                     sided,
                                     alpha,
                                     targetNumberOfEvents,
                                     actualNumberOfEvents,
                                     allocationRatioPlanned,
                                     hazardRate,
                                     ...,
                                     hazardRatio = NULL,
                                     alternative = NULL) {
  if (targetPowerType == "conditionalPower") {
    # calculate the parameter being tested by the hypothesis
    theta <- ifelse(is.null(hazardRatio), log(hazardRate[2] / hazardRate[1]), log(hazardRatio))

    # calculate the actual information level
    actualInformationLevel <- actualNumberOfEvents *
      allocationRatioPlanned / (allocationRatioPlanned + 1) * 1 / (allocationRatioPlanned + 1)
    # calculate the target information level
    targetInformationLevel <- targetNumberOfEvents *
      allocationRatioPlanned / (allocationRatioPlanned + 1) * 1 / (allocationRatioPlanned + 1)

    if (sided == 1) {
      if (is.null(alternative)) {
        testStatistic <- 
          (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(targetInformationLevel) +
             theta * (targetInformationLevel - actualInformationLevel) +
             sqrt(targetInformationLevel - actualInformationLevel) *
             qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(actualInformationLevel))
        testStatistic <- 
          cbind(testStatistic,
                (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(targetInformationLevel) -
                   theta * (targetInformationLevel - actualInformationLevel) +
                   sqrt(targetInformationLevel - actualInformationLevel) *
                   qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(actualInformationLevel),
                deparse.level = 0)
      } else {
        if (alternative == "less") {
          testStatistic <- 
            (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(targetInformationLevel) +
               theta * (targetInformationLevel - actualInformationLevel) +
               sqrt(targetInformationLevel - actualInformationLevel) *
               qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(actualInformationLevel))
        } else {
          testStatistic <- 
            (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(targetInformationLevel) -
               theta * (targetInformationLevel - actualInformationLevel) +
               sqrt(targetInformationLevel - actualInformationLevel) *
               qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(actualInformationLevel)
        }
      }
    } else {
      testStatistic <- 
        (qnorm(p = 1 - alpha / 2, mean = 0, sd = 1) * sqrt(targetInformationLevel) +
           theta * (targetInformationLevel - actualInformationLevel) +
           sqrt(targetInformationLevel - actualInformationLevel) *
           qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(actualInformationLevel))
      testStatistic <- 
        cbind(testStatistic,
              (qnorm(p = 1 - alpha / 2, mean = 0, sd = 1) * sqrt(targetInformationLevel) -
                 theta * (targetInformationLevel - actualInformationLevel) +
                 sqrt(targetInformationLevel - actualInformationLevel) *
                 qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(actualInformationLevel),
              deparse.level = 0)
      testStatistic <- seq(from = min(testStatistic), to = max(testStatistic), length.out = 1e+06)

      conditionalPower <- getConditionalPowerSurvival(
        testStatistic = testStatistic,
        sided = sided,
        alpha = alpha,
        targetNumberOfEvents = targetNumberOfEvents,
        actualNumberOfEvents = actualNumberOfEvents,
        allocationRatioPlanned = allocationRatioPlanned,
        hazardRatio = hazardRatio,
        alternative = alternative)

      testStatistic <- testStatistic[which(conditionalPower >= targetPower)]
      testStatistic <- c(min(testStatistic), max(testStatistic))
    }
  } else {
    # calculate the actual information level
    actualInformationLevel <- actualNumberOfEvents *
      allocationRatioPlanned / (allocationRatioPlanned + 1) * 1 / (allocationRatioPlanned + 1)
    # calculate the target information level
    targetInformationLevel <- targetNumberOfEvents *
      allocationRatioPlanned / (allocationRatioPlanned + 1) * 1 / (allocationRatioPlanned + 1)

    if (sided == 1) {
      if (is.null(alternative)) {
        testStatistic <- 
          (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
             sqrt(targetInformationLevel - actualInformationLevel) *
             qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(targetInformationLevel))
        testStatistic <- 
          cbind(testStatistic,
                (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
                   sqrt(targetInformationLevel - actualInformationLevel) *
                   qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(targetInformationLevel),
                deparse.level = 0)
      } else {
        if (alternative == "less") {
          testStatistic <- 
            (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
               sqrt(targetInformationLevel - actualInformationLevel) *
               qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(targetInformationLevel))
        } else {
          testStatistic <- 
            (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
               sqrt(targetInformationLevel - actualInformationLevel) *
               qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(targetInformationLevel)
        }
      }
    } else {
      testStatistic <- 
        (qnorm(p = 1 - alpha / 2, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
           sqrt(targetInformationLevel - actualInformationLevel) *
           qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(targetInformationLevel))
      testStatistic <- 
        cbind(testStatistic,
              (qnorm(p = 1 - alpha / 2, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
                 sqrt(targetInformationLevel - actualInformationLevel) *
                 qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(targetInformationLevel),
              deparse.level = 0)
      testStatistic <- seq(from = min(testStatistic), to = max(testStatistic), length.out = 1e+06)

      predictivePower <- getPredictivePowerSurvival(
        testStatistic = testStatistic,
        sided = sided,
        alpha = alpha,
        actualNumberOfEvents = actualNumberOfEvents,
        targetNumberOfEvents = targetNumberOfEvents,
        allocationRatioPlanned = allocationRatioPlanned,
        alternative = alternative)

      testStatistic <- testStatistic[which(predictivePower >= targetPower)]
      testStatistic <- c(min(testStatistic), max(testStatistic))
    }
  }

  return(testStatistic)
}

# # validation with Example 2 from PASS Chapter 701
# targetPower <- 0.63454
# targetPowerType <- "conditionalPower"
# sided <- 1
# alpha <- 0.025
# targetNumberOfEvents <- 200
# actualNumberOfEvents <- 100
# allocationRatioPlanned <- 1
# hazardRatio <- 0.8
# alternative <- "less"
# 
# getTestStatisticSurvival(
#   targetPower = targetPower,
#   targetPowerType = targetPowerType,
#   testStatistic = testStatistic,
#   sided = sided,
#   alpha = alpha,
#   targetNumberOfEvents = targetNumberOfEvents,
#   actualNumberOfEvents = actualNumberOfEvents,
#   allocationRatioPlanned = allocationRatioPlanned,
#   hazardRatio = hazardRatio,
#   alternative = alternative)
# 
# targetPower <- 0.80743
# targetPowerType <- "predictivePower"
# sided <- 1
# alpha <- 0.025
# targetNumberOfEvents <- 200
# actualNumberOfEvents <- 100
# allocationRatioPlanned <- 1
# hazardRatio <- 0.8
# alternative <- "less"
# 
# getTestStatisticSurvival(
#   targetPower = targetPower,
#   targetPowerType = targetPowerType,
#   sided = sided,
#   alpha = alpha,
#   targetNumberOfEvents = targetNumberOfEvents,
#   actualNumberOfEvents = actualNumberOfEvents,
#   allocationRatioPlanned = allocationRatioPlanned,
#   hazardRatio = hazardRatio,
#   alternative = alternative)

################################################################################
