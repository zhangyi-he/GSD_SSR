#' @title Get the test statistic
#' @author Zhangyi He

#' version 0.1

# # call the getConditionalPower functions
# source("./getConditionalPower.R")
# # call the getPreditivePower functions
# source("./getPreditivePower.R")

################################################################################

#' return the test statistic given the target power and maximum sample size for testing means in two samples
getTestStatisticMeans <- function(targetPower,
                                  targetPowerType = c("conditional", "predictive"),
                                  alternative = c("two.sided", "less", "greater"),
                                  alpha,
                                  actualNumberOfSubjects,
                                  maxNumberOfSubjects,
                                  allocationRatio,
                                  mu,
                                  sigma) {
  if (targetPowerType == "conditional") {
    # calculate the parameter being tested by the hypothesis
    theta <- mu[2] - mu[1]

    # calculate the actual information level
    actualInformationLevel <- 1 /
      (sigma[1]^2 / (actualNumberOfSubjects * (allocationRatio[1] / sum(allocationRatio))) +
         sigma[2]^2 / (actualNumberOfSubjects * (allocationRatio[2] / sum(allocationRatio))))
    # calculate the maximum information level
    maxInformationLevel <- 1 /
      (sigma[1]^2 / (maxNumberOfSubjects * (allocationRatio[1] / sum(allocationRatio))) +
         sigma[2]^2 / (maxNumberOfSubjects * (allocationRatio[2] / sum(allocationRatio))))

    if (alternative == "less") {
      testStatistic <-
        (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(maxInformationLevel) +
           theta * (maxInformationLevel - actualInformationLevel) +
           sqrt(maxInformationLevel - actualInformationLevel) *
           qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(actualInformationLevel))
    } else if (alternative == "greater") {
      testStatistic <-
        (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(maxInformationLevel) -
           theta * (maxInformationLevel - actualInformationLevel) +
           sqrt(maxInformationLevel - actualInformationLevel) *
           qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(actualInformationLevel)
    } else {
      testStatistic <-
        (qnorm(p = 1 - alpha / 2, mean = 0, sd = 1) * sqrt(maxInformationLevel) +
           theta * (maxInformationLevel - actualInformationLevel) +
           sqrt(maxInformationLevel - actualInformationLevel) *
           qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(actualInformationLevel))
      testStatistic <-
        cbind(testStatistic,
              (qnorm(p = 1 - alpha / 2, mean = 0, sd = 1) * sqrt(maxInformationLevel) -
                 theta * (maxInformationLevel - actualInformationLevel) +
                 sqrt(maxInformationLevel - actualInformationLevel) *
                 qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(actualInformationLevel),
              deparse.level = 0)
      testStatistic <- seq(from = min(testStatistic), to = max(testStatistic), length.out = 1e+06)

      conditionalPower <- getConditionalPowerMeans(
        testStatistic = testStatistic,
        alternative = alternative,
        alpha = alpha,
        actualNumberOfSubjects = actualNumberOfSubjects,
        maxNumberOfSubjects = maxNumberOfSubjects,
        allocationRatio = allocationRatio,
        mu = mu,
        sigma = sigma)

      testStatistic <- testStatistic[which(conditionalPower >= targetPower)]
      testStatistic <- c(min(testStatistic), max(testStatistic))
    }
  } else {
    # calculate the actual information level
    actualInformationLevel <- 1 /
      (sigma[1]^2 / (actualNumberOfSubjects * (allocationRatio[1] / sum(allocationRatio))) +
         sigma[2]^2 / (actualNumberOfSubjects * (allocationRatio[2] / sum(allocationRatio))))
    # calculate the maximum information level
    maxInformationLevel <- 1 /
      (sigma[1]^2 / (maxNumberOfSubjects * (allocationRatio[1] / sum(allocationRatio))) +
         sigma[2]^2 / (maxNumberOfSubjects * (allocationRatio[2] / sum(allocationRatio))))

    if (alternative == "less") {
      testStatistic <-
        (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
           sqrt(maxInformationLevel - actualInformationLevel) *
           qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(maxInformationLevel))
    } else if (alternative == "greater") {
      testStatistic <-
        (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
           sqrt(maxInformationLevel - actualInformationLevel) *
           qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(maxInformationLevel)
    } else {
      testStatistic <-
        (qnorm(p = 1 - alpha / 2, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
           sqrt(maxInformationLevel - actualInformationLevel) *
           qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(maxInformationLevel))
      testStatistic <-
        cbind(testStatistic,
              (qnorm(p = 1 - alpha / 2, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
                 sqrt(maxInformationLevel - actualInformationLevel) *
                 qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(maxInformationLevel),
              deparse.level = 0)
      testStatistic <- seq(from = min(testStatistic), to = max(testStatistic), length.out = 1e+06)

      predictivePower <- getPredictivePowerMeans(
        testStatistic = testStatistic,
        alternative = alternative,
        alpha = alpha,
        actualNumberOfSubjects = actualNumberOfSubjects,
        maxNumberOfSubjects = maxNumberOfSubjects,
        allocationRatio = allocationRatio,
        sigma = sigma)

      testStatistic <- testStatistic[which(predictivePower >= targetPower)]
      testStatistic <- c(min(testStatistic), max(testStatistic))
    }
  }

  return(testStatistic)
}

# # validation with Example 2 from PASS Chapter 433
# targetPower <- 0.43342
# targetPowerType <- "conditional"
# alternative <- "two.sided"
# alpha <- 0.05
# actualNumberOfSubjects <- 60
# maxNumberOfSubjects <- 120
# allocationRatio <- c(1, 1)
# mu <- c(0, 0.5)
# sigma <- c(4, 4)
#
# getTestStatisticMeans(
#   targetPower = targetPower,
#   targetPowerType = targetPowerType,
#   alternative = alternative,
#   alpha = alpha,
#   actualNumberOfSubjects = actualNumberOfSubjects,
#   maxNumberOfSubjects = maxNumberOfSubjects,
#   allocationRatio = allocationRatio,
#   mu = mu,
#   sigma = sigma)
#
# targetPower <- 0.8504
# targetPowerType <- "predictivePower"
# alternative <- "two.sided"
# alpha <- 0.05
# actualNumberOfSubjects <- 60
# maxNumberOfSubjects <- 120
# allocationRatio <- c(1, 1)
# mu <- c(0, 0.5)
# sigma <- c(4, 4)
#
# getTestStatisticMeans(
#   targetPower = targetPower,
#   targetPowerType = targetPowerType,
#   alternative = alternative,
#   alpha = alpha,
#   actualNumberOfSubjects = actualNumberOfSubjects,
#   maxNumberOfSubjects = maxNumberOfSubjects,
#   allocationRatio = allocationRatio,
#   mu = mu,
#   sigma = sigma)

########################################

#' return the test statistic given the target power and maximum sample size for testing rates in two samples
getTestStatisticRates <- function(targetPower,
                                  targetPowerType = c("conditional", "predictive"),
                                  alternative = c("two.sided", "less", "greater"),
                                  alpha,
                                  actualNumberOfSubjects,
                                  maxNumberOfSubjects,
                                  allocationRatio,
                                  pi) {
  if (targetPowerType == "conditional") {
    # calculate the parameter being tested by the hypothesis
    theta <- pi[2] - pi[1]

    # calculate the actual information level
    actualInformationLevel <- 1 / (((pi[2] + pi[1]) / 2) * (1 - ((pi[2] + pi[1]) / 2))) * 1 /
      (1 / (actualNumberOfSubjects * (allocationRatio[1] / sum(allocationRatio))) +
         1 / (actualNumberOfSubjects * (allocationRatio[2] / sum(allocationRatio))))
    # calculate the maximum information level
    maxInformationLevel <- 1 / (((pi[2] + pi[1]) / 2) * (1 - ((pi[2] + pi[1]) / 2))) * 1 /
      (1 / (maxNumberOfSubjects * (allocationRatio[1] / sum(allocationRatio))) +
         1 / (maxNumberOfSubjects * (allocationRatio[2] / sum(allocationRatio))))

    if (alternative == "less") {
      testStatistic <-
        (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(maxInformationLevel) +
           theta * (maxInformationLevel - actualInformationLevel) +
           sqrt(maxInformationLevel - actualInformationLevel) *
           qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(actualInformationLevel))
    } else if (alternative == "greater") {
      testStatistic <-
        (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(maxInformationLevel) -
           theta * (maxInformationLevel - actualInformationLevel) +
           sqrt(maxInformationLevel - actualInformationLevel) *
           qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(actualInformationLevel)
    } else {
      testStatistic <-
        (qnorm(p = 1 - alpha / 2, mean = 0, sd = 1) * sqrt(maxInformationLevel) +
           theta * (maxInformationLevel - actualInformationLevel) +
           sqrt(maxInformationLevel - actualInformationLevel) *
           qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(actualInformationLevel))
      testStatistic <-
        cbind(testStatistic,
              (qnorm(p = 1 - alpha / 2, mean = 0, sd = 1) * sqrt(maxInformationLevel) -
                 theta * (maxInformationLevel - actualInformationLevel) +
                 sqrt(maxInformationLevel - actualInformationLevel) *
                 qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(actualInformationLevel),
              deparse.level = 0)
      testStatistic <- seq(from = min(testStatistic), to = max(testStatistic), length.out = 1e+06)

      conditionalPower <- getConditionalPowerRates(
        testStatistic = testStatistic,
        alternative = alternative,
        alpha = alpha,
        actualNumberOfSubjects = actualNumberOfSubjects,
        maxNumberOfSubjects = maxNumberOfSubjects,
        allocationRatio = allocationRatio,
        pi = pi)

      testStatistic <- testStatistic[which(conditionalPower >= targetPower)]
      testStatistic <- c(min(testStatistic), max(testStatistic))
    }
  } else {
    # calculate the actual information level
    actualInformationLevel <- 1 / (((pi[2] + pi[1]) / 2) * (1 - ((pi[2] + pi[1]) / 2))) * 1 /
      (1 / (actualNumberOfSubjects * (allocationRatio[1] / sum(allocationRatio))) +
         1 / (actualNumberOfSubjects * (allocationRatio[2] / sum(allocationRatio))))
    # calculate the maximum information level
    maxInformationLevel <- 1 / (((pi[2] + pi[1]) / 2) * (1 - ((pi[2] + pi[1]) / 2))) * 1 /
      (1 / (maxNumberOfSubjects * (allocationRatio[1] / sum(allocationRatio))) +
         1 / (maxNumberOfSubjects * (allocationRatio[2] / sum(allocationRatio))))

    if (alternative == "less") {
      testStatistic <-
        (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
           sqrt(maxInformationLevel - actualInformationLevel) *
           qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(maxInformationLevel))
    } else if (alternative == "greater") {
      testStatistic <-
        (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
           sqrt(maxInformationLevel - actualInformationLevel) *
           qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(maxInformationLevel)
    } else {
      testStatistic <-
        (qnorm(p = 1 - alpha / 2, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
           sqrt(maxInformationLevel - actualInformationLevel) *
           qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(maxInformationLevel))
      testStatistic <-
        cbind(testStatistic,
              (qnorm(p = 1 - alpha / 2, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
                 sqrt(maxInformationLevel - actualInformationLevel) *
                 qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(maxInformationLevel),
              deparse.level = 0)
      testStatistic <- seq(from = min(testStatistic), to = max(testStatistic), length.out = 1e+06)

      predictivePower <- getPredictivePowerRates(
        testStatistic = testStatistic,
        alternative = alternative,
        alpha = alpha,
        actualNumberOfSubjects = actualNumberOfSubjects,
        maxNumberOfSubjects = maxNumberOfSubjects,
        allocationRatio = allocationRatio,
        pi = pi)

      testStatistic <- testStatistic[which(predictivePower >= targetPower)]
      testStatistic <- c(min(testStatistic), max(testStatistic))
    }
  }

  return(testStatistic)
}

# # validation with Example 2 from PASS Chapter 202
# targetPower <- 0.16858
# targetPowerType <- "conditional"
# alternative <- "greater"
# alpha <- 0.025
# actualNumberOfSubjects <- 60
# maxNumberOfSubjects <- 120
# allocationRatio <- c(1, 1)
# pi <- c(0.6, 0.7)
#
# getTestStatisticRates(
#   targetPower = targetPower,
#   targetPowerType = targetPowerType,
#   alternative = alternative,
#   alpha = alpha,
#   actualNumberOfSubjects = actualNumberOfSubjects,
#   maxNumberOfSubjects = maxNumberOfSubjects,
#   allocationRatio = allocationRatio,
#   pi = pi)
#
# targetPower <- 0.29262
# targetPowerType <- "predictivePower"
# alternative <- "greater"
# alpha <- 0.025
# actualNumberOfSubjects <- 60
# maxNumberOfSubjects <- 120
# allocationRatio <- c(1, 1)
# pi <- c(0.6, 0.7)
#
# getTestStatisticRates(
#   targetPower = targetPower,
#   targetPowerType = targetPowerType,
#   alternative = alternative,
#   alpha = alpha,
#   actualNumberOfSubjects = actualNumberOfSubjects,
#   maxNumberOfSubjects = maxNumberOfSubjects,
#   allocationRatio = allocationRatio,
#   pi = pi)

########################################

#' return the test statistic given the target power and maximum sample size for testing the hazard ratio in two samples
getTestStatisticSurvival <- function(targetPower,
                                     targetPowerType = c("conditional", "predictive"),
                                     alternative = c("two.sided", "less", "greater"),
                                     alpha,
                                     actualNumberOfEvents,
                                     maxNumberOfEvents,
                                     allocationRatio,
                                     lambda,
                                     ...,
                                     hazardRatio = NULL) {
  if (targetPowerType == "conditional") {
    hazardRatio <- ifelse(is.null(hazardRatio), lambda[2] / lambda[1], hazardRatio)

    # calculate the parameter being tested by the hypothesis
    theta <- log(hazardRatio)

    # calculate the actual information level
    actualInformationLevel <- actualNumberOfEvents *
      allocationRatio[1] / sum(allocationRatio) * allocationRatio[2] / sum(allocationRatio)
    # calculate the maximum information level
    maxInformationLevel <- maxNumberOfEvents *
      allocationRatio[1] / sum(allocationRatio) * allocationRatio[2] / sum(allocationRatio)

    if (alternative == "less") {
      testStatistic <-
        (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(maxInformationLevel) +
           theta * (maxInformationLevel - actualInformationLevel) +
           sqrt(maxInformationLevel - actualInformationLevel) *
           qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(actualInformationLevel))
    } else if (alternative == "greater") {
      testStatistic <-
        (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(maxInformationLevel) -
           theta * (maxInformationLevel - actualInformationLevel) +
           sqrt(maxInformationLevel - actualInformationLevel) *
           qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(actualInformationLevel)
    } else {
      testStatistic <-
        (qnorm(p = 1 - alpha / 2, mean = 0, sd = 1) * sqrt(maxInformationLevel) +
           theta * (maxInformationLevel - actualInformationLevel) +
           sqrt(maxInformationLevel - actualInformationLevel) *
           qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(actualInformationLevel))
      testStatistic <-
        cbind(testStatistic,
              (qnorm(p = 1 - alpha / 2, mean = 0, sd = 1) * sqrt(maxInformationLevel) -
                 theta * (maxInformationLevel - actualInformationLevel) +
                 sqrt(maxInformationLevel - actualInformationLevel) *
                 qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(actualInformationLevel),
              deparse.level = 0)
      testStatistic <- seq(from = min(testStatistic), to = max(testStatistic), length.out = 1e+06)

      conditionalPower <- getConditionalPowerSurvival(
        testStatistic = testStatistic,
        alternative = alternative,
        alpha = alpha,
        actualNumberOfEvents = actualNumberOfEvents,
        maxNumberOfEvents = maxNumberOfEvents,
        allocationRatio = allocationRatio,
        hazardRatio = hazardRatio)

      testStatistic <- testStatistic[which(conditionalPower >= targetPower)]
      testStatistic <- c(min(testStatistic), max(testStatistic))
    }
  } else {
    # calculate the actual information level
    actualInformationLevel <- actualNumberOfEvents *
      allocationRatio[1] / sum(allocationRatio) * allocationRatio[2] / sum(allocationRatio)
    # calculate the maximum information level
    maxInformationLevel <- maxNumberOfEvents *
      allocationRatio[1] / sum(allocationRatio) * allocationRatio[2] / sum(allocationRatio)

    if (alternative == "less") {
      testStatistic <-
        (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
           sqrt(maxInformationLevel - actualInformationLevel) *
           qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(maxInformationLevel))
    } else if (alternative == "greater") {
      testStatistic <-
        (qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
           sqrt(maxInformationLevel - actualInformationLevel) *
           qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(maxInformationLevel)
    } else {
      testStatistic <-
        (qnorm(p = 1 - alpha / 2, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
           sqrt(maxInformationLevel - actualInformationLevel) *
           qnorm(p = targetPower, mean = 0, sd = 1)) / (-sqrt(maxInformationLevel))
      testStatistic <-
        cbind(testStatistic,
              (qnorm(p = 1 - alpha / 2, mean = 0, sd = 1) * sqrt(actualInformationLevel) +
                 sqrt(maxInformationLevel - actualInformationLevel) *
                 qnorm(p = targetPower, mean = 0, sd = 1)) / sqrt(maxInformationLevel),
              deparse.level = 0)
      testStatistic <- seq(from = min(testStatistic), to = max(testStatistic), length.out = 1e+06)

      predictivePower <- getPredictivePowerSurvival(
        testStatistic = testStatistic,
        alternative = alternative,
        alpha = alpha,
        actualNumberOfEvents = actualNumberOfEvents,
        maxNumberOfEvents = maxNumberOfEvents,
        allocationRatio = allocationRatio)

      testStatistic <- testStatistic[which(predictivePower >= targetPower)]
      testStatistic <- c(min(testStatistic), max(testStatistic))
    }
  }

  return(testStatistic)
}

# # validation with Example 2 from PASS Chapter 701
# targetPower <- 0.63454
# targetPowerType <- "conditional"
# alternative <- "less"
# alpha <- 0.025
# actualNumberOfEvents <- 100
# maxNumberOfEvents <- 200
# allocationRatio <- c(1, 1)
# hazardRatio <- 0.8
#
# getTestStatisticSurvival(
#   targetPower = targetPower,
#   targetPowerType = targetPowerType,
#   alternative = alternative,
#   alpha = alpha,
#   actualNumberOfEvents = actualNumberOfEvents,
#   maxNumberOfEvents = maxNumberOfEvents,
#   allocationRatio = allocationRatio,
#   hazardRatio = hazardRatio)
#
# targetPower <- 0.80743
# targetPowerType <- "predictivePower"
# alternative <- "less"
# alpha <- 0.025
# actualNumberOfEvents <- 100
# maxNumberOfEvents <- 200
# allocationRatio <- c(1, 1)
# hazardRatio <- 0.8
#
# getTestStatisticSurvival(
#   targetPower = targetPower,
#   targetPowerType = targetPowerType,
#   alternative = alternative,
#   alpha = alpha,
#   actualNumberOfEvents = actualNumberOfEvents,
#   maxNumberOfEvents = maxNumberOfEvents,
#   allocationRatio = allocationRatio,
#   hazardRatio = hazardRatio)

################################################################################
