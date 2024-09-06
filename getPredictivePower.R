#' @title Get predictive power
#' @author Zhangyi He

################################################################################

#' return the predictive power for testing means in two samples
getPredictivePowerMeans <- function(testStatistic,
                                    sided,
                                    alpha,
                                    actualNumberOfSubjects,
                                    targetNumberOfSubjects,
                                    allocationRatioPlanned,
                                    sigma,
                                    ...,
                                    alternative = NULL) {
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
      predictivePower <- 
        pnorm(q = (-testStatistic * sqrt(targetInformationLevel) -
                     qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel)) /
                sqrt(targetInformationLevel - actualInformationLevel),
              mean = 0, sd = 1)
      predictivePower <- 
        cbind(predictivePower,
              pnorm(q = (testStatistic * sqrt(targetInformationLevel) -
                           qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel)) /
                      sqrt(targetInformationLevel - actualInformationLevel),
                    mean = 0, sd = 1), 
              deparse.level = 0)
    } else {
      if (alternative == "less") {
        predictivePower <- 
          pnorm(q = (-testStatistic * sqrt(targetInformationLevel) -
                       qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel)) /
                  sqrt(targetInformationLevel - actualInformationLevel),
                mean = 0, sd = 1)
      } else {
        predictivePower <- 
          pnorm(q = (testStatistic * sqrt(targetInformationLevel) -
                       qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel)) /
                  sqrt(targetInformationLevel - actualInformationLevel),
                mean = 0, sd = 1)
      }
    }
  } else {
    predictivePower <-
      pnorm(q = (-abs(testStatistic) * sqrt(targetInformationLevel) -
                   qnorm(1 - alpha / 2) * sqrt(actualInformationLevel)) /
              sqrt(targetInformationLevel - actualInformationLevel),
            mean = 0, sd = 1) +
      pnorm(q = (abs(testStatistic) * sqrt(targetInformationLevel) -
                   qnorm(1 - alpha / 2) * sqrt(actualInformationLevel)) /
              sqrt(targetInformationLevel - actualInformationLevel),
            mean = 0, sd = 1)
  }

  return(predictivePower)
}

# # validation with Example 2 from PASS Chapter 433
# testStatistic <- 2.12
# sided <- 2
# alpha <- 0.05
# actualNumberOfSubjects <- 60
# targetNumberOfSubjects <- 120
# allocationRatioPlanned <- 1
# mu <- c(0, 0.5)
# sigma <- c(4, 4)
# 
# getPredictivePowerMeans(
#   testStatistic = testStatistic,
#   sided = sided,
#   alpha = alpha,
#   actualNumberOfSubjects = actualNumberOfSubjects,
#   targetNumberOfSubjects = targetNumberOfSubjects,
#   allocationRatioPlanned = allocationRatioPlanned,
#   sigma = sigma)

########################################

#' return the predictive power for testing rates in two samples
getPredictivePowerRates <- function(testStatistic,
                                    sided,
                                    alpha,
                                    actualNumberOfSubjects,
                                    targetNumberOfSubjects,
                                    allocationRatioPlanned,
                                    rate,
                                    ...,
                                    alternative = NULL) {
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
      predictivePower <- 
        pnorm(q = (-testStatistic * sqrt(targetInformationLevel) -
                     qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel)) /
                sqrt(targetInformationLevel - actualInformationLevel),
              mean = 0, sd = 1)
      predictivePower <- 
        cbind(predictivePower,
              pnorm(q = (testStatistic * sqrt(targetInformationLevel) -
                           qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel)) /
                      sqrt(targetInformationLevel - actualInformationLevel),
                    mean = 0, sd = 1), 
              deparse.level = 0)
    } else {
      if (alternative == "less") {
        predictivePower <- 
          pnorm(q = (-testStatistic * sqrt(targetInformationLevel) -
                       qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel)) /
                  sqrt(targetInformationLevel - actualInformationLevel),
                mean = 0, sd = 1)
      } else {
        predictivePower <- 
          pnorm(q = (testStatistic * sqrt(targetInformationLevel) -
                       qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel)) /
                  sqrt(targetInformationLevel - actualInformationLevel),
                mean = 0, sd = 1)
      }
    }
  } else {
    predictivePower <-
      pnorm(q = (-abs(testStatistic) * sqrt(targetInformationLevel) -
                   qnorm(1 - alpha / 2) * sqrt(actualInformationLevel)) /
              sqrt(targetInformationLevel - actualInformationLevel),
            mean = 0, sd = 1) +
      pnorm(q = (abs(testStatistic) * sqrt(targetInformationLevel) -
                   qnorm(1 - alpha / 2) * sqrt(actualInformationLevel)) /
              sqrt(targetInformationLevel - actualInformationLevel),
            mean = 0, sd = 1)
  }

  return(predictivePower)
}

# # validation with Example 1 from PASS Chapter 202
# testStatistic <- c(0, 0.5, 1, 1.5, 2, 2.5)
# sided <- 1
# alpha <- 0.025
# actualNumberOfSubjects <- 60
# targetNumberOfSubjects <- 120
# allocationRatioPlanned <- 1
# rate <- c(0.6, 0.7)
# alternative <- "greater"
# 
# getPredictivePowerRates(
#   testStatistic = testStatistic,
#   sided = sided,
#   alpha = alpha,
#   actualNumberOfSubjects = actualNumberOfSubjects,
#   targetNumberOfSubjects = targetNumberOfSubjects,
#   allocationRatioPlanned = allocationRatioPlanned,
#   rate = rate,
#   alternative = alternative)

########################################

#' return the predictive power for testing the hazard ratio in two samples
getPredictivePowerSurvival <- function(testStatistic,
                                       sided,
                                       alpha,
                                       targetNumberOfEvents,
                                       actualNumberOfEvents,
                                       allocationRatioPlanned,
                                       ...,
                                       alternative = NULL) {
  # calculate the actual information level
  actualInformationLevel <- actualNumberOfEvents *
    allocationRatioPlanned / (allocationRatioPlanned + 1) * 1 / (allocationRatioPlanned + 1)
  # calculate the target information level
  targetInformationLevel <- targetNumberOfEvents *
    allocationRatioPlanned / (allocationRatioPlanned + 1) * 1 / (allocationRatioPlanned + 1)

  if (sided == 1) {
    if (is.null(alternative)) {
      predictivePower <- 
        pnorm(q = (-testStatistic * sqrt(targetInformationLevel) -
                     qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel)) /
                sqrt(targetInformationLevel - actualInformationLevel),
              mean = 0, sd = 1)
      predictivePower <- 
        cbind(predictivePower,
              pnorm(q = (testStatistic * sqrt(targetInformationLevel) -
                           qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel)) /
                      sqrt(targetInformationLevel - actualInformationLevel),
                    mean = 0, sd = 1), 
              deparse.level = 0)
    } else {
      if (alternative == "less") {
        predictivePower <- 
          pnorm(q = (-testStatistic * sqrt(targetInformationLevel) -
                       qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel)) /
                  sqrt(targetInformationLevel - actualInformationLevel),
                mean = 0, sd = 1)
      } else {
        predictivePower <- 
          pnorm(q = (testStatistic * sqrt(targetInformationLevel) -
                       qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel)) /
                  sqrt(targetInformationLevel - actualInformationLevel),
                mean = 0, sd = 1)
      }
    }
  } else {
    predictivePower <-
      pnorm(q = (-abs(testStatistic) * sqrt(targetInformationLevel) -
                   qnorm(1 - alpha / 2) * sqrt(actualInformationLevel)) /
              sqrt(targetInformationLevel - actualInformationLevel),
            mean = 0, sd = 1) +
      pnorm(q = (abs(testStatistic) * sqrt(targetInformationLevel) -
                   qnorm(1 - alpha / 2) * sqrt(actualInformationLevel)) /
              sqrt(targetInformationLevel - actualInformationLevel),
            mean = 0, sd = 1)
  }

  return(predictivePower)
}

# # validation with Example 1 from PASS Chapter 701
# testStatistic <- c(-3, -2.5, -2, -1.5, -1)
# sided <- 1
# alpha <- 0.025
# targetNumberOfEvents <- 200
# actualNumberOfEvents <- 100
# allocationRatioPlanned <- 1
# alternative <- "less"
# 
# getPredictivePowerSurvival(
#   testStatistic = testStatistic,
#   sided = sided,
#   alpha = alpha,
#   targetNumberOfEvents = targetNumberOfEvents,
#   actualNumberOfEvents = actualNumberOfEvents,
#   allocationRatioPlanned = allocationRatioPlanned,
#   alternative = alternative)

################################################################################
