#' @title Get conditional power
#' @author Zhangyi He

################################################################################

#' return the conditional power for testing means in two samples
getConditionalPowerMeans <- function(testStatistic,
                                     sided,
                                     alpha,
                                     actualNumberOfSubjects,
                                     targetNumberOfSubjects,
                                     allocationRatioPlanned,
                                     mu,
                                     sigma,
                                     ...,
                                     alternative = NULL) {
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
      conditionalPower <- 
        pnorm(q = (-testStatistic * sqrt(actualInformationLevel) -
                     qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(targetInformationLevel) -
                     theta * (targetInformationLevel - actualInformationLevel)) /
                sqrt(targetInformationLevel - actualInformationLevel),
              mean = 0, sd = 1)
      conditionalPower <- 
        cbind(conditionalPower,
              pnorm(q = (testStatistic * sqrt(actualInformationLevel) -
                           qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(targetInformationLevel) +
                           theta * (targetInformationLevel - actualInformationLevel)) /
                      sqrt(targetInformationLevel - actualInformationLevel),
                    mean = 0, sd = 1), 
              deparse.level = 0)
    } else {
      if (alternative == "less") {
        conditionalPower <- 
          pnorm(q = (-testStatistic * sqrt(actualInformationLevel) -
                       qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(targetInformationLevel) -
                       theta * (targetInformationLevel - actualInformationLevel)) /
                  sqrt(targetInformationLevel - actualInformationLevel),
                mean = 0, sd = 1)
      } else {
        conditionalPower <- 
          pnorm(q = (testStatistic * sqrt(actualInformationLevel) -
                       qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(targetInformationLevel) +
                       theta * (targetInformationLevel - actualInformationLevel)) /
                  sqrt(targetInformationLevel - actualInformationLevel),
                mean = 0, sd = 1)
      }
    }
  } else {
    conditionalPower <-
      pnorm(q = (-testStatistic * sqrt(actualInformationLevel) -
                   qnorm(1 - alpha / 2) * sqrt(targetInformationLevel) -
                   theta * (targetInformationLevel - actualInformationLevel)) /
              sqrt(targetInformationLevel - actualInformationLevel),
            mean = 0, sd = 1) +
      pnorm(q = (testStatistic * sqrt(actualInformationLevel) -
                   qnorm(1 - alpha / 2) * sqrt(targetInformationLevel) +
                   theta * (targetInformationLevel - actualInformationLevel)) /
              sqrt(targetInformationLevel - actualInformationLevel),
            mean = 0, sd = 1)
  }

  return(conditionalPower)
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
# getConditionalPowerMeans(
#   testStatistic = testStatistic,
#   sided = sided,
#   alpha = alpha,
#   actualNumberOfSubjects = actualNumberOfSubjects,
#   targetNumberOfSubjects = targetNumberOfSubjects,
#   allocationRatioPlanned = allocationRatioPlanned,
#   mu = mu,
#   sigma = sigma)

########################################

#' return the conditional power for testing rates in two samples
getConditionalPowerRates <- function(testStatistic,
                                     sided,
                                     alpha,
                                     actualNumberOfSubjects,
                                     targetNumberOfSubjects,
                                     allocationRatioPlanned,
                                     rate,
                                     ...,
                                     alternative = NULL) {
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
      conditionalPower <- 
        pnorm(q = (-testStatistic * sqrt(actualInformationLevel) -
                     qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(targetInformationLevel) -
                     theta * (targetInformationLevel - actualInformationLevel)) /
                sqrt(targetInformationLevel - actualInformationLevel),
              mean = 0, sd = 1)
      conditionalPower <- 
        cbind(conditionalPower,
              pnorm(q = (testStatistic * sqrt(actualInformationLevel) -
                           qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(targetInformationLevel) +
                           theta * (targetInformationLevel - actualInformationLevel)) /
                      sqrt(targetInformationLevel - actualInformationLevel),
                    mean = 0, sd = 1), 
              deparse.level = 0)
    } else {
      if (alternative == "less") {
        conditionalPower <- 
          pnorm(q = (-testStatistic * sqrt(actualInformationLevel) -
                       qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(targetInformationLevel) -
                       theta * (targetInformationLevel - actualInformationLevel)) /
                  sqrt(targetInformationLevel - actualInformationLevel),
                mean = 0, sd = 1)
      } else {
        conditionalPower <- 
          pnorm(q = (testStatistic * sqrt(actualInformationLevel) -
                       qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(targetInformationLevel) +
                       theta * (targetInformationLevel - actualInformationLevel)) /
                  sqrt(targetInformationLevel - actualInformationLevel),
                mean = 0, sd = 1)
      }
    }
  } else {
    conditionalPower <-
      pnorm(q = (-testStatistic * sqrt(actualInformationLevel) -
                   qnorm(1 - alpha / 2) * sqrt(targetInformationLevel) -
                   theta * (targetInformationLevel - actualInformationLevel)) /
              sqrt(targetInformationLevel - actualInformationLevel),
            mean = 0, sd = 1) +
      pnorm(q = (testStatistic * sqrt(actualInformationLevel) -
                   qnorm(1 - alpha / 2) * sqrt(targetInformationLevel) +
                   theta * (targetInformationLevel - actualInformationLevel)) /
              sqrt(targetInformationLevel - actualInformationLevel),
            mean = 0, sd = 1)
  }

  return(conditionalPower)
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
# getConditionalPowerRates(
#   testStatistic = testStatistic,
#   sided = sided,
#   alpha = alpha,
#   actualNumberOfSubjects = actualNumberOfSubjects,
#   targetNumberOfSubjects = targetNumberOfSubjects,
#   allocationRatioPlanned = allocationRatioPlanned,
#   rate = rate,
#   alternative = alternative)

########################################

#' return the conditional power for testing the hazard ratio in two samples
getConditionalPowerSurvival <- function(testStatistic,
                                        sided,
                                        alpha,
                                        targetNumberOfEvents,
                                        actualNumberOfEvents,
                                        allocationRatioPlanned,
                                        hazardRate,
                                        ...,
                                        hazardRatio = NULL,
                                        alternative = NULL) {
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
      conditionalPower <- 
        pnorm(q = (-testStatistic * sqrt(actualInformationLevel) -
                     qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(targetInformationLevel) -
                     theta * (targetInformationLevel - actualInformationLevel)) /
                sqrt(targetInformationLevel - actualInformationLevel),
              mean = 0, sd = 1)
      conditionalPower <- 
        cbind(conditionalPower,
              pnorm(q = (testStatistic * sqrt(actualInformationLevel) -
                           qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(targetInformationLevel) +
                           theta * (targetInformationLevel - actualInformationLevel)) /
                      sqrt(targetInformationLevel - actualInformationLevel),
                    mean = 0, sd = 1), 
              deparse.level = 0)
    } else {
      if (alternative == "less") {
        conditionalPower <- 
          pnorm(q = (-testStatistic * sqrt(actualInformationLevel) -
                       qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(targetInformationLevel) -
                       theta * (targetInformationLevel - actualInformationLevel)) /
                  sqrt(targetInformationLevel - actualInformationLevel),
                mean = 0, sd = 1)
      } else {
        conditionalPower <- 
          pnorm(q = (testStatistic * sqrt(actualInformationLevel) -
                       qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(targetInformationLevel) +
                       theta * (targetInformationLevel - actualInformationLevel)) /
                  sqrt(targetInformationLevel - actualInformationLevel),
                mean = 0, sd = 1)
      }
    }
  } else {
    conditionalPower <-
      pnorm(q = (-testStatistic * sqrt(actualInformationLevel) -
                   qnorm(1 - alpha / 2) * sqrt(targetInformationLevel) -
                   theta * (targetInformationLevel - actualInformationLevel)) /
              sqrt(targetInformationLevel - actualInformationLevel),
            mean = 0, sd = 1) +
      pnorm(q = (testStatistic * sqrt(actualInformationLevel) -
                   qnorm(1 - alpha / 2) * sqrt(targetInformationLevel) +
                   theta * (targetInformationLevel - actualInformationLevel)) /
              sqrt(targetInformationLevel - actualInformationLevel),
            mean = 0, sd = 1)
  }

  return(conditionalPower)
}

# # validation with Example 1 from PASS Chapter 701
# testStatistic <- c(-3, -2.5, -2, -1.5, -1)
# sided <- 1
# alpha <- 0.025
# targetNumberOfEvents <- 200
# actualNumberOfEvents <- 100
# allocationRatioPlanned <- 1
# hazardRatio <- 0.8
# alternative <- "less"
# 
# getConditionalPowerSurvival(
#   testStatistic = testStatistic,
#   sided = sided,
#   alpha = alpha,
#   targetNumberOfEvents = targetNumberOfEvents,
#   actualNumberOfEvents = actualNumberOfEvents,
#   allocationRatioPlanned = allocationRatioPlanned,
#   hazardRatio = hazardRatio,
#   alternative = alternative)

################################################################################
