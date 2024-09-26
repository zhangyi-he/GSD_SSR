#' @title Get the predictive power
#' @author Zhangyi He

#' version 0.1

################################################################################

#' return the predictive power for testing means in two samples
getPredictivePowerMeans <- function(testStatistic,
                                    alternative = c("two.sided", "less", "greater"),
                                    alpha,
                                    actualNumberOfSubjects,
                                    maxNumberOfSubjects,
                                    allocationRatio,
                                    sigma) {
  # calculate the actual information level
  actualInformationLevel <- 1 /
    (sigma[1]^2 / (actualNumberOfSubjects * (allocationRatio[1] / sum(allocationRatio))) +
       sigma[2]^2 / (actualNumberOfSubjects * (allocationRatio[2] / sum(allocationRatio))))
  # calculate the maximum information level
  maxInformationLevel <- 1 /
    (sigma[1]^2 / (maxNumberOfSubjects * (allocationRatio[1] / sum(allocationRatio))) +
       sigma[2]^2 / (maxNumberOfSubjects * (allocationRatio[2] / sum(allocationRatio))))

  if (alternative == "less") {
    predictivePower <-
      pnorm(q = (-testStatistic * sqrt(maxInformationLevel) -
                   qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel)) /
              sqrt(maxInformationLevel - actualInformationLevel),
            mean = 0, sd = 1)
  } else if (alternative == "greater") {
    predictivePower <-
      pnorm(q = (testStatistic * sqrt(maxInformationLevel) -
                   qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel)) /
              sqrt(maxInformationLevel - actualInformationLevel),
            mean = 0, sd = 1)
  } else {
    predictivePower <-
      pnorm(q = (-abs(testStatistic) * sqrt(maxInformationLevel) -
                   qnorm(1 - alpha / 2) * sqrt(actualInformationLevel)) /
              sqrt(maxInformationLevel - actualInformationLevel),
            mean = 0, sd = 1) +
      pnorm(q = (abs(testStatistic) * sqrt(maxInformationLevel) -
                   qnorm(1 - alpha / 2) * sqrt(actualInformationLevel)) /
              sqrt(maxInformationLevel - actualInformationLevel),
            mean = 0, sd = 1)
  }

  return(predictivePower)
}

# # validation with Example 2 from PASS Chapter 433
# testStatistic <- 2.12
# alternative <- "two.sided"
# alpha <- 0.05
# actualNumberOfSubjects <- 60
# maxNumberOfSubjects <- 120
# allocationRatio <- c(1, 1)
# mu <- c(0, 0.5)
# sigma <- c(4, 4)
#
# getPredictivePowerMeans(
#   testStatistic = testStatistic,
#   alternative = alternative,
#   alpha = alpha,
#   actualNumberOfSubjects = actualNumberOfSubjects,
#   maxNumberOfSubjects = maxNumberOfSubjects,
#   allocationRatio = allocationRatio,
#   sigma = sigma)

########################################

#' return the predictive power for testing rates in two samples
getPredictivePowerRates <- function(testStatistic,
                                    alternative = c("two.sided", "less", "greater"),
                                    alpha,
                                    actualNumberOfSubjects,
                                    maxNumberOfSubjects,
                                    allocationRatio,
                                    pi) {
  # calculate the actual information level
  actualInformationLevel <- 1 / (((pi[2] + pi[1]) / 2) * (1 - ((pi[2] + pi[1]) / 2))) * 1 /
    (1 / (actualNumberOfSubjects * (allocationRatio[1] / sum(allocationRatio))) +
       1 / (actualNumberOfSubjects * (allocationRatio[2] / sum(allocationRatio))))
  # calculate the maximum information level
  maxInformationLevel <- 1 / (((pi[2] + pi[1]) / 2) * (1 - ((pi[2] + pi[1]) / 2))) * 1 /
    (1 / (maxNumberOfSubjects * (allocationRatio[1] / sum(allocationRatio))) +
       1 / (maxNumberOfSubjects * (allocationRatio[2] / sum(allocationRatio))))

  if (alternative == "less") {
    predictivePower <-
      pnorm(q = (-testStatistic * sqrt(maxInformationLevel) -
                   qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel)) /
              sqrt(maxInformationLevel - actualInformationLevel),
            mean = 0, sd = 1)
  } else if (alternative == "greater") {
    predictivePower <-
      pnorm(q = (testStatistic * sqrt(maxInformationLevel) -
                   qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel)) /
              sqrt(maxInformationLevel - actualInformationLevel),
            mean = 0, sd = 1)
  } else {
    predictivePower <-
      pnorm(q = (-abs(testStatistic) * sqrt(maxInformationLevel) -
                   qnorm(1 - alpha / 2) * sqrt(actualInformationLevel)) /
              sqrt(maxInformationLevel - actualInformationLevel),
            mean = 0, sd = 1) +
      pnorm(q = (abs(testStatistic) * sqrt(maxInformationLevel) -
                   qnorm(1 - alpha / 2) * sqrt(actualInformationLevel)) /
              sqrt(maxInformationLevel - actualInformationLevel),
            mean = 0, sd = 1)
  }

  return(predictivePower)
}

# # validation with Example 1 from PASS Chapter 202
# testStatistic <- c(0, 0.5, 1, 1.5, 2, 2.5)
# alternative <- "greater"
# alpha <- 0.025
# actualNumberOfSubjects <- 60
# maxNumberOfSubjects <- 120
# allocationRatio <- c(1, 1)
# pi <- c(0.6, 0.7)
#
# getPredictivePowerRates(
#   testStatistic = testStatistic,
#   alternative = alternative,
#   alpha = alpha,
#   actualNumberOfSubjects = actualNumberOfSubjects,
#   maxNumberOfSubjects = maxNumberOfSubjects,
#   allocationRatio = allocationRatio,
#   pi = pi)

########################################

#' return the predictive power for testing the hazard ratio in two samples
getPredictivePowerSurvival <- function(testStatistic,
                                       alternative = c("two.sided", "less", "greater"),
                                       alpha,
                                       actualNumberOfEvents,
                                       maxNumberOfEvents,
                                       allocationRatio) {
  # calculate the actual information level
  actualInformationLevel <- actualNumberOfEvents *
    allocationRatio[1] / sum(allocationRatio) * allocationRatio[2] / sum(allocationRatio)
  # calculate the maximum information level
  maxInformationLevel <- maxNumberOfEvents *
    allocationRatio[1] / sum(allocationRatio) * allocationRatio[2] / sum(allocationRatio)

  if (alternative == "less") {
    predictivePower <-
      pnorm(q = (-testStatistic * sqrt(maxInformationLevel) -
                   qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel)) /
              sqrt(maxInformationLevel - actualInformationLevel),
            mean = 0, sd = 1)
  } else if (alternative == "greater") {
    predictivePower <-
      pnorm(q = (testStatistic * sqrt(maxInformationLevel) -
                   qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(actualInformationLevel)) /
              sqrt(maxInformationLevel - actualInformationLevel),
            mean = 0, sd = 1)
  } else {
    predictivePower <-
      pnorm(q = (-abs(testStatistic) * sqrt(maxInformationLevel) -
                   qnorm(1 - alpha / 2) * sqrt(actualInformationLevel)) /
              sqrt(maxInformationLevel - actualInformationLevel),
            mean = 0, sd = 1) +
      pnorm(q = (abs(testStatistic) * sqrt(maxInformationLevel) -
                   qnorm(1 - alpha / 2) * sqrt(actualInformationLevel)) /
              sqrt(maxInformationLevel - actualInformationLevel),
            mean = 0, sd = 1)
  }

  return(predictivePower)
}

# # validation with Example 1 from PASS Chapter 701
# testStatistic <- c(-3, -2.5, -2, -1.5, -1)
# alternative <- "less"
# alpha <- 0.025
# actualNumberOfEvents <- 100
# maxNumberOfEvents <- 200
# allocationRatio <- c(1, 1)
#
# getPredictivePowerSurvival(
#   testStatistic = testStatistic,
#   alternative = alternative,
#   alpha = alpha,
#   actualNumberOfEvents = actualNumberOfEvents,
#   maxNumberOfEvents = maxNumberOfEvents,
#   allocationRatio = allocationRatio)

################################################################################
