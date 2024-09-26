#' @title Get the conditional power
#' @author Zhangyi He

#' version 0.1

################################################################################

#' return the conditional power for testing means in two samples
getConditionalPowerMeans <- function(testStatistic,
                                     alternative = c("two.sided", "less", "greater"),
                                     alpha,
                                     actualNumberOfSubjects,
                                     maxNumberOfSubjects,
                                     allocationRatio,
                                     mu,
                                     sigma) {
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
    conditionalPower <-
      pnorm(q = (-testStatistic * sqrt(actualInformationLevel) -
                   qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(maxInformationLevel) -
                   theta * (maxInformationLevel - actualInformationLevel)) /
              sqrt(maxInformationLevel - actualInformationLevel),
            mean = 0, sd = 1)
  } else if (alternative == "greater") {
    conditionalPower <-
      pnorm(q = (testStatistic * sqrt(actualInformationLevel) -
                   qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(maxInformationLevel) +
                   theta * (maxInformationLevel - actualInformationLevel)) /
              sqrt(maxInformationLevel - actualInformationLevel),
            mean = 0, sd = 1)
  } else {
    conditionalPower <-
      pnorm(q = (-testStatistic * sqrt(actualInformationLevel) -
                   qnorm(1 - alpha / 2) * sqrt(maxInformationLevel) -
                   theta * (maxInformationLevel - actualInformationLevel)) /
              sqrt(maxInformationLevel - actualInformationLevel),
            mean = 0, sd = 1) +
      pnorm(q = (testStatistic * sqrt(actualInformationLevel) -
                   qnorm(1 - alpha / 2) * sqrt(maxInformationLevel) +
                   theta * (maxInformationLevel - actualInformationLevel)) /
              sqrt(maxInformationLevel - actualInformationLevel),
            mean = 0, sd = 1)
  }

  return(conditionalPower)
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
# getConditionalPowerMeans(
#   testStatistic = testStatistic,
#   alternative = alternative,
#   alpha = alpha,
#   actualNumberOfSubjects = actualNumberOfSubjects,
#   maxNumberOfSubjects = maxNumberOfSubjects,
#   allocationRatio = allocationRatio,
#   mu = mu,
#   sigma = sigma)

########################################

#' return the conditional power for testing rates in two samples
getConditionalPowerRates <- function(testStatistic,
                                     alternative = c("two.sided", "less", "greater"),
                                     alpha,
                                     actualNumberOfSubjects,
                                     maxNumberOfSubjects,
                                     allocationRatio,
                                     pi) {
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
    conditionalPower <-
      pnorm(q = (-testStatistic * sqrt(actualInformationLevel) -
                   qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(maxInformationLevel) -
                   theta * (maxInformationLevel - actualInformationLevel)) /
              sqrt(maxInformationLevel - actualInformationLevel),
            mean = 0, sd = 1)
  } else if (alternative == "greater") {
    conditionalPower <-
      pnorm(q = (testStatistic * sqrt(actualInformationLevel) -
                   qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(maxInformationLevel) +
                   theta * (maxInformationLevel - actualInformationLevel)) /
              sqrt(maxInformationLevel - actualInformationLevel),
            mean = 0, sd = 1)
  } else {
    conditionalPower <-
      pnorm(q = (-testStatistic * sqrt(actualInformationLevel) -
                   qnorm(1 - alpha / 2) * sqrt(maxInformationLevel) -
                   theta * (maxInformationLevel - actualInformationLevel)) /
              sqrt(maxInformationLevel - actualInformationLevel),
            mean = 0, sd = 1) +
      pnorm(q = (testStatistic * sqrt(actualInformationLevel) -
                   qnorm(1 - alpha / 2) * sqrt(maxInformationLevel) +
                   theta * (maxInformationLevel - actualInformationLevel)) /
              sqrt(maxInformationLevel - actualInformationLevel),
            mean = 0, sd = 1)
  }

  return(conditionalPower)
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
# getConditionalPowerRates(
#   testStatistic = testStatistic,
#   alternative = alternative,
#   alpha = alpha,
#   actualNumberOfSubjects = actualNumberOfSubjects,
#   maxNumberOfSubjects = maxNumberOfSubjects,
#   allocationRatio = allocationRatio,
#   pi = pi)

########################################

#' return the conditional power for testing the hazard ratio in two samples
getConditionalPowerSurvival <- function(testStatistic,
                                        alternative = c("two.sided", "less", "greater"),
                                        alpha,
                                        actualNumberOfEvents,
                                        maxNumberOfEvents,
                                        allocationRatio,
                                        lambda,
                                        ...,
                                        hazardRatio = NULL) {
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
    conditionalPower <-
      pnorm(q = (-testStatistic * sqrt(actualInformationLevel) -
                   qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(maxInformationLevel) -
                   theta * (maxInformationLevel - actualInformationLevel)) /
              sqrt(maxInformationLevel - actualInformationLevel),
            mean = 0, sd = 1)
  } else if (alternative == "greater") {
    conditionalPower <-
      pnorm(q = (testStatistic * sqrt(actualInformationLevel) -
                   qnorm(p = 1 - alpha, mean = 0, sd = 1) * sqrt(maxInformationLevel) +
                   theta * (maxInformationLevel - actualInformationLevel)) /
              sqrt(maxInformationLevel - actualInformationLevel),
            mean = 0, sd = 1)
  } else {
    conditionalPower <-
      pnorm(q = (-testStatistic * sqrt(actualInformationLevel) -
                   qnorm(1 - alpha / 2) * sqrt(maxInformationLevel) -
                   theta * (maxInformationLevel - actualInformationLevel)) /
              sqrt(maxInformationLevel - actualInformationLevel),
            mean = 0, sd = 1) +
      pnorm(q = (testStatistic * sqrt(actualInformationLevel) -
                   qnorm(1 - alpha / 2) * sqrt(maxInformationLevel) +
                   theta * (maxInformationLevel - actualInformationLevel)) /
              sqrt(maxInformationLevel - actualInformationLevel),
            mean = 0, sd = 1)
  }

  return(conditionalPower)
}

# # validation with Example 1 from PASS Chapter 701
# testStatistic <- c(-3, -2.5, -2, -1.5, -1)
# alternative <- "less"
# alpha <- 0.025
# actualNumberOfEvents <- 100
# maxNumberOfEvents <- 200
# allocationRatio <- c(1, 1)
# hazardRatio <- 0.8
#
# getConditionalPowerSurvival(
#   testStatistic = testStatistic,
#   alternative = alternative,
#   alpha = alpha,
#   actualNumberOfEvents = actualNumberOfEvents,
#   maxNumberOfEvents = maxNumberOfEvents,
#   allocationRatio = allocationRatio,
#   hazardRatio = hazardRatio)

################################################################################
