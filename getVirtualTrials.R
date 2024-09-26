#' @title Get virtual trials
#' @author Zhangyi He

#' version 0.1

# install and load required packages
# if (!require("randomizeR")) {
#   install.packages("randomizeR")
#   library("randomizeR")
# }
# if (!require("DescTools")) {
#   install.packages("DescTools")
#   library("DescTools")
# }

################################################################################

#' return virtual trials with subject information (except for responses) and simulation settings
getSubjectInformation <- function(numberOfTrials,
                                  numberOfSubjects,
                                  centreLabels,
                                  treatmentLabels,
                                  allocationRatio,
                                  blockSizes,
                                  followupTime,
                                  dropoutRates,
                                  accrualTimes,
                                  accrualIntensities,
                                  seed) {
  set.seed(seed = ifelse(is.null(seed), 21, seed))

  numberOfCentres <- length(centreLabels)

  # simulate enrolment
  enrolmentStatus <- vector(mode = "list", length = numberOfCentres)
  cat("running enrolment... \n")
  progressBar <- txtProgressBar(min = 1, max = numberOfCentres, style = 3)
  for (j in 1:numberOfCentres) {
    setTxtProgressBar(pb = progressBar, value = j)
    enrolmentStatus[[j]] <-
      matrix(data = rpois(n = (accrualTimes[j, 2] - accrualTimes[j, 1]) *
                            numberOfTrials, lambda = accrualIntensities[j]),
             nrow = accrualTimes[j, 2] - accrualTimes[j, 1], ncol = numberOfTrials)
  }
  cat("\n")
  cat("Done! \n")

  # simulate randomisation
  subjectInformation <- vector(mode = "list", length = numberOfTrials)
  cat("running randomisation... \n")
  progressBar <- txtProgressBar(min = 1, max = numberOfTrials, style = 3)
  for (i in 1:numberOfTrials) {
    setTxtProgressBar(pb = progressBar, value = i)
    subjectInfo <- NULL
    for (j in 1:numberOfCentres) {
      randomisationTime <- accrualTimes[j, 1] +
        rep(x = 1:length(enrolmentStatus[[j]][, i]), times = enrolmentStatus[[j]][, i])
      observationTime <- randomisationTime + followupTime
      treatmentAllocation <- getRandList(genSeq(rpbrPar(
        N = ceiling((ceiling(length(randomisationTime) / LCM(blockSizes))) *
                      LCM(blockSizes) / sum(allocationRatio)) * sum(allocationRatio),
        rb = blockSizes,
        K = length(treatmentLabels),
        ratio = allocationRatio,
        groups = treatmentLabels), r = 1, seed = seed))[1:length(randomisationTime)]
      subjectInfo <-
        rbind(subjectInfo,
              data.frame(centre = centreLabels[j],
                         treatment = treatmentAllocation,
                         start = randomisationTime,
                         stop = observationTime))
    }
    subjectInfo <- subjectInfo[order(subjectInfo$start), ]
    subjectInfo <- subjectInfo[1:min(numberOfSubjects, nrow(subjectInfo)), ]
    subjectInfo <- cbind(id = 1:nrow(subjectInfo), subjectInfo)
    subjectInfo$stop[which(subjectInfo$treatment == treatmentLabels[1] &
                             rbinom(n = nrow(subjectInfo), size = 1, prob = dropoutRates[1]) == 1)] <- NA
    subjectInfo$stop[which(subjectInfo$treatment == treatmentLabels[2] &
                             rbinom(n = nrow(subjectInfo), size = 1, prob = dropoutRates[2]) == 1)] <- NA
    subjectInformation[[i]] <- as.data.frame(subjectInfo)
    rownames(subjectInformation[[i]]) <- NULL
    colnames(subjectInformation[[i]]) <- c("id", "centre", "treatment", "randomisationTime", "observationTime")
  }
  cat("\n")
  cat("Done! \n")

  return(list(subjectInformation = subjectInformation,
              simulationSettings = list(
                numberOfTrials = numberOfTrials,
                numberOfSubjects = numberOfSubjects,
                centreLabels = centreLabels,
                treatmentLabels = treatmentLabels,
                allocationRatio = allocationRatio,
                blockSizes = blockSizes,
                followupTime = followupTime,
                dropoutRates = dropoutRates,
                accrualTimes = accrualTimes,
                accrualIntensities = accrualIntensities,
                seed = ifelse(is.null(seed), 21, seed))))
}

# # test with the OSCAR trial (Young et al., 2013)
# numberOfTrials <- 1e+05
# numberOfSubjects <- 1006
# centreLabels <- 1:29
# treatmentLabels <- c("SOC", "HFOV")
# allocationRatio <- c(1, 1)
# blockSizes <- c(4, 6)
# followupTime <- 30
# dropoutRates <- c(0.03, 0.03)
# accrualTimes <- matrix(data = c(rep(x = 0, times = length(centreLabels)),
#                                 rep(x = 30 * 12 * 5, times = length(centreLabels))),
#                        nrow = length(centreLabels), ncol = 2)
# accrualIntensities <- rep(x = 5.5 / 7 / 29, times = length(centreLabels))
# seed <- 21
#
# system.time(trials <- getSubjectInformation(
#   numberOfTrials = numberOfTrials,
#   numberOfSubjects = numberOfSubjects,
#   centreLabels = centreLabels,
#   treatmentLabels = treatmentLabels,
#   allocationRatio = allocationRatio,
#   blockSizes = blockSizes,
#   followupTime = followupTime,
#   dropoutRates = dropoutRates,
#   accrualTimes = accrualTimes,
#   accrualIntensities = accrualIntensities,
#   seed = seed))

########################################

#' return virtual trials with subject information and simulation settings for testing means in two samples
addSubjectOutcomesMeans <- function(trials,
                                    mu,
                                    sigma) {
  set.seed(seed = trials$simulationSettings$seed)

  # update simulation settings
  trials$simulationSettings$mu <- mu
  trials$simulationSettings$sigma <- sigma

  numberOfTrials <- trials$simulationSettings$numberOfTrials
  numberOfSubjects <- trials$simulationSettings$numberOfSubjects
  treatmentLabels <- trials$simulationSettings$treatmentLabels

  # simulate outcomes
  cat("adding outcomes... \n")
  progressBar <- txtProgressBar(min = 1, max = numberOfTrials, style = 3)
  subjectOutcomes <-
    list(matrix(data = rnorm(n = numberOfSubjects * numberOfTrials,
                             mean = mu[1], sd = sigma[1]),
                nrow = numberOfSubjects, ncol = numberOfTrials),
         matrix(data = rnorm(n = numberOfSubjects * numberOfTrials,
                             mean = mu[2], sd = sigma[2]),
                nrow = numberOfSubjects, ncol = numberOfTrials))
  for (i in 1:numberOfTrials) {
    setTxtProgressBar(pb = progressBar, value = i)
    subjectInfo <- trials$subjectInformation[[i]]
    subjectInfo$response <- NA
    index <- which(subjectInfo$treatment == treatmentLabels[1] &
                     !is.na(subjectInfo$observationTime))
    subjectInfo$response[index] <- subjectOutcomes[[1]][1:length(index), i]
    index <- which(subjectInfo$treatment == treatmentLabels[2] &
                     !is.na(subjectInfo$observationTime))
    subjectInfo$response[index] <- subjectOutcomes[[2]][1:length(index), i]
    trials$subjectInformation[[i]] <- subjectInfo
  }
  cat("\n")
  cat("Done! \n")

  return(trials)
}

# # test
# numberOfTrials <- 1e+05
# numberOfSubjects <- 1006
# centreLabels <- 1:29
# treatmentLabels <- c("SOC", "HFOV")
# allocationRatio <- c(1, 1)
# blockSizes <- c(4, 6)
# followupTime <- 30
# dropoutRates <- c(0.03, 0.03)
# accrualTimes <- matrix(data = c(rep(x = 0, times = length(centreLabels)),
#                                 rep(x = 30 * 12 * 5, times = length(centreLabels))),
#                        nrow = length(centreLabels), ncol = 2)
# accrualIntensities <- rep(x = 5.5 / 7 / 29, times = length(centreLabels))
# seed <- 21
#
# system.time(trials <- getSubjectInformation(
#   numberOfTrials = numberOfTrials,
#   numberOfSubjects = numberOfSubjects,
#   centreLabels = centreLabels,
#   treatmentLabels = treatmentLabels,
#   allocationRatio = allocationRatio,
#   blockSizes = blockSizes,
#   followupTime = followupTime,
#   dropoutRates = dropoutRates,
#   accrualTimes = accrualTimes,
#   accrualIntensities = accrualIntensities,
#   seed = seed))
#
# mu <- c(0, 1)
# sigma <- c(1, 1)
# system.time(trials <- addSubjectOutcomesMeans(
#   trials = trials,
#   mu = mu,
#   sigma = sigma))

########################################

#' return virtual trials with subject information and simulation settings for testing rates in two samples
addSubjectOutcomesRates <- function(trials,
                                    pi) {
  set.seed(seed = trials$simulationSettings$seed)

  # update simulation settings
  trials$simulationSettings$pi <- pi

  numberOfTrials <- trials$simulationSettings$numberOfTrials
  numberOfSubjects <- trials$simulationSettings$numberOfSubjects
  treatmentLabels <- trials$simulationSettings$treatmentLabels

  # simulate outcomes
  cat("adding outcomes... \n")
  progressBar <- txtProgressBar(min = 1, max = numberOfTrials, style = 3)
  subjectOutcomes <-
    list(matrix(data = rbinom(n = numberOfSubjects * numberOfTrials,
                              size = 1, prob = pi[1]),
                nrow = numberOfSubjects, ncol = numberOfTrials),
         matrix(data = rbinom(n = numberOfSubjects * numberOfTrials,
                              size = 1, prob = pi[2]),
                nrow = numberOfSubjects, ncol = numberOfTrials))
  for (i in 1:numberOfTrials) {
    setTxtProgressBar(pb = progressBar, value = i)
    subjectInfo <- trials$subjectInformation[[i]]
    subjectInfo$response <- NA
    index <- which(subjectInfo$treatment == treatmentLabels[1] &
                     !is.na(subjectInfo$observationTime))
    subjectInfo$response[index] <- subjectOutcomes[[1]][1:length(index), i]
    index <- which(subjectInfo$treatment == treatmentLabels[2] &
                     !is.na(subjectInfo$observationTime))
    subjectInfo$response[index] <- subjectOutcomes[[2]][1:length(index), i]
    trials$subjectInformation[[i]] <- subjectInfo
  }
  cat("\n")
  cat("Done! \n")

  return(trials)
}

# # test
# numberOfTrials <- 1e+05
# numberOfSubjects <- 1006
# centreLabels <- 1:29
# treatmentLabels <- c("SOC", "HFOV")
# allocationRatio <- c(1, 1)
# blockSizes <- c(4, 6)
# followupTime <- 30
# dropoutRates <- c(0.03, 0.03)
# accrualTimes <- matrix(data = c(rep(x = 0, times = length(centreLabels)),
#                                 rep(x = 30 * 12 * 5, times = length(centreLabels))),
#                        nrow = length(centreLabels), ncol = 2)
# accrualIntensities <- rep(x = 5.5 / 7 / 29, times = length(centreLabels))
# seed <- 21
#
# system.time(trials <- getSubjectInformation(
#   numberOfTrials = numberOfTrials,
#   numberOfSubjects = numberOfSubjects,
#   centreLabels = centreLabels,
#   treatmentLabels = treatmentLabels,
#   allocationRatio = allocationRatio,
#   blockSizes = blockSizes,
#   followupTime = followupTime,
#   dropoutRates = dropoutRates,
#   accrualTimes = accrualTimes,
#   accrualIntensities = accrualIntensities,
#   seed = seed))
#
# pi <- c(0.45, 0.36)
# system.time(trials <- addSubjectOutcomesRates(
#   trials = trials,
#   pi = pi))

########################################

#' return virtual trials with subject information and simulation settings for testing the hazard ratio in two samples
addSubjectOutcomesSurvival <- function(trials,
                                       lambda) {
  set.seed(seed = trials$simulationSettings$seed)

  # update simulation settings
  trials$simulationSettings$lambda <- lambda

  numberOfTrials <- trials$simulationSettings$numberOfTrials
  numberOfSubjects <- trials$simulationSettings$numberOfSubjects
  treatmentLabels <- trials$simulationSettings$treatmentLabels
  followupTime <- trials$simulationSettings$followupTime

  # simulate outcomes
  cat("adding outcomes... \n")
  progressBar <- txtProgressBar(min = 1, max = numberOfTrials, style = 3)
  subjectOutcomes <-
    list(matrix(data = ceiling(rexp(n = numberOfSubjects * numberOfTrials,
                                    rate = lambda[1])),
                nrow = numberOfSubjects, ncol = numberOfTrials),
         matrix(data = ceiling(rexp(n = numberOfSubjects * numberOfTrials,
                                    rate = lambda[2])),
                nrow = numberOfSubjects, ncol = numberOfTrials))
  for (i in 1:numberOfTrials) {
    setTxtProgressBar(pb = progressBar, value = i)
    subjectInfo <- trials$subjectInformation[[i]]
    subjectInfo$time <- NA
    index <- which(subjectInfo$treatment == treatmentLabels[1])
    subjectInfo$time[index] <- subjectOutcomes[[1]][1:length(index), i]
    index <- which(subjectInfo$treatment == treatmentLabels[2])
    subjectInfo$time[index] <- subjectOutcomes[[2]][1:length(index), i]
    # censor
    subjectInfo$status <- ifelse(subjectInfo$time <= followupTime, 1, 0)
    subjectInfo$time <- pmin(subjectInfo$time, followupTime)
    # dropout
    index <- which(is.na(subjectInfo$observationTime))
    subjectInfo$status[index] <- 0
    subjectInfo$time[index] <-
      ceiling(runif(n = length(index), min = 0, max = subjectInfo$time[index]))
    subjectInfo$observationTime <- subjectInfo$randomisationTime + subjectInfo$time

    trials$subjectInformation[[i]] <- subjectInfo
  }
  cat("\n")
  cat("Done! \n")

  return(trials)
}

# # test
# numberOfTrials <- 1e+05
# numberOfSubjects <- 1006
# centreLabels <- 1:29
# treatmentLabels <- c("SOC", "HFOV")
# allocationRatio <- c(1, 1)
# blockSizes <- c(4, 6)
# followupTime <- 30
# dropoutRates <- c(0.03, 0.03)
# accrualTimes <- matrix(data = c(rep(x = 0, times = length(centreLabels)),
#                                 rep(x = 30 * 12 * 5, times = length(centreLabels))),
#                        nrow = length(centreLabels), ncol = 2)
# accrualIntensities <- rep(x = 5.5 / 7 / 29, times = length(centreLabels))
# seed <- 21
#
# system.time(trials <- getSubjectInformation(
#   numberOfTrials = numberOfTrials,
#   numberOfSubjects = numberOfSubjects,
#   centreLabels = centreLabels,
#   treatmentLabels = treatmentLabels,
#   allocationRatio = allocationRatio,
#   blockSizes = blockSizes,
#   followupTime = followupTime,
#   dropoutRates = dropoutRates,
#   accrualTimes = accrualTimes,
#   accrualIntensities = accrualIntensities,
#   seed = seed))
#
# lambda <- c(0.05, 0.03)
# system.time(trials <- addSubjectOutcomesSurvival(
#   trials = trials,
#   lambda = lambda))

################################################################################
