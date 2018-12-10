#' Estimate Gremlin's Model - Hierarchical MNP
#'
#' The function estimates the model described in "Gremlin's in the Data: Identifying the
#' Information Content of Research Subjects" using a hierarchical multinomial probit model
#'
#'
#' @param data A matrix containing the raw data.  The format of this data is a N row matrix.
#'   With the first column the respondnet identifier, followed by the design number, then
#'   Number of Observations indexes for the alternative chosen for each of the tasks.  The
#'   index is an integer with 1 representing the first alternative and n representing the
#'   final alternative.
#' @param design A Matrix representing the coded design file.  The design file should be
#'   formatted as a matrix with number of versions \* number of tasks \* number of alternatives
#'   rows and number of paramters + 3 columns.  The first column contains the version number,
#'   The second columns contains the task number, the third number contains the alternative,
#'   and the remaining columns contain the coded design.  A generic Sawtooth Software design
#'   file can be converted to this format using the \code{codeSawtoothDesign} function.
#'  @param Priors A datastructure that contains the priors for to the model.  Can be null indicating
#'   the use of default priors or must contain a full prior specification.
#'  @param numSegments (Default = 2) The number of segments for the scale factor
#'  @param R The number of repitions in the chain
#'  @param covariates (Optional) A matrix of covariates for the model.  One row per respondent with the
#'   respondent identified in the first column.
#'  @param constraints (Optional) a vector of length n-param sepecifying the constraints
#'    to impose on the parameters or NULL.  a 1 indicates the parameter is constrained to be positive
#'    a -1 constrains to be  negative, and a 0 indicates no constraint.
#'  @return A data structure containing the draws from the complete MCMC chain
#'
#'  @seealso \code{\link{codeSawtoothDesign}} to create the design from a Sawtooth file,
#'  \code{\link{extractDraws}} for examples on how to process the output
#'
estimateGremlinsMNPModel <- function(data, design, R, Priors = NULL, nSegments = 2, covariates = NULL, constraints = NULL, segmentCovariates = NULL, startingValues = NULL) {

  gremlinsEnv$nParams <- ncol(design) - 3
  gremlinsEnv$nUnits <- nrow(data)
  gremlinsEnv$nTasks <- max(design[,2])
  gremlinsEnv$nConcepts <- max(design[,3])
  gremlinsEnv$nSegments <- nSegments

  design <- as.matrix(design)
  data <- as.matrix(data)
  #Check for covariates on segment Variables

  # TODO: Not implemented
  if(!is.null(segmentCovariates)) {
    gremlinsEnv$useDelta = TRUE
    gremlinsEnv$nDeltaParams = ncol(segmentCovariates)
  } else {
    gremlinsEnv$useDelta = FALSE
  }

  #Set up covariates to just an intercept if not supplied
  if(is.null(covariates)) {
    covariates = as.matrix(rep(1, gremlinsEnv$nUnits))
    gremlinsEnv$nCovariates = 1
  } else {
    if(!is.matrix(covariates) || nrow(covariates) != gremlinsEnv$nUnits) {
      stop("covariates must be a matrix with
           nrows equal to the number of units")
    }
    gremlinsEnv$nCovariates = ncol(covariates)
    }

  if(is.null(constraints)) {
    constraints = double(gremlinsEnv$nParams)
  } else {
    if(length(constraints) != gremlinsEnv$nParams || !all(constraints %in% c(-1, 0, 1))) {
      stop("constraints must contain a value fo reach parameter and
           can only be -1, 0, or 1 for negative, no constraint, or positive")
    }
    }
  #validate and Set priors
  if(is.null(Priors)) {
    Priors = list(
      mu_not = matrix(0, ncol = gremlinsEnv$nParams, nrow = gremlinsEnv$nCovariates),
      Ainv = solve(100*diag(ncol(covariates))),

      nu_not = ncol(covariates) + 3,
      V_not = (ncol(covariates) + 3)*diag(gremlinsEnv$nParams),

      lambdaScale = 5,
      lambdaShape = 5,

      psi_k = rep(5, gremlinsEnv$nSegments)

    )
  } else {
    validatePriors(Priors, useDelta = useDelta)
  }

  if(is.null(startingValues)) {
    startingValues = list(slope = matrix(0, nrow = gremlinsEnv$nUnits, ncol = gremlinsEnv$nParams),
                          slopeBar = double(gremlinsEnv$nParams),
                          slopeCov =  diag(gremlinsEnv$nParams),
                          lambda = seq(1, by=25, length.out = gremlinsEnv$nSegments),
                          segMembership = sample(1:gremlinsEnv$nSegments, size = gremlinsEnv$nUnits, replace = TRUE),
                          phi_lambda = rep(1/gremlinsEnv$nSegments, gremlinsEnv$nSegments))
  } else {
    # TODO: Validate Startings values
  }


  # Initialize Storage
  gremlinsEnv$acceptanceRate <- matrix(0, nrow=3, ncol=2) # Three different jump sizes

  slope <- array(0, dim=c(R, gremlinsEnv$nUnits, gremlinsEnv$nParams))
  slope[1,,] <- startingValues$slope

  slopeBar <- array(double(R * gremlinsEnv$nParams), c(R, gremlinsEnv$nParams))
  slopeBar[1,] <- startingValues$slopeBar


  slopeCov <- array(0, dim=c(R, gremlinsEnv$nBetaSegments, gremlinsEnv$nParams, gremlinsEnv$nParams))
  slopeCov[1,,] <- startingValues$slopeCov

  lambda <- matrix(0, nrow = R, ncol = gremlinsEnv$nSegments)
  lambda[1,] <- startingValues$lambda

  K <- matrix(0, nrow = R, ncol = gremlinsEnv$nUnits)
  K[1,] <- startingValues$segMembership

  phi_lambda <- matrix(0, ncol = gremlinsEnv$nSegments, nrow = R)
  phi_lambda[1, ] <- startingValues$phi_lambda

  cat("Beginning MCMC Routine\n")

  # Run MCMC
  for(rep in 2:R) {
    #; Provide status update on code
    if(rep %% 10 == 0) {
      cat(paste("Completing iteration: ", rep, "\n"))
    }


    for(ind in 1:gremlinsEnv$nUnits) {
      respData <- data[ind, -c(1:2)] # Drop the RespID and Design Version
      respDesign <- design[design[,1] == data[ind, 2], -c(1:3)] # Drop the Version, Task, and Concept Columns
      slope[rep, ind, ] <- generateIndividualSlopeProbit(respData, respDesign,
                                                   slope[rep - 1, ind, ], lambda[rep-1, K[rep-1,ind]],
                                                   slopeBar[rep-1, ], slopeCov[rep-1,,], constraints)
      K[rep, ind] <- generateSegmentMembershipProbit(respData, respDesign, slope[rep, ind, ], lambda[rep - 1,], phi_lambda[rep - 1,], Priors)
    }
    multiregOut <- multivariateRegression(slope[rep,,], as.matrix(covariates), Priors)
    slopeBar[rep,] <- multiregOut$Beta
    slopeCov[rep,,] <- multiregOut$Sigma

    lambda[rep,] <- generateLambdaProbit(lambda[rep - 1,], data, design, slope[rep,,], K[rep-1,], Priors)
    phi_lambda[rep,] <- drawPhi(K[rep,], Priors$psi_k)
  }

  #; Save the output

  return(list(slope = slope, slopeBar = slopeBar, slopeCov = slopeCov, lambda = lambda, segmentMembership = K, lambdaSegProb = phi_lambda))

}

generateIndividualSlopeProbit <- function(data, design, currentSlope, lambda, slopeBar, slopeCov, constraints) {

  jumpSizeIndex <- sample(seq_along(gremlinsEnv$jumpSizeProbs), 1, prob = gremlinsEnv$jumpSizeProbs)
  jumpSize <- gremlinsEnv$jumpSizes[jumpSizeIndex]

  proposedSlope <- double(length(currentSlope))
  constraintProbability <- 0

  for(i in seq_along(currentSlope)) {
    proposedSlope[i] = currentSlope[i] + rnorm(1, 0, 1)*jumpSize
    if(constraints[i] != 0) {
      numTries <- 0
      while(sign(proposedSlope[i]) != constraints[i]) {
        proposedSlope[i] = currentSlope[i] + rnorm(1, 0, 1)*jumpSize
        numTries <- numTries + 1
        if(numTries == gremlinsEnv$totalConstraintTries) {
          warning("Constraint could not be satisfied for slope of paramter ", i, ".")
          proposedSlope[i] <- currentSlope[i]
        } else {
          constraintProbability <- constraintProbability + pnorm(0, currentSlope[i], jumpSize, lower.tail = (constraints[i] < 0), log.p = TRUE)
          - pnorm(0, proposedSlope[i], jumpSize, lower.tail = (constraints[i] < 0), log.p = TRUE)
        }
      }
    }
  }

  currentLL <- gremlinsLogLikelihoodProbit(data, design, currentSlope, lambda)
  proposedLL <- gremlinsLogLikelihoodProbit(data, design, proposedSlope, lambda)

  currentPrior <- -0.5 * t(currentSlope - slopeBar)%*%solve(slopeCov)%*%(currentSlope - slopeBar)
  proposedPrior <- -0.5 * t(proposedSlope - slopeBar)%*%solve(slopeCov)%*%(proposedSlope - slopeBar)

  logAcceptProb <- proposedLL + proposedPrior - currentLL - currentPrior + constraintProbability

  alpha <- log(runif(1))

  if(alpha < logAcceptProb) {
    currentSlope <- proposedSlope
    gremlinsEnv$acceptanceRate[jumpSizeIndex, 1] <- gremlinsEnv$acceptanceRate[jumpSizeIndex, 1] + 1 # Column one is the acceptance, columns 2 is the number of tries
  }

  gremlinsEnv$acceptanceRate[jumpSizeIndex, 2] <- gremlinsEnv$acceptanceRate[jumpSizeIndex, 2] + 1

  return(currentSlope)
}

generateLambdaProbit <- function(currentLambda, data, design, currentSlope, K, priors) {
  #; This is a Metropolis-Hastings step with one little twist.  The lambdas are constrained to be in
  #; increasing magnitude.  This is accomplished using rejection sampling for the truncated distributions.
  #; This is equivalent to applying a truncated prior, but is more efficient due to the potentially high rejection
  #; rate on the draws.  Because of this it is necessary to correct for the non-symetric proposal distribution in
  #; the acceptance probability step.

  #; Draw the lambdas.  (Lambda[1] is always 1)
  lambda <- currentLambda

  for(k in 2:gremlinsEnv$nSegments) {
    if(k < gremlinsEnv$nSegments) {
      proposedLambda <- lambda[k] + ((lambda[k+1] - lambda[k-1])/6) * rnorm(1, 0, 1)
      while(proposedLambda < lambda[k-1] || proposedLambda > lambda[k+1]) {
        proposedLambda <- lambda[k] + ((lambda[k+1] - lambda[k-1])/6) * rnorm(1, 0, 1)
      }
    } else {
      proposedLambda <- lambda[k] + ((lambda[k] - lambda[1])/(6*gremlinsEnv$nSegments)) * rnorm(1, 0, 1)
      while(proposedLambda < lambda[k-1]) {
        proposedLambda <- lambda[k] + ((lambda[k] - lambda[1])/(6*gremlinsEnv$nSegments)) * rnorm(1, 0, 1)
      }
    }

    cIndLL <- double(gremlinsEnv$nUnits)
    pIndLL <- double(gremlinsEnv$nUnits)

    for(ind in 1:gremlinsEnv$nUnits) {
      if(K[ind] == k) {

        cIndLL[ind] <- gremlinsLogLikelihoodProbit(data[ind, -c(1:2)], design[design[,1] == data[ind, 2], -c(1:3)], currentSlope[ind,], lambda[k])
        pIndLL[ind] <- gremlinsLogLikelihoodProbit(data[ind, -c(1:2)], design[design[,1] == data[ind, 2], -c(1:3)], currentSlope[ind,], proposedLambda)


      }
    }

    cLL <- sum(cIndLL)
    pLL <- sum(pIndLL)
    cPrior <- (priors$lambdaShape + 1) * log(lambda[k]) - lambda[k]/priors$lambdaScale
    pPrior <- (priors$lambdaShape + 1) * log(proposedLambda) - proposedLambda/priors$lambdaScale

    lap <- pLL + pPrior - cLL - cPrior

    if(log(runif(1)) < lap) {
      lambda[k] <- proposedLambda
    }

  }


  return(lambda)

}

generateSegmentMembershipProbit <- function(data, design, slope, lambda, phi_lambda, priors) {
  prob <- double(gremlinsEnv$nSegments)
  for(k in 1:gremlinsEnv$nSegments) {
    prob[k] <- gremlinsLogLikelihoodProbit(data, design, slope, lambda[k])
    if(k != 1) {
      prob[k] <- prob[k] + (priors$lambdaShape - 1) * log(lambda[k]) - lambda[k] / priors$lambdaScale
    }
    prob[k] <- prob[k] + log(phi_lambda[k])
  }

  prob <- exp(prob)
  return(which.max(rmultinom(1, 1, prob)))
}

gremlinsLogLikelihoodProbit <- function(data, design, slope, lambda) {
  logLike <- double(gremlinsEnv$nTasks)
  for(i in 1:gremlinsEnv$nTasks) {
    taskDesign <- design[((i-1)*gremlinsEnv$nConcepts + 1):(i*gremlinsEnv$nConcepts), ]
    logLike[i] <- llmnp(slope, lambda*diag(gremlinsEnv$nConcepts), taskDesign, data[i], 30) #30 should be sufficient
  }
  return(sum(logLike))
}
