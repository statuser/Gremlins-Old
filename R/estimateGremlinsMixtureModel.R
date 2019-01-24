#' Estimate Gremlin's Model - Hierarchical MNL
#'
#' The function estimates the model described in "Gremlin's in the Data: Identifying the
#' Information Content of Research Subjects" using a hierarchical multinomial logit model
#'
#'
#' @param data A matrix containing the raw data.  The format of this data is a N row matrix.
#'   With the first column the respondnet identifier, followed by the design number, then
#'   Number of Observations indexes for the alternative chosen for each of the tasks.  The
#'   index is an integer with 1 representing the first alternative and n representing the
#'   final alternative.
#' @param design A Matrix representing the coded design file.  The design file should be
#'   formatted as a matrix with number of versions * number of tasks * number of alternatives
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
estimateGremlinsMixtureModel <- function(data, design, R, Priors = NULL, nSegments = 2, nBetaSegments = 2, covariates = NULL, constraints = NULL, segmentCovariates = NULL, startingValues = NULL) {

  gremlinsEnv$nParams <- ncol(design) - 3
  gremlinsEnv$nUnits <- nrow(data)
  gremlinsEnv$nTasks <- max(design[,2])
  gremlinsEnv$nConcepts <- max(design[,3])
  gremlinsEnv$nSegments <- nSegments
  gremlinsEnv$nBetaSegments <- nBetaSegments

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
        Ainv = solve(1*diag(ncol(covariates))),

        nu_not = gremlinsEnv$nParams + 3,
        V_not = (gremlinsEnv$nParams + 3)*diag(gremlinsEnv$nParams),

        lambdaScale = 5,
        lambdaShape = 5,

        psi_k = rep(5, gremlinsEnv$nSegments),
        psi_l = rep(5, gremlinsEnv$nBetaSegments)

      )
  } else {
    validatePriors(Priors, useDelta = useDelta)
  }

  if(is.null(startingValues)) {
    # This is a problem.  We need a better way  to find starting values:
    # Can we run a standard MNL model and then cluster separately on the fit for lambda segments
    # and the cluster on the betas for the regular model.

    #startingValues = findStartingValues()


    slopeCov <- array(double(gremlinsEnv$nBetaSegments*gremlinsEnv$nParams*gremlinsEnv$nParams),
                c(gremlinsEnv$nBetaSegments, gremlinsEnv$nParams, gremlinsEnv$nParams))
    for(l in 1:gremlinsEnv$nBetaSegments) {
      slopeCov[l,,] <- diag(gremlinsEnv$nParams)
    }

    startingValues = list(slope = matrix(0, nrow = gremlinsEnv$nUnits, ncol = gremlinsEnv$nParams),
      slopeBar = matrix(double(gremlinsEnv$nParams*gremlinsEnv$nBetaSegments), nrow = gremlinsEnv$nBetaSegments),
      slopeCov = slopeCov,
      lambda = seq(1, by=25, length.out = gremlinsEnv$nSegments),
      segMembership = sample(1:gremlinsEnv$nSegments, size = gremlinsEnv$nUnits, replace = TRUE),
      betaSegMembership = sample(1:gremlinsEnv$nBetaSegments, size = gremlinsEnv$nUnits, replace = TRUE),
      phi_lambda = rep(1/gremlinsEnv$nSegments, gremlinsEnv$nSegments),
      phi_beta = rep(1/gremlinsEnv$nBetaSegments, gremlinsEnv$nBetaSegments))
  } else {
    # TODO: Validate Startings values
  }


  # Initialize Storage
  gremlinsEnv$acceptanceRate <- matrix(0, nrow=3, ncol=2) # Three different jump sizes

  slope <- array(0, dim=c(R, gremlinsEnv$nUnits, gremlinsEnv$nParams))
  slope[1,,] <- startingValues$slope

  slopeBar <- array(double(R * gremlinsEnv$nParams * gremlinsEnv$nBetaSegments), c(R, gremlinsEnv$nBetaSegments, gremlinsEnv$nParams))
  for(segment in 1:gremlinsEnv$nBetaSegments) {
    slopeBar[1,segment,] <- startingValues$slopeBar[segment,]
  }


  slopeCov <- array(0, dim=c(R, gremlinsEnv$nBetaSegments, gremlinsEnv$nParams, gremlinsEnv$nParams))
  slopeCov[1,,,] <- startingValues$slopeCov

  lambda <- matrix(0, nrow = R, ncol = gremlinsEnv$nSegments)
  lambda[1,] <- startingValues$lambda

  K <- matrix(0, nrow = R, ncol = gremlinsEnv$nUnits)
  K[1,] <- startingValues$segMembership

  L <- matrix(0, nrow=R, ncol = gremlinsEnv$nUnits)
  L[1,] <- startingValues$betaSegMembership

  phi_lambda <- matrix(0, ncol = gremlinsEnv$nSegments, nrow = R)
  phi_lambda[1, ] <- startingValues$phi_lambda

  phi_beta <- matrix(0, ncol = gremlinsEnv$nBetaSegments, nrow = R)
  phi_beta[1, ] <- startingValues$phi_beta

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
      slope[rep, ind, ] <- generateIndividualSlope(respData, respDesign,
         slope[rep - 1, ind, ], lambda[rep-1, K[rep-1,ind]],
         slopeBar[rep-1, L[rep-1, ind],], slopeCov[rep-1,L[rep - 1, ind],,], constraints)
      	 K[rep, ind] <- generateSegmentMembership(respData, respDesign, slope[rep, ind, ], lambda[rep - 1,], phi_lambda[rep - 1,], Priors)
      	 L[rep, ind] <- generateBetaSegmentMembership(slope[rep,ind, ], slopeBar[rep - 1,,], slopeCov[rep - 1,,,], phi_beta[rep - 1,])
    }
    for(l in 1:gremlinsEnv$nBetaSegments) {
      multiregOut <- multivariateRegression(slope[rep,L[rep,] == l,], as.matrix(covariates[L[rep,] == l]), Priors)
      slopeBar[rep,l,] <- multiregOut$Beta
      slopeCov[rep,l,,] <- multiregOut$Sigma
    }
    lambda[rep,] <- generateLambda(lambda[rep - 1,], data, design, slope[rep,,], K[rep-1,], Priors)
    phi_lambda[rep,] <- drawPhi(K[rep,], Priors$psi_k)
    phi_beta[rep,] <- drawPhi(L[rep,], Priors$psi_l)
  }

  #; Save the output

  return(list(slope = slope, slopeBar = slopeBar, slopeCov = slopeCov, lambda = lambda, segmentMembership = K, betaSegmentMembership = L, lambdaSegProb = phi_lambda, betaSegProb = phi_beta))

}


generateIndividualSlope <- function(data, design, currentSlope, lambda, slopeBar, slopeCov, constraints) {

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

  currentLL <- gremlinsLogLikelihood(data, design, currentSlope, lambda)
  proposedLL <- gremlinsLogLikelihood(data, design, proposedSlope, lambda)

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

#' Utility function to comute a multivariate regression.  Based on information from "Bayesian Statistics and Marketing (2005)"
#' @importFrom MASS mvrnorm
multivariateRegression <- function(data, design, priors) {
  U = chol(priors$Ainv)
  R = rbind(design, U)
  Q = rbind(data, U %*% priors$mu_not)
  BTilde = chol2inv(chol(crossprod(R))) %*% crossprod(R, Q)

  S = crossprod(Q - R%*%BTilde)

  Sigma <- chol2inv(chol(drop(rWishart(1, priors$nu_not + nrow(design), chol2inv(chol(priors$V_not + S))))))

  mean <- as.vector(BTilde)
  variance <- Sigma %x% chol2inv(chol(crossprod(design) + priors$Ainv))

  Beta <- mean + rnorm(length(mean), 0, 1) %*% chol(variance)
  Beta <- matrix(Beta, nrow = ncol(design))

  return(list(Beta = Beta, Sigma = Sigma))

}

mvn <- function(mean, sigma) {
  drop(mean + rnorm(length(mean)) %*% chol(sigma))
}

generateLambda <- function(currentLambda, data, design, currentSlope, K, priors) {
  #; This is a Metropolis-Hastings step with one little twist.  The lambdas are constrained to be in
  #; increasing magnitude.  This is accomplished using rejection sampling for the truncated distributions.
  #; This is equivalent to applying a truncated prior, but is more efficient due to the potentially high rejection
  #; rate on the draws.  Because of this it is necessary to correct for the non-symetric proposal distribution in
  #; the acceptance probability step.

  #; Draw the lambdas.  (Lambda[1] is always 1)
  if(gremlinsEnv$nSegments == 1) {
    return(1)
  }

  lambda <- currentLambda

  for(k in 2:gremlinsEnv$nSegments) {
    if(k < gremlinsEnv$nSegments) {
      proposedLambda <- lambda[k] + ((lambda[k+1] - lambda[k-1])/6) * rnorm(1, 0, 1)
      while(proposedLambda < lambda[k-1] || proposedLambda > lambda[k+1]) {
        proposedLambda <- lambda[k] + ((lambda[k+1] - lambda[k-1])/6) * rnorm(1, 0, 1)
        cat("Lambda Draw Violated Contraints")
      }
    } else {
      proposedLambda <- lambda[k] + ((lambda[k] - lambda[1])/(6*gremlinsEnv$nSegments)) * rnorm(1, 0, 1)
      # while(proposedLambda < lambda[k-1]) {
      while(proposedLambda < lambda[k-1]) {
        proposedLambda <- lambda[k] + ((lambda[k] - lambda[1])/(6*gremlinsEnv$nSegments)) * rnorm(1, 0, 1)
        cat("Lambda Draw Violated Contraints")
      }
    }

    cIndLL <- double(gremlinsEnv$nUnits)
    pIndLL <- double(gremlinsEnv$nUnits)

    for(ind in 1:gremlinsEnv$nUnits) {
      if(K[ind] == k) {

        cIndLL[ind] <- gremlinsLogLikelihood(data[ind, -c(1:2)], design[design[,1] == data[ind, 2], -c(1:3)], currentSlope[ind,], lambda[k])
        pIndLL[ind] <- gremlinsLogLikelihood(data[ind, -c(1:2)], design[design[,1] == data[ind, 2], -c(1:3)], currentSlope[ind,], proposedLambda)
      }
    }

    cLL <- sum(cIndLL)
    pLL <- sum(pIndLL)
    cPrior <- (priors$lambdaShape - 1) * log(lambda[k]) - lambda[k]/priors$lambdaScale # * sum(K[i == k])
    pPrior <- (priors$lambdaShape - 1) * log(proposedLambda) - proposedLambda/priors$lambdaScale # * sum(K[i == k])

    lap <- pLL + pPrior - cLL - cPrior

    if(log(runif(1)) < lap) {
      lambda[k] <- proposedLambda
    }

  }


  return(lambda)

}

generateSegmentMembership <- function(data, design, slope, lambda, phi_lambda, priors) {
  prob <- double(gremlinsEnv$nSegments)
  for(k in 1:gremlinsEnv$nSegments) {
    prob[k] <- gremlinsLogLikelihood(data, design, slope, lambda[k])
    prob[k] <- prob[k] + log(phi_lambda[k])
  }

  prob <- exp(prob)
  return(which.max(rmultinom(1, 1, prob)))
}

generateBetaSegmentMembership <- function(slope, slopeBar, slopeCov, phi_beta) {
  if(gremlinsEnv$nBetaSegments > 1) {
    prob <- double(gremlinsEnv$nBetaSegments)
    for(l in 1:gremlinsEnv$nBetaSegments) {
      U = chol(solve(slopeCov[l,,]))
      slopeDiff <- U %*% (slope - slopeBar[l,])
      prob[l] <- prod(diag(U)) * exp(-0.5 * crossprod(slopeDiff)) * phi_beta[l]
    }

    return(which.max(rmultinom(1, 1, prob)))
  } else {
    return(1)
  }
}

gremlinsLogLikelihood <- function(data, design, slope, lambda) {
  utils <- design %*% slope / lambda
  utils <- matrix(utils, ncol=gremlinsEnv$nConcepts, byrow=TRUE)
  return(sum(utils[cbind(seq_along(data), data)] - log(rowSums(exp(utils)))))
}

validatePriors <- function(Priors, useDelta = FALSE) {

  if(is.null(Priors$mu_not) ||
     ncol(Priors$mu_not) != gremlinsEnv$nParams ||
     nrow(Priors$mu_not) != gremlinsEnv$nCovariates ||
     typeof(Priors$mu_not) != "double") {
    stop("mu_not should be a matrix of numerics values
         with columns equal to the number of parameters
         and rows equal to the number of covariates.")
  }

  if(is.null(Priors$Ainv) ||
     !is.matrix(Priors$Ainv) ||
     ncol(Priors$Ainv) != nrow(Priors$Ainv) ||
     ncol(Priors$Ainv) != gremlinsEnv$nCovariates ||
     typeof(Priors$Ainv) != "double") {
    stop("Ainv should be a square matrix of numeric values
         with dimensions equal to the number of covariates")
  }



  #TODO Finish validating the priors

}

drawPhi <- function(K, priors) {
  N_k <- double(length(priors))
  for(i in 1:length(priors)) {
    N_k[i] <- sum(K == i)
  }
  return(rdirichlet(N_k + priors))
}

rdirichlet <- function(alpha) {
  y <- rgamma(length(alpha), alpha)
  return(y/sum(y))
}
