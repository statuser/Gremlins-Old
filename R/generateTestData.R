#' Generate Test data give an sepcific design
#'
#' @param nRespondents The number of respondents
#' @param design The design file to use for the simulation
#' @param covariates The covariates to be used from the upper level of heterogeniety
#' @param nSegments (Default = 2) The number of segments for the scale factor.  Must be an integer
#' @return The matrix of choices
#'
#' @examples
#'
#'
generateTestData <- function(nRespondents, design, covariates = NULL, constraints = NULL, nSegments = 2) {

    nParam <- ncol(design) - 3
    nTasks <- length(unique(design[,2]))
    nConcepts <- length(unique(design[,3]))

    if(is.null(constraints)) {
      useConstraints <- FALSE
    } else {
      # - Check to make sure that the constraints are valid
      if(!all(constraints %in% c(-1, 0, 1))) {
        stop("The specified constraints are not valid.\nConstraints must be of the -1, 0, or 1")
      }
      if(length(constraints) != nParam) {
        stop("The specified constraints are not valid.\nThe constraints vector have
             one value of each parameter")
      }
    }

    if(is.null(covariates)) {
      covariates <- rep(1, nRespondents)
    } else {
      stop("Covariates are not yet implemented")
    }

    is.wholenumber <-
      function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    stopifnot(is.wholenumber(nSegments))


    # We need to generate parameters for the following parameters:
    # - beta_n - The individual betas - Should be draws from a MVN with mean beta_bar and variariance Sigma
    # - beta_bar - a P length vector of overall means
    # - sigma - a number of segments vector.  The first element is 1
    # - K - probability of belonging to each segment of lambda


    beta_bar <- rnorm(nParam, 0, 1)
    if(useConstraints) {
      # Swap the sign of any parameters that violate the constraints
      beta_bar[beta_bar * constraints < 0] <- beta_bar[beta_bar * constraints < 0] * -1
    }

    Sigma <- diag(1/rgamma(nParam, shape = 10, scale = 1/10)) # All variances must be positive a rate of 1/2 seems to give reasonable values

   phi_lambda <- rdirichlet(rep(1/nSegments, nSegments))

   beta <- matrix(0, ncol=nParam, nrow=nRespondents)
  K <- double(nRespondents)
    for(i in 1:nRespondents) {
      K[i] <- sample(1:nSegments, 1, prob = phi_lambda)
      beta[i,] <- beta_bar + rnorm(nParam, 0, 1)%*%chol(Sigma)
      if(useConstraints) {
        # Make sure that the constraints aren't violated.  We can't just switch the sign in this case
        badConstraints <- which(beta[i,] * constraints < 0)
        for(j in badConstraints) {
          while(beta[i, j] * constraints[j] < 0) {
            beta[i, j] = beta_bar[j] + rnorm(1, 0, 1) * sqrt(Sigma[j])
          }
        }

      }
    }
  lambda <- double(nSegments)

  lambda[1] <- 1
  if(nSegments > 1) {
    for(i in 2:nSegments) {
      lambda[i] <- rgamma(1, shape = (i-1)*5, scale = (i-1)*5)
      while(lambda[i] <= lambda[i-1]) {
        lambda[i] <- rgamma(1, shape = (i-1)*5, scale = (i-1)*5)
      }
    }
  }

    cat("Generating synthetic data with the following settings\n")
    cat("beta-bar:\n")
    print(beta_bar)
    cat("Sigma:\n")
    print(Sigma)
    cat("lambda:\n")
    print(lambda)

    data <- matrix(0, nrow=nRespondents, ncol = nTasks + 2)
    parameters <- matrix(0, nrow=nRespondents, ncol = nParam)
    segmentMembership <- double(nRespondents)

    for(i in 1:nRespondents) {
      data[i,1] <- i
      data[i,2] <- version <- sample(1:max(design[,1]), 1)


      respDesign <- design[design[,1] %in% version,]

      respDesign <- respDesign[,-c(1:3)] #Drop the columns for the control variables

      xbeta <- exp(as.matrix(respDesign) %*% (beta[i,]/lambda[K[i]]))
      xbeta <- matrix(xbeta, ncol=nConcepts, byrow=TRUE)
      sumXBeta <- apply(xbeta, 1, sum) %x% t(rep(1, nConcepts))
      probs <- xbeta/sumXBeta
      y <- apply(probs, 1, function(x) {sample(1:nConcepts, 1, prob = x)})
      data[i, 3:ncol(data)] <- y


    }

    return(list(data = data, trueValues = list(slope = beta, slopeBar = beta_bar, slopeCov = Sigma, lambda = lambda, segMembership = K,  phi_lambda = phi_lambda)))

}
