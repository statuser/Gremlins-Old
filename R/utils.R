library(data.table)

getCameraDataOutput = function(fname){
  # Get camera data from file
  mcmcData = read.csv(fname, header=FALSE, sep = ' ')

  mcmcHeaders = c("iteration",'individual','sig2','segment', 'Canon',
                  'Sony', 'Nikon', 'Panasonic', 'MegaPixels', 'Zoom', 'HDVideo', 'Swivel', 'Wifi', 'Price')

  covariates = c('Canon', 'Sony', 'Nikon', 'Panasonic', 'MegaPixels', 'Zoom', 'HDVideo', 'Swivel', 'Wifi', 'Price')
  names(mcmcData) = mcmcHeaders
  return(mcmcData)
}

flattenOutput = function(output){
  n_ind = dim(output[['slope']])[2]
  n_betas = dim(output[['slope']])[3]
  n_draw = dim(output[['slope']])[1]

  p = 4 + n_betas
  n_row = n_draw*n_ind
  mcmc_ind_data = array(0, dim=c(n_row, p))


  prev_row_idx = 0
  for(i in 1:n_draw){
    row_idx = i*n_ind
    prev_row_idx = prev_row_idx + 1

    # Iteration number
    mcmc_ind_data[prev_row_idx:row_idx,1] = i
    # Respondent number
    mcmc_ind_data[prev_row_idx:row_idx,2] = 1:n_ind

    # Segment number
    mcmc_ind_data[prev_row_idx:row_idx,4] = output[['segmentMembership']][i,]

    mcmc_ind_data[prev_row_idx:row_idx,3] = output[['lambda']][i,2]


    # Betas
    mcmc_ind_data[prev_row_idx:row_idx,5:p] = output[['slope']][i,,]
    prev_row_idx = row_idx

  }

  mcmc_ind_data[,3] = ifelse(mcmc_ind_data[,4] == 1, 1, mcmc_ind_data[,3])
  return(mcmc_ind_data)
}


getParameterSummary = function(mcmcData, covariates){
  # Function to retrieve summary table of MCMC results
  mcmcData = as.data.frame(flattenOutput(mcmcData))
  mcmcHeaders = c("iteration",'individual','sig2','segment')
  mcmcHeaders = c(mcmcHeaders, covariates)
  names(mcmcData) = mcmcHeaders


  mcmcDataScaled = mcmcData[,covariates]*1/mcmcData$sig2
  mcmcDataScaled = cbind(mcmcData[,c("iteration",'individual','sig2','segment')],mcmcDataScaled)
  mcmcDataScaled = data.table(mcmcDataScaled)


  #means = mcmcDataScaled[ ,.(mean(Canon),mean(Sony),mean(Nikon), mean(Panasonic), mean(MegaPixels), mean(Zoom), mean(HDVideo), mean(Swivel), mean(Wifi), mean(Price)), by = segment]
  means = mcmcDataScaled[, lapply(.SD, mean), by = segment]
  means = as.data.frame(means)
  means = means[,c('segment',covariates)]



  #sd = mcmcDataScaled[ ,.(sd(Canon),sd(Sony),sd(Nikon), sd(Panasonic), sd(MegaPixels), sd(Zoom), sd(HDVideo), sd(Swivel), sd(Wifi), sd(Price)), by = segment]
  sd = mcmcDataScaled[, lapply(.SD, sd), by = segment]
  sd = as.data.frame(sd)
  sd = sd[,c('segment',covariates)]



  allMeans = colMeans(as.data.frame(mcmcDataScaled)[,covariates])
  allSd = apply(as.data.frame(mcmcDataScaled)[,covariates],2,sd)


  nStates = max(means$segment)
  paramSummary = matrix(0,nrow=length(covariates),ncol=2*nStates + 2)


  for(b in 1:length(covariates)){
    beta = covariates[b]


    for(n in 1:(nStates*2)){
      state = ceiling(n/2)

      if(n %% 2 == 0){
        paramSummary[b,n] = sd[sd$segment == state,beta]
      }
      else{
        #print(state)
        #print(n)
        paramSummary[b,n] = means[means$segment == state,beta]
      }

      ncols = dim(paramSummary)[2]

      paramSummary[b,ncols] = allSd[beta]
      paramSummary[b,ncols-1] = allMeans[beta]


    }

  }


  paramSummary = as.data.frame(paramSummary)
  col_names = array("",nStates*2 + 2)
  for(i in 1:(nStates*2)){
    segment = ceiling(i/2)
    if(i %% 2 == 1){
      header = paste("mean_segment_",segment,sep='')
    }
    else{
      header = paste("sd_segment_",segment,sep='')
    }
    col_names[i] = header
  }

  col_names[nStates*2 + 1] = "meanAll"
  col_names[nStates*2 + 2] = "sdAll"



  names(paramSummary) = col_names
  rownames(paramSummary) = covariates
  return(paramSummary)
}

simulateTask <- function(task, betas) {

  xbeta = exp(task%*%t(betas))
  sum_logit = sum(xbeta[,1])
  probs = xbeta / sum_logit


  return(probs)
}



calcExpectedShares = function(task,mcmcBetas,numBetas){
  # Function to calculate expected revenue
  # Removed the niter argument
  # Function expects avg beta_i's for each indiviual
  # Best when used with slopeOverSig2 estimates
  # Added segmentProbs for weighting by probability of being in high info. group

  unique_iter = unique(mcmcBetas[,1])

  out = matrix(0,nrow=length(unique_iter),ncol=nrow(task))
  rowCtn = 1
  for(j in unique_iter){
    {cat("\r",paste("Simulating choices at MCMC iteration: ", j, sep=''))}

    betas_iter = mcmcBetas[mcmcBetas[,1] == as.numeric(j),]

    iter_results = matrix(0,nrow=nrow(betas_iter),ncol=nrow(task))
    for(i in 1:nrow(betas_iter) ){

      beta = as.vector(betas_iter[i,])
      # Assume that mcmcBetas has four columns before slopes:
      #   - Iteration, repondent ID, lambda, segment
      colStart = 5
      beta = beta[colStart:length(beta)]

      iter_results[i,] = t(simulateTask(task,t(beta)) )

    }
    totalMkt = colSums(iter_results)
    mktShare = totalMkt/sum(totalMkt)
    out[rowCtn,] = mktShare
    rowCtn = rowCtn + 1
  }
  return(out)

}


marketShareSimulation = function(mcmcData, choice_names, holdoutTask){
  # Perform market share simulation using MCMC from beta_i's
  # Arguments:
  #   mcmcData: (DataFrame) Individual-level output from MCMC algorithm
  #               - row i: draws for individual ii at iteration j
  #   choice_names: (array) Names of choices in task:
  #               - e.g., c('Panasonic', 'Nikon', 'Sony', 'Canon', 'None')
  #   holdOutTask: (matrix?) A matrix encoding the choice task for the simulation
  #               - row i: denotes the choice_names of each covariate (in column j)
  #                 e.g.,
  #                      Canon Sony Nikon Panasonic MegaPixels Zoom HDVideo Swivel Wifi Price
  #                      [1,]     0    0     0         1          0    1       0      1    0  1.20
  #                      [2,]     0    0     1         0          0    1       0      0    1  1.25
  #                      [3,]     0    1     0         0          1    0       0      1    0  1.30
  #                      [4,]     1    0     0         0          1    1       1      1    0  1.40
  #                      [5,]     0    0     0         0          0    0       0      0    0  0.00

  mcmcData = flattenOutput(mcmcData)

  if(!is.matrix(holdoutTask)){
    holdoutTask = as.matrix(holdoutTask)
  }

  mktShareSim = calcExpectedShares(holdoutTask,mcmcData,10)
  mktShare = as.data.frame(cbind(colMeans(mktShareSim),apply(mktShareSim,2,sd)))
  names(mktShare) = c("Mean","Standard Deviation")
  rownames(mktShare) = choice_names

  return(mktShare)
}


calcExpectedRevenue = function(task,betas,priceCol){
  # Function to calculate expected revenue
  # Removed the niter argument
  # Function expects avg beta_i's for each indiviual
  # Best when used with slopeOverSig2 estimates


  out = matrix(0,nrow=nrow(betas),ncol=nrow(task))

  for(i in 1:nrow(betas) ){
    beta = as.vector(betas[i,])
    sims = simulateTask(task,t(beta))
    out[i,] = sims * task[,priceCol]

  }
  return(colSums(out))

}


priceOptimizationSimulation = function(mcmcDta,simTask,price_grid,optimBrand,numCovariates,price_col,intervalRange = c(.15,.85)){

  iterations = unique(mcmcDta[,1])
  optimalPriceArray = rep(0,length(iterations))

  start = proc.time()

  brandOptimizationMatrix =  matrix(0,nrow=length(price_grid),ncol=length(iterations))
  max_iter = iterations*length(price_grid)
  prog_cntr = 1
  iterationCnt = 1
  for(i in iterations){


    draws_i = mcmcDta[mcmcDta[,1] == i,]
    nColsDraws = dim(draws_i)[2]

    betas_i = draws_i[,(nColsDraws-numCovariates + 1):nColsDraws]



    revenue_i = matrix(0,nrow=length(price_grid),ncol=nrow(simTask))
    idx = 1
    for(cantd_price in price_grid)
    {
      {cat("\r",paste("MCMC Iteration:", i, "Candidate Price:", cantd_price, sep=' '))}
      # Assuming price column is last!!
      simTask[optimBrand,numCovariates] = cantd_price
      simulation_j = calcExpectedRevenue(simTask,betas_i,price_col)/nrow(betas_i)

      revenue_i[idx,] = simulation_j

      idx = idx + 1
      prog_cntr = prog_cntr + 1
    }


    brandOptimizationMatrix[,iterationCnt] = revenue_i[,optimBrand]

    # Insert which.max code here
    optimalPrice_i = price_grid[which.max(revenue_i[,optimBrand])]*100
    optimalPriceArray[iterationCnt] = optimalPrice_i
    iterationCnt = iterationCnt + 1

  }

  end = proc.time()

  end - start

  priceMeans = rowMeans(brandOptimizationMatrix)
  priceSd = apply(brandOptimizationMatrix,1,sd)


  lowError = priceMeans - priceSd
  upperError = priceMeans + priceSd

  optimalPriceMean = mean(optimalPriceArray)
  optimalPriceInterval = quantile(optimalPriceArray,intervalRange,names=F)

  optimalPriceOut = c(optimalPriceMean,optimalPriceInterval)

  out = list(cbind(priceMeans,lowError,upperError),optimalPriceOut)

  return(out)

}






runPriceOptimization = function(mcmcDataList, holdoutTask, minPrice, maxPrice, priceStep, featureToOptimize, price_col){

  mcmcData = flattenOutput(mcmcDataList)


  if(!is.matrix(holdoutTask)){
    holdoutTask = as.matrix(holdoutTask)
  }

  numCovariates = dim(holdoutTask)[2]
  price_grid = seq(minPrice, maxPrice, priceStep)

  optim1 = priceOptimizationSimulation(mcmcDta=mcmcData, holdoutTask,price_grid = price_grid,
                                       optimBrand=featureToOptimize,numCovariates=numCovariates,price_col=price_col)



  return(optim1)
}
