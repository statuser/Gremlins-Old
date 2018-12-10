source("utils.R")

load("data-raw/output.Rdata")
tmp = flattenOutput(output)

camera_covariates  = c('Canon', 'Sony', 'Nikon', 'Panasonic', 'MegaPixels', 'Zoom', 'HDVideo', 'Swivel', 'Wifi', 'Price','NA')
# Default holdout task
holdoutTask <- as.data.frame(rbind(c(0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1.2),
                                  c(0, 0, 0, 1, 0, 0, 1, 0, 0, 1,1.25),
                                  c(0, 0, 1, 0, 0, 1, 0, 0, 1, 0,1.30),
                                  c(0, 1, 0, 0, 0, 1, 1, 1, 1, 0,1.4 ),
                                  c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)))

names(holdoutTask) = camera_covariates
param_summary = getParameterSummary(output,camera_covariates)
market_share = marketShareSimulation(output,brands,holdoutTask=holdoutTask)
runPriceOptimization(output, holdoutTask, 0, 1.75, .25, 4, 11)


