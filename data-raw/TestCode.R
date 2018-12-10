set.seed(1234561234)

designFile <- read.csv("../Data Sets/AGPL Data/CBCExport/Q23_design.csv")
designFile$Input.Planning.Assistance. <- as.factor(designFile$Input.Planning.Assistance.)
designFile$Advertising. <- as.factor(designFile$Advertising.)
designFile$Proximity.to.data.points. <- as.factor(designFile$Proximity.to.data.points.)
designFile$Mobile.App.Availability. <- as.factor(designFile$Mobile.App.Availability.)
designFile$Price <- as.factor(designFile$Price)
designFile <- model.matrix(~-1 + ., designFile)
designFile <- as.data.frame(designFile)
# designFile <- do.call(rbind, replicate(50, designFile, simplify=FALSE))
names(designFile)[1] <- "version"
names(designFile)[2] <- "scenario"

testDemand <- generateTestData(500, design=designFile)

nParams <- ncol(designFile) - 3

testDemand$data <- as.data.frame(testDemand$data)
names(testDemand$data) <- c("respId", "version", "Task1", "Task2", "Task3", "Task4", "Task5", "Task6", "Task7", "Task8", "Task9", "Task10", "Task11", "Task12")
output = estimateGremlinsModel(testDemand$data, designFile, R = 5000, startingValues = testDemand$trueValues)

output = estimateGremlinsMNPModel(testDemand$data, designFile, R = 5000, startingValues = testDemand$trueValues)

cor(colMeans(output$segmentMembership - 1), testDemand$trueValues$segMembership - 1)

matplot(output$slope[,1,], type="l")
matplot(output$slope[,100,], type="l")
matplot(output$slope[,200,], type="l")
matplot(output$slope[,300,], type="l")

plot(colMeans(output$slope[,1,]), testDemand$trueValues$slope[1,])
abline(0,1)

plot(colMeans(output$slope[,100,]), testDemand$trueValues$slope[100,])
abline(0,1)

plot(colMeans(output$slope[,200,]), testDemand$trueValues$slope[200,])
abline(0,1)

plot(colMeans(output$slope[,300,]), testDemand$trueValues$slope[300,])
abline(0,1)




plot(colMeans(output$segmentMembership), testDemand$trueValues$segMembership)
plot(output$lambda[,2], type="l")
lines(output$lambda[,2], col="red")

matplot(output$slopeBar, type = "l")
plot(colMeans(output$slopeBar), testDemand$trueValues$slopeBar)
abline(0,1)

matplot(output$slopeCov, type="l")
plot(colMeans(output$slopeCov), testDemand$trueValues$slopeCov)
abline(0,1)

plot(density(rgamma(1000, shape=10, scale = 4)))


indBetas <- simplify2array(replicate(5000, testDemand$trueValues$slope, simplify=FALSE))
indBetas <- aperm(indBetas, c(3,1,2))
baseHitRates <- sapply(1:500, inSampleSimulation, design=designFile, data=testDemand$data)

indBetas <- output$slope
inSampleHitRates <- sapply(1:500, inSampleSimulation, design=designFile, data=testDemand$data)

plot(baseHitRates, inSampleHitRates)
cor(baseHitRates, inSampleHitRates)
plot(colMeans(output$segmentMembership), baseHitRates)
plot(colMeans(output$segmentMembership), inSampleHitRates)

#' Read in the camera dataset and run the model
#'
library(tidyverse)

cameraData <- read.csv(system.file("extdata", "CameraFullData.csv", package = "gremlins"))
cameraDesign <- read.csv(system.file("extdata", "CameraDesign.csv", package = "gremlins"))


#' We need to code of the dataset for use in the model.  The big challenge is figuring out
#' whether we want to generate a general function to do this with a general dataset.
#'

cameraDesign$brand <- as.factor(cameraDesign$brand)
cameraDesign$pixels <- as.factor(cameraDesign$pixels)
cameraDesign$zoom <- as.factor(cameraDesign$zoom)
cameraDesign$video <- as.factor(cameraDesign$video)
cameraDesign$swivel <- as.factor(cameraDesign$swivel)
cameraDesign$wifi <- as.factor(cameraDesign$wifi)

cameraDesign <- model.matrix(~-1 + ., cameraDesign)
cameraDesign <- as.tibble(cameraDesign)
#' Include the none option
#' Always concept 5 for all tasks for all version

versions <- unique(cameraDesign$version)
tasks <- unique(cameraDesign$scenario)
concepts <- unique(cameraDesign$concept)

none <- do.call(rbind, replicate(length(versions)*length(tasks), c(5, rep(0, 10), 1),  simplify = FALSE))
none <- as.tibble(cbind(rep(versions, each = length(tasks)), rep(tasks, length(versions)), none))
cameraDesign <- cameraDesign %>%
  add_column(none = 0)
names(none) <- names(cameraDesign)


cameraDesign <- cameraDesign %>%
  bind_rows(none) %>%
  arrange(version, scenario, concept)

# cameraDesign %>% print(n = Inf, width = Inf)

#' Run the model
nParams <- ncol(cameraDesign) - 3

output = estimateGremlinsModel(cameraData, cameraDesign, R = 10000)


startingValues = list(slope = output$slope[10000,,]*(pi/sqrt(6)),
                      slopeBar = output$slopeBar[10000,],
                      slopeCov = output$slopeCov[10000,,],
                      lambda = output$lambda[10000,],
                      segMembership = output$segmentMembership[10000,])

output2 <- estimateGremlinsMNPModel(cameraData, cameraDesign, R = 5000)

output2a <- output2
startingValues = list(slope = output2$slope[5000,,],
                      slopeBar = output2$slopeBar[5000,],
                      slopeCov = output2$slopeCov[5000,,],
                      lambda = output2$lambda[5000,],
                      segMembership = output2$segmentMembership[5000,])


output2 <- estimateGremlinsMNPModel(cameraData, cameraDesign, prior, startingValues = startingValues, R = 5000)

load("data-raw/output.RData")
probs <- colMeans(output2$segmentMembership[2500:5000,]) - 1
probs.logit <- colMeans(output$segmentMembership[7500:10000,]) - 1

jagsCode <- read_delim("data-raw/gremlin_results_JAGS.txt", ";")

plot(probs, jagsCode$probit)
cor(probs, jagsCode$probit)
matplot(output2$segmentMembership[,2], type="l")
output$segmentMembership[,2]

names(jagsCode)

#' Run the code with the 2 tasks held out.
#'
#'
set.seed(07282018)

holdoutTasks <- sample(1:16, 2, replace=FALSE)
estimationTasks <- (1:16)[-holdoutTasks]

holdOutTaskNames <- paste("Task", holdoutTasks, sep="")


cameraValidationData <- cameraData %>%
  select(xsurvnum, version, !!holdOutTaskNames)

cameraTrainingData <- cameraData %>%
  select(-!!holdOutTaskNames)

cameraValidationDesign <- cameraDesign %>%
  filter(scenario %in% holdoutTasks)

cameraTrainingDesign <- cameraDesign %>%
  filter(!(scenario %in% holdoutTasks))

output = estimateGremlinsModel(cameraTrainingData, cameraTrainingDesign, prior, R = 50000)

indBetas <- output$slope

calculateHitRate <- function(betas, lambdas, design, actualChoice) {
  xbeta <- (betas/lambdas) %*% t(design)
  exbeta <- exp(xbeta)
  probChoice <- exbeta/rowSums(exbeta)
  results <- apply(probChoice, 1, function(x) { sample(seq_along(x), 1, prob = x) })
  mean(results == actualChoice)
}

savedDraws <- 25001:50000

simulationValidation <- function(i, task) {
  lambdas <- output$lambda[cbind(savedDraws, output$segmentMembership[savedDraws, i])]
  tempDesign <- as.matrix(cameraValidationDesign[cameraValidationDesign$version == cameraValidationData$version[i] & cameraValidationDesign$scenario == task,-c(1:3)])
  calculateHitRate(indBetas[savedDraws, i, ], lambdas, tempDesign, cameraValidationData[i, paste0("Task", task)])
}


inSampleSimulation <- function(i, design = cameraTrainingDesign, data = cameraTrainingData) {
  hitrate <- double(14)
  lambdas <- output$lambda[cbind(savedDraws, output$segmentMembership[savedDraws, i])]
  for(j in estimationTasks) {
    tempDesign <- as.matrix(design[design[,1] == data[i,2] & design[,2] == j,-c(1:3)])
    hitrate[j] <- calculateHitRate(indBetas[savedDraws, i, ], lambdas, tempDesign, data[i,paste0("Task", j)])
  }
  return(mean(hitrate))
}


results <- sapply(1:501, simulationValidation, task = holdoutTasks[1])
results2 <- sapply(1:501, simulationValidation, task = holdoutTasks[2])


outOfSampleHitRates <- rowMeans(cbind(results, results2))

inSampleHitRates <- sapply(1:501, inSampleSimulation)

segMem <- colMeans(output$segmentMembership[savedDraws, ] - 1)
segLabel <- ifelse(segMem < 0.5, "Info Rich", "Info Poor")

packageOut <- cbind(segMem, outOfSampleHitRates, inSampleHitRates, segLabel)
write.csv(packageOut, "data-raw/packageOutput.csv")


plot(results, colMeans(output$segmentMembership)-1)
points(results2, colMeans(output$segmentMembership)-1, col="blue")
points(outOfSampleHitRates, colMeans(output$segmentMembership)-1, col="red")

plot(inSampleHitRates, outOfSampleHitRates)
abline(0,1)

plot(inSampleHitRates, colMeans(output$segmentMembership)-1)

plot(jagsCode$pm_outhit_logit, jagsCode$logit, col="green")
 cor(jagsCode$pm_outhit_probit, outOfSampleHitRates)
cor(jagsCode$pm_outhit_probit, jagsCode$pm_outhit_logit)

evalPoints <- seq(.1, 25, by = .1)
llike2 <- matrix(0, ncol = length(evalPoints), nrow = 501)
for(ind in 1:501) {
  llike2[ind, ] <- sapply(evalPoints, function(x) {gremlinsLogLikelihood(cameraTrainingData[ind,-c(1,2)], cameraTrainingDesign[cameraTrainingDesign[,1] == cameraTrainingData[ind, 2], -c(1:3)], colMeans(output$slope[savedDraws, ind, ]), x)})
}
