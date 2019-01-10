# Attributes:
#   
#   -- Input Planning Assistance --
# 1. Provides planning assistance for<br>seed variety/trait selection
# 2. Provides planning assistance for<br>seed variety/trait selection<br> <u>and</u> fungicides
# 
# -- Advertising --
# 1. No ad sponsor
# 2. Image 1
# 3. Image 2
# 4. Image 3
# 
# -- Proximity to Data Points --
# 1.Recommendations based on farms within 15-100 mile radius of yours 
# 2. Recommendations based on farms within 30-100 mile radius of yours
# 3. Recommendations based on farms within 60-100 mile radius of yours
# 4. Recommendations based on farms within 100 mile radius of yours

# -- Mobile App availability --
# 1. Access system via PC
# 2. Access system via pC and mobile app

# -- Price --
# 1. 1200 per year
# 2. 1800 per year
# 3. 2400 per year
# 4. 3000 per year
# 5. 4800 per year






# Read in Data #
inputData = read.csv(file = "./CBCExport/Q23_responses.csv", header = TRUE, stringsAsFactors = FALSE)

inputDesign = read.csv(file = "./CBCExport/Q23_design.csv", header = TRUE, stringsAsFactors = FALSE)

# Number of units, tasks, concepts #
numUnits = nrow(inputData)
numTasks = length(unique(inputDesign$Task))
numConcepts = length(unique(inputDesign$Concept))


# Rename column headers #
names(inputDesign) = c('Version','Task','Concept','input_planning_assistance','advertising','proximity_to_data',
                       'mobile_app_availability','price')

names(inputData)[1:2] = c('sys_SequentialRespondentNumber','sys_CBCVersion_CBC')

# Drop Fixed Columns #

dropCols = c('Q23_Fixed1','Q23_Fixed1_none','Q23_Fixed2','Q23_Fixed2_none') ## Delete non-random columns?
inputData = inputData[,!names(inputData) %in% dropCols]



levels(as.factor(inputDesign$input_planning_assistance))
levels(as.factor(inputDesign$advertising))
levels(as.factor(inputDesign$proximity_to_data))
levels(as.factor(inputDesign$mobile_app_availability))


## Code cost per assessment ##
price = c(1200,1800,2400,3000,4800)



#Linear Variables #
inputDesign$price[inputDesign$price != 0] = price[inputDesign$price]


# Lengths for matrix # 
l1 = length(levels(as.factor(inputDesign$input_planning_assistance)))
l2 = length(levels(as.factor(inputDesign$advertising)))
l3 = length(levels(as.factor(inputDesign$proximity_to_data)))
l4 = length(levels(as.factor(inputDesign$mobile_app_availability)))



# Factor Variables - Code the Design using dummy#
inputDesign$input_planning_assistance = factor(inputDesign$input_planning_assistance, levels=as.character(1:max(inputDesign$input_planning_assistance)))
inputDesign$advertising = factor(inputDesign$advertising, levels=as.character(1:max(inputDesign$advertising)))
inputDesign$proximity_to_data = factor(inputDesign$proximity_to_data, levels=as.character(1:max(inputDesign$proximity_to_data)))
inputDesign$mobile_app_availability = factor(inputDesign$mobile_app_availability, levels=as.character(1:max(inputDesign$mobile_app_availability)))


holdOutTasks = 2
p = 4 + l1-1 + l2-1 + l3-1 + l4-1 + 1 + 1
# 1 + 1 represents the price column and the none option

#p = 15 #; Number of columns in the coded design + the number of control columns
trainingSet = matrix(0, nrow=numUnits*(numTasks - holdOutTasks)*(numConcepts+1), ncol=p)
validationSet = matrix(0, nrow=numUnits*holdOutTasks*(numConcepts+1), ncol=p)






for(unit in 1:numUnits) {
  respNum = inputData[unit, ]$sys_SequentialRespondentNumber
  version = inputData[unit, ]$sys_CBCVersion_CBC
  
  if(is.null(version)) {
    error(resp, "version = NULL")
  }
  

  respDesign = inputDesign[inputDesign$Version %in% version, -c(1, 2, 3)]  # remove "Version", 'Task', 'Concept' 
  tempDesign = model.matrix(~., model.frame(~., respDesign, na.action=na.pass))
  tempDesign[is.na(tempDesign)] = 0
  tempDesign = tempDesign[,-1]
  
  #; Include the none option
  concepts = rep(1:numConcepts, numTasks)
  tasks = rep(1:numTasks, 1, each = numConcepts)
  tempDesign = cbind(tasks, concepts, tempDesign, 0)
  noneOption = matrix(0, nrow=numTasks, ncol=ncol(tempDesign))
  noneOption[,ncol(noneOption)] = 1
  noneOption[,1] = 1:numTasks
  noneOption[,2] = numConcepts + 1
  
  tempDesign = rbind(tempDesign, noneOption)
  
  tempDesign = tempDesign[order(tempDesign[,1], tempDesign[,2]),]
  
  #; These tasks are a duak response none meaning that the outside option is captured in a second queston
  responses = inputData[unit, -c(1,2)]
  responses = matrix(responses, ncol=2, byrow=TRUE)
  #responses[responses[,2] == 2,1] = 5
  responses[responses[,2] == 2,1] = numConcepts + 1
  responses = unlist(responses[,1])
  
  indexes = cbind(1:numTasks, responses)
  y = matrix(0, nrow=numTasks, ncol=numConcepts+1)
  y[indexes] = 1
  y = matrix(t(y), ncol=1)
  
  tempDesign = cbind(respNum, y, tempDesign)
  tempDesign = tempDesign[,c(1, 3, 4, 2, 5:ncol(tempDesign))]
  
  tasksToHoldOut = sample(1:numTasks, holdOutTasks)
  
  trainingSample = tempDesign[!tempDesign[,2] %in% tasksToHoldOut,]
  holdOutSample = tempDesign[tempDesign[,2] %in% tasksToHoldOut,]
  
  tIndex = (unit - 1)*nrow(trainingSample) + 1
  vIndex = (unit - 1)*nrow(holdOutSample) + 1
  
  trainingSet[tIndex:(tIndex+nrow(trainingSample)-1), ] = trainingSample
  validationSet[vIndex:(vIndex+nrow(holdOutSample)-1), ] = holdOutSample
  
}


write.csv(trainingSet, "AGPL-Training.csv")
write.csv(validationSet, "AGPL-Validation.csv")

