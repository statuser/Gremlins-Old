# Conjoint response: dual response
# Des of attributes: backward
# Number: "just right"
# Price Range: 79-279

data <- read.csv(file = "RawCameraData.csv", header = TRUE, stringsAsFactors = FALSE)

# For respondents that picked the none option, there are two tasks for that choice.
# The first is X.0 and the second is X.1.  The X.0 task will always be a none option.
# The X.1 task with be the preference choice not accounting for the none.  These two choices have
# a correlated error term and should not be analyzed as an extra piece of information.
# We avoid this problem by dropping the preference data and analyzing it as if the only
# choice made as the none option.


dataSingleChoice <- subset(data, !grepl("\\.1$", data$scenario, perl=TRUE))

write.csv(dataSingleChoice, "CameraData.csv")


responseTime <- read.csv(file = "cellXTime.csv", header = TRUE, stringsAsFactors = FALSE)

cameraResponseTimes <- subset(responseTime, responseTime$ID %in% unique(data$xsurvnum))
dim(cameraResponseTimes)
length(unique(data$xsurvnum))
write.csv(cameraResponseTimes, "cameraTime.csv")


#Re-frame the data into Sawtooth Style Format.

nTotalTasks <- nrow(dataSingleChoice)/5 #5 choices per set including the none.
concepts <- rep(1:5, nTotalTasks)

dataSingleChoice <- as.tibble(dataSingleChoice)
dataSingleChoice <- dataSingleChoice %>% add_column(concept = concepts, .before = "y1")

cameraResponses <- dataSingleChoice %>%
  select(xsurvnum, version, scenario, concept, y1) %>%
  mutate(scenario, scenario = paste("Task", scenario, sep ="")) %>%
  filter(y1 == 1) %>%
  select(-y1) %>%
  spread(scenario, concept) %>%
  select(xsurvnum, version, Task1, Task2, Task3, Task4, Task5, Task6, Task7,
         Task8, Task9, Task10, Task11, Task12, Task13, Task14, Task15, Task16)

cameraDesign <- dataSingleChoice %>%
  select(-y1) %>%
  mutate(brand = brand_1 + 2*brand_2 + 3*brand_3 + 4*brand_4) %>%
  select(-brand_1, -brand_2, -brand_3, -brand_4) %>%
  mutate(pixels = pixels_2 + 1, zoom = zoom_2 + 1, video = video_2 + 1, swivel=swivel_2 + 1, wifi=wifi_2 + 1) %>%
  select(-pixels_2, -zoom_2, -video_2, -swivel_2, -wifi_2) %>%
  filter(concept != 5) %>%
  select(version, scenario, concept, brand, pixels, zoom, video, swivel, wifi, price_lin) %>%
  mutate(price_lin = round(exp(price_lin)/100, digits=2)) %>%
  filter(!duplicated(.)) %>%
  arrange(version, scenario, concept)


write_csv(cameraDesign, "CameraDesign.csv")

write_csv(cameraResponses, "CameraFullData.csv")


#' Split the data into a test and training data set.
#'
#' This is basically just peeling off a set number of tasks
#' for each respondent.  Because there are ten blocks to the design
#' It should be OK to just peel off 4 of the 16 random tasks in a fixed
#' way.  The previous model use a random set of tasks for each respondent
#' This makes the data files much more complicated.


set.seed(07282018)

holdoutTasks <- sample(1:16, 4, replace=FALSE)

holdOutTaskNames <- paste("Task", holdoutTasks, sep="")


cameraValidationData <- cameraResponses %>%
  select(xsurvnum, version, !!holdOutTaskNames)

cameraTrainingData <- cameraResponses %>%
  select(-!!holdOutTaskNames)

cameraValidationDesign <- cameraDesign %>%
  filter(scenario %in% holdoutTasks)

cameraTrainingDesign <- cameraDesign %>%
  filter(!(scenario %in% holdoutTasks))

write_csv(cameraValidationData, "CameraValidationData.csv")
write_csv(cameraValidationDesign, "CameraValidationDesign.csv")
write_csv(cameraTrainingData, "CameraTrainingData.csv")
write_csv(cameraTrainingDesign, "CameraTrainingDesign.csv")
