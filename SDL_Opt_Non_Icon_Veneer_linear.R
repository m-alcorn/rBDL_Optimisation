

library(doSNOW)
library(foreach)
library(SWTools)
library(lubridate)
library(tidyverse) 
library(animation)
library(zoo)


# Set Variables -----------------------------------------------------------
# getwd()
scenVars <- list( 
  target = 5.28,
  pFolder = "E:/ACT_Uplift/Models/PopulationV3",# Main working folder
  numScenarios = as.integer(16),                            # Number of times to run model
  scenFile = "Population_Function_Setup.csv",  # Contains full list of scenario values
  numParallels = as.integer(4),                             # Number of cores to run on
  startServer =  as.integer(9876) ,                         # Initial server port to run model on
  SourceVersion = "E:/Source/550"   ,           # Location of RiverSystems.Forms.exe
  sourceProjectName = "UMB_ACT_POPV3.rsproj",   # Source model to run
  ScenarioInputSet = "BDL.PopAndCrop"   ,                  # Name of scenario input set
  iScenario = "Upper_Murrumbidgee_River_System_Model", # Source scenario
  outputFolder = "outputs"   ,                  # dir to store model outputs (under pFolder)
  inputSetFile = "ACT_PopulationV3_InputSet.txt", # Name of Scenario Input Set File (Relative Path/Reload-on-run)
  runFlag = as.logical("TRUE"),
  summaryFilename = "BDL_Non_Icon_Take_OptV3_Veneer_AnnualMeans.csv",
  minPop = 1100000,
  maxPop = 1500000,
  min_iFactor = 1.0,
  max_iFactor = 3.0,
  startDate = "01/07/1890",
  endDate = "30/06/2009",
  startWindow = "01/07/1895",
  endWindow = "30/06/2009",
  baseURL = "http://localhost:9876"
  
)

source("E:/R/shSDL/SDL_Opt_Functions.R")

# Setup files and folders-----------------------------------------------------
# Check if all files/folders exist
# Flag success/failure
fileSetup(scenVars)

#-----------------------------------------------------------------------------
# Open the "Scenario File"
if(file.exists(paste0(scenVars[["pFolder"]][1],"/",scenVars[["scenFile"]][1]))){
  sf <-  read.csv(paste0(scenVars[["pFolder"]][1],"/",scenVars[["scenFile"]][1]))
} else {
  print(paste0("The base scenario file - ", 
               paste0(scenVars[["pFolder"]][1],"/",scenVars[["scenFile"]][1]), 
               " - does not exist. Please check for file before continuing..")
  )
}

# Step 1.------------------------------------------------------------------- 
# start measuring time taken for grid searching
start_time <- Sys.time()
#Run model 4 times (#TODO in parallel) to get starting curve

#Set up population points
# Using min, min*1.43, min*1.86, max
popRange <- scenVars[["max_iFactor"]][1]-scenVars[["min_iFactor"]][1]

startVals <- c(scenVars[["min_iFactor"]][1],
               scenVars[["min_iFactor"]][1]+popRange/4,
               scenVars[["max_iFactor"]][1]-popRange/4,
               scenVars[["max_iFactor"]][1]
)
print("Running model via Veneer Plugin 2 times.....")
# for(i in 1:length(startVals)){
# i <- 1

for(i in 1:2){
  # Set Function Variable (or write Input Set)
  
  writeInputSetFile3(iFactor = startVals[i],scenVars,sf)
  
  #Run Source
  VeneerRunSource(StartDate = scenVars[["startDate"]][1],
                  EndDate = scenVars[["endDate"]][1],
                  # InputSet = scenVars[["ScenarioInputSet"]][1],
                  baseURL = scenVars[["baseURL"]][1]
                  
  )
  #Get Result
  thisrow <- getResRowVeneer(scenVars,i,iFactor = startVals[i])
  if(i==1){
    resDF <- thisrow
  } else{
    resDF <- rbind(resDF,thisrow)
  }
  
  
}
print("Finished initialising residual function.")


# Step 2.-------------------------------------------------------------------

# Get the minimum value of the phi from the intial 4 runs

lastRoot <- getLinearRoot(scenVars,i)
x <- 3
i <- 3

#   Loop in
phi <- resDF$Residual[i-1]
print(paste0("    Current Residual Function (Phi)= ",phi))
# print(paste0("    Current Residual Function (Phi)= ",phi))
print(paste0("    Current Irrigation Factor Estimate = ",lastRoot))

while( phi > 0.00001){
  print(paste0("Starting run: ",i))
  print(paste0("    Current Residual Function (Phi)= ",phi))
  print(paste0("    Current Irrigation Factor Estimate = ",lastRoot))
  
  writeInputSetFile3(iFactor= lastRoot,scenVars,sf)
  
  #Run Source
  VeneerRunSource(StartDate = scenVars[["startDate"]][1],
                  EndDate = scenVars[["endDate"]][1],
                  # InputSet = scenVars[["ScenarioInputSet"]][1],
                  baseURL = scenVars[["baseURL"]][1]
                  
  )
  #Get Result
  thisrow <- getResRowVeneer(scenVars,i,iFactor = lastRoot)
  if(i==1){
    resDF <- thisrow
  } else{
    resDF <- rbind(resDF,thisrow)
  }
  lastRoot <- getLinearRoot(scenVars,x)  #getRoot(scenVars,x)
  if(lastRoot<scenVars[["min_iFactor"]][1]||lastRoot>scenVars[["max_iFactor"]][1]){
    lastRoot <- sample(scenVars[["min_iFactor"]][1]:scenVars[["max_iFactor"]][1],1)
    print("Warning.... New solution out of bounds, generating random solution.")
  }
    
  
  # 1. Open SUmmary CSV
  phi <- resDF$Residual[i]
  print(phi)

  i <- i + 1
  x <- x + 1
  
}
#Get the last result
thisrow <- getResRowVeneer(scenVars,i,iFactor = lastRoot)
if(i==1){
  resDF <- thisrow
} else{
  resDF <- rbind(resDF,thisrow)
}
# Timer ends here
end_time <- Sys.time()  # stop measuring time taken for grid searching

# Calculate and print the time taken
execution_time <- round(end_time - start_time,2)  # Calculate elapsed time

print(paste0("-------------SDL Non Icon Take Optimisation Complete-------------"))
print(paste0("Optimum Irrigation Factor = ",lastRoot))
print(paste0("Number of runs required = ",i))
print(paste("Execution Time: ", format(execution_time)))  # Print grid search execution time
print(paste0("Model Run Period: ",scenVars[["startDate"]][1],"-",scenVars[["endDate"]]))
print(as.matrix(scenVars))
print(paste0("---------------------------------------------------"))



d <- data.frame(Population=resDF$Population,Residual=resDF$Residual)

d$Population2 <- d$Population^2


qm <- lm( Residual ~ Population + Population2, data = d)

qmDat <- data.frame(x=seq(min(resDF$Population),max(resDF$Population),length.out = 100), y =1)
a <- qm[["coefficients"]][["Population2"]]
b <- qm[["coefficients"]][["Population"]]
c <- qm[["coefficients"]][["(Intercept)"]]

qmDat$y <- a*(qmDat$x^2) + b*qmDat$x + c
# summary(qm)
# 3. Predict the root (phi=0)
root <- -(qm[["coefficients"]][["Population"]])/(2*qm[["coefficients"]][["Population2"]])

# Plot
f <- function(x) a*(x^2)+b*x+c
gp <- d %>%
  ggplot( aes(x=Population, y=Residual)) +
  geom_function(fun = f)+
  geom_point(aes(alpha = 0.7),shape=21, color="black", fill="#69b3a2", size=4) +
  ggtitle("Phi")+
  geom_point(aes(x=lastRoot,y=phi,alpha=0.7),shape=21, color="black", fill="orange", size=6) +
  theme_minimal()+
  theme(legend.position = "none")
plot(gp)
ggsave(file = "ACT_SDL_Population.png",
       plot = gp,
       path = "E:/ACT_Uplift/Models/PopulationV3/outputs")





