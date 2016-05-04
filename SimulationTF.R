
########################################################################
## Simulation file for scenarios presented in                         ##
## "Martingale residual-based method to control for confounders       ##
##  measured only in a validation sample in time-to-event analysis",  ##
##  Burne R.M., Abrahamowicz M., Stat-in-Med 2016                     ##
##                                                                    ##
## Code: R.M. Burne                                                   ##
## Date edited: 03/05/2016                                            ##
########################################################################

# This code uses functions and parameters in DataGenerationTF.R, FunctionsTF.R, ParametersTF.R, 
#  runs simulations and obtains results

## Simulation set-up
# Source functions and parameters
source("ParametersTF.R")
source("DataGenerationTF.R")
source("FunctionsTF.R")

# Load necessary packages:
library(survival)
library(doParallel)
library(doRNG)

# Use multiple cores (optional, but reduces time)
max.cores <- 4
cl <- makeCluster(max.cores)
registerDoParallel(cl)

n.reps <- 1000
simname <- "TimeFixedSim"

for(index in 1:length(parameters)){
  datalist <- foreach(k=1:n.reps) %dorng% {DataFunction(k,parameters[[index]]) }
  save(datalist,file=paste0("DataTF/Datalist",index,simname,".RData"))
  
  system.time(output <- foreach(k=1:n.reps) %dorng% {SimFunction(dat=datalist[[k]],var.obs=c("C1","C2"),var.unmeas=c("U1","U2"),
                                                     exposure="X",time="y",event="e",p=parameters[[index]]) })
  save(output,file=paste0("OutputTF/Output",index,simname,".RData"))
} 


## Combine results (Results for Table 1)
bias.tab <- matrix(nrow=length(parameters),ncol=4)
rel.SD <- matrix(nrow=length(parameters),ncol=3)
rel.RMSE <- matrix(nrow=length(parameters),ncol=3)

for(index in 1:length(parameters)){
  load(paste0("OutputTF/Output",index,simname,".RData"))
  models <- names(output[[1]]$models)
  coefficients <- lapply(models,function(name){unlist(lapply(lapply(lapply(lapply(output,"[[","models"),"[[",name),"[[","coefficients"),"[","X","coef"))})
  names(coefficients) <- models
  
  bias <- lapply(coefficients,function(x){mean(x - log(parameters[[index]]$HR.X))})
  bias.tab[index,] <- unlist(bias)
  SD <- lapply(coefficients,sd)
  rel.SD[index,] <- unlist(lapply(SD,function(x){x/SD$MRImp}))[1:3]
  RMSE <- lapply(models,function(name){sqrt(bias[[name]]^2 + SD[[name]]^2)}); names(RMSE) <- models
  rel.RMSE[index,] <- unlist(lapply(RMSE,function(x){x/RMSE$MRImp}))[1:3]
}




