
########################################################################
## Simulation file for time-varying scenarios presented in            ##
## "Martingale residual-based method to control for confounders       ##
##  measured only in a validation sample in time-to-event analysis",  ##
##  Burne R.M., Abrahamowicz M., Stat-in-Med 2016                     ##
##                                                                    ##
## Code: R.M. Burne                                                   ##
## Date edited: 03/05/2016                                            ##
########################################################################

# This code uses functions and parameters in DataGenerationTD.R, FunctionsTD.R, ParametersTD.R, 
#  runs simulations and obtains results


## Simulation set-up
# Source functions and parameters
source("ParametersTD.R")
source("DataGenerationTD.R")
source("FunctionsTD.R")

# Load necessary packages:
library(survival)
library(doParallel)
library(doRNG)

# Use multiple cores (optional, but reduces time)
max.cores <- 4
cl <- makeCluster(max.cores)
registerDoParallel(cl)

n.reps <- 1000
simname <- "TimeVaryingSim"

for(index in 1:length(parameters)){
  system.time(datalist <- foreach(k=1:n.reps) %dorng% {DataGeneration(k,parameters[[index]]) })
  save(datalist,file=paste0("DataTD/Datalist",index,simname,".RData"))
  # Note: if applying in parallel, this may use a large amount of memory. 
  # It may be necessary to split the datalist if using many cores / a computer with limited memory.
  system.time(output <- foreach(k=1:n.reps) %dorng% {SimFunction(dat=datalist[[k]],var.obs=c("C1","C2"),var.unmeas=c("U1","U2"),
                                                                 exposure="X",start="Start",stop="Stop",event="Event") })
  save(output,file=paste0("OutputTD/Output",index,simname,".RData"))
} 


## Combine results (Results for Table 3)
bias.tab <- matrix(nrow=length(parameters),ncol=3)
rel.SD <- matrix(nrow=length(parameters),ncol=2)
rel.RMSE <- matrix(nrow=length(parameters),ncol=2)

for(index in 1:length(parameters)){
  load(paste0("OutputTD/Output",index,simname,".RData"))
  models <- names(output[[1]])
  coefficients <- lapply(models,function(name){unlist(lapply(lapply(lapply(output,"[[",name),"[[","coefficients"),"[","X","coef"))})
  names(coefficients) <- models
  
  bias <- lapply(coefficients,function(x){mean(x - log(parameters[[index]]$HR.X))})
  bias.tab[index,] <- unlist(bias)
  SD <- lapply(coefficients,sd)
  rel.SD[index,] <- unlist(lapply(SD,function(x){x/SD$MR.Imp}))[1:2]
  RMSE <- lapply(models,function(name){sqrt(bias[[name]]^2 + SD[[name]]^2)}); names(RMSE) <- models
  rel.RMSE[index,] <- unlist(lapply(RMSE,function(x){x/RMSE$MR.Imp}))[1:2]
}