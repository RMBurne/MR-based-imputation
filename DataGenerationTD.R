
########################################################################
## Data generation function for time-dependent simulation scenarios   ##
##  presented in                                                      ##
## "Martingale residual-based method to control for confounders       ##
##  measured only in a validation sample in time-to-event analysis",  ##
##  Burne R.M., Abrahamowicz M., Stat-in-Med 2016                     ##
##                                                                    ##
## Code: R.M. Burne                                                   ##
## Date edited: 03/05/2016                                            ##
########################################################################

# This code contains the function DataFunctionTD(), called in Simulation.R, 
#  and ConfoundersFunction() called in DataFunction

####################
## DATAGENERATION ##
####################
## Arguments: 
#  k: index to apply over
#  p: a list of parameters (contained in Parameters.R)

DataGeneration <- function(k, p){
  require(MASS);require(PermAlgo)
  system.time(datamain <- DataFunctionTD(p, n.subjects=10000, incidence.rate=0.1))
  system.time(datavs <- DataFunctionTD(p, n.subjects=1000,incidence.rate=0.2))
  datamain$U1 <- NA; datamain$U2 <- NA
  data <- rbind(datamain,datavs)
  rownames(data) <- 1:nrow(data)
  return(data)
}

####################
## DATAFUNCTIONTD ##
####################
## Arguments: 
#  p: a list of parameters (contained in Parameters.R)
#  n.subjects: number of subjects
#  incidence.rate: incidence rate for (exponentially) generated event times

DataFunctionTD <- function(p, n.subjects, incidence.rate){
  
  ## Generate confounders (using ConfoundersFunction())
  conf.mat <- do.call("rbind",lapply(1:n.subjects, ConfoundersFunction, p))
  
  ## Generate time of event & censoring
  # Assume that the marginal distribution of event times is the exponential distribution
  # Formula for survival at time t given an exponential dist: S(t)=exp{-lambda*t}, so lambda=-log{S(t)}/t,
  # where S(t)=1-incidence.rate at t
  
  if(p$InformativeCensoring==FALSE){
    eventtimes <- ceiling(rexp(n.subjects, -log(1-incidence.rate)/p$f.up))
    
    # Only administrative censoring, i.e. at the end on follow-up, if no event occurred
    censortimes <- rep(p$f.up, times=n.subjects) # Censoring times = end of follow-up
    
    covmat <- conf.mat[,colnames(conf.mat) %in% c("X","C1","C2","U1","U2")]
    
    dimnames(covmat) <- NULL
    
    gendata <- permalgorithm(numSubjects=n.subjects, maxTime=p$f.up, Xmat=covmat,
                             eventRandom=eventtimes, censorRandom=censortimes, 
                             betas=c(log(p$HR.X),log(p$HR.C1),log(p$HR.C2),log(p$HR.U1),log(p$HR.U2)),
                             XmatNames=c("X", "C1","C2","U1","U2"))
  }
  
  
  if(p$InformativeCensoring==TRUE){
    # Will separately generate event times (|X,C1,C2,U1,U2,R) and censoring times (|R)
    eventtimes <- ceiling(rexp(n.subjects, -log(1-incidence.rate)/p$f.up))
    censortimes <- ceiling(rexp(n.subjects, -log(1-incidence.rate)/p$f.up))
    endstudy <- rep(p$f.up, times=n.subjects) 
    
    covmat <- conf.mat[,colnames(conf.mat) %in% c("patid","X","C1","C2","U1","U2","R")]
    
    dimnames(covmat) <- NULL
    
    gendata <- permalgorithm(numSubjects=n.subjects, maxTime=p$f.up, Xmat=covmat,
                             eventRandom=eventtimes, censorRandom=endstudy, 
                             betas=c(0,log(p$HR.X),log(p$HR.C1),log(p$HR.C2),log(p$HR.U1),log(p$HR.U2),log(2)),
                             XmatNames=c("patid","X", "C1","C2","U1","U2","R"))
    
    censdata <- permalgorithm(numSubjects=n.subjects, maxTime=p$f.up, Xmat=covmat,
                              eventRandom=censortimes, censorRandom=endstudy, 
                              betas=c(0,0,0,0,0,0,log(2)),
                              XmatNames=c("patid","X", "C1","C2","U1","U2","R"))
    
    maxsurv <- gendata$Fup[c(diff(gendata$Id)!=0,TRUE)]
    maxcens <- censdata$Fup[c(diff(censdata$Id)!=0,TRUE)]
    maxtime <- apply(cbind(maxsurv,maxcens),1,min)
    
    gendata <- do.call("rbind",lapply(unique(gendata$patid),function(x){gendata[gendata$patid==x,][1:maxtime[x],]}))
    temp <- do.call("rbind",lapply(unique(censdata$patid),function(x){censdata[censdata$patid==x,][1:maxtime[x],]}))
    gendata$cens <- temp$Event
  }
  
  rownames(gendata) <- 1:nrow(gendata)
  return(gendata)
}


#########################
## CONFOUNDERSFUNCTION ##
#########################

## Arguments: 
#  k: index to apply over
#  p: a list of parameters (contained in Parameters.R)

ConfoundersFunction <- function(k, p){
  patid <- k
  rho <- 0.8
  # 1. generate covariates for 10 time intervals for each individual
  
  # R is a (not common) risk factor for the event, which will also be related to censoring time
  # (used for informative censoring)
  R <- rep(rbinom(n=1,p=0.3,size=1),p$f.up)
  
  # C1 is an observed confounder with population mean 0 and sd between 0.5 and 2.5
  sigmaC1 <- runif(n=1,0.5,2.5)
  muC1 <- rep(rnorm(n=1,0,1),p$f.up)
  H <- abs(outer(1:p$f.up, 1:p$f.up, "-"))
  V1 <-  sigmaC1^2* rho^H
  C1 <- mvrnorm(n=1,muC1,V1)
  
  # C2 is an observed confounder with population mean 0 and sd between 0.5 and 2.5
  sigmaC2 <- runif(n=1,0.5,2.5)
  muC2 <- rep(rnorm(n=1,0,1),p$f.up)
  V2 <-  sigmaC2^2* rho^H
  C2 <- mvrnorm(n=1,muC2,V2)
  
  # U1 is an unobserved binomial confounder (rare)
  pU1 <- runif(n=1,0.1,0.3)
  U1 <- rbinom(n=p$f.up,p=pU1,size=1)
  
  # U2 is an unobserved confounder with mean 0 and sd between 0.5 and 2.5
  sigmaU2 <- runif(n=1,0.5,2.5)
  muU2 <- rep(rnorm(n=1,0,1),p$f.up)
  VU2 <-  sigmaU2^2* rho^H
  U2 <- mvrnorm(n=1,muU2,VU2)
  
  # Odds ratios (moderate)
  mednorm <- 1.3
  medbin <- 2
  
  # Generate potential exposure for each time:
  PSgs <- exp(log(p$OR.int)+log(p$OR.C1)*C1+log(p$OR.C2)*C2+log(p$OR.U1)*U1+log(p$OR.U2)*U2)/(1+exp(log(p$OR.int)+log(p$OR.C1)*C1+log(p$OR.C2)*C2+log(p$OR.U1)*U1+log(p$OR.U2)*U2))
  xvec <- unlist(lapply(PSgs,function(PSgs){rbinom(1,1,PSgs)}))
  xvec <- as.vector(xvec)
  
  # For each exposure status, generate the number of times it is repeated
  X <- NULL
  j <- 1
  Texp <- rep(NA,p$f.up)
  while(j <= 10){
    exp <- xvec[j]
    t.exp <- as.integer(round(runif(1,1,5)*1*(exp==1) + runif(1,2,8)*1*(exp==0)))
    Texp[j] <- t.exp
    X <- c(X,rep(exp,t.exp))
    j <- j+t.exp
  }
  X <- X[1:10]
  
  cov <- cbind(rep(patid,p$f.up),X,C1,C2,U1,U2,R,1:p$f.up)
  colnames(cov) <- c("patid","X","C1","C2","U1","U2","R","time")
  return(cov)
}



