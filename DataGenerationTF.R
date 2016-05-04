
########################################################################
## Data generation function for simulation scenarios presented in     ##
## "Martingale residual-based method to control for confounders       ##
##  measured only in a validation sample in time-to-event analysis",  ##
##  Burne R.M., Abrahamowicz M., Stat-in-Med 2016                     ##
##                                                                    ##
## Code: R.M. Burne                                                   ##
## Date edited: 03/05/2016                                            ##
########################################################################

# This code contains the function DataFunction(), called in Simulation.R 

##################
## DATAFUNCTION ##
##################
## Arguments: 
#  k: index to apply over
#  p: a list of parameters (contained in Parameters.R)

DataFunction <- function(k,p){
  C1 <- rnorm(p$N,0,1)
  C2 <- rnorm(p$N,0,1)
  U1 <- rnorm(p$N,0,1)
  
  if(p$Scenario=="MAR"){
    # Induce correlation between C1 (standardized age) and U2
    beta0 <- log(0.2/0.8); beta1 <- log(1.3)
    p.U2 <- exp(beta0+beta1*C1)/(1+exp(beta0+beta1*C1))
    U2 <- rbinom(p$N,1,p.U2)
  } else {
    U2 <- rbinom(p$N,1,0.6)
  }
  
  # Exposure
  PSgs <- (1 + exp(-log(p$OR.int) - log(p$OR.C1)*C1 - log(p$OR.C2)*C2 - log(p$OR.U1)*U1 - log(p$OR.U2)*U2))^(-1)
  X <- unlist(lapply(PSgs,function(PSgs){rbinom(1,1,PSgs)}))
  X <- as.vector(X)
  
  # Event time t
  if(p$Scenario=="ExternalVS"){
    unif <- runif(p$N,0,1)
    t <- rep(NA,p$N)
    index.main <- 1:(p$N - p$val.size)
    index.vs <- (p$N - p$val.size + 1):p$N
    
    t[index.main] <- apply(as.matrix(cbind(unif,X,C1,C2,U1,U2))[index.main,],1,function(x){
      unif <- x[1];X <- x[2];C1<- x[3];C2 <- x[4]; U1 <- x[5]; U2 <- x[6]
      -log(unif)/(p$rate*exp(log(p$HR.X)*X + log(p$HR.C1)*C1 + log(p$HR.C2)*C2 + log(p$HR.U1)*U1 + log(p$HR.U2)*U2))})
    
    t[index.vs] <- apply(as.matrix(cbind(unif,X,C1,C2,U1,U2))[index.vs,],1,function(x){
      unif <- x[1];X <- x[2];C1<- x[3];C2 <- x[4]; U1 <- x[5]; U2 <- x[6]
      -log(unif)/(2*p$rate*exp(log(p$HR.X)*X + log(p$HR.C1)*C1 + log(p$HR.C2)*C2 + log(p$HR.U1)*U1 + log(p$HR.U2)*U2))})
  } else {
    unif <- runif(p$N,0,1)
    t <- apply(as.matrix(cbind(unif,X,C1,C2,U1,U2)),1,function(x){
      unif <- x[1];X <- x[2];C1<- x[3];C2 <- x[4]; U1 <- x[5]; U2 <- x[6]
      -log(unif)/(p$rate*exp(log(p$HR.X)*X + log(p$HR.C1)*C1 + log(p$HR.C2)*C2 + log(p$HR.U1)*U1 + log(p$HR.U2)*U2))})
  }
  
  ## Censoring y: dependent on scenario
  if(p$Scenario == "RandomAndAdminCens"){
    randcens <- rexp(length(t),1/(mean(t)))
    y1 <- apply(cbind(t,randcens),1,function(x){min(x[1],x[2])})
    e1 <- 1*(y1==t)
    # Administrative censoring time after 10% events have been *observed*
    admincens <- quantile(y1[e1==1],.1/mean(e1))
    y <- sapply(y1,function(x){min(x,admincens)})
  } else {
    if(p$Scenario=="RandomCens"){
      # Fully random censoring
      randcens <- rweibull(length(t),scale=(1/p$rate),shape=p$shape)
      y <- apply(cbind(t,randcens),1,function(x){min(x[1],x[2])}) 
      } else {
        # Everything else - admin censoring at 10th percentile
        cutoff <- quantile(t,.1)
        y <- sapply(t,function(x){min(x,cutoff)})
      }
    }
  
  ## Event e
  e <- 1*(y==t)
  
  ## Combine data
  data <- data.frame(C1, C2, U1, U2, X, t, y, e)
  
  # Inclusion in VS
  if(p$Scenario=="MAR"){
    # for MAR scenarios P(incl in VS) depends directly on C1
    p.vs <- exp(data$C1)/(1+exp(data$C1))
    VS.ind <- sample(1:nrow(data),p$val.size,replace=FALSE,prob=p.vs)
    VS <- data[VS.ind,]
    main <- data[-VS.ind,]
  }
  if(p$Scenario=="MNAR"){
    # for MNAR scenarios P(incl in VS) depends directly on U1
    p.vs <- exp(data$U1)/(1+exp(data$U1))
    VS.ind <- sample(1:nrow(data),p$val.size,replace=FALSE,prob=p.vs)
    VS <- data[VS.ind,]
    main <- data[-VS.ind,]
  }
  if(p$Scenario=="ExternalVS"){
    VS <- data[index.vs,]
    main <- data[index.main,]
  }
  if(p$Scenario %in% c("Base","RandomAndAdminCens","RandomCens")){
    sampleVS <- sample(1:p$N,p$val.size,replace=FALSE)
    VS <- data[sampleVS,]
    main <- data[-sampleVS,]
  }
  
  main$U1 <- rep(NA,nrow(main))
  main$U2 <- rep(NA,nrow(main))
  data.miss <- rbind(main,VS)
  rownames(data.miss) <- 1:p$N
  return(data.miss)
}
