
########################################################################
## Functions to perform simulation scenarios presented in             ##
## "Martingale residual-based method to control for confounders       ##
##  measured only in a validation sample in time-to-event analysis",  ##
##  Burne R.M., Abrahamowicz M., Stat-in-Med 2016                     ##
##                                                                    ##
## Code: R.M. Burne                                                   ##
## Date edited: 03/05/2016                                            ##
########################################################################

# This code contains the functions SimFunction() and ImputeFunction(), called in Simulations.R

###############################
## MARTINGALE MODEL FUNCTION ##
###############################
## Arguments:
# dat: data frame with combined main data and VS
# var.obs: character vector with names of measured confounders eg c("C1","C2")
# var.unmeas: character vector with names of unmeasured confounders eg c("U1","U2")
# exposure: character vector length 1 with name of exposure eg "X"
# time: character vector length 1 with name of time variable eg "y"
# event: character vector length 1 with name of censoring indicator variable eg "e"
# models: takes a character vector argument, with options c("MRImp", "LogTImp", "PSC", "Standard")
# p: parameter set

SimFunction <- function(dat, var.obs, var.unmeas, exposure, time, event, p){
  require(survival)
  dat <- data.frame(dat)
  
  
  # Get log(y) for imputation
  dat$logT <- log(dat[,time])
  
  # Get martingale residual and PS, based on whether validation sample is external or not
  if(p$Scenario=="ExternalVS"){
    # Split data, since models are run separately for External VS
    datamain <- dat[!(complete.cases(dat)),]
    dataVS <- dat[complete.cases(dat),]
    rm(dat);gc()
    
    # Get martingale residuals in data
    formula <- paste("Surv(",time,",",event,")~",exposure,sep="")
    for (i in 1:length(var.obs)){
      formula<-paste(formula,"+",var.obs[i],sep="")
    }
      
    naivemodelVS <- coxph(eval(parse(text=formula)),data=dataVS)
    naivemodelmain <- coxph(eval(parse(text=formula)),data=datamain)
      
    # Obtain martingale residuals #
    dataVS$MR <- residuals(naivemodelVS,type="martingale")
    datamain$MR <- residuals(naivemodelmain,type="martingale")
    
    # Error-prone PS for PSC:
    formula<-paste(exposure,"~",sep="")
    for (i in 1:length(var.obs)){
      formula<-paste(formula,"+",var.obs[i],sep="")
    }
    PSep.modelVS <- glm(eval(parse(text=formula)), data=dataVS, family=binomial)
    PSep.modelmain <- glm(eval(parse(text=formula)), data=datamain, family=binomial)
    dataVS$PSep <- fitted.values(PSep.modelVS)
    datamain$PSep <- fitted.values(PSep.modelmain)
    
  } else { # (If scenario not ExternalVS)
    # Get martingale residuals in data
    formula <- paste("Surv(",time,",",event,")~",exposure,sep="")
    for (i in 1:length(var.obs)){
      formula<-paste(formula,"+",var.obs[i],sep="")
    }
      
    naivemodel <- coxph(eval(parse(text=formula)),data=dat)
      
    # Obtain martingale residuals #
    dat$MR <- residuals(naivemodel,type="martingale")
    
    # Get error-prone PS for PSC:
    formula<-paste(exposure,"~",sep="")
    for (i in 1:length(var.obs)){
      formula<-paste(formula,"+",var.obs[i],sep="")
    }
    temp.model <- glm(eval(parse(text=formula)), data=dat, family=binomial)
    dat$PSep <- fitted.values(temp.model);rm(temp.model)
    
    # Get main and VS
    datamain <- dat[!(complete.cases(dat)),]
    dataVS <- dat[complete.cases(dat),]
    rm(dat);gc()
  }
  
  # PSC: get corrected PS in datamain
  formula<-paste(exposure,"~",sep="")
  for (i in 1:length(var.obs)){
    formula<-paste(formula,"+",var.obs[i],sep="")
  }
  for (i in 1:length(var.unmeas)){
    formula<-paste(formula,"+",var.unmeas[i],sep="")
  }
  PSgs.model <- glm(eval(parse(text=formula)), data=dataVS, family=binomial)
  dataVS$PSgs <- fitted.values(PSgs.model)
    
  formula <- paste("PSgs ~ PSep + ", exposure,sep="")
  temp.model <- lm(eval(parse(text=formula)),data=dataVS)
  datamain$PSgs <- predict(temp.model,newdata=datamain);rm(temp.model)
  data.PSC <- rbind(datamain,dataVS)
  
  # Imputation steps
  temp <- datamain[,!names(datamain) %in% var.unmeas]
  U.MR <- sapply(var.unmeas,ImputeFunction,var.obs=var.obs,exposure=exposure,VS=dataVS,main=datamain,timevar=NA,event=event,imputation="martingale")
  data.MR <- rbind(cbind(temp,U.MR),dataVS)

  U.logT <- sapply(var.unmeas,ImputeFunction,var.obs=var.obs,exposure=exposure,VS=dataVS,main=datamain,timevar="logT",event=event,imputation="timevar")
  data.logT <- rbind(cbind(temp,U.logT),dataVS)
  
  ## FINAL MODELS AND OUTPUT
  to.return <- list()
  to.return$models <- list()
  
  formula.full<-paste("Surv(",time,",",event,")~",exposure,sep="")
  for (i in 1:length(var.obs)){
    formula.full<-paste(formula.full,"+",var.obs[i],sep="")
  }
  for (i in 1:length(var.unmeas)){
    formula.full<-paste(formula.full,"+",var.unmeas[i],sep="")
  }
  
  # Standard model:
  formula<-paste("Surv(",time,",",event,")~",exposure,sep="")
  for (i in 1:length(var.obs)){
    formula<-paste(formula,"+",var.obs[i],sep="")
  }
  model.Standard <- coxph(eval(parse(text=formula)), data=datamain)
  to.return$models$Standard <- summary(model.Standard)
  
  #PSC:
  formula<-paste("Surv(",time,",",event,") ~ ",exposure,"+ PSgs",sep="")
  model.PSC <- coxph(eval(parse(text=formula)),data=data.PSC)
  to.return$models$PSC <- summary(model.PSC)
  
  # Imputation with log T:
  model.LogTImp <- coxph(eval(parse(text=formula.full)),data=data.logT)
  to.return$models$LogTImp <- summary(model.LogTImp)
  
  # Imputation with MR:
  model.MRImp <- coxph(eval(parse(text=formula.full)),data=data.MR)
  to.return$models$MRImp <- summary(model.MRImp)

  # Assess surrogacy:
  formula <- paste("Surv(",time,",",event,")~",exposure,"+PSgs",sep="")
  m1 <- coxph(eval(parse(text=formula)),data=data.PSC)
  formula <- paste(formula,"+PSep",sep="")
  m2 <- coxph(eval(parse(text=formula)),data=data.PSC)
    
  surr.test <- anova(m1,m2)
  surr <- exp(m1$loglik[2] - m2$loglik[2])*100
  surr.test$SurrogacyPct <- c(NA,surr)
  to.return$Surrogacy <- surr.test
  
  return(to.return)
}


#######################
##IMPUTATION FUNCTION##
#######################
## Arguments:
# U.name: character vector of length 1 with name of the unmeasured confounder to impute
# var.obs: character vector with observed variables (to be used in imputation model)
# exposure: character vector length 1 with name of exposure eg "X"
# VS: validation sample (data.frame)
# main: main data (data.frame)
# timevar: If time variable is used in imputation model (used for LogT imputation)
# event: character vector length 1 with name of censoring indicator variable eg "e"
# imputation: options "martingale" or "timevar" - for either MR-based imputation or LogT imputation


ImputeFunction <- function(U.name,var.obs,exposure,VS,main,timevar,event,imputation){
  
  ## Create variable "type" - binomial / gaussian (continuous) ##
  type <- if(all(VS[,names(VS)==U.name] %in% c(0,1))){"binomial"} else {"gaussian"}
  
  ## Fit model dependent on type ## 
  # Define the formula: #
  formula <- switch(imputation, 
                    martingale = paste0(exposure, "+MR") ,          
                    timevar = paste0(exposure, "+",timevar,"+",event))
  
  for (i in 1:length(var.obs)){
    formula<-paste(formula,"+",var.obs[i],sep="")
  }
  
  
  imp.model <- switch(type,
                      gaussian = lm(eval(parse(text=paste(U.name,"~",formula))),data=VS),
                      binomial = glm(eval(parse(text=paste(U.name,"~",formula))),data=VS,family=binomial))
  
  coef <- imp.model$coef
  
  temp <- names(coef)[2:length(coef)]
  XmatVS <- as.matrix(cbind(1,subset(VS,select=temp)))
  Xmatmain <- as.matrix(cbind(1,subset(main,select=temp)))
  
  ## Impute step. ##
  u.imp <- switch(type,
                  binomial = apply(Xmatmain, 1,
                                   function(x){
                                     mu <- sum(coef*x)
                                     p <- exp(mu)/(1+exp(mu))
                                     new.u <- rbinom(1,1,p)
                                   }),
                  gaussian = apply(Xmatmain, 1,
                                   function(x){
                                     mu <- sum(coef*x)
                                     var <- summary(imp.model)$sigma^2*
                                       (1+t(as.matrix(x))%*%solve(t(XmatVS)%*%XmatVS)%*%as.matrix(x))    
                                     #Prediction variance.
                                     new.u <- rnorm(1,mean=mu,sd=sqrt(var))}))
  
  rm(Xmatmain);rm(XmatVS) 
  gc()
  return(u.imp)  
}
