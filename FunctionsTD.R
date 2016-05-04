
########################################################################
## Functions to perform time-varying simulations presented in         ##
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

SimFunction <- function(dat, var.obs, var.unmeas, exposure, start, stop, event){
  require(survival)
  datamain <- dat[!(complete.cases(dat)),]
  dataVS <- dat[complete.cases(dat),]
  rm(dat);gc()
  
  # Get martingale residual, and imputed U
  formula <- paste("Surv(",start,",",stop,",",event,")~",exposure,sep="")
  for (i in 1:length(var.obs)){
    formula<-paste(formula,"+",var.obs[i],sep="")
  }
  
  naivemodelVS <- coxph(eval(parse(text=formula)),data=dataVS)
  naivemodelmain <- coxph(eval(parse(text=formula)),data=datamain)
  
  dataVS$MR <- residuals(naivemodelVS, type="martingale")
  datamain$MR <- residuals(naivemodelmain, type="martingale")
  
  U.MI <- sapply(var.unmeas, ImputeFunction, var.obs=var.obs, exposure=exposure, VS=dataVS, main=datamain)
  
  # PSC: get gold standard and error-prone PS, and impute gold standard PS in datamain
  # PSep:
  formula <- paste(exposure,"~",sep=""); 
  for (i in 1:length(var.obs)){ 
    formula<-paste(formula,"+",var.obs[i],sep="")  
  }
  PSep.model.VS <- glm(eval(parse(text=formula)), data=dataVS, family=binomial)
  dataVS$PSep <- fitted.values(PSep.model.VS)
  PSep.model.main <- glm(eval(parse(text=formula)), data=datamain, family=binomial)
  datamain$PSep <- fitted.values(PSep.model.main)
  
  # PSgs:
  for (i in 1:length(var.unmeas)){ 
    formula <- paste(formula,"+",var.unmeas[i],sep="")  
  }
  PSgs.model.VS <- glm(eval(parse(text=formula)), data=dataVS, family=binomial)
  dataVS$PSgs <- fitted.values(PSgs.model.VS)
  
  # impute PS in main:
  formula <- paste("PSgs ~ PSep + ", exposure,sep="")
  pred.PSgs.model <- lm(eval(parse(text=formula)),data=dataVS)
  datamain$PSgs <- predict(pred.PSgs.model, type="response", newdata=datamain)
  
  ## Combine imputed values and create dataset
  mainimp <- cbind(datamain[,!names(datamain) %in% var.unmeas],U.MI)
  rm(U.MI)  
  rm(datamain);gc()
  
  ## FINAL MODELS
  # formulae:
  biased.formula <- paste0("Surv(",start,",",stop,",",event,")~",exposure); for (i in 1:length(var.obs)){biased.formula<-paste0(biased.formula,"+",var.obs[i])}
  full.formula <- biased.formula; for(i in 1:length(var.unmeas)){full.formula<-paste0(full.formula,"+",var.unmeas[i])}
  PSC.formula <- paste0("Surv(",start,",",stop,",",event,") ~ ",exposure,"+PSgs")
  
  to.return <- list()
  
  # Models to return:
  model.Standard <- coxph(eval(parse(text=biased.formula)), data=mainimp) 
  to.return$Standard <- summary(model.Standard)
  
  model.PSC <- coxph(eval(parse(text=PSC.formula)),data=mainimp)
  to.return$PSC <- summary(model.PSC)
  
  model.MRImp <- coxph(eval(parse(text=full.formula)), data=mainimp) 
  to.return$MR.Imp <- summary(model.MRImp)  
  
  rm(dataVS);rm(mainimp)
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

ImputeFunction <- function(U.name, var.obs, exposure, VS, main){
  
  ## Create variable "type" - binomial / gaussian (continuous) ##
  #### At present we only deal with binary/normal data
  type <- if(all(VS[,names(VS)==U.name] %in% c(0,1))){"binomial"} else {"gaussian"}
  
  ## Fit model dependent on type ##
  
  # Define the formula: #
  formula <- paste(exposure, "+MR", sep="")
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
