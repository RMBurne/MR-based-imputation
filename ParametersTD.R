
########################################################################
## Parameter sets for time-varying simulations presented in paper     ##
## "Martingale residual-based method to control for confounders       ##
##  measured only in a validation sample in time-to-event analysis",  ##
##  Burne R.M., Abrahamowicz M., Stat-in-Med 2016                     ##
##                                                                    ##
## Code: R.M. Burne                                                   ##
## Date edited: 02/05/2016                                            ##
########################################################################

parameters <- list()

parameters[[1]] <- list("HR.X"=1, "HR.C1"=1.3, "HR.C2"=1.3, "HR.U1"=2, "HR.U2"=1.3, 
                        "OR.int"=0.35, "OR.C1"=1.3, "OR.C2"=1.3, "OR.U1"=2, "OR.U2"=1.3, 
                        "f.up"=10, "InformativeCensoring"=FALSE)

parameters[[2]] <- list("HR.X"=1, "HR.C1"=1.3, "HR.C2"=1.3, "HR.U1"=2, "HR.U2"=1.3, 
                        "OR.int"=0.35, "OR.C1"=1.3, "OR.C2"=1.3, "OR.U1"=0.5, "OR.U2"=0.75, 
                        "f.up"=10, "InformativeCensoring"=FALSE)

parameters[[3]] <- list("HR.X"=1, "HR.C1"=1.3, "HR.C2"=1.3, "HR.U1"=3, "HR.U2"=1.5, 
                        "OR.int"=0.35, "OR.C1"=1.3, "OR.C2"=1.3, "OR.U1"=0.3, "OR.U2"=0.7, 
                        "f.up"=10, "InformativeCensoring"=FALSE)

parameters[[4]] <- list("HR.X"=1.5, "HR.C1"=1.3, "HR.C2"=1.3, "HR.U1"=2, "HR.U2"=1.3, 
                        "OR.int"=0.35, "OR.C1"=1.3, "OR.C2"=1.3, "OR.U1"=2, "OR.U2"=1.3, 
                        "f.up"=10, "InformativeCensoring"=FALSE)

parameters[[5]] <- list("HR.X"=1, "HR.C1"=1.3, "HR.C2"=1.3, "HR.U1"=2, "HR.U2"=1.3, 
                        "OR.int"=0.35, "OR.C1"=1.3, "OR.C2"=1.3, "OR.U1"=2, "OR.U2"=1.3, 
                        "f.up"=10, "InformativeCensoring"=TRUE)

parameters[[6]] <- list("HR.X"=1.5, "HR.C1"=1.3, "HR.C2"=1.3, "HR.U1"=2, "HR.U2"=1.3, 
                        "OR.int"=0.35, "OR.C1"=1.3, "OR.C2"=1.3, "OR.U1"=2, "OR.U2"=1.3, 
                        "f.up"=10, "InformativeCensoring"=TRUE)