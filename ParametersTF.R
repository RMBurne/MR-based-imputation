
########################################################################
## Parameter sets for simulation scenarios presented in paper         ##
## "Martingale residual-based method to control for confounders       ##
##  measured only in a validation sample in time-to-event analysis",  ##
##  Burne R.M., Abrahamowicz M., Stat-in-Med 2016                     ##
##                                                                    ##
## Code: R.M. Burne                                                   ##
## Date edited: 02/05/2016                                            ##
########################################################################

parameters <- list()

parameters[[1]] <- list("N"=10000, "val.size"=1000, "HR.X"=1, "HR.C1"=1.3, "HR.C2"=1.3, "HR.U1"=1.3, "HR.U2"=2,
                        "OR.int"=0.25, "OR.C1"=1.3, "OR.C2"=1.3, "OR.U1"=1.3, "OR.U2"=2, "rate"=0.01, "Scenario"="Base")

# Change True HR
parameters[[2]] <- list("N"=10000, "val.size"=1000, "HR.X"=1.2, "HR.C1"=1.3, "HR.C2"=1.3, "HR.U1"=1.3, "HR.U2"=2,
                        "OR.int"=0.25, "OR.C1"=1.3, "OR.C2"=1.3, "OR.U1"=1.3, "OR.U2"=2, "rate"=0.01, "Scenario"="Base")
parameters[[3]] <- list("N"=10000, "val.size"=1000, "HR.X"=1.5, "HR.C1"=1.3, "HR.C2"=1.3, "HR.U1"=1.3, "HR.U2"=2,
                        "OR.int"=0.25, "OR.C1"=1.3, "OR.C2"=1.3, "OR.U1"=1.3, "OR.U2"=2, "rate"=0.01, "Scenario"="Base")

# Change strength of unmeasured confounding
parameters[[4]] <- list("N"=10000, "val.size"=1000, "HR.X"=1, "HR.C1"=1.3, "HR.C2"=1.3, "HR.U1"=1.5, "HR.U2"=3,
                        "OR.int"=0.25, "OR.C1"=1.3, "OR.C2"=1.3, "OR.U1"=1.5, "OR.U2"=3, "rate"=0.01, "Scenario"="Base")
parameters[[5]] <- list("N"=10000, "val.size"=1000, "HR.X"=1, "HR.C1"=1.3, "HR.C2"=1.3, "HR.U1"=1.1, "HR.U2"=1.5,
                        "OR.int"=0.25, "OR.C1"=1.3, "OR.C2"=1.3, "OR.U1"=1.1, "OR.U2"=1.5, "rate"=0.01, "Scenario"="Base")

# Violation of surrogacy
parameters[[6]] <- list("N"=10000, "val.size"=1000, "HR.X"=1, "HR.C1"=1.3, "HR.C2"=1.3, "HR.U1"=1.3, "HR.U2"=2,
                        "OR.int"=0.25, "OR.C1"=1.3, "OR.C2"=1.3, "OR.U1"=0.75, "OR.U2"=0.5, "rate"=0.01, "Scenario"="Base")
parameters[[7]] <- list("N"=10000, "val.size"=1000, "HR.X"=1, "HR.C1"=1.3, "HR.C2"=1.3, "HR.U1"=1.1, "HR.U2"=1.5,
                        "OR.int"=0.25, "OR.C1"=1.3, "OR.C2"=1.3, "OR.U1"=0.9, "OR.U2"=0.7, "rate"=0.01, "Scenario"="Base")
parameters[[8]] <- list("N"=10000, "val.size"=1000, "HR.X"=1, "HR.C1"=1.3, "HR.C2"=1.3, "HR.U1"=1.3, "HR.U2"=2,
                        "OR.int"=0.25, "OR.C1"=1.3, "OR.C2"=0.75, "OR.U1"=1.3, "OR.U2"=0.5, "rate"=0.01, "Scenario"="Base")
parameters[[9]] <- list("N"=10000, "val.size"=1000, "HR.X"=1, "HR.C1"=1.3, "HR.C2"=1.1, "HR.U1"=1.3, "HR.U2"=1.5,
                        "OR.int"=0.25, "OR.C1"=1.3, "OR.C2"=0.9, "OR.U1"=1.3, "OR.U2"=0.7, "rate"=0.01, "Scenario"="Base")
parameters[[10]] <- list("N"=10000, "val.size"=1000, "HR.X"=1, "HR.C1"=1.3, "HR.C2"=1.3, "HR.U1"=1.3, "HR.U2"=2,
                        "OR.int"=0.25, "OR.C1"=1.3, "OR.C2"=0.75, "OR.U1"=1.3, "OR.U2"=2, "rate"=0.01, "Scenario"="Base")
parameters[[11]] <- list("N"=10000, "val.size"=1000, "HR.X"=1, "HR.C1"=1.3, "HR.C2"=1.1, "HR.U1"=1.3, "HR.U2"=2,
                        "OR.int"=0.25, "OR.C1"=1.3, "OR.C2"=0.9, "OR.U1"=1.3, "OR.U2"=2, "rate"=0.01, "Scenario"="Base")

# Change size of validation sample
parameters[[12]] <- list("N"=10000, "val.size"=500, "HR.X"=1, "HR.C1"=1.3, "HR.C2"=1.3, "HR.U1"=1.3, "HR.U2"=2,
                        "OR.int"=0.25, "OR.C1"=1.3, "OR.C2"=1.3, "OR.U1"=1.3, "OR.U2"=2, "rate"=0.01, "Scenario"="Base")
parameters[[13]] <- list("N"=10000, "val.size"=500, "HR.X"=1, "HR.C1"=1.3, "HR.C2"=1.3, "HR.U1"=1.3, "HR.U2"=2,
                        "OR.int"=0.25, "OR.C1"=1.3, "OR.C2"=1.3, "OR.U1"=0.75, "OR.U2"=0.5, "rate"=0.01, "Scenario"="Base")
parameters[[14]] <- list("N"=10000, "val.size"=250, "HR.X"=1, "HR.C1"=1.3, "HR.C2"=1.3, "HR.U1"=1.3, "HR.U2"=2,
                        "OR.int"=0.25, "OR.C1"=1.3, "OR.C2"=1.3, "OR.U1"=1.3, "OR.U2"=2, "rate"=0.01, "Scenario"="Base")
parameters[[15]] <- list("N"=10000, "val.size"=2500, "HR.X"=1, "HR.C1"=1.3, "HR.C2"=1.3, "HR.U1"=1.3, "HR.U2"=2,
                        "OR.int"=0.25, "OR.C1"=1.3, "OR.C2"=1.3, "OR.U1"=0.75, "OR.U2"=0.5, "rate"=0.01, "Scenario"="Base")

# Selection into validation sample
parameters[[16]] <- list("N"=10000, "val.size"=1000, "HR.X"=1, "HR.C1"=1.3, "HR.C2"=1.3, "HR.U1"=1.3, "HR.U2"=2,
                        "OR.int"=0.25, "OR.C1"=1.3, "OR.C2"=1.3, "OR.U1"=1.3, "OR.U2"=2, "rate"=0.01, "Scenario"="MAR")
parameters[[17]] <- list("N"=10000, "val.size"=1000, "HR.X"=1, "HR.C1"=1.3, "HR.C2"=1.3, "HR.U1"=1.3, "HR.U2"=2,
                        "OR.int"=0.25, "OR.C1"=1.3, "OR.C2"=1.3, "OR.U1"=1.3, "OR.U2"=2, "rate"=0.01, "Scenario"="MNAR")

# External validation sample (different event rate)
parameters[[18]] <- list("N"=10000, "val.size"=1000, "HR.X"=1, "HR.C1"=1.3, "HR.C2"=1.3, "HR.U1"=1.3, "HR.U2"=2,
                        "OR.int"=0.25, "OR.C1"=1.3, "OR.C2"=1.3, "OR.U1"=1.3, "OR.U2"=2, "rate"=0.01, "Scenario"="ExternalVS")

# Underlying hazard, censoring
parameters[[19]] <- list("N"=10000, "val.size"=1000, "HR.X"=1, "HR.C1"=1.3, "HR.C2"=1.3, "HR.U1"=1.3, "HR.U2"=2,
                        "OR.int"=0.25, "OR.C1"=1.3, "OR.C2"=1.3, "OR.U1"=1.3, "OR.U2"=2, "rate"=0.01, 
                        "Scenario"="RandomAndAdminCens")
parameters[[20]] <- list("N"=10000, "val.size"=1000, "HR.X"=1, "HR.C1"=1.3, "HR.C2"=1.3, "HR.U1"=1.3, "HR.U2"=2,
                         "OR.int"=0.25, "OR.C1"=1.3, "OR.C2"=1.3, "OR.U1"=1.3, "OR.U2"=2, "rate"=0.01, 
                         "Scenario"="RandomCens", "shape"=1)
parameters[[21]] <- list("N"=10000, "val.size"=1000, "HR.X"=1, "HR.C1"=1.3, "HR.C2"=1.3, "HR.U1"=1.3, "HR.U2"=2,
                         "OR.int"=0.25, "OR.C1"=1.3, "OR.C2"=1.3, "OR.U1"=1.3, "OR.U2"=2, "rate"=0.01, 
                         "Scenario"="RandomCens", "shape"=2)


