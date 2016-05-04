# Simulation code for martingale residual-based imputation method

The code in this repository was used to perform the simulations presented in the paper: Burne and Abrahamowicz (2016) "Martingale residual-based method to control for confounders measured only in a validation sample in time-to-event analysis"

## Time-fixed simulations (Section 3.1)

The simulations in this section are run by running `SimulationsTF.R`, which sources `ParametersTF.R`, `DataGenerationTF.R` and `FunctionsTF.R`.
*  `ParametersTF.R` contains a list of parameter sets, with indices aligning with rows in Table 1
*  `DataGenerationTF.R` contains a function to generate the data, and
*  `FunctionsTF.R` contains functions which are used to perform the analyses and obtain results given in Table 1.

## Time-varying simulations (Section 3.2)

The simulations in this section are run by running `SimulationsTD.R`, which sources `ParametersTD.R`, `DataGenerationTD.R` and `FunctionsTD.R`.
*  `ParametersTD.R` contains a list of parameter sets, with indices aligning with rows in Table 3
*  `DataGenerationTD.R` contains a function to generate the data, and
*  `FunctionsTD.R` contains functions which are used to perform the analyses and obtain results given in Table 3.

