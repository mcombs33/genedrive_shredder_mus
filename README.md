# genedrive_shredder_mus
Functions and scripts for running X-Shredder and Y-Shredder gene drive simulations for an island population of Mus musculus

# Functions: all functions found in GDShredMus_Function.R
## Yshred
Simulates demographic, dispersal, and Y-Shredder gene drive processes across time steps. Requires input parameters and starting abundances for each sex-age-genotype combination.

## Yshred.v2
Same as Yshred except juvenile dispersal rate is set to same value as adult dispersal rate. Used for initial set of simulations in manuscript.

## Xshred
Simulates demographic, dispersal, and Y-Shredder gene drive processes across time steps. Requires input parameters and starting abundances for each sex-age-genotype combination.

## Xshred.v2
Same as Xshred except juvenile dispersal rate is set to same value as adult dispersal rate. Used for initial set of simulations in manuscript.

## calculateCentroids
Calculates XY position of patch centroids within square patch lattice. This function is implemented within the natMigration function.

## natMigration
Determines destination patch for each migrant and final abundance in each patch following migration. 
Requires distance threshold of permissible difference between dispersal distance and distance between starting location and centroids of possible destination patches.

## randBirthYS
Determines discrete number of births of each potential genotype in the Y-Shredder system using weighted sampling.

## randBirthXS
Determines discrete number of births of each potential genotype in the X-Shredder system using weighted sampling.

## randpois
Determines discrete numbers of individuals experiencing dispersal or demographic processes for each patch, given relevant parameter and patch abundances.

## randpoisLitter
Determines discrete number of individuals in each litter, capped at minimum and maximum values for house mice.

## MatedFemByPoly
Determines the number of mated females in each patch given male availability and mating system characteristics (i.e. max number of male mates, probabilities of female mating with 1 - n mates).

## minCheck
Determines minimum value. Used to determine number of mated females.

# Scripts:
## GDShredMus_XInitialSim.R
Script to setup and run X-shredder gene drive simulations for house mouse metapopulation. Run parameters established within script.

## GDShredMus_YInitialSim.R
Script to setup and run Y-shredder gene drive simulations for house mouse metapopulation. Run parameters established within script.

## GDShredMus_XSensAnLimited_LHS.R
Script to generate parameters for limited sensitivity analysis of X-Shredder system using latin hypercube sampling. Note that parameters must be combined with static parameters in correct column order before use.

## GDShredMus_XSensAnLimited.R
Script to run limited sensitivity analysis of X-Shredder gene drive.

## GDShredMus_YSensAnLimited_LHS.R
Script to generate parameters for limited sensitivity analysis of Y-Shredder system using latin hypercube sampling. Note that parameters must be combined with static parameters in correct column order before use.

## GDShredMus_YSensAnLimited.R
Script to run limited sensitivity analysis of X-Shredder gene drive.

## GDShredMus_XSensAnGlobal_LHS.R
Script to generate parameters for global sensitivity analysis of X-Shredder system using latin hypercube sampling. Note that parameters must be combined with static parameters in correct column order before use.

## GDShredMus_XSensAnGlobal.R
Script to run global sensitivity analysis of X-Shredder gene drive.

## GDShredMus_YSensAnGlobal_LHS.R
Script to generate parameters for global sensitivity analysis of Y-Shredder system using latin hypercube sampling. Note that parameters must be combined with static parameters in correct column order before use.

## GDShredMus_YSensAnGlobal.R
Script to run global sensitivity analysis of Y-Shredder gene drive.

