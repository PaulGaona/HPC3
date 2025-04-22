# functions and packages
source(".//funs.R")
source("./packages.R")
# script for running analysis
source("./med_multi_sims.R")

# results from multi_sims
save(large.cov.large.n, file = "large.cov.large.n2.rda")
