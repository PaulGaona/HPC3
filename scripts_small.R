# functions and packages
source(".//funs.R")
source("./packages.R")
# script for running analysis
source("./multi_sims.R")

# results from multi_sims
save(small.cov.small.n, file = "small.cov.small.n2.rda")
