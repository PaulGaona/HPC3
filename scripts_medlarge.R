# functions and packages
source(".//funs.R")
source("./packages.R")
# script for running analysis
source("./medlarge_multi_sims.R")

# results from multi_sims
save(medlarge.cov.medlarge.n, file = "medlarge.cov.medlarge.n2.rda")
