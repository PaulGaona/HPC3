# functions and packages
source(".//funs.R")
source("./packages.R")
# script for running analysis
source("./med_multi_sims.R")

# results from multi_sims
save(med.cov.med.n, file = "med.cov.med.n2.rda")
