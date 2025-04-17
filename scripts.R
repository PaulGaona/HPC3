# functions and packages
source(".//funs.R")
source("./packages.R")
# script for running analysis
# change when needed to
#source("./static_sim_and_setting.R")
source("./multi_sims.R")
# save mse as rda file
#save(mse.all, file = "mse.all.rda")
# save lasso covariates as rda file
#save(cov.kept.list, file = "lasso.cov.kept.rda")

# results from multi_sims
#save(small.cov.small.n, file = "small.cov.small.n.rda")
save(large.cov.large.n, file = "large.cov.large.n.rda")
