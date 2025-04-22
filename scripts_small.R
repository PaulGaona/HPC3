# functions and packages
source(".//funs.R")
source("./packages.R")
# script for running analysis
# change when needed to
#source("./static_sim_and_setting.R")
#source("./large_multi_sims.R")
source("./multi_sims.R")
# save mse as rda file
#save(mse.all, file = "mse.all.rda")
# save lasso covariates as rda file
#save(cov.kept.list, file = "lasso.cov.kept.rda")

# results from multi_sims
#save(large.cov.large.n, file = "large.cov.large.n.rda")
save(small.cov.small.n, file = "small.cov.small.n2.rda")
