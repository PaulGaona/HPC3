# functions and packages
source("./funs.R")
source("./packages.R")
# script for running analysis

source("./static_sim_and_setting.R")

# save mse as rda file
save(mse.all, file = "mse.all.rda")
# save lasso covariates as rda file
save(lasso.cov.kept, file = "lasso.cov.kept.rda")

