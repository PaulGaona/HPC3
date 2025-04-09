# parameter combination
p.vals <- c(10, 25, 50)  # covs
n.obs <- c(50, 100, 250, 500, 1000)  # obs
settings <- 1:9
nsims <- 100

# parallel backend (get cpu cores)
n_cores <- length( parallelly::availableWorkers() )/2
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# necessary functions and packages to cluster
clusterEvalQ(cl, {
  library(glmnet)
  library(rpart)
  library(randomForest)
  library(RFPlus)
})

# df to store all parameter combinations
sim.combos <- expand.grid(
  p = p.vals,
  n = n.obs,
  setting = settings,
  sim = 1:nsims
)

# initialize results data frame
sim.res <- data.frame()

# parallel sims
set.seed(92617)
sim.res <- foreach(i = 1:nrow(sim.combos), .errorhandling = 'pass', .packages = c("glmnet", "rpart", "randomForest", "RFPlus")) %dopar% {
  # extract parameters
  n <- sim.combos$n[i]
  p <- sim.combos$p[i]
  setting <- sim.combos$setting[i]
  sim_num <- sim.combos$sim[i]
  sd <- 1

  # run sims
  tryCatch({
    ###### PREP
    # DGP (sims)
    sim.data <- sparse.sims(n = n, p = p, sd = sd, setting = setting)

    og.train <- as.data.frame(sim.data[[1]]$train)
    og.test <- as.data.frame(sim.data[[1]]$test)

    # Oracle model
    oracle.train.model.list <- oracle.mods(data = og.train, setting = setting)
    oracle.train.model <- oracle.train.model.list[[1]]

    # MDI framework
    mdi.def <- rfplus(x = subset(og.train, select = -y), y = og.train$y,
                      normalize_stumps = TRUE, normalize_raw = TRUE,
                      ntree = 1, replace = FALSE, sample.fraction = 1, mtry = ncol(og.train) - 1,
                      min.bucket = round(20/3), min.node.size = 20)

    psis.def.og.train <- mdi.def$psis_train[[1]]
    psis.def.og.train.df <- cbind(y = og.train$y, psis.def.og.train)

    psis.def.og.test <- get_train_test_psis(mdi.def, data.frame(subset(og.test, select = -y)))$test[[1]]
    og.and.psis.og.test <- cbind(og.test, psis.def.og.test)

    # Augmented Data frames
    num.nodes.depth <- c(1,3,7,15,31) # up to depth 4

    for (j in 1:length(num.nodes.depth)) {
      if (num.nodes.depth[j] < ncol(psis.def.og.train)) {
        assign(paste0("og.and.d", j-1, ".train.df"),
               cbind(og.train, subset(psis.def.og.train, select = c(1:num.nodes.depth[j]))))
      } else {
        assign(paste0("og.and.d", j-1, ".train.df"),
               cbind(og.train, psis.def.og.train))
      }
    }
    og.and.psis.og.train <- cbind(og.train, psis.def.og.train)

    ##### MODELS
    # Linear models
    lin.og.train <- lm(y ~ ., data = og.train)
    lin.psis.og.train <- lm(y ~ ., data = psis.def.og.train.df)
    lin.d0.train <- lm(y ~ ., data = og.and.d0.train.df)
    lin.d1.train <- lm(y ~ ., data = og.and.d1.train.df)
    lin.d2.train <- lm(y ~ ., data = og.and.d2.train.df)
    lin.d3.train <- lm(y ~ ., data = og.and.d3.train.df)
    lin.d4.train <- lm(y ~ ., data = og.and.d4.train.df)
    lin.all.train <- lm(y ~ ., data = og.and.psis.og.train)

    # Lasso models
    las.og.train <- cv.glmnet(as.matrix(og.train[, -1]), og.train$y, alpha = 1)
    las.psis.og.train <- cv.glmnet(as.matrix(psis.def.og.train.df[, -1]), psis.def.og.train.df$y, alpha = 1)
    las.d0.train <- cv.glmnet(as.matrix(og.and.d0.train.df[, -1]), og.and.d0.train.df$y, alpha = 1)
    las.d1.train <- cv.glmnet(as.matrix(og.and.d1.train.df[, -1]), og.and.d1.train.df$y, alpha = 1)
    las.d2.train <- cv.glmnet(as.matrix(og.and.d2.train.df[, -1]), og.and.d2.train.df$y, alpha = 1)
    las.d3.train <- cv.glmnet(as.matrix(og.and.d3.train.df[, -1]), og.and.d3.train.df$y, alpha = 1)
    las.d4.train <- cv.glmnet(as.matrix(og.and.d4.train.df[, -1]), og.and.d4.train.df$y, alpha = 1)
    las.all.train <- cv.glmnet(as.matrix(og.and.psis.og.train[, -1]), og.and.psis.og.train$y, alpha = 1)

    # Ridge models
    rid.og.train <- cv.glmnet(as.matrix(og.train[, -1]), og.train$y, alpha = 0)
    rid.psis.og.train <- cv.glmnet(as.matrix(psis.def.og.train.df[, -1]), psis.def.og.train.df$y, alpha = 0)
    rid.d0.train <- cv.glmnet(as.matrix(og.and.d0.train.df[, -1]), og.and.d0.train.df$y, alpha = 0)
    rid.d1.train <- cv.glmnet(as.matrix(og.and.d1.train.df[, -1]), og.and.d1.train.df$y, alpha = 0)
    rid.d2.train <- cv.glmnet(as.matrix(og.and.d2.train.df[, -1]), og.and.d2.train.df$y, alpha = 0)
    rid.d3.train <- cv.glmnet(as.matrix(og.and.d3.train.df[, -1]), og.and.d3.train.df$y, alpha = 0)
    rid.d4.train <- cv.glmnet(as.matrix(og.and.d4.train.df[, -1]), og.and.d4.train.df$y, alpha = 0)
    rid.all.train <- cv.glmnet(as.matrix(og.and.psis.og.train[, -1]), og.and.psis.og.train$y, alpha = 0)

    # Tree-based methods
    tree.train <- rpart(y ~ ., data = og.train)
    rf.train <- randomForest(x = og.train[, -1], y = og.train$y)

    ##### PREDS
    oracle.test <- predict(oracle.train.model, newdata = og.test)

    # Linear predictions
    lin.test <- predict(lin.og.train, newdata = og.test)
    lin.psis.test <- predict(lin.psis.og.train, newdata = psis.def.og.test)
    lin.d0.test <- predict(lin.d0.train, newdata = og.and.psis.og.test)
    lin.d1.test <- predict(lin.d1.train, newdata = og.and.psis.og.test)
    lin.d2.test <- predict(lin.d2.train, newdata = og.and.psis.og.test)
    lin.d3.test <- predict(lin.d3.train, newdata = og.and.psis.og.test)
    lin.d4.test <- predict(lin.d4.train, newdata = og.and.psis.og.test)
    lin.all.test <- predict(lin.all.train, newdata = og.and.psis.og.test)

    # Lasso predictions
    newx.start <- as.matrix(og.and.psis.og.test[, -1])
    las.test <- predict(las.og.train, newx = as.matrix(og.test[, -1]), s = "lambda.1se")
    las.psis.test <- predict(las.psis.og.train, newx = as.matrix(psis.def.og.test), s = "lambda.1se")
    las.d0.test <- predict(las.d0.train, newx = newx.start[,colnames(og.and.d0.train.df[-1])], s = "lambda.1se")
    las.d1.test <- predict(las.d1.train, newx = newx.start[,colnames(og.and.d1.train.df[-1])], s = "lambda.1se")
    las.d2.test <- predict(las.d2.train, newx = newx.start[,colnames(og.and.d2.train.df[-1])], s = "lambda.1se")
    las.d3.test <- predict(las.d3.train, newx = newx.start[,colnames(og.and.d3.train.df[-1])], s = "lambda.1se")
    las.d4.test <- predict(las.d4.train, newx = newx.start[,colnames(og.and.d4.train.df[-1])], s = "lambda.1se")
    las.all.test <- predict(las.all.train, newx = newx.start, s = "lambda.1se")

    # Ridge predictions
    rid.test <- predict(rid.og.train, newx = as.matrix(og.test[, -1]), s = "lambda.1se")
    rid.psis.test <- predict(rid.psis.og.train, newx = as.matrix(psis.def.og.test), s = "lambda.1se")
    rid.d0.test <- predict(rid.d0.train, newx = newx.start[,colnames(og.and.d0.train.df[-1])], s = "lambda.1se")
    rid.d1.test <- predict(rid.d1.train, newx = newx.start[,colnames(og.and.d1.train.df[-1])], s = "lambda.1se")
    rid.d2.test <- predict(rid.d2.train, newx = newx.start[,colnames(og.and.d2.train.df[-1])], s = "lambda.1se")
    rid.d3.test <- predict(rid.d3.train, newx = newx.start[,colnames(og.and.d3.train.df[-1])], s = "lambda.1se")
    rid.d4.test <- predict(rid.d4.train, newx = newx.start[,colnames(og.and.d4.train.df[-1])], s = "lambda.1se")
    rid.all.test <- predict(rid.all.train, newx = newx.start, s = "lambda.1se")

    # Tree predictions
    tree.def.test <- predict(tree.train, newdata = og.test)
    rf.test <- predict(rf.train, newdata = og.test)

    ###### MSE
    mse.1sim <- c(
      oracle = mean((oracle.test - og.test$y)^2),
      lin = mean((lin.test - og.test$y)^2),
      lin.psis = mean((lin.psis.test - og.test$y)^2),
      lin.d0 = mean((lin.d0.test - og.test$y)^2),
      lin.d1 = mean((lin.d1.test - og.test$y)^2),
      lin.d2 = mean((lin.d2.test - og.test$y)^2),
      lin.d3 = mean((lin.d3.test - og.test$y)^2),
      lin.d4 = mean((lin.d4.test - og.test$y)^2),
      lin.all = mean((lin.all.test - og.test$y)^2),
      lasso = mean((las.test - og.test$y)^2),
      lasso.psis = mean((las.psis.test - og.test$y)^2),
      lasso.d0 = mean((las.d0.test - og.test$y)^2),
      lasso.d1 = mean((las.d1.test - og.test$y)^2),
      lasso.d2 = mean((las.d2.test - og.test$y)^2),
      lasso.d3 = mean((las.d3.test - og.test$y)^2),
      lasso.d4 = mean((las.d4.test - og.test$y)^2),
      lasso.all = mean((las.all.test - og.test$y)^2),
      ridge = mean((rid.test - og.test$y)^2),
      ridge.psis = mean((rid.psis.test - og.test$y)^2),
      ridge.d0 = mean((rid.d0.test - og.test$y)^2),
      ridge.d1 = mean((rid.d1.test - og.test$y)^2),
      ridge.d2 = mean((rid.d2.test - og.test$y)^2),
      ridge.d3 = mean((rid.d3.test - og.test$y)^2),
      ridge.d4 = mean((rid.d4.test - og.test$y)^2),
      ridge.all = mean((rid.all.test - og.test$y)^2),
      tree.def = mean((tree.def.test - og.test$y)^2),
      rf = mean((rf.test - og.test$y)^2)
    )

    # lasso non-zero coefficients
    cov.kept.list <- list(
      all_og = rownames(coef(las.og.train, s = las.og.train$lambda.1se))[which(coef(las.og.train, s = las.og.train$lambda.1se) != 0)],
      all_nodes = rownames(coef(las.psis.og.train, s = las.psis.og.train$lambda.1se))[which(coef(las.psis.og.train, s = las.psis.og.train$lambda.1se) != 0)],
      all_og_d0 = rownames(coef(las.d0.train, s = las.d0.train$lambda.1se))[which(coef(las.d0.train, s = las.d0.train$lambda.1se) != 0)],
      all_og_d1 = rownames(coef(las.d1.train, s = las.d1.train$lambda.1se))[which(coef(las.d1.train, s = las.d1.train$lambda.1se) != 0)],
      all_og_d2 = rownames(coef(las.d2.train, s = las.d2.train$lambda.1se))[which(coef(las.d2.train, s = las.d2.train$lambda.1se) != 0)],
      all_og_d3 = rownames(coef(las.d3.train, s = las.d3.train$lambda.1se))[which(coef(las.d3.train, s = las.d3.train$lambda.1se) != 0)],
      all_og_d4 = rownames(coef(las.d4.train, s = las.d4.train$lambda.1se))[which(coef(las.d4.train, s = las.d4.train$lambda.1se) != 0)],
      all_og_all_nodes = rownames(coef(las.all.train, s = las.all.train$lambda.1se))[which(coef(las.all.train, s = las.all.train$lambda.1se) != 0)]
    )

    # sim.res as a data frame row
    res.1sim.row <- data.frame(
      n = n,
      p = p,
      setting = setting,
      sim = sim_num,
      t(mse.1sim),
      stringsAsFactors = FALSE
    )

    # information about non-zero coefficients
    res.1sim.row$cov_kept_all_og <- list(cov.kept.list$all_og)
    res.1sim.row$cov_kept_all_nodes <- list(cov.kept.list$all_nodes)
    res.1sim.row$cov_kept_all_og_d0 <- list(cov.kept.list$all_og_d0)
    res.1sim.row$cov_kept_all_og_d1 <- list(cov.kept.list$all_og_d1)
    res.1sim.row$cov_kept_all_og_d2 <- list(cov.kept.list$all_og_d2)
    res.1sim.row$cov_kept_all_og_d3 <- list(cov.kept.list$all_og_d3)
    res.1sim.row$cov_kept_all_og_d4 <- list(cov.kept.list$all_og_d4)
    res.1sim.row$cov_kept_all_og_all_nodes <- list(cov.kept.list$all_og_all_nodes)

    return(res.1sim.row)
  }, error = function(e) {
    # error check return (prob should delete)
    res.1sim.row <- data.frame(
      n = n,
      p = p,
      setting = setting,
      sim = sim_num,
      error = TRUE,
      error_message = as.character(e),
      stringsAsFactors = FALSE
    )
    return(res.1sim.row)
  })
}

# end
stopCluster(cl)

small.cov.small.n <- sim.res



