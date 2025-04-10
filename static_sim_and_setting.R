n = 50
p = 25
sd = 1
setting = 1

#set.seed(92617)
sim.data <- sparse.sims(n = n, p = p, sd = sd, setting = setting)

og.train <- as.data.frame(sim.data[[1]]$train)
og.test <- as.data.frame(sim.data[[1]]$test)

##### Oracle
### oracle model
oracle.train.model.list <- oracle.mods(data = og.train, setting = setting)
oracle.train.model <- oracle.train.model.list[[1]]
##### MDI framework
# NODES
# poor performance if wanted to can uncomment code
#mdi <- rfplus(x = subset(og.train, select = -y), y = og.train$y,
#              normalize_stumps = TRUE, normalize_raw = TRUE,
#              ntree = 1, mtry = ncol(og.train) - 1, replace = FALSE,
#              sample.fraction = 1, min.bucket = 1, min.node.size = 1)
#psis.og.train <- mdi$psis_train[[1]]
#nodes.df <- cbind(y = og.train$y, psis.og.train)

### mdi+ default (unable to prune possible to prune via depth)
mdi.def <- rfplus(x = subset(og.train, select = -y), y = og.train$y,
                  normalize_stumps = TRUE, normalize_raw = TRUE,
                  ntree = 1, replace = FALSE, sample.fraction = 1, mtry = ncol(og.train) - 1,
                  min.bucket = round(20/3), min.node.size = 20)

psis.def.og.train <- mdi.def$psis_train[[1]]
psis.def.og.train.df <- cbind(y = og.train$y, psis.def.og.train)
og.and.psis.og.train <- cbind(og.train, psis.def.og.train)

psis.def.og.test <- get_train_test_psis(mdi.def, data.frame(subset(og.test, select = -y)))$test[[1]]
og.and.psis.og.test <- cbind(og.test, psis.def.og.test)
# Augmented Data frames
# cross validaiton set up
num.nodes.depth <- c(1,3,7,15,31) # up to depth 4
#set.seed(92617)
# 5 fold
tot.n <- nrow(og.train)
folds <- sample(rep(1:5, length.out = tot.n))

cv.res <- data.frame(depth = character(), mse = numeric(), stringsAsFactors = FALSE)
lasso.cv.res <- data.frame(depth = character(), nodes = numeric(), mse = numeric(), stringsAsFactors = FALSE)
ridge.cv.res <- data.frame(depth = character(), nodes = numeric(), mse = numeric(), stringsAsFactors = FALSE)

for (i in 1:length(num.nodes.depth)) {
  depth.name <- paste0("d", i-1)

  # Create dataset for current depth only (don't store all in memory)
  if (num.nodes.depth[i] < ncol(psis.def.og.train)) {
    curr.depth.df <- cbind(og.train, subset(psis.def.og.train, select = c(1:num.nodes.depth[i])))
  } else {
    curr.depth.df <- cbind(og.train, psis.def.og.train)
  }

  fold.mse <- numeric(5)
  # for linear model
  # For each fold
  for (k in 1:5) {
    # training and validation indeces
    cv.train.indices <- which(folds != k)
    cv.valid.indices <- which(folds == k)

    cv.train.data <- curr.depth.df[cv.train.indices, ]
    cv.valid.data <- curr.depth.df[cv.valid.indices, ]

    cv.mod <- lm(y ~ ., data = cv.train.data)
    cv.pred <- predict(cv.mod, newdata = cv.valid.data)

    # Calculate MSE for this fold
    fold.mse[k] <- mean((cv.valid.data$y - cv.pred)^2)

    # Clean up to free memory
    rm(cv.train.data, cv.valid.data, cv.mod, cv.pred)
    gc(verbose = FALSE)
  }
  # Average MSE across folds
  cv.avg.mse <- mean(fold.mse)

  cv.res <- rbind(cv.res, data.frame(
    depth = depth.name,
    nodes = num.nodes.depth[i],
    mse = cv.avg.mse
  ))
  # regularization x and y
  x <- as.matrix(curr.depth.df[, -1])
  y <- curr.depth.df$y

  # for lasso
  cv.lasso <- cv.glmnet(x, y, alpha = 1, nfolds = 5)

  min.mse.lasso <- min(cv.lasso$cvm)

  lasso.cv.res <- rbind(lasso.cv.res, data.frame(
    depth = depth.name,
    nodes = num.nodes.depth[i],
    mse = min.mse.lasso,
    min.lambda = cv.lasso$lambda.min,
    min.lambda.1se = cv.lasso$lambda.1se
  ))

  # for ridge
  cv.ridge <- cv.glmnet(x, y, alpha = 0, nfolds = 5)
  min.mse.ridge <- min(cv.ridge$cvm)

  ridge.cv.res <- rbind(ridge.cv.res, data.frame(
    depth = depth.name,
    nodes = num.nodes.depth[i],
    mse = min.mse.ridge,
    min.lambda = cv.ridge$lambda.min,
    min.lambda.1se = cv.ridge$lambda.1se
  ))

  # Clean up to free memory
  rm(curr.depth.df, fold.mse,cv.lasso, x, y, cv.ridge, min.mse.lasso, min.mse.ridge)

  gc(verbose = FALSE)
}
# get dataframe for best depth from cv.res
best.depth <- cv.res[which.min(cv.res$mse), "depth"]
best.depth.ind <- as.numeric(gsub("d", "", best.depth)) + 1
# lasso
best.depth.lasso <- lasso.cv.res[which.min(lasso.cv.res$mse), "depth"]
best.depth.lasso.ind <- as.numeric(gsub("d", "", best.depth.lasso)) + 1
# ridge
best.depth.ridge <- ridge.cv.res[which.min(ridge.cv.res$mse), "depth"]
best.depth.ridge.ind <- as.numeric(gsub("d", "", best.depth.ridge)) + 1

if (num.nodes.depth[best.depth.ind] < ncol(psis.def.og.train)) {
  cv.og.and.node.train.df <- cbind(og.train, subset(psis.def.og.train, select = c(1:num.nodes.depth[best.depth.ind])))
} else {
  cv.og.and.node.train.df <- cbind(og.train, psis.def.og.train)
}

if (num.nodes.depth[best.depth.lasso.ind] < ncol(psis.def.og.train)) {
  cv.og.and.node.train.df <- cbind(og.train, subset(psis.def.og.train, select = c(1:num.nodes.depth[best.depth.lasso.ind])))
} else {
  cv.og.and.node.train.df <- cbind(og.train, psis.def.og.train)
}

if (num.nodes.depth[best.depth.ridge.ind] < ncol(psis.def.og.train)) {
  cv.og.and.node.train.df <- cbind(og.train, subset(psis.def.og.train, select = c(1:num.nodes.depth[best.depth.ridge.ind])))
} else {
  cv.og.and.node.train.df <- cbind(og.train, psis.def.og.train)
}


##### Statistical Models
# linear
lin.og.train <- lm(y ~ ., data = og.train) # linear - all og covariates
lin.psis.og.train <- lm(y ~ ., data = psis.def.og.train.df) # linear - all nodes
lin.dcv.train <- lm(y ~ ., data = cv.og.and.node.train.df) # linear - best model from cv
lin.all.train <- lm(y ~ ., data = og.and.psis.og.train) # linear - all og covariates and all nodes
# lasso
las.og.train <- cv.glmnet(as.matrix(og.train[, -1]), og.train$y, alpha = 1) # lasso - all og covariates
las.psis.og.train <- cv.glmnet(as.matrix(psis.def.og.train.df[, -1]), psis.def.og.train.df$y, alpha = 1) # lasso - all nodes
las.dcv.train <- cv.glmnet(as.matrix(cv.og.and.node.train.df[, -1]), cv.og.and.node.train.df$y, alpha = 1) # lasso - best model from cv
las.all.train <- cv.glmnet(as.matrix(og.and.psis.og.train[, -1]), og.and.psis.og.train$y, alpha = 1) # lasso - all og covariates and all nodes
# return non-0 coefficients
# ridge
rid.og.train <- cv.glmnet(as.matrix(og.train[, -1]), og.train$y, alpha = 0) # ridge - all og covariates
rid.psis.og.train <- cv.glmnet(as.matrix(psis.def.og.train.df[, -1]), psis.def.og.train.df$y, alpha = 0) # ridge - all nodes
rid.dcv.train <- cv.glmnet(as.matrix(cv.og.and.node.train.df[, -1]), cv.og.and.node.train.df$y, alpha = 0) # ridge - best model from cv
rid.all.train <- cv.glmnet(as.matrix(og.and.psis.og.train[, -1]), og.and.psis.og.train$y, alpha = 0) # ridge - all og covariates and all nodes
#### tree based methods
# tree
tree.train <- rpart(y ~ ., data = og.train)
# random forest
rf.train <- randomForest(x = og.train[, -1], y = og.train$y)
################# Predictions
oracle.test <- predict(oracle.train.model, newdata = og.test) # oracle
# linear
lin.test <- predict(lin.og.train, newdata = og.test) # linear - all og covariates
lin.psis.test <- predict(lin.psis.og.train, newdata = psis.def.og.test) # linear - all nodes
lin.dcv.test <- predict(lin.dcv.train, newdata = og.and.psis.og.test) # linear - best model from cv
lin.all.test <- predict(lin.all.train, newdata = og.and.psis.og.test) # linear - all og covariates and all nodes
# lasso
newx.start <- as.matrix(og.and.psis.og.test[, -1])

las.test <- predict(las.og.train, newx = as.matrix(og.test[, -1]), s = "lambda.1se") # lasso - all og covariates
las.psis.test <- predict(las.psis.og.train, newx = as.matrix(psis.def.og.test), s = "lambda.1se") # lasso - all nodes
las.dcv.test <- predict(las.dcv.train, newx = newx.start[,colnames(cv.og.and.node.train.df[, -1])], s = "lambda.1se") # lasso - best model from cv
las.all.test <- predict(las.all.train, newx = newx.start, s = "lambda.1se") # lasso - all og covariates and all nodes
# ridge
rid.test <- predict(rid.og.train, newx = as.matrix(og.test[, -1]), s = "lambda.1se") # ridge - all og covariates
rid.psis.test <- predict(rid.psis.og.train, newx = as.matrix(psis.def.og.test), s = "lambda.1se") # ridge - all nodes
rid.dcv.test <- predict(rid.dcv.train, newx = newx.start[,colnames(cv.og.and.node.train.df[, -1])], s = "lambda.1se") # ridge - best model from cv
rid.all.test <- predict(rid.all.train, newx = newx.start, s = "lambda.1se") # ridge - all og covariates and all nodes
# tree
tree.def.test <- predict(tree.train, newdata = og.test) # def tree
# random forest
rf.test <- predict(rf.train, newdata = og.test) # random forest
################# MSE
mse.all <- rbind(
  oracle = mean((og.test$y - oracle.test)^2),
  lin.og = mean((og.test$y - lin.test)^2),
  lin.psis = mean((og.test$y - lin.psis.test)^2),
  lin.dcv = mean((og.test$y - lin.dcv.test)^2),
  lin.all = mean((og.test$y - lin.all.test)^2),
  las.og = mean((og.test$y - las.test)^2),
  las.psis = mean((og.test$y - las.psis.test)^2),
  las.dcv = mean((og.test$y - las.dcv.test)^2),
  las.all = mean((og.test$y - las.all.test)^2),
  rid.og = mean((og.test$y - rid.test)^2),
  rid.psis = mean((og.test$y - rid.psis.test)^2),
  rid.dcv = mean((og.test$y - rid.dcv.test)^2),
  rid.all = mean((og.test$y - rid.all.test)^2),
  tree.def = mean((og.test$y - tree.def.test)^2),
  rf = mean((og.test$y - rf.test)^2)
)
mod.names <- c("oracle",
             "lin.og","lin.psis","lin.dcv","lin.all",
             "las.og","las.psis","las.dcv","las.all",
             "rid.og","rid.psis","rid.dcv","rid.all",
             "tree.def", "rf")

rownames(mse.all) <- mod.names

mse.all
# list of vectors with non-zero coefficient names
cov.kept.list <- list(  rownames(coef(las.og.train, s = las.og.train$lambda.1se))[which(coef(las.og.train, s = las.og.train$lambda.1se) != 0)],
                        rownames(coef(las.psis.og.train, s = las.psis.og.train$lambda.1se))[which(coef(las.psis.og.train, s = las.psis.og.train$lambda.1se) != 0)],
                        rownames(coef(las.dcv.train, s = las.dcv.train$lambda.1se))[which(coef(las.dcv.train, s = las.dcv.train$lambda.1se) != 0)],
                        rownames(coef(las.all.train, s = las.all.train$lambda.1se))[which(coef(las.all.train, s = las.all.train$lambda.1se) != 0)]
)

names(cov.kept.list) <- c("og", "nodes", "cv", "all")

cov.kept.list
