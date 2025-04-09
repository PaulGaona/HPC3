n = 100
p = 10
sd = 1
setting = 2

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

psis.def.og.test <- get_train_test_psis(mdi.def, data.frame(subset(og.test, select = -y)))$test[[1]]
og.and.psis.og.test <- cbind(og.test, psis.def.og.test)
# Augmented Data frames
num.nodes.depth <- c(1,3,7,15,31) # up to depth 4

for (i in 1:length(num.nodes.depth)) {
  if (num.nodes.depth[i] < ncol(psis.def.og.train)) {
    assign( paste0( "og.and.d", i-1 ,".train.df"), cbind(og.train, subset(psis.def.og.train, select = c(1:num.nodes.depth[i]) ) )
    ) } else {
      assign(paste0("og.and.d", i-1 ,".train.df"), cbind(og.train, psis.def.og.train))
    }
}
og.and.psis.og.train <- cbind(og.train, psis.def.og.train)
##### Statistical Models
# linear
lin.og.train <- lm(y ~ ., data = og.train) # linear - all og covariates
lin.psis.og.train <- lm(y ~ ., data = psis.def.og.train.df) # linear - all nodes
lin.d0.train <- lm(y ~ ., data = og.and.d0.train.df) # linear - all og covariates and [depth 0 ,nodes 1 ]
lin.d1.train <- lm(y ~ ., data = og.and.d1.train.df) # linear - all og covariates and [depth 1 ,nodes 3 ]
lin.d2.train <- lm(y ~ ., data = og.and.d2.train.df) # linear - all og covariates and [depth 2 ,nodes 7 ]
lin.d3.train <- lm(y ~ ., data = og.and.d3.train.df) # linear - all og covariates and [depth 3 ,nodes 15 ]
lin.d4.train <- lm(y ~ ., data = og.and.d4.train.df) # linear - all og covariates and [depth 4 ,nodes 31 ]
lin.all.train <- lm(y ~ ., data = og.and.psis.og.train) # linear - all og covariates and all nodes
# lasso
las.og.train <- cv.glmnet(as.matrix(og.train[, -1]), og.train$y, alpha = 1) # lasso - all og covariates
las.psis.og.train <- cv.glmnet(as.matrix(psis.def.og.train.df[, -1]), psis.def.og.train.df$y, alpha = 1) # lasso - all nodes
las.d0.train <- cv.glmnet(as.matrix(og.and.d0.train.df[, -1]), og.and.d0.train.df$y, alpha = 1) # lasso - all og covariates and [depth 0 ,nodes 1 ]
las.d1.train <- cv.glmnet(as.matrix(og.and.d1.train.df[, -1]), og.and.d1.train.df$y, alpha = 1) # lasso - all og covariates and [depth 1 ,nodes 3 ]
las.d2.train <- cv.glmnet(as.matrix(og.and.d2.train.df[, -1]), og.and.d2.train.df$y, alpha = 1) # lasso - all og covariates and [depth 2 ,nodes 7 ]
las.d3.train <- cv.glmnet(as.matrix(og.and.d3.train.df[, -1]), og.and.d3.train.df$y, alpha = 1) # lasso - all og covariates and [depth 3 ,nodes 15 ]
las.d4.train <- cv.glmnet(as.matrix(og.and.d4.train.df[, -1]), og.and.d4.train.df$y, alpha = 1) # lasso - all og covariates and [depth 4 ,nodes 31 ]
las.all.train <- cv.glmnet(as.matrix(og.and.psis.og.train[, -1]), og.and.psis.og.train$y, alpha = 1) # lasso - all og covariates and all nodes
# return non-0 coefficients
# ridge
rid.og.train <- cv.glmnet(as.matrix(og.train[, -1]), og.train$y, alpha = 0) # ridge - all og covariates
rid.psis.og.train <- cv.glmnet(as.matrix(psis.def.og.train.df[, -1]), psis.def.og.train.df$y, alpha = 0) # ridge - all nodes
rid.d0.train <- cv.glmnet(as.matrix(og.and.d0.train.df[, -1]), og.and.d0.train.df$y, alpha = 0) # ridge - all og covariates and [depth 0 ,nodes 1 ]
rid.d1.train <- cv.glmnet(as.matrix(og.and.d1.train.df[, -1]), og.and.d1.train.df$y, alpha = 0) # ridge - all og covariates and [depth 1 ,nodes 3 ]
rid.d2.train <- cv.glmnet(as.matrix(og.and.d2.train.df[, -1]), og.and.d2.train.df$y, alpha = 0) # ridge - all og covariates and [depth 2 ,nodes 7 ]
rid.d3.train <- cv.glmnet(as.matrix(og.and.d3.train.df[, -1]), og.and.d3.train.df$y, alpha = 0) # ridge - all og covariates and [depth 3 ,nodes 15 ]
rid.d4.train <- cv.glmnet(as.matrix(og.and.d4.train.df[, -1]), og.and.d4.train.df$y, alpha = 0) # ridge - all og covariates and [depth 4 ,nodes 31 ]
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
lin.d0.test <- predict(lin.d0.train, newdata = og.and.psis.og.test) # linear - all og covariates and [depth 0 ,nodes 1 ]
lin.d1.test <- predict(lin.d1.train, newdata = og.and.psis.og.test) # linear - all og covariates and [depth 1 ,nodes 3 ]
lin.d2.test <- predict(lin.d2.train, newdata = og.and.psis.og.test) # linear - all og covariates and [depth 2 ,nodes 7 ]
lin.d3.test <- predict(lin.d3.train, newdata = og.and.psis.og.test) # linear - all og covariates and [depth 3 ,nodes 15 ]
lin.d4.test <- predict(lin.d4.train, newdata = og.and.psis.og.test) # linear - all og covariates and [depth 4 ,nodes 31 ]
lin.all.test <- predict(lin.all.train, newdata = og.and.psis.og.test) # linear - all og covariates and all nodes
# lasso
newx.start <- as.matrix(og.and.psis.og.test[, -1])

las.test <- predict(las.og.train, newx = as.matrix(og.test[, -1]), s = "lambda.1se") # lasso - all og covariates
las.psis.test <- predict(las.psis.og.train, newx = as.matrix(psis.def.og.test), s = "lambda.1se") # lasso - all nodes
las.d0.test <- predict(las.d0.train, newx = newx.start[,colnames(og.and.d0.train.df[-1])], s = "lambda.1se") # lasso - all og covariates and [depth 0 ,nodes 1 ]
las.d1.test <- predict(las.d1.train, newx = newx.start[,colnames(og.and.d1.train.df[-1])], s = "lambda.1se") # lasso - all og covariates and [depth 1 ,nodes 3 ]
las.d2.test <- predict(las.d2.train, newx = newx.start[,colnames(og.and.d2.train.df[-1])], s = "lambda.1se") # lasso - all og covariates and [depth 2 ,nodes 7 ]
las.d3.test <- predict(las.d3.train, newx = newx.start[,colnames(og.and.d3.train.df[-1])], s = "lambda.1se") # lasso - all og covariates and [depth 3 ,nodes 15 ]
las.d4.test <- predict(las.d4.train, newx = newx.start[,colnames(og.and.d4.train.df[-1])], s = "lambda.1se") # lasso - all og covariates and [depth 4 ,nodes 31 ]
las.all.test <- predict(las.all.train, newx = newx.start, s = "lambda.1se") # lasso - all og covariates and all nodes
# ridge
rid.test <- predict(rid.og.train, newx = as.matrix(og.test[, -1]), s = "lambda.1se") # ridge - all og covariates
rid.psis.test <- predict(rid.psis.og.train, newx = as.matrix(psis.def.og.test), s = "lambda.1se") # ridge - all nodes
rid.d0.test <- predict(rid.d0.train, newx = newx.start[,colnames(og.and.d0.train.df[-1])], s = "lambda.1se") # ridge - all og covariates and [depth 0 ,nodes 1 ]
rid.d1.test <- predict(rid.d1.train, newx = newx.start[,colnames(og.and.d1.train.df[-1])], s = "lambda.1se") # ridge - all og covariates and [depth 1 ,nodes 3 ]
rid.d2.test <- predict(rid.d2.train, newx = newx.start[,colnames(og.and.d2.train.df[-1])], s = "lambda.1se") # ridge - all og covariates and [depth 2 ,nodes 7 ]
rid.d3.test <- predict(rid.d3.train, newx = newx.start[,colnames(og.and.d3.train.df[-1])], s = "lambda.1se") # ridge - all og covariates and [depth 3 ,nodes 15 ]
rid.d4.test <- predict(rid.d4.train, newx = newx.start[,colnames(og.and.d4.train.df[-1])], s = "lambda.1se") # ridge - all og covariates and [depth 4 ,nodes 31 ]
rid.all.test <- predict(rid.all.train, newx = newx.start, s = "lambda.1se") # ridge - all og covariates and all nodes
# tree
tree.def.test <- predict(tree.train, newdata = og.test) # def tree
# random forest
rf.test <- predict(rf.train, newdata = og.test) # random forest
################# MSE
mse.all <- rbind(
  mean((oracle.test - og.test$y)^2), # oracle
  mean((lin.test - og.test$y)^2), # linear - all og covariates
  mean((lin.psis.test - og.test$y)^2), # linear - all nodes
  mean((lin.d0.test - og.test$y)^2), # linear - all og covariates and [depth 0 ,nodes 1 ]
  mean((lin.d1.test - og.test$y)^2), # linear - all og covariates and [depth 1 ,nodes 3 ]
  mean((lin.d2.test - og.test$y)^2), # linear - all og covariates and [depth 2 ,nodes 7 ]
  mean((lin.d3.test - og.test$y)^2), # linear - all og covariates and [depth 3 ,nodes 15 ]
  mean((lin.d4.test - og.test$y)^2), # linear - all og covariates and [depth 4 ,nodes 31 ]
  mean((lin.all.test - og.test$y)^2), # linear - all og covariates and all nodes
  mean((las.test - og.test$y)^2), # lasso - all og covariates
  mean((las.psis.test - og.test$y)^2), # lasso - all nodes
  mean((las.d0.test - og.test$y)^2), # lasso - all og covariates and [depth 0 ,nodes 1 ]
  mean((las.d1.test - og.test$y)^2), # lasso - all og covariates and [depth 1 ,nodes 3 ]
  mean((las.d2.test - og.test$y)^2), # lasso - all og covariates and [depth 2 ,nodes 7 ]
  mean((las.d3.test - og.test$y)^2), # lasso - all og covariates and [depth 3 ,nodes 15 ]
  mean((las.d4.test - og.test$y)^2), # lasso - all og covariates and [depth 4 ,nodes 31 ]
  mean((las.all.test - og.test$y)^2), # lasso - all og covariates and all nodes
  mean((rid.test - og.test$y)^2), # ridge - all og covariates
  mean((rid.psis.test - og.test$y)^2), # ridge - all nodes
  mean((rid.d0.test - og.test$y)^2), # ridge - all og covariates and [depth 0 ,nodes 1 ]
  mean((rid.d1.test - og.test$y)^2), # ridge - all og covariates and [depth 1 ,nodes 3 ]
  mean((rid.d2.test - og.test$y)^2), # ridge - all og covariates and [depth 2 ,nodes 7 ]
  mean((rid.d3.test - og.test$y)^2), # ridge - all og covariates and [depth 3 ,nodes 15 ]
  mean((rid.d4.test - og.test$y)^2), # ridge - all og covariates and [depth 4 ,nodes 31 ]
  mean((rid.all.test - og.test$y)^2), # ridge - all og covariates and all nodes
  mean((tree.def.test - og.test$y)^2), # def tree
  mean((rf.test - og.test$y)^2) # random forest
)
mod.names <- c("oracle",
               "lin","lin.psis","lin.d0","lin.d1","lin.d2","lin.d3","lin.d4","lin.all",
               "lasso","lasso.psis","lasso.d0","lasso.d1","lasso.d2","lasso.d3","lasso.d4","lasso.all",
               "ridge","ridge.psis","ridge.d0","ridge.d1","ridge.d2","ridge.d3","ridge.d4","ridge.all",
               "tree.def", "rf"
)
rownames(mse.all) <- mod.names


# list of vectors with non-zero coefficient names
cov.kept.list <- list(
  rownames(coef(las.og.train, s = las.og.train$lambda.1se))[which(coef(las.og.train, s = las.og.train$lambda.1se) != 0)],
  rownames(coef(las.psis.og.train, s = las.psis.og.train$lambda.1se))[which(coef(las.psis.og.train, s = las.psis.og.train$lambda.1se) != 0)],
  rownames(coef(las.d0.train, s = las.d0.train$lambda.1se))[which(coef(las.d0.train, s = las.d0.train$lambda.1se) != 0)],
  rownames(coef(las.d1.train, s = las.d1.train$lambda.1se))[which(coef(las.d1.train, s = las.d1.train$lambda.1se) != 0)],
  rownames(coef(las.d2.train, s = las.d2.train$lambda.1se))[which(coef(las.d2.train, s = las.d2.train$lambda.1se) != 0)],
  rownames(coef(las.d3.train, s = las.d3.train$lambda.1se))[which(coef(las.d3.train, s = las.d3.train$lambda.1se) != 0)],
  rownames(coef(las.d4.train, s = las.d4.train$lambda.1se))[which(coef(las.d4.train, s = las.d4.train$lambda.1se) != 0)],
  rownames(coef(las.all.train, s = las.all.train$lambda.1se))[which(coef(las.all.train, s = las.all.train$lambda.1se) != 0)]
)

names(cov.kept.list) <- c("all og", "all nodes",
                              "all og + 0 depth nodes", "all og + 1 depth nodes",
                              "all og + 2 depth nodes", "all og + 3 depth nodes",
                              "all og + 4 depth nodes", "all og + all nodes")
