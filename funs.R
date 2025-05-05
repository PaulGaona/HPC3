#
##
### delete.NULLs , a function to delete null/empty entries in a list
##
#
deleteNULL.Lists  <-  function(x.list){
  x.list[unlist(lapply(x.list, length) != 0)]
}

################################################################################
################################################################################
############################# HAVE TO FORK RFPlus ##############################
################################################################################
################################################################################
# all sparse simulation functions:
# this is a collection of functions that simulate data from different sparse models
{

  ###### 4 mdi+ examples

  # 1. Linear model: E[Y|X] = sum_{j=1}^5 X_j
  lin.model.sparse <- function(n, p, sd = 1) {
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    y <- rowSums(X[, 1:5]) + rnorm(n, mean = 0, sd = sd)

    data <- data.frame(y,X)
    colnames(data)[-1] <- paste0("X", 1:p)
    return(data)
  }
  # 2. Locally-spiky-sparse (LSS) model: E[Y|X] = sum_{m=1}^3 1(X_{2m-1} > 0) * 1(X_{2m} > 0)
  lss.model.sparse <- function(n, p, sd = 1) {

    X <- matrix(rnorm(n * p), nrow = n, ncol = p)

    y <- 0
    for (m in 1:3) {
      y <- y + (X[, 2*m-1] > 0) * (X[, 2*m] > 0)
    }
    y <- y + rnorm(n, mean = 0, sd = sd)

    data <- data.frame(y,X)
    colnames(data)[-1] <- paste0("X", 1:p)
    return(data)
  }
  # 3. Polynomial interaction model: E[Y|X] = sum_{m=1}^3 X_{2m-1} + sum_{m=1}^3 X_{2m-1} * X_{2m}
  poly.inter.sparse <- function(n, p, sd = 1) {

    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    y <- 0
    for (m in 1:3) {
      y <- y + X[, 2*m-1] + X[, 2*m-1] * X[, 2*m]
    }
    y <- y + rnorm(n, mean = 0, sd = sd)

    data <- data.frame(y,X)
    colnames(data)[-1] <- paste0("X", 1:p)
    return(data)
  }
  # 4. Linear + LSS model: E[Y|X] = sum_{m=1}^3 X_{2m-1} + sum_{m=1}^3 1(X_{2m-1} > 0) * 1(X_{2m} > 0)
  lin.lss.sparse <- function(n, p, sd = 1) {

    X <- matrix(rnorm(n * p), nrow = n, ncol = p)

    y <- 0
    for (m in 1:3) {
      y <- y + X[, 2*m-1] + (X[, 2*m-1] > 0) * (X[, 2*m] > 0)
    }
    y <- y + rnorm(n, mean = 0, sd = sd)

    # Return data frame
    data <- data.frame(y,X)
    colnames(data)[-1] <- paste0("X", 1:p)
    return(data)
  }

  ###### 2 breiman examples

  # 5. Breiman CART Simulation: E[Y|X] = X1 + 2*X2 + X3 if X1 = 1, X5 + 2*X6 + X7 if X1 = -1
  Breiman.CART.sim.b <- function(n, p, sd = 1) {

    X <- matrix(0, nrow = n, ncol = p)
    X[, 1] <- sample(c(-1, 1), size = n, replace = TRUE)
    for (j in 2:10) {
      X[, j] <- sample(c(-1, 0, 1), size = n, replace = TRUE, prob = c(1/3, 1/3, 1/3))
    }


    if (p > 10) {
      X[, 11:p] <- matrix(rnorm(n * (p - 10)), nrow = n, ncol = (p - 10))
    }

    Z <- rnorm(n, mean = 0, sd = sd)
    y <- ifelse(X[, 1] == 1,
                X[, 2] + 2 * X[, 3] + X[, 4] + Z,
                X[, 5] + 2 * X[, 6] + X[, 7] + Z)


    data <- data.frame(y,X)
    colnames(data)[-1] <- paste0("X", 1:p)
    return(data)
  }
  # 6. Breiman CART Simulation: E[Y|X] = 4*X2 + 3*X3 + 2*X4 if X1 = 1, 4*X5 + 3*X6 + 2*X7*X8 if X1 = -1
  Breiman.CART.sim2b <- function(n, p, sd = 1) {
    X <- matrix(0, nrow = n, ncol = p)

    X[, 1] <- sample(c(0, 1), size = n, replace = TRUE)
    for (j in 2:10) {
      X[, j] <- rnorm(n)
    }

    if (p > 10) {
      X[, 11:p] <- matrix(rnorm(n * (p - 10)), nrow = n, ncol = (p - 10))
    }

    Z <- rnorm(n, mean = 0, sd = sd)
    y <- ifelse(X[, 1] == 1,
                4 * X[, 2] + 3 * X[, 3] + 2 * X[, 4] + Z,
                4 * X[, 5] + 3 * X[, 6] + 2 * X[, 7] * X[, 8] + Z)

    data <- data.frame(y,X)
    colnames(data)[-1] <- paste0("X", 1:p)
    return(data)
  }

  ###### 3 freidman examples

  # 7. Friedman 1: E[Y|X] = 10*sin(pi*X1*X2) + 20*(X3-0.5)^2 + 10*X4 + 5*X5
  friedman1.sparse <- function(n, p, sd = 1) {
    X <- matrix(runif(p * n), ncol = p)
    y <- 10 * sin(pi * X[, 1] * X[, 2]) +
      20 * (X[, 3] - 0.5)^2 +
      10 * X[, 4] +
      5 * X[, 5] +
      rnorm(n, mean = 0, sd = sd)
    data <- data.frame(y,X)
    colnames(data)[-1] <- paste0("X", 1:p)
    return(data)
  }
  # 8. Friedman 2: E[Y|X] = sqrt(X1^2 + (X2*X3 - 1/(X2*X4))^2)
  friedman2.sparse <- function(n, p, sd = 125) {
    X <- matrix(runif(p * n), ncol = p)
    x.sparse <- cbind(
      runif(n, min = 0, max = 100),
      runif(n, min = 40 * pi, max = 560 * pi),
      runif(n, min = 0, max = 1),
      runif(n, min = 1, max = 11)
    )
    y <- sqrt(
      x.sparse[, 1]^2 +
        (x.sparse[, 2] * x.sparse[, 3] - 1 / (x.sparse[, 2] * x.sparse[, 4]))^2
    ) + rnorm(n, mean = 0, sd = sd)
    X[, 1:4] <- x.sparse
    data <- data.frame(y,X)
    colnames(data)[-1] <- paste0("X", 1:p)
    return(data)
  }
  # 9. Friedman 3: E[Y|X] = atan((X2*X3 - 1/(X2*X4))/X1)
  friedman3.sparse <- function(n, p, sd = 0.5) {
    X <- matrix(runif(p * n), ncol = p)
    x.sparse <- cbind(
      runif(n, min = 0, max = 100),
      runif(n, min = 40 * pi, max = 560 * pi),
      runif(n, min = 0, max = 1),
      runif(n, min = 1, max = 11)
    )
    y <- atan(
      (x.sparse[, 2] * x.sparse[, 3] - 1 / (x.sparse[, 2] * x.sparse[, 4])) / x.sparse[, 1]
    ) + rnorm(n, mean = 0, sd = sd)
    X[, 1:4] <- x.sparse
    data <- data.frame(y,X)
    colnames(data)[-1] <- paste0("X", 1:p)
    return(data)
  }
}

################################################################################
############################## oracles
# this is a collection of functions that define the oracle models from the sparse
# models above
{
  # mdi+ oracles
  # 1a. Linear model
  lin.oracle <- function(data) {
    lm(y ~ 0+X1 + X2 + X3 + X4 + X5, data = data)
  }
  # 1b. Locally-spiky-sparse (LSS) model
  lss.oracle <- function(data) {
    lm(y ~ 0+I(as.numeric(X1 > 0 & X2 > 0)) +
         I(as.numeric(X3 > 0 & X4 > 0)) +
         I(as.numeric(X5 > 0 & X6 > 0)),
       data = data)
  }
  # 1c. Polynomial interaction
  poly.inter.oracle <- function(data) {
    lm(y ~ 0+X1 + X3 + X5 + X1:X2 + X3:X4 + X5:X6, data = data)
  }
  # 1d. Linear + LSS model
  lin.lss.oracle <- function(data) {
    lm(y ~ 0+X1 + X3 + X5 +
         I(as.numeric(X1 > 0 & X2 > 0)) +
         I(as.numeric(X3 > 0 & X4 > 0)) +
         I(as.numeric(X5 > 0 & X6 > 0)),
       data = data)
  }
  # Breiman oracles:
  # 2b. Breiman CART Simulation
  breiman.oracle <- function(data) {
    lm(y ~ 0+I(as.numeric(X1 == 1) * X2) +
         I(as.numeric(X1 == 1) * X3) +
         I(as.numeric(X1 == 1) * X4) +
         I(as.numeric(X1 == -1) * X5) +
         I(as.numeric(X1 == -1) * X6) +
         I(as.numeric(X1 == -1) * X7),
       data = data)
  }
  # 2c. Breiman CART Simulation (continuous and interaction)
  breiman.oracle2 <- function(data) {
    lm(y ~ 0+I(as.numeric(X1 == 1) * X2) +
         I(as.numeric(X1 == 1) * X3) +
         I(as.numeric(X1 == 1) * X4) +
         I(as.numeric(X1 == 0) * X5) +
         I(as.numeric(X1 == 0) * X6) +
         I(as.numeric(X1 == 0) * X7 * X8),
       data = data)
  }

  # freidman oracles
  # 3a. Friedman 1
  friedman1.oracle <- function(data) {
    lm(y ~ 0+I(sin(pi * X1 * X2)) +
         I((X3 - 0.5)^2) +
         X4 + X5,
       data = data)
  }
  # 3b. Friedman 2
  friedman2.oracle <- function(data) {
    lm(y ~ 0+I(sqrt(X1^2 + (X2 * X3 - 1/(X2 * X4))^2)),
       data = data)
  }
  # 3c. Friedman 3
  friedman3.oracle <- function(data) {
    lm(y ~ 0+I(atan((X2 * X3 - 1/(X2 * X4)) / X1)),
       data = data)
  }
}

################################################################################
# this is a collection of functions that simulate data from different sparse models

sparse.sims <- function(n = n, p = p, sd = sd,setting = setting,
                        prop = c(0.8, 0.2), validation.set = FALSE){
  # propertion check
  if (sum(prop) != 1) {
    stop("The proportions (prop) must sum to 1.")
  }

  # sample sizes for training, testing, and validation sets
  n.train <- floor(prop[1] * n)
  n.test <- floor(prop[2] * n)
  n.valid <- n - n.train - n.test

  # indice assignment: training, testing, and validation sets
  all.indices <- sample(1:n, size = n, replace = FALSE)
  train.indices <- all.indices[1:n.train]
  test.indices <- all.indices[(n.train + 1):(n.train + n.test)]
  val.indices <- all.indices[(n.train + n.test + 1):n]

  results <- list()

  # 1. Linear model: E[Y|X] = sum_{j=1}^5 X_j
  if (1 %in% setting) {
    data <- lin.model.sparse(n, p, sd = sd)

    train.data <- data[train.indices, ]
    test.data <- data[test.indices, ]
    val.data <- if (validation.set) data[val.indices, ] else NULL

    results$lin <- list(train = train.data, test = test.data, val = val.data)
  }
  # 2. Locally-spiky-sparse (LSS) model: E[Y|X] = sum_{m=1}^3 1(X_{2m-1} > 0) * 1(X_{2m} > 0)
  if (2 %in% setting) {
    data <- lss.model.sparse(n, p, sd = sd)

    train.data <- data[train.indices, ]
    test.data <- data[test.indices, ]
    val.data <- if (validation.set) data[val.indices, ] else NULL

    results$lss <- list(train = train.data, test = test.data, val = val.data)
  }
  # 3. Polynomial interaction model: E[Y|X] = sum_{m=1}^3 X_{2m-1} + sum_{m=1}^3 X_{2m-1} * X_{2m}
  if (3 %in% setting) {
    data <- poly.inter.sparse(n, p, sd = sd)

    train.data <- data[train.indices, ]
    test.data <- data[test.indices, ]
    val.data <- if (validation.set) data[val.indices, ] else NULL

    results$poly.inter <- list(train = train.data, test = test.data, val = val.data)
  }
  # 4. Linear + LSS model: E[Y|X] = sum_{m=1}^3 X_{2m-1} + sum_{m=1}^3 1(X_{2m-1} > 0) * 1(X_{2m} > 0)
  if (4 %in% setting) {
    data <- lin.lss.sparse(n, p, sd = sd)

    train.data <- data[train.indices, ]
    test.data <- data[test.indices, ]
    val.data <- if (validation.set) data[val.indices, ] else NULL

    results$lin.lss <- list(train = train.data, test = test.data, val = val.data)
  }
  # 5. Breiman CART Simulation: E[Y|X] = X1 + 2*X2 + X3 if X1 = 1, X5 + 2*X6 + X7 if X1 = -1
  if (5 %in% setting) {
    data <- Breiman.CART.sim.b(n, p, sd = sd)

    train.data <- data[train.indices, ]
    test.data <- data[test.indices, ]
    val.data <- if (validation.set) data[val.indices, ] else NULL

    results$Breiman.CART.b <- list(train = train.data, test = test.data, val = val.data)
  }
  # 6. Breiman CART Simulation: E[Y|X] = 4*X2 + 3*X3 + 2*X4 if X1 = 1, 4*X5 + 3*X6 + 2*X7*X8 if X1 = -1
  if (6 %in% setting) {
    data <- Breiman.CART.sim2b(n, p, sd = sd)

    train.data <- data[train.indices, ]
    test.data <- data[test.indices, ]
    val.data <- if (validation.set) data[val.indices, ] else NULL

    results$Breiman.CART.2b <- list(train = train.data, test = test.data, val = val.data)
  }
  # 7. Friedman 1: E[Y|X] = 10*sin(pi*X1*X2) + 20*(X3-0.5)^2 + 10*X4 + 5*X5
  if (7 %in% setting) {
    data <- friedman1.sparse(n, p, sd = sd)

    train.data <- data[train.indices, ]
    test.data <- data[test.indices, ]
    val.data <- if (validation.set) data[val.indices, ] else NULL

    results$friedman1 <- list(train = train.data, test = test.data, val = val.data)
  }
  # 8. Friedman 2: E[Y|X] = sqrt(X1^2 + (X2*X3 - 1/(X2*X4))^2)
  if (8 %in% setting) {
    data <- friedman2.sparse(n, p, sd = sd)

    train.data <- data[train.indices, ]
    test.data <- data[test.indices, ]
    val.data <- if (validation.set) data[val.indices, ] else NULL

    results$friedman2 <- list(train = train.data, test = test.data, val = val.data)
  }
  # 9. Friedman 3: E[Y|X] = atan((X2*X3 - 1/(X2*X4))/X1)
  if (9 %in% setting) {
    data <- friedman3.sparse(n, p, sd = sd)

    train.data <- data[train.indices, ]
    test.data <- data[test.indices, ]
    val.data <- if (validation.set) data[val.indices, ] else NULL

    results$friedman3 <- list(train = train.data, test = test.data, val = val.data)
  }

  results <- deleteNULL.Lists(results)
  return(results)
}
# this is a collection of functions that simulate the different oracle
# models
oracle.mods <- function(data = train,setting){
  oracles <- list()

  # 1. Linear model
  if (1 %in% setting) {
    oracles$lin <- lin.oracle(data)
  }
  # 2. Locally-spiky-sparse (LSS) model
  if (2 %in% setting) {
    oracles$lss <- lss.oracle(data)
  }
  # 3. Polynomial interaction model
  if (3 %in% setting) {
    oracles$poly.inter <- poly.inter.oracle(data)
  }
  # 4. Linear + LSS model
  if (4 %in% setting) {
    oracles$lin.lss <- lin.lss.oracle(data)
  }
  # 5. Breiman CART Simulation
  if (5 %in% setting) {
    oracles$Breiman.CART.b <- breiman.oracle(data)
  }
  # 6. Breiman CART Simulation
  if (6 %in% setting) {
    oracles$Breiman.CART.2b <- breiman.oracle2(data)
  }
  # 7. Friedman 1
  if (7 %in% setting) {
    oracles$friedman1 <- friedman1.oracle(data)
  }
  # 8. Friedman 2
  if (8 %in% setting) {
    oracles$friedman2 <- friedman2.oracle(data)
  }
  # 9. Friedman 3
  if (9 %in% setting) {
    oracles$friedman3 <- friedman3.oracle(data)
  }

  oracles <- deleteNULL.Lists(oracles)
  return(oracles)
}

################################################################################

# this is a cross-validated function getting a modified MSE

cv.model.fun <- function(og.train = og.train,
                         psis.def.og.train = psis.def.og.train,
                         num.nodes.depth = num.nodes.depth,
                         folds = folds) {
  # Initialize result dataframe
  cv.res <- data.frame(depth = character(), nodes = numeric(), mse = numeric(), sd.mse = numeric())

  # Extract node numbers from column names
  node.nums <- as.numeric(gsub(".node", "", colnames(psis.def.og.train)))

  for (i in 1:length(num.nodes.depth)) {
    depth.name <- paste0("d", i-1)

    # Select columns where node number is less than the current depth threshold
    sltd.cols <- which(node.nums < num.nodes.depth[i])
    # Create dataset for current depth only
    if (length(sltd.cols) > 0) {
      curr.depth.df <- cbind(og.train, psis.def.og.train[, sltd.cols, drop = FALSE])
    } else {
      curr.depth.df <- og.train
    }

    fold.mse <- numeric(max(folds))
    # For each fold
    for (k in seq(max(folds))) {
      # Training and validation indices
      cv.train.indices <- which(folds != k)
      cv.valid.indices <- which(folds == k)

      cv.train.data <- curr.depth.df[cv.train.indices, ]
      cv.valid.data <- curr.depth.df[cv.valid.indices, ]

      cv.mod <- lm(y ~ ., data = cv.train.data)
      cv.pred <- predict(cv.mod, newdata = cv.valid.data)

      # Calculate MSE for this fold
      num <- length(sltd.cols)
      den <- num.nodes.depth[i]

      fold.mse[k] <- (sum((cv.valid.data$y - cv.pred)^2)) / (length(cv.pred)*(num/den))

      # Clean up to free memory
      rm(cv.train.data, cv.valid.data, cv.mod, cv.pred)
      gc(verbose = FALSE)
    }

    # Average MSE across folds
    cv.avg.mse <- mean(fold.mse)
    cv.sd.mse <- sd(fold.mse)

    cv.res <- rbind(cv.res, data.frame(
      depth = depth.name,
      nodes = length(sltd.cols),
      mse = cv.avg.mse,
      sd.mse = cv.sd.mse
    ))

    # Clean up to free memory
    rm(curr.depth.df, fold.mse)
    gc(verbose = FALSE)
  }

  # Return results
  return(list(lin.cv.res = cv.res))
}
