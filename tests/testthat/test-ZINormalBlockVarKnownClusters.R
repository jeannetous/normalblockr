###############################################################################
## ZINormalBlockVarKnownClusters with intercept-only zero-inflation (formula's
## right-hand side is `| 1`, unlike test-covar-ZINormalBlockVarKnownClusters.R's
## X0-dependent one): diagonal and spherical residual covariance, with and
## without sparsity, and the moment-based heuristic fit.
testdata <- readRDS("testdata/testdata_normal_zi.RDS")
Y <- testdata$Y
X <- testdata$X ; X <- cbind(X, rnorm(nrow(X))) ; colnames(X) <- c("X1", "X2")
C <- testdata$parameters$C
data <- NormalBlockData$new(Y, X, formula = ~ 0 + X1 | 1)

test_that("zero inflated normal block with diagonal residual covariance and known clusters", {
  ## Diagonal model
  model <- ZINormalBlockVarKnownClusters$new(data, C)
  model$optimize()
  expect_lt(model$BIC, 5600)
  expect_gt(model$loglik, -2700)
  expect_lt(Metrics::rmse(model$fitted, Y), 0.9)

  model <- ZINormalBlockVarKnownClusters$new(data, C, sparsity = 2)
  model$optimize()
  expect_gt(model$loglik, -2700)
  expect_lt(Metrics::rmse(model$fitted, Y), 0.9)
})

test_that("zero inflated normal block with spherical residual covariance and known clusters", {
  ## Spherical model
  ctrl <- NB_control(noise_covariance = "spherical")
  model <- ZINormalBlockVarKnownClusters$new(data, C, control = ctrl)
  model$optimize()
  expect_gt(model$loglik, -2700)
  expect_lt(Metrics::rmse(model$fitted, Y), 0.9)

  model <- ZINormalBlockVarKnownClusters$new(data, C, sparsity = 2, control = ctrl)
  model$optimize()
  expect_gt(model$loglik, -2700)
  expect_lt(Metrics::rmse(model$fitted, Y), 0.9)
})

test_that("zero inflated normal block with known clusters, heuristic", {
  model <- normalblockr:::ZINormalBlockVarKnownClusters$new(data, C, sparsity = 0.05,
                               control = NB_control(heuristic = TRUE))
  model$optimize()
  expect_lt(Metrics::rmse(model$fitted, Y), 2)
})
