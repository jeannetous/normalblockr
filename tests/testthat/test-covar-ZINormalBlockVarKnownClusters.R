###############################################################################
## ZINormalBlockVarKnownClusters with a covariate-dependent zero-inflation
## layer, i.e. X0 != NULL (R/NormalBlockData.R): the excess-zero probability
## is regressed on its own design matrix rather than fit as one intercept per
## variable. Known clusters, diagonal covariance, with and without sparsity.
testdata <- readRDS("testdata/testdata_normal_covardep_zi.RDS")
Y  <- testdata$Y
X  <- testdata$X
X0 <- testdata$X0
C  <- testdata$parameters$C
data  <- NormalBlockData$new(Y, X, X0 = X0)

test_that("zero inflated normal block with diagonal residual covariance and known clusters", {
  ## Diagonal model
  model <- ZINormalBlockVarKnownClusters$new(data, C)
  model$optimize()
  expect_lt(model$BIC, 4300)
  expect_gt(model$loglik, -2000)
  expect_lt(Metrics::rmse(model$fitted, Y), 0.75)

  model <- ZINormalBlockVarKnownClusters$new(data, C, sparsity = 2)
  model$optimize()
  expect_gt(model$loglik, -2000)
  expect_lt(Metrics::rmse(model$fitted, Y), 0.75)
})
