###############################################################################
## NormalBlockVarCollectionClustersSparsity (R/NormalBlockVarCollectionClustersSparsity.R):
## a collection crossing several q with a sparsity path at each, both the
## plain VEM path and the moment-based heuristic one.
## Use pre-save testdata (seed are hard to handle in testhat)
testdata <- readRDS("testdata/testdata_normal.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C ; q <- ncol(C)

data  <- normalblockr:::NormalBlockData$new(Y, X)

test_that("normal block with spherical residual covariance and unknown clusters", {
  model <- normalblockr:::NormalBlockVarCollectionClustersSparsity$new(data, c(2,3,4))
  model$optimize()
  model_BIC <- model$get_best_model("BIC")
  expect_lt(model_BIC$BIC, 5509)
  expect_gt(model_BIC$loglik, -2656)
})

test_that("normal block with spherical residual covariance and unknown clusters heuristic", {
  data <- NormalBlockData$new(Y, X)
  model <- normalblockr:::NormalBlockVarCollectionClustersSparsity$new(data, c(2,3,4),
                                              control = NB_control(heuristic = TRUE))
  model$optimize()
  model_1 <- model$get_model(3, 0.1)
  expect_lt(Metrics::rmse(model_1$fitted, Y), 2.4)
})
