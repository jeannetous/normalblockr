###############################################################################
## NormalBlockVarUnknownClusters (R/NormalBlockVarUnknownClusters.R): the
## variational (VEM) counterpart of test-NormalBlockVarKnownClusters.R --
## diagonal and spherical covariance, with and without sparsity, the
## heuristic fit, and that the recovered clustering matches ground truth.
## Use pre-save testdata (seed are hard to handle in testhat)
testdata <- readRDS("testdata/testdata_normal.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C ; q <- ncol(C)

test_that("normal block with diagonal residual covariance and unknown clusters", {
  data <- NormalBlockData$new(Y, X)
  model <- NormalBlockVarUnknownClusters$new(data, q)
  model$optimize()
  expect_gt(model$loglik, -2620)
  expect_lt(Metrics::rmse(model$fitted, Y), 1.1)
  expect_equal(model$sparsity, 0)
  expect_equal(matching_group_scores(model$clustering, get_clusters(C)), 1)

  model_sparse <- NormalBlockVarUnknownClusters$new(data, q, sparsity = 0.05)
  model_sparse$optimize()
  expect_equal(model_sparse$sparsity, 0.05)

  model$sparsity <- 0.05
  expect_equal(model$sparsity, 0.05)
  model$optimize()
  expect_equal(model_sparse$loglik, model$loglik, tolerance = 1e-2)
})


test_that("normal block with spherical residual covariance and unknown clusters", {
  data <- NormalBlockData$new(Y, X)

  ctrl <- NB_control(noise_covariance = "spherical")
  model <- NormalBlockVarUnknownClusters$new(data, q, control = ctrl)
  model$optimize()
  expect_gt(model$loglik, -2650)
  expect_lt(Metrics::rmse(model$fitted, Y), 1.1)
  expect_equal(matching_group_scores(model$clustering, get_clusters(C)), 1)

  model <- NormalBlockVarUnknownClusters$new(data, q, sparsity = 0.05, control = ctrl)
  model$optimize()

})

test_that("normal block with unknown clusters, heuristic", {

  data <- NormalBlockData$new(Y, X)
  model <- NormalBlockVarUnknownClusters$new(data, q, sparsity = 0.05,
                          control = NB_control(heuristic = TRUE))
  model$optimize()
  expect_lt(Metrics::rmse(model$fitted, Y), 2.4)
  expect_equal(matching_group_scores(model$clustering, get_clusters(C)), 1)
})
