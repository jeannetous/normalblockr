###############################################################################
###############################################################################
## Use pre-save testdata (seed are hard to handle in testhat)
testdata <- readRDS("testdata/testdata_normal.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C

###############################################################################

test_that("normal block with diagonal residual covariance and known clusters", {
  data <- normal_data$new(Y, X)

  model <- normalblockr:::NB_fixed_blocks$new(data, C)
  model$optimize()
  expect_gt(model$loglik, -2600)
  expect_lt(Metrics::rmse(model$fitted, Y), 1)

  model_sparse <- normalblockr:::NB_fixed_blocks$new(data, C, sparsity = 0.05)
  model_sparse$optimize()
  expect_gt(model_sparse$loglik, -2600)
  expect_lt(Metrics::rmse(model_sparse$fitted, Y), 1)

  model$sparsity <- 0.05
  model$optimize()
  expect_equal(model_sparse$loglik, model$loglik, tolerance = 1e-2)
  expect_equal(model_sparse$model_par$OmegaQ, model$model_par$OmegaQ, tolerance = 1e-2)

})

test_that("normal block with spherical residual covariance and known clusters", {
  data <- normal_data$new(Y, X)

  ## Spherical model
  ctrl <- NB_control(noise_covariance = "spherical")
  model <- NB_fixed_blocks$new(data, C, control = ctrl)
  model$optimize()
  expect_gt(model$loglik, -2630)
  expect_lt(Metrics::rmse(model$fitted, Y), 1)

  model <- NB_fixed_blocks$new(data, C, sparsity = 0.01, control = ctrl)
  model$optimize()
  expect_gt(model$loglik, -2630)
  expect_lt(Metrics::rmse(model$fitted, Y), 1)
})

test_that("normal block with known clusters, heuristic", {
  data <- normal_data$new(Y, X)
  model <- NB_fixed_blocks$new(data, C, sparsity = 0.05,
                                              control = NB_control(heuristic = TRUE))
  model$optimize()
  expect_lt(Metrics::rmse(model$fitted, Y), 2.4)
})
