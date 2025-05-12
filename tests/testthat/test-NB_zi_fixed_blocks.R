###############################################################################
###############################################################################
testdata <- readRDS("testdata/testdata_normal_zi.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C

test_that("zero inflated normal block with diagonal residual covariance and known clusters", {
  data <- normal_data$new(Y, X)

  model <- normalblockr:::NB_zi_fixed_blocks$new(data, C)
  model$optimize()
  expect_lt(model$BIC, 5410)
  expect_gt(model$loglik, -2600)
  expect_lt(Metrics::rmse(model$fitted, Y), 3)

  model <- normalblockr:::NB_zi_fixed_blocks$new(data, C, sparsity = 2)
  model$optimize()
  expect_gt(model$loglik, -2600)
  expect_lt(Metrics::rmse(model$fitted, Y), 3)
})

test_that("zero inflated normal block with spherical residual covariance and known clusters", {
  data <- normal_data$new(Y, X)

  ## Spherical model
  ctrl <- normalblockr:::NB_control(noise_covariance = "spherical")
  model <- normalblockr:::NB_zi_fixed_blocks$new(data, C, control = ctrl)
  model$optimize()
  expect_gt(model$loglik, -2600)
  expect_lt(Metrics::rmse(model$fitted, Y), 3)

  model <- normalblockr:::NB_zi_fixed_blocks$new(data, C, sparsity = 2, control = ctrl)
  model$optimize()
  expect_gt(model$loglik, -2620)
  expect_lt(Metrics::rmse(model$fitted, Y), 3)
})

test_that("zero inflated normal block with known clusters, heuristic", {
  data <- normal_data$new(Y, X)
  model <- normalblockr:::NB_zi_fixed_blocks$new(data, C, sparsity = 0.05,
                               control = NB_control(heuristic = TRUE))
  model$optimize()
  expect_lt(Metrics::rmse(model$fitted, Y), 4)
})
