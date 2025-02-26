###############################################################################
###############################################################################
testdata <- readRDS("testdata/testdata_normal_zi.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C

###############################################################################
###############################################################################


test_that("zero inflated normal block with diagonal residual covariance and known clusters", {
  data <- normal_data$new(Y, X)

  model <- NB_zi_fixed_blocks$new(data, C)
  model$optimize()
  expect_lt(model$BIC, 5400)
  expect_gt(model$loglik, -2600)
  expect_lt(Metrics::rmse(model$fitted, Y), 0.7)

  model <- NB_zi_fixed_blocks$new(data, C, penalty = 0.1)
  model$optimize()
  expect_gt(model$loglik, -2600)
  expect_lt(Metrics::rmse(model$fitted, Y), 1)

})

test_that("zero inflated normal block with spherical residual covariance and known clusters", {
  data <- normal_data$new(Y, X)

  ## Spherical model
  ctrl <- NB_control(noise_covariance = "spherical")
  model <- NB_zi_fixed_blocks$new(data, C, control = ctrl)
  model$optimize()
  expect_gt(model$loglik, -2600)
  expect_lt(Metrics::rmse(model$fitted, Y), 0.7)

  model <- NB_zi_fixed_blocks$new(data, C, penalty = 0.1, control = ctrl)
  model$optimize()
  expect_gt(model$loglik, -2600)
  expect_lt(Metrics::rmse(model$fitted, Y), 0.7)
})

test_that("zero inflated normal block with known clusters, heuristic", {
  data <- normal_data$new(Y, X)
  model <- NB_zi_fixed_blocks$new(data, C, penalty = 0.05,
                               control = NB_control(heuristic = TRUE))
  model$optimize()
  expect_lt(Metrics::rmse(model$fitted, Y), 2)
})
