###############################################################################
###############################################################################
testdata <- readRDS("testdata/testdata_normal_zi.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C
data <- NBData$new(Y, X)

test_that("zero inflated normal block with diagonal residual covariance and known clusters", {
  ## Diagonal model
  model <- ZINB_fixed_blocks$new(data, C)
  model$optimize()
  expect_lt(model$BIC, 3500)
  expect_gt(model$loglik, -1600)
  expect_lt(Metrics::rmse(model$fitted, Y), 0.7)

  model <- ZINB_fixed_blocks$new(data, C, sparsity = 2)
  model$optimize()
  expect_gt(model$loglik, -1650)
  expect_lt(Metrics::rmse(model$fitted, Y), 0.7)
})

test_that("zero inflated normal block with spherical residual covariance and known clusters", {
  ## Spherical model
  ctrl <- NB_control(noise_covariance = "spherical")
  model <- ZINB_fixed_blocks$new(data, C, control = ctrl)
  model$optimize()
  expect_gt(model$loglik, -1610)
  expect_lt(Metrics::rmse(model$fitted, Y), 0.7)

  model <- ZINB_fixed_blocks$new(data, C, sparsity = 2, control = ctrl)
  model$optimize()
  expect_gt(model$loglik, -1650)
  expect_lt(Metrics::rmse(model$fitted, Y), 0.7)
})

test_that("zero inflated normal block with known clusters, heuristic", {
  model <- normalblockr:::ZINB_fixed_blocks$new(data, C, sparsity = 0.05,
                               control = NB_control(heuristic = TRUE))
  model$optimize()
  expect_lt(Metrics::rmse(model$fitted, Y), 2)
})
