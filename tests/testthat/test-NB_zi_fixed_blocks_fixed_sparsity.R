###############################################################################
###############################################################################
testdata <- readRDS("testdata/testdata_normal_zi.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C

###############################################################################
###############################################################################
test_that("NB_zi_fixed_blocks: check dimensions, optimization and field access", {
  data  <- normal_data$new(Y, X)
  model <- NB_zi_fixed_blocks_fixed_sparsity_diagonal$new(data, C, penalty = 0.2)
  model$optimize()
  expect_lt(model$BIC, 6460)
  expect_gt(model$loglik,  -3100)
  model <- NB_zi_fixed_blocks_fixed_sparsity_diagonal$new(data, C, penalty = 0.2,
                                                          control = NB_control(inference_method = "heuristic"))
  model$optimize()
  expect_lt(Metrics::rmse(model$fitted, Y), 1.9)
})
