###############################################################################
###############################################################################
## Use pre-save testdata (seed are hard to handle in testhat)
testdata <- readRDS("testdata/testdata_normal.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C ; Q <- ncol(C)


###############################################################################
###############################################################################

test_that("normal: check dimensions and field access", {
  data <- normal_data$new(Y, X)
  model <- NB_fixed_Q_diagonal$new(data, Q)
  model$optimize()
  expect_gt(model$loglik, -2700)
  model <- NB_fixed_Q_diagonal$new(data, Q, penalty = 0.05,
                                   control = NB_control(inference_method = "heuristic"))
  model$optimize()
  expect_equal(matching_group_scores(model$clustering, get_clusters(C)), 1)
  expect_lt(Metrics::rmse(model$fitted, Y), 2.4)
})
