###############################################################################
###############################################################################
## Use pre-save testdata (seed are hard to handle in testhat)
testdata <- readRDS("testdata/testdata_normal.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C ; Q <- ncol(C)


###############################################################################
###############################################################################

data  <- normalblockr:::NBData$new(Y, X)

test_that("normal block with spherical residual covariance and unknown clusters", {
  model <- normalblockr:::NB_unknown_Q_changing_sparsity$new(data, c(2,3,4))
  model$optimize()
  model_BIC <- model$get_best_model("BIC")
  expect_lt(model_BIC$BIC, 5509)
  expect_gt(model_BIC$loglik, -2656)
})

test_that("normal block with spherical residual covariance and unknown clusters heuristic", {
  data <- NBData$new(Y, X)
  model <- normalblockr:::NB_unknown_Q_changing_sparsity$new(data, c(2,3,4),
                                              control = NB_control(heuristic = TRUE))
  model$optimize()
  model_1 <- model$get_model(3, 0.1)
  expect_lt(Metrics::rmse(model_1$fitted, Y), 2.4)
})
