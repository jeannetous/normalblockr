###############################################################################
###############################################################################
## Use pre-save testdata (seed are hard to handle in testhat)
testdata <- readRDS("testdata/testdata_normal.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$param$C
Q <- ncol(C)

###############################################################################
###############################################################################

test_that("NB_unknown: check dimensions, optimization and field access", {
  model <- NB_unknown$new(Y, X, c(3, 6, 4, 5))
  model$optimize()
  best_model <- model$get_best_model("BIC")
  true_model <- model$get_model(Q)
  expect_equal(true_model$Q, Q)
  expect_equal(model$n, nrow(Y))
  expect_equal(model$p, ncol(Y))
  expect_equal(model$d, ncol(X))
  expect_lt(best_model$BIC, 5450)
  expect_gt(true_model$loglik, -2616)
  model_sparse <- NB_unknown$new(Y, X, c(3, 6, 4, 5), c(0.01, 0.04, 0.02, 0.03))
  model_sparse$optimize()
})
