###############################################################################
###############################################################################
## Use pre-save testdata (seed are hard to handle in testhat)
testdata <- readRDS("testdata/testdata_normal.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C
Q <- ncol(C)

###############################################################################
###############################################################################

test_that("NB_fixed_blocks_sparse: check dimensions, optimization and field access", {
  model <- normal_block(Y, X, c(3,4,5), sparsity = TRUE)
  model_BIC <- model$get_best_model("BIC")
  params <- model_BIC$model_par
  expect_equal(model_BIC$n, nrow(Y))
  expect_equal(model_BIC$p, ncol(Y))
  expect_equal(model_BIC$d, ncol(X))
  expect_lt(model_BIC$BIC, 5510)
  expect_gt(model_BIC$loglik, -2656)
  model_Q <- model$get_model(Q)
  model_Q_StARS <- model_Q$get_best_model("StARS")
  expect_lt(model_Q_StARS$penalty_term, 0.45)
})

