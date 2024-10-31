###############################################################################
###############################################################################
## Use pre-save testdata (seed are hard to handle in testhat)
testdata <- readRDS("testdata/testdata_sparse_normal.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$C
Q <- ncol(C)

###############################################################################
###############################################################################

test_that("NB_fixed_blocks_sparse: check dimensions, optimization and field access", {
  model <- normal_block(Y, X, C, sparsity = TRUE)
  model_BIC <- model$get_best_model("BIC")
  params <- model_BIC$model_par
  expect_equal(model_BIC$n, nrow(Y))
  expect_equal(model_BIC$p, ncol(Y))
  expect_equal(model_BIC$d, ncol(X))
  expect_lt(model_BIC$BIC, 100650)
  expect_gt(model_BIC$loglik, -50000)
  model_StARS <- model$get_best_model("StARS")
  expect_lt(model_StARS$penalty_term, 0.05)
})

