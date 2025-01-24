###############################################################################
###############################################################################
testdata <- readRDS("testdata/testdata_sparse_ZInormal.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$C
Q <- ncol(C)

###############################################################################
###############################################################################

test_that("NB_fixed_Q_zi_sparse: check dimensions, optimization and field access", {
  model <- normal_block(Y, X, Q, sparsity = TRUE, zero_inflation = TRUE)
  model_BIC <- model$get_best_model("BIC")
  params <- model_BIC$model_par
  expect_equal(model_BIC$n, nrow(Y))
  expect_equal(model_BIC$n, nrow(model_BIC$var_par$rho))
  expect_equal(model_BIC$p, ncol(Y))
  expect_equal(model_BIC$d, ncol(X))
  expect_equal(model_BIC$n_edges, 1)
  expect_lt(model_BIC$BIC, 100500)
  expect_gt(model_BIC$loglik, -50000)
})

