###############################################################################
###############################################################################

## Use pre-save testdata (seed are hard to handle in testhat)
testdata <- readRDS("testdata/testdata_normal_zi.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C
Q <- ncol(C)

###############################################################################
###############################################################################

test_that("NB_fixed_blocks_zi: check dimensions, optimization and field access", {
  model <- NB_fixed_blocks_zi_diagonal$new(Y, X, C)
  model$optimize(niter = 60)
  params <- model$model_par
  expect_equal(model$n, nrow(Y))
  expect_equal(model$p, ncol(Y))
  expect_equal(model$d, ncol(X))
  expect_lt(model$BIC,  5394)
  expect_gt(model$loglik, -2553)
  model_sparse <- NB_fixed_blocks_zi_spherical$new(Y, X, C, penalty = 0.05)
  model_sparse$optimize(threshold = 1e-5)
  model_sparse <- normal_block(Y, X, C, sparsity = 10, zero_inflation = TRUE, control = normal_block_param(sparsity_weights = matrix(rep(1, Q * Q), nrow = Q)))
  expect_lt(sum(diag(model_sparse$model_par$omegaQ)), 0.5)
})
