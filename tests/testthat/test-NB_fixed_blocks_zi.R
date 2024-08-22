###############################################################################
###############################################################################

## Use pre-save testdata (seed are hard to handle in testhat)
testdata <- readRDS("testdata/testdata_ZInormal.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$C
Q <- ncol(C)

###############################################################################
###############################################################################

test_that("NB_fixed_blocks_zi: check dimensions, optimization and field access", {
  model <- NB_fixed_blocks_zi$new(Y, X, C)
  model$optimize(niter = 60)
  params <- model$model_par
  expect_equal(model$n, nrow(Y))
  expect_equal(model$p, ncol(Y))
  expect_equal(model$d, ncol(X))
  expect_lt(model$BIC,  100150)
  expect_gt(model$loglik, -48950)
  model_sparse <- NB_fixed_blocks_zi$new(Y, X, C, sparsity = 0.05)
  model_sparse$optimize(threshold = 1e-5)
})
