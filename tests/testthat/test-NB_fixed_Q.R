###############################################################################
###############################################################################

## Use pre-save testdata (seed are hard to handle in testhat)
testdata <- readRDS("testdata/testdata_normal.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$C
Q <- ncol(C)

###############################################################################
###############################################################################

test_that("NB_fixed_Q: check dimensions, optimization and field access", {
  expect_true(inherits(model <- NB_fixed_Q$new(Y, X, Q), "NB_fixed_Q"))
  expect_silent(model$optimize(niter = 60))
  params <- model$model_par
  expect_equal(model$n, nrow(Y))
  expect_equal(model$p, ncol(Y))
  expect_equal(model$d, ncol(X))
  expect_lt(model$BIC, 36250)
  expect_gt(model$loglik, -17690)
})

test_that("NB_fixed_Q: sparsity works", {
  expect_true(inherits(model_sparse <- NB_fixed_Q$new(Y, X, Q, penalty = 0.05), "NB_fixed_Q"))
  expect_silent(model_sparse$optimize(niter = 60))
})

test_that("NB_fixed_Q: works with Q=1", {
  expect_true(inherits(model <- NB_fixed_Q$new(Y, X, 1), "NB_fixed_Q"))
  expect_silent(model$optimize(niter = 60))
})
