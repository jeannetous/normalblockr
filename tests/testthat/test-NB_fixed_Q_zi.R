###############################################################################
###############################################################################
testdata <- readRDS("testdata/testdata_ZInormal.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$C
Q <- ncol(C)

###############################################################################
###############################################################################

test_that("NB_fixed_Q_zi: check dimensions, optimization and field access", {
  expect_true(inherits(model <- NB_fixed_Q_zi$new(Y, X, Q), "NB_fixed_Q_zi"))
  expect_silent(model$optimize())
  params <- model$model_par
  expect_equal(model$n, nrow(Y))
  expect_equal(model$p, ncol(Y))
  expect_equal(model$d, ncol(X))
  expect_lt(model$BIC, 100500)
  expect_gt(model$loglik, -49050)
})

test_that("NB_fixed_Q_zi: sparsity works", {
  expect_true(inherits(model_sparse <- NB_fixed_Q_zi$new(Y, X, Q, penalty = 0.05), "NB_fixed_Q_zi"))
  expect_silent(model_sparse$optimize(niter = 60))
})

test_that("NB_fixed_Q_zi: works with Q=1", {
  expect_true(inherits(model <- NB_fixed_Q_zi$new(Y, X, 1), "NB_fixed_Q_zi"))
  expect_silent(model$optimize(niter = 60))
})
