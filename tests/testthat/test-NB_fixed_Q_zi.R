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
  model <- NB_fixed_Q_zi$new(Y, X, Q)
  model$optimize(niter = 60)
  params <- model$model_par
  expect_equal(model$n, nrow(Y))
  expect_equal(model$p, ncol(Y))
  expect_equal(model$d, ncol(X))
  expect_lt(model$BIC, 295100)
  expect_gt(model$loglik,  -52800)
  model_sparse <- NB_fixed_Q_zi$new(Y, X, Q, sparsity = 0.05)
  model_sparse$optimize(niter = 60)
})
