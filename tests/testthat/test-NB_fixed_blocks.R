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

test_that("NB_fixed_blocks: check dimensions, optimization and field access", {
  model <- NB_fixed_blocks$new(Y, X, C)
  model$optimize()
  params <- model$model_par
  expect_equal(model$n, nrow(Y))
  expect_equal(model$p, ncol(Y))
  expect_equal(model$d, ncol(X))
  expect_lt(model$BIC, 36100)
  expect_gt(model$loglik, -17700)
  model_sparse <- NB_fixed_blocks$new(Y, X, C, sparsity <- 0.2)
  model_sparse$optimize()
})
