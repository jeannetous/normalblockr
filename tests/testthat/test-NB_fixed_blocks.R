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

test_that("NB_fixed_blocks_diagonal: check dimensions, optimization and field access", {
  model <- normal_block(Y, X, C, noise_cov = "diagonal")
  params <- model$model_par
  expect_equal(model$n, nrow(Y))
  expect_equal(model$p, ncol(Y))
  expect_equal(model$d, ncol(X))
  expect_lt(model$BIC, 5391)
  expect_gt(model$loglik, -2595)
  model_sparse <- normal_block(Y, X, C, sparsity = 0.1, noise_cov = "diagonal")
})

test_that("NB_fixed_blocks_spherical: check dimensions, optimization and field access", {
  model <- normal_block(Y, X, C, noise_cov = "spherical")
  params <- model$model_par
  expect_equal(model$n, nrow(Y))
  expect_equal(model$p, ncol(Y))
  expect_equal(model$d, ncol(X))
  expect_lt(model$BIC, 5460)
  expect_gt(model$loglik, -2639)
  model_sparse <- normal_block(Y, X, C, sparsity = 0.1, noise_cov = "spherical")
})
