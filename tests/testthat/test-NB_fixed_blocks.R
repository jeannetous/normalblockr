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

test_that("NB_fixed_blocks_diagonal: check dimensions, optimization and field access", {
  model <- normal_block(Y, X, C, noise_cov = "diagonal")
  params <- model$model_par
  expect_equal(model$n, nrow(Y))
  expect_equal(model$p, ncol(Y))
  expect_equal(model$d, ncol(X))
  expect_lt(model$BIC, 36100)
  expect_gt(model$loglik, -17700)
  model_sparse <- normal_block(Y, X, C, sparsity = 0.1, noise_cov = "diagonal")
})

test_that("NB_fixed_blocks_spherical: check dimensions, optimization and field access", {
  model <- normal_block(Y, X, C, noise_cov = "spherical")
  params <- model$model_par
  expect_equal(model$n, nrow(Y))
  expect_equal(model$p, ncol(Y))
  expect_equal(model$d, ncol(X))
  expect_lt(model$BIC, 37400)
  expect_gt(model$loglik, -18300)
  model_sparse <- normal_block(Y, X, C, sparsity = 0.1, noise_cov = "spherical")
})
