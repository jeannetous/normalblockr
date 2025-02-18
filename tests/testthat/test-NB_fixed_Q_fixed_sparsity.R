###############################################################################
###############################################################################
## Use pre-save testdata (seed are hard to handle in testhat)
testdata <- readRDS("testdata/testdata_normal.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C


###############################################################################
###############################################################################

test_that("normal: check dimensions and field access", {
  data <- normal_data$new(Y, X)
  model <- NB_fixed_blocks_fixed_sparsity_diagonal$new(data, C)
  model$optimize()
  expect_gt(model$loglik, -2600)
})
