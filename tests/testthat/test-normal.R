###############################################################################
###############################################################################
## Use pre-save testdata (seed are hard to handle in testhat)
testdata <- readRDS("testdata/testdata_normal.RDS")
Y <- testdata$Y
X <- testdata$X


###############################################################################
###############################################################################

test_that("normal: check dimensions and field access", {
  data <- normal_data$new(Y, X)
  model <- normal$new(data)
  expect_equal(model$n, nrow(Y))
  expect_equal(model$p, ncol(Y))
  expect_equal(model$d, ncol(X))
})
