###############################################################################
###############################################################################
testdata <- readRDS("testdata/testdata_normal_zi.RDS")
Y <- testdata$Y
X <- testdata$X

###############################################################################
###############################################################################
test_that("normal_diag_zi: check dimensions, optimization and field access", {
  data <- normalblockr:::normal_data$new(Y, X)
  model <- normalblockr:::normal_diag_zi$new(data)
  model$optimize()
  expect_equal(model$n, nrow(Y))
  expect_equal(model$p, ncol(Y))
  expect_equal(model$d, ncol(X))
  expect_lt(model$BIC, 6700)
  expect_gt(model$loglik,  -3200)
})
