###############################################################################
###############################################################################

testdata <- readRDS("testdata/testdata_ZInormal.RDS")
Y <- testdata$Y
X <- testdata$X

###############################################################################
###############################################################################

test_that("normal_zi: check dimensions, optimization and field access", {
  model <- normal_zi$new(Y, X)
  model$optimize()
  expect_equal(model$n, nrow(Y))
  expect_equal(model$p, ncol(Y))
  expect_equal(model$d, ncol(X))
  expect_lt(model$BIC, 107409)
  expect_gt(model$loglik,  -52564)
})
