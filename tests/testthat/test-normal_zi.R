###############################################################################
###############################################################################

testdata <- readRDS("testdata/testdata_ZInormal.RDS")
Y <- testdata$Y
X <- testdata$X
Q <- ncol(C)

###############################################################################
###############################################################################

test_that("normal_zi: check dimensions, optimization and field access", {
  model <- normal_zi$new(Y, X)
  model$optimize(niter = 60)
  params <- model$get_model_parameters()
  expect_equal(params$n, nrow(Y))
  expect_equal(params$p, ncol(Y))
  expect_equal(params$d, ncol(X))
  expect_lt(model$BIC, 82510)
  expect_gt(model$loglik,  -40400)
})
