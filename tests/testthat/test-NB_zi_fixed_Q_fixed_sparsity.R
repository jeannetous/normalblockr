###############################################################################
###############################################################################
testdata <- readRDS("testdata/testdata_normal_zi.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C ; Q <- ncol(C)

###############################################################################
###############################################################################
test_that("NB_zi_fixed_Q: check dimensions, optimization and field access", {
  data  <- normal_data$new(Y, X)
  model <- NB_zi_fixed_Q_fixed_sparsity_diagonal$new(data, Q, penalty = 0.2)
  model$optimize()
})
