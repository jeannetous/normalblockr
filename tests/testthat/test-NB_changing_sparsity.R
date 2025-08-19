###############################################################################
###############################################################################
## Use pre-save testdata (seed are hard to handle in testhat)
testdata <- readRDS("testdata/testdata_normal.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C ; Q <- ncol(C)


###############################################################################
###############################################################################
data <- NB_data$new(Y, X)

test_that("normal block with changing sparsity, integrated inference", {
  model <- normalblockr:::NB_changing_sparsity$new(data, Q, control = NB_control(verbose = FALSE))
  model$optimize()
  model$stability_selection()
  model_StARS <- model$get_best_model("StARS")
  expect_gt(model_StARS$loglik, -2635)
})


test_that("normal block with fixed clusters, spherical residual covariance and heuristic inference", {
  model <- normalblockr:::NB_changing_sparsity$new(data, C,
                                                   control = NB_control(noise_covariance = "spherical",
                                                                        heuristic = TRUE,
                                                                        verbose = FALSE))
  model$optimize()
  model$stability_selection()
  model_StARS <- model$get_best_model("StARS")
  expect_lt(Metrics::rmse(model_StARS$fitted, Y), 2.4)
})
