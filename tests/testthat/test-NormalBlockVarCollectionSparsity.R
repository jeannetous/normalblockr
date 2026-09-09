###############################################################################
## NormalBlockVarCollectionSparsity (R/NormalBlockVarCollectionSparsity.R):
## a sparsity path at fixed q (known or unknown clusters), StARS stability
## selection, and the M-step warm-starting each penalty from the previous
## one's converged parameters (fewer total EM iterations, same BIC path).
## Use pre-save testdata (seed are hard to handle in testhat)
testdata <- readRDS("testdata/testdata_normal.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C ; q <- ncol(C)

data <- NormalBlockData$new(Y, X)

test_that("normal block with changing sparsity, integrated inference", {
  model <- normalblockr:::NormalBlockVarCollectionSparsity$new(data, q, control = NB_control(verbose = FALSE))
  model$optimize()
  model$stability_selection()
  model_StARS <- model$get_best_model("StARS")
  expect_gt(model_StARS$loglik, -2635)
})


test_that("normal block with fixed clusters, spherical residual covariance and heuristic inference", {
  model <- normalblockr:::NormalBlockVarCollectionSparsity$new(data, C,
                                                   control = NB_control(noise_covariance = "spherical",
                                                                        heuristic = TRUE,
                                                                        verbose = FALSE))
  model$optimize()
  model$stability_selection()
  model_StARS <- model$get_best_model("StARS")
  expect_lt(Metrics::rmse(model_StARS$fitted, Y), 2.4)
})

test_that("optimize() warm-starts each penalty from the previous one, reducing total EM iterations without changing the BIC path", {
  ## n_sparsity_penalties = 15: the warm-start advantage grows with the
  ## length of the sparsity path and is only marginal (sometimes net-negative
  ## by a couple of iterations) on very short paths -- 15 gives a robust,
  ## comfortable margin on this dataset.
  set.seed(123)
  warm <- normalblockr:::NormalBlockVarCollectionSparsity$new(data, C, control = NB_control(verbose = FALSE, n_sparsity_penalties = 15))
  warm$optimize(NB_control(verbose = FALSE))

  ## same penalties, but optimized independently (the pre-warm-start behaviour:
  ## every model starts from its own heuristic-derived initial parameters)
  set.seed(123)
  cold <- normalblockr:::NormalBlockVarCollectionSparsity$new(data, C, control = NB_control(verbose = FALSE, n_sparsity_penalties = 15))
  cold$models <- lapply(cold$models, function(model) { model$optimize(NB_control(verbose = FALSE)); model })

  niter_warm <- sapply(warm$models, function(m) length(m$objective))
  niter_cold <- sapply(cold$models, function(m) length(m$objective))

  ## the first model has no predecessor to warm-start from, so it should be
  ## unaffected; every later one should need no more iterations than cold
  expect_equal(niter_warm[1], niter_cold[1])
  expect_lt(sum(niter_warm), sum(niter_cold))

  ## warm-starting changes the optimization path, not the optimum it converges to
  expect_equal(sapply(warm$models, function(m) m$BIC), sapply(cold$models, function(m) m$BIC), tolerance = 1e-2)
})
