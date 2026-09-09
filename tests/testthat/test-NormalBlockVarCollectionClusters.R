###############################################################################
## NormalBlockVarCollectionClusters (R/NormalBlockVarCollectionClusters.R,
## R/NormalBlockCollectionClusters.R): fitting a
## collection over several q at once, including a per-q clustering_init list
## and the heuristic path. The second half of the file exercises refine()
## (candidates_split()/candidates_merge()) on the collection.
## Use pre-save testdata (seed are hard to handle in testhat)
testdata <- readRDS("testdata/testdata_normal_zi.RDS")
Y <- testdata$Y ; X <- testdata$X
C <- testdata$parameters$C ; q <- ncol(C)

test_that("normal block with diagonal residual covariance and unknown q", {
  clustering_init_3 <- normalblockr:::get_clusters(C)
  clustering_init_2 <- clustering_init_3 ; clustering_init_2[which(clustering_init_2 == 3) ] = 2
  clustering_init_1 <- clustering_init_2 ; clustering_init_1[which(clustering_init_1 == 2) ] = 1
  clustering_init   <- list(clustering_init_1, clustering_init_2, clustering_init_3)

  data  <- normalblockr:::NormalBlockData$new(Y, X)
  model <- normalblockr:::NormalBlockVarCollectionClusters$new(data, c(1, 2, 3), zero_inflation = FALSE,
                            control = NB_control(clustering_init = clustering_init))
  model$optimize()
  model_BIC <- model$get_best_model("BIC")
  expect_lt(model_BIC$BIC, 7455)
})

test_that("normal block with unknown q, heuristic", {
  data  <- NormalBlockData$new(Y, X)
  model <- normalblockr:::NormalBlockVarCollectionClusters$new(data, c(2, 3, 4), sparsity = 0.05,
                            control = NB_control(heuristic = TRUE))
  model$optimize()
  model_3 <- model$get_model(3)
  expect_lt(Metrics::rmse(model_3$fitted, Y), 2.9)
})

###############################################################################
## refine(): tries a short split from each model's smaller-q neighbor and/or a
## short merge from its larger-q neighbor, keeping a candidate only if it's
## strictly better -- see its roxygen for the empirical rationale
## (inst/clustering_initialization_benchmark and the brca_rppa/university
## comparison that motivated it, and why "merge" was added alongside "split").
set.seed(7)
ex2 <- generate_normal_block_var_data(n = 60, p = 30, d = 1, q = 5)
data2 <- NormalBlockData$new(ex2$Y, ex2$X)

test_that("NB_control()'s refine option defaults to FALSE", {
  expect_false(NB_control()$refine)
})

test_that("refine() never makes a model's deviance worse, and skips non-contiguous q gaps", {
  coll <- NormalBlockVarCollectionClusters$new(data2, c(2, 3, 5), control = NB_control(verbose = FALSE))
  coll$optimize(control = list(niter = 100, threshold = 1e-4, verbose = FALSE))
  deviance_before <- sapply(coll$models, function(m) m$deviance)

  expect_no_error(coll$refine())
  deviance_after <- sapply(coll$models, function(m) m$deviance)
  expect_true(all(deviance_after <= deviance_before + 1e-8))
})

test_that("refine(directions = 'split') only tries split candidates, never merge", {
  coll <- NormalBlockVarCollectionClusters$new(data2, 2:5, control = NB_control(verbose = FALSE))
  coll$optimize(control = list(niter = 100, threshold = 1e-4, verbose = FALSE))
  expect_no_error(coll$refine(directions = "split"))
})

test_that("refine(directions = 'merge') only tries merge candidates, never split", {
  coll <- NormalBlockVarCollectionClusters$new(data2, 2:5, control = NB_control(verbose = FALSE))
  coll$optimize(control = list(niter = 100, threshold = 1e-4, verbose = FALSE))
  expect_no_error(coll$refine(directions = "merge"))
})

test_that("refine() with both directions is never worse than either alone", {
  coll_split <- NormalBlockVarCollectionClusters$new(data2, 2:5, control = NB_control(verbose = FALSE))
  coll_split$optimize(control = list(niter = 100, threshold = 1e-4, verbose = FALSE))
  coll_merge <- coll_split$clone(deep = TRUE)
  coll_both  <- coll_split$clone(deep = TRUE)

  coll_split$refine(directions = "split")
  coll_merge$refine(directions = "merge")
  coll_both$refine(directions = c("split", "merge"))

  dev_split <- sapply(coll_split$models, function(m) m$deviance)
  dev_merge <- sapply(coll_merge$models, function(m) m$deviance)
  dev_both  <- sapply(coll_both$models, function(m) m$deviance)

  expect_true(all(dev_both <= dev_split + 1e-8))
  expect_true(all(dev_both <= dev_merge + 1e-8))
})

test_that("NB_control(refine = TRUE) makes optimize() call refine() automatically", {
  coll_auto <- NormalBlockVarCollectionClusters$new(data2, 2:5, control = NB_control(verbose = FALSE, refine = TRUE))
  coll_auto$optimize(control = list(niter = 100, threshold = 1e-4, verbose = FALSE))

  coll_manual <- NormalBlockVarCollectionClusters$new(data2, 2:5, control = NB_control(verbose = FALSE, refine = FALSE))
  coll_manual$optimize(control = list(niter = 100, threshold = 1e-4, verbose = FALSE))
  coll_manual$refine()

  expect_equal(sapply(coll_auto$models, function(m) m$deviance),
               sapply(coll_manual$models, function(m) m$deviance))
})

