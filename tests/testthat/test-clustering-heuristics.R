###############################################################################
## The clustering-heuristic registry (private$clustering_methods, in
## R/NormalBlockBase.R): the four user-selectable algorithms
## (kmeans/ward2/sbm/spectral, via NB_control(clustering_init = )), the ward2
## fallback used whenever the chosen one collapses to fewer than q clusters,
## and how comparable the optima they lead to are.
##
## How a collection reuses the q-independent part of these is a separate
## concern, in test-shared-initialization.R.
###############################################################################
set.seed(123)
ex   <- generate_normal_block_var_data(n = 30, p = 10, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)

test_that("every user-selectable clustering heuristic name produces a valid q-cluster init", {
  for (approx in c("kmeans", "ward2", "sbm", "spectral")) {
    model <- NormalBlockVarUnknownClusters$new(data, 3, control = NB_control(verbose = FALSE, clustering_init = approx))
    expect_no_error(model$optimize(control = list(niter = 2, threshold = -1)))
    expect_equal(model$q, 3)
  }
})

test_that("heuristic_clustering() falls back to ward2 when the chosen method collapses clusters", {
  model <- NormalBlockVarUnknownClusters$new(data, 3, control = NB_control(verbose = FALSE, clustering_init = "kmeans"))
  priv  <- model$.__enclos_env__$private

  ## Force the chosen heuristic to always return a single cluster, regardless of R
  priv$clustering_methods$kmeans <- function(R, q) rep(1, ncol(R))

  R <- matrix(rnorm(data$n * data$p), data$n, data$p)
  C <- suppressWarnings(priv$heuristic_clustering(R))

  expect_equal(dim(C), c(data$p, 3))
  expect_equal(length(unique(get_clusters(C))), 3)
})

test_that("ward2 does not error on a (near-)constant column (cor() = NA, e.g. a rare ZI variable)", {
  R <- matrix(rnorm(data$n * data$p), data$n, data$p)
  R[, 1] <- 0; R[1, 1] <- 5 # a single non-zero residual for variable 1
  priv <- NormalBlockVarUnknownClusters$new(data, 3, control = NB_control(verbose = FALSE))$.__enclos_env__$private
  expect_no_warning(cl <- priv$clustering_methods$ward2(R, 3))
  expect_length(cl, data$p)
  expect_equal(length(unique(cl)), 3)
})

test_that("the heuristics reach comparable optima on a variance-block model", {
  td   <- readRDS("testdata/testdata_normal.RDS")
  q    <- ncol(td$parameters$C)
  data <- NormalBlockData$new(td$Y, td$X)

  model_kmeans <- NormalBlockVarUnknownClusters$new(data, q, control = NB_control(clustering_init = "kmeans"))
  model_kmeans$optimize()

  model_ward2 <- NormalBlockVarUnknownClusters$new(data, q, control = NB_control(clustering_init = "ward2"))
  model_ward2$optimize()

  model_sbm <- NormalBlockVarUnknownClusters$new(data, q, control = NB_control(clustering_init = "sbm"))
  model_sbm$optimize()

  expect_gt(model_kmeans$loglik, -2650)
  expect_gt(model_ward2$loglik, -2650)
  expect_gt(model_sbm$loglik, -2650)

  # expect_gte(model_kmeans$loglik, model_sbm$loglik)
  # expect_gte(model_ward2$loglik, model_sbm$loglik)

  expect_gt(aricode::ARI(model_kmeans$clustering, model_sbm$clustering), 0.75)
  expect_gt(aricode::ARI(model_kmeans$clustering, model_ward2$clustering), 0.75)
})


testdata <- readRDS("testdata/testdata_normal_zi.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C ; q <- ncol(C)
data <- NormalBlockData$new(Y, X)

test_that("the heuristics reach comparable optima on a zero-inflated model", {
  td   <- readRDS("testdata/testdata_normal_zi.RDS")
  q    <- ncol(td$parameters$C)
  data <- NormalBlockData$new(td$Y, td$X)

  model_ward2 <- ZINormalBlockVarUnknownClusters$new(data, q, control = NB_control(clustering_init = "ward2"))
  model_ward2$optimize()

  model_kmeans <- ZINormalBlockVarUnknownClusters$new(data, q, control = NB_control(clustering_init = "kmeans"))
  model_kmeans$optimize()

  model_sbm <- ZINormalBlockVarUnknownClusters$new(data, q, control = NB_control(clustering_init = "sbm"))
  model_sbm$optimize()

  expect_gt(model_kmeans$loglik, -2750)
  expect_gt(model_ward2$loglik, -2750)
  expect_gt(model_sbm$loglik, -2750)

  expect_gte(model_kmeans$loglik, model_sbm$loglik)
  expect_gte(model_ward2$loglik, model_sbm$loglik)

  expect_gte(aricode::ARI(model_kmeans$clustering, model_sbm$clustering)  , 0.49)
  expect_gte(aricode::ARI(model_kmeans$clustering, model_ward2$clustering), 0.99)
})

test_that("spectral() caps its eigenvector count at cov(R)'s numerical rank", {
  ## R = X %*% B (the mean-block family's own clustering input) is rank <= d
  ## by construction; eigen() still returns q vectors when asked for more
  ## than that, completing the null space arbitrarily -- not from the data at
  ## all. Using them pollutes the embedding kmeans clusters on. This is the
  ## regime the mean-block family is routinely in (small d, q well above it).
  set.seed(1)
  ex   <- generate_normal_block_mean_data(n = 100, p = 30, d = 2, q = 5)
  data <- NormalBlockData$new(ex$Y, ex$X)
  R    <- data$X %*% data$ols_fit()$B

  expect_equal(sum(eigen(cov(R), symmetric = TRUE, only.values = TRUE)$values > 1e-8), 2)

  priv <- NormalBlockMeanUnknownClusters$new(data, 5, control = NB_control(verbose = FALSE))$.__enclos_env__$private
  cl   <- priv$clustering_methods$spectral(R, 5)
  expect_length(unique(cl), 5) # still produces the requested number of clusters

  ## a full-rank input (the variance-block family's residuals) is untouched:
  ## the cap never binds there, so behaviour for that family doesn't change
  exv   <- generate_normal_block_var_data(n = 100, p = 30, d = 1, q = 4)
  datav <- NormalBlockData$new(exv$Y, exv$X)
  Rv    <- ols_residuals(datav)
  expect_equal(sum(eigen(cov(Rv), symmetric = TRUE, only.values = TRUE)$values > 1e-8 * max(eigen(cov(Rv), symmetric = TRUE, only.values = TRUE)$values)),
               ncol(Rv))
})
