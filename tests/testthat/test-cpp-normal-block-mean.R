###############################################################################
## These tests check that the Rcpp/RcppArmadillo (V)EM core
## (NormalBlockMeanKnownClusters_fit / NormalBlockMeanUnknownClusters_fit,
## src/exports.cpp) reproduces *exactly* (up to numerical precision) the
## EM/VEM recursion implemented in R (private$EM_optimize_R, kept in
## R/NormalBlockMean*Clusters.R while the port is validated), starting from
## the same initial parameters, both unpenalized (sparsity = 0, plain
## inversion) and penalized (sparsity > 0, in-package graphical lasso,
## called back from C++, see src/omega_estimation.h).
##
## The comparison runs with accelerate = FALSE: the SQUAREM extrapolation
## the C++ core applies by default deliberately changes the trajectory (it is
## checked separately below), so only the plain recursion can be compared
## trace-for-trace.
##
## Both recursions are driven through a fresh model object each time so that
## private$optim_initialize() returns the same starting point: the clustering
## is passed explicitly (rather than re-derived by a randomized heuristic) for
## the unknown-clusters case.
## Use pre-save testdata (seed are hard to handle in testhat)
testdata <- readRDS("testdata/testdata_normal_mean_block.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C ; q <- ncol(C)

niter     <- 6
threshold <- -1 # never trigger early stopping: forces exactly `niter` iterations

test_that("NormalBlockMeanKnownClusters_fit matches the R recursion (unpenalized/sparse)", {
  data <- NormalBlockData$new(Y, X)

  for (sparsity in c(0, 0.05)) {
    nc <- "full"
    model <- NormalBlockMeanKnownClusters$new(data, C, sparsity = sparsity,
                                              control = NB_control(verbose = FALSE,
                                                                   noise_covariance = nc))
    init <- model$.__enclos_env__$private$optim_initialize()
    ref  <- model$.__enclos_env__$private$EM_optimize_R(
      list(niter = niter, threshold = threshold))

    res <- NormalBlockMeanKnownClusters_fit(
      Y = data$Y, X = data$X, C = C, B0 = init$B, Omega0 = init$Omega,
      sparsity = sparsity, sparsity_weights = model$sparsity_weights,
      noise_covariance = nc, niter = niter, threshold = threshold, accelerate = FALSE)

    expect_equal(res$B,     ref$B,     tolerance = 1e-8)
    expect_equal(res$Omega, ref$Omega, tolerance = 1e-8)
    expect_equal(res$objective, ref$ll_list, tolerance = 1e-8)
  }
})

test_that("NormalBlockMeanUnknownClusters_fit matches the R recursion (unpenalized/sparse)", {
  data <- NormalBlockData$new(Y, X)
  fixed_point_niter <- 5

  for (sparsity in c(0, 0.05)) {
    nc <- "full"
    ctrl <- NB_control(verbose = FALSE, clustering_init = get_clusters(C),
                       noise_covariance = nc)
    model_R   <- NormalBlockMeanUnknownClusters$new(data, q, sparsity = sparsity, control = ctrl)
    model_cpp <- NormalBlockMeanUnknownClusters$new(data, q, sparsity = sparsity, control = ctrl)

    init <- model_cpp$.__enclos_env__$private$optim_initialize()
    ref  <- model_R$.__enclos_env__$private$EM_optimize_R(
      list(niter = niter, threshold = threshold, fixed_point_niter = fixed_point_niter))

    res <- NormalBlockMeanUnknownClusters_fit(
      Y = data$Y, X = data$X, B0 = init$B, Omega0 = init$Omega, tau0 = init$tau,
      sparsity = sparsity, sparsity_weights = model_R$sparsity_weights,
      noise_covariance = nc, fixed_point_niter = fixed_point_niter,
      niter = niter, threshold = threshold, accelerate = FALSE)

    expect_equal(res$B,      ref$B,      tolerance = 1e-8)
    expect_equal(res$Omega,  ref$Omega,  tolerance = 1e-8)
    expect_equal(res$C,      ref$C,      tolerance = 1e-8)
    expect_equal(res$alpha,  ref$alpha,  tolerance = 1e-8)
    expect_equal(res$Psi,    ref$Psi,    tolerance = 1e-8)
    expect_equal(res$Phi,    ref$Phi,    tolerance = 1e-8)
    expect_equal(res$Lambda, ref$Lambda, tolerance = 1e-8)
    expect_equal(res$objective, ref$ll_list, tolerance = 1e-8)
  }
})

test_that("the ELBO trace of the C++ core is non-decreasing", {
  data  <- NormalBlockData$new(Y, X)
  model <- NormalBlockMeanUnknownClusters$new(data, q, control = NB_control(verbose = FALSE))
  model$optimize(control = NB_control(verbose = FALSE))
  expect_true(all(diff(model$objective) >= -1e-6))
})

test_that("the SQUAREM extrapolation reaches at least the plain recursion's ELBO", {
  data <- NormalBlockData$new(Y, X)
  ctrl <- NB_control(verbose = FALSE, clustering_init = get_clusters(C),
                     noise_covariance = "full")
  init <- NormalBlockMeanUnknownClusters$new(data, q, control = ctrl
            )$.__enclos_env__$private$optim_initialize()

  fit <- function(accelerate) {
    NormalBlockMeanUnknownClusters_fit(
      Y = data$Y, X = data$X, B0 = init$B, Omega0 = init$Omega, tau0 = init$tau,
      sparsity = 0, sparsity_weights = matrix(1, ncol(Y), ncol(Y)) - diag(ncol(Y)),
      noise_covariance = "full", fixed_point_niter = 5, niter = 50, threshold = 1e-6,
      accelerate = accelerate)
  }
  plain <- fit(FALSE)
  accel <- fit(TRUE)

  expect_gte(tail(accel$objective, 1), tail(plain$objective, 1) - 1e-6)
  expect_true(all(diff(accel$objective) >= -1e-6)) # still monotone
})

test_that("the diagonal/spherical Sigma variants also match the R recursion", {
  data <- NormalBlockData$new(Y, X)
  fixed_point_niter <- 5

  for (nc in c("diagonal", "spherical")) {
    ctrl <- NB_control(verbose = FALSE, clustering_init = get_clusters(C), noise_covariance = nc)
    model_R   <- NormalBlockMeanUnknownClusters$new(data, q, control = ctrl)
    model_cpp <- NormalBlockMeanUnknownClusters$new(data, q, control = ctrl)

    init <- model_cpp$.__enclos_env__$private$optim_initialize()
    ref  <- model_R$.__enclos_env__$private$EM_optimize_R(
      list(niter = niter, threshold = threshold, fixed_point_niter = fixed_point_niter))

    res <- NormalBlockMeanUnknownClusters_fit(
      Y = data$Y, X = data$X, B0 = init$B, Omega0 = init$Omega, tau0 = init$tau,
      sparsity = 0, sparsity_weights = model_R$sparsity_weights, noise_covariance = nc,
      fixed_point_niter = fixed_point_niter, niter = niter, threshold = threshold,
      accelerate = FALSE)

    expect_equal(res$B,     ref$B,     tolerance = 1e-8)
    expect_equal(res$Omega, ref$Omega, tolerance = 1e-8)
    expect_equal(res$C,     ref$C,     tolerance = 1e-8)
    expect_equal(res$objective, ref$ll_list, tolerance = 1e-8)
    ## the whole point of these variants: Omega stays diagonal
    expect_equal(res$Omega, diag(diag(res$Omega)))
  }
})

test_that("a diagonal Sigma lifts the n > p requirement", {
  narrow <- NormalBlockData$new(Y[1:(ncol(Y) - 10), , drop = FALSE],
                                X[1:(ncol(Y) - 10), , drop = FALSE])
  ## only a full Sigma has to be inverted, and it is no longer the default
  expect_error(normal_block(narrow, blocks = q, model = "mean",
                            control = NB_control(verbose = FALSE, noise_covariance = "full")),
               "need n > p", fixed = TRUE)
  expect_no_error(normal_block(narrow, blocks = q, model = "mean",
                               control = NB_control(verbose = FALSE)))
  for (nc in c("diagonal", "spherical"))
    expect_no_error(normal_block(narrow, blocks = q, model = "mean",
                                 control = NB_control(verbose = FALSE, noise_covariance = nc)))
})
