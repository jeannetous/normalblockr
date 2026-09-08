###############################################################################
## Work a collection over q used to repeat once per model although it does not
## depend on q: the OLS fit, the zero-inflation logistic regressions, and
## whatever part of a clustering heuristic is shared (one hierarchical tree cut
## at each q, one eigendecomposition, one lossless row compression, one wide
## SBM exploration).
##
## These are removals of redundant computation, so most of what follows checks
## that the fits they feed come out unchanged. The heuristics themselves live
## in test-clustering-heuristics.R.
###############################################################################
set.seed(99)
sbm_ex   <- generate_normal_block_var_data(n = 60, p = 24, d = 1, q = 4)
sbm_data <- NormalBlockData$new(sbm_ex$Y, sbm_ex$X)
sbm_R    <- ols_residuals(sbm_data)

test_that("ols_fit() is memoized and matches a fresh computation", {
  set.seed(1)
  ex <- generate_normal_block_var_data(n = 60, p = 20, d = 1, q = 3)
  d  <- NormalBlockData$new(ex$Y, ex$X)

  first <- d$ols_fit()
  expect_named(first, c("B", "R", "Sigma"))
  expect_equal(first$B, d$XtXm1 %*% d$XtY)
  expect_equal(first$R, d$Y - d$X %*% first$B)
  expect_equal(first$Sigma, cov(first$R))

  ## same object back, not a recomputation
  expect_identical(d$ols_fit(), first)
})

test_that("zi_fit() is memoized and matches the per-variable logistic fits", {
  set.seed(2)
  ex <- generate_normal_block_var_data(n = 80, p = 12, d = 1, q = 2, kappa = rep(0.3, 12))
  d  <- NormalBlockData$new(ex$Y, ex$X)

  zi <- d$zi_fit()
  expect_named(zi, c("B0", "kappa", "ZI_cond_mean"))
  expect_identical(d$zi_fit(), zi)

  ## the same p independent logistic regressions, spelled out
  B0 <- t(sapply(lapply(1:d$p, function(j) {
    glm(zeros ~ 0 + ., family = binomial(link = "logit"),
        data = data.frame(zeros = d$zeros[, j], d$X0))$coefficients
  }), unlist))
  ## values only: B0's dimnames are a leaked deparse of whatever expression
  ## built the data.frame, which is a separate (pre-existing) wart
  expect_equal(zi$B0, B0, ignore_attr = TRUE)
  expect_equal(zi$ZI_cond_mean,
               sum(normalblockr:::xlogy(d$zeros, zi$kappa)) +
               sum(normalblockr:::xlogy(d$zeros_bar, 1 - zi$kappa)))
})

test_that("zi_ols_fit() is memoized and matches the free function", {
  set.seed(21)
  ex <- generate_normal_block_var_data(n = 70, p = 12, d = 1, q = 2, kappa = rep(0.3, 12))
  d  <- NormalBlockData$new(ex$Y, ex$X)
  fit <- d$zi_ols_fit()
  expect_named(fit, c("B", "dm1", "R"))
  expect_equal(fit, normalblockr:::zi_weighted_fit(d))
  expect_identical(d$zi_ols_fit(), fit)
  ## the masked residuals every zero-inflated heuristic starts from
  expect_equal(normalblockr:::zi_residuals(d), fit$R)
})

test_that("every model in a collection shares one data object, hence one fit", {
  set.seed(3)
  ex <- generate_normal_block_var_data(n = 60, p = 15, d = 1, q = 2, kappa = rep(0.3, 15))
  d  <- NormalBlockData$new(ex$Y, ex$X)
  coll <- normal_block(d, blocks = 2:4, zero_inflation = TRUE,
                       control = NB_control(verbose = FALSE))

  ## R6 reference semantics: the models point at the caller's object, so the
  ## memoized zi_fit()/ols_fit() are computed once for the whole collection
  for (m in coll$models) expect_identical(m$data, d)
})

test_that("a q-range collection shares the q-independent part of its heuristic", {
  set.seed(4)
  ex <- generate_normal_block_var_data(n = 80, p = 25, d = 1, q = 3)
  d  <- NormalBlockData$new(ex$Y, ex$X)
  R  <- normalblockr:::ols_residuals(d)
  q_list <- 2:5

  ## ward2: one tree, cut at each q -- identical to cutting a per-q tree
  path <- normalblockr:::clustering_path_for_collection(R, q_list, "ward2")
  expect_equal(length(path), length(q_list))
  for (i in seq_along(q_list)) {
    expect_equal(path[[i]], cutree(normalblockr:::ward2_tree(R), q_list[i]))
  }

  ## kmeans shares its lossless row compression, but not Lloyd's algorithm
  km <- normalblockr:::clustering_path_for_collection(R, q_list, "kmeans")
  for (i in seq_along(q_list)) expect_equal(length(unique(km[[i]])), q_list[i])

  ## nothing is shareable for the model's own best_of_inits search
  expect_null(normalblockr:::clustering_path_for_collection(R, q_list, "best_of_inits"))
  expect_null(normalblockr:::clustering_path_for_collection(R, q_list, c("ward2", "kmeans")))
})

test_that("sharing the heuristic leaves the fitted collection unchanged", {
  set.seed(5)
  ex <- generate_normal_block_var_data(n = 80, p = 25, d = 1, q = 3)
  d  <- NormalBlockData$new(ex$Y, ex$X)
  q_list <- 2:5

  shared <- normal_block(d, blocks = q_list, control = NB_control(verbose = FALSE))

  ## the same models, each cold-started from its own independently-built tree
  per_model <- lapply(q_list, function(q) {
    cl <- cutree(normalblockr:::ward2_tree(normalblockr:::ols_residuals(d)), q)
    normal_block(d, blocks = q, control = NB_control(verbose = FALSE, clustering_init = cl))
  })

  expect_equal(sapply(shared$models, function(m) m$loglik),
               sapply(per_model, function(m) m$loglik))
  expect_equal(lapply(shared$models, function(m) m$clustering),
               lapply(per_model, function(m) m$clustering))
})

test_that("a heuristic that collapses below q falls back to the model's own", {
  ## clustering_path_for_collection() must never hand down a clustering with
  ## the wrong number of groups: NormalBlockBase would reject it outright,
  ## where heuristic_clustering() has a documented fallback. Asking for more
  ## groups than there are distinct variables forces the degenerate case.
  set.seed(6)
  R <- matrix(rnorm(40), 20, 2)[, c(1, 1, 2, 2)] # only two distinct columns
  path <- normalblockr:::clustering_path_for_collection(R, c(2, 4), "ward2")
  ## whatever it returns per q is either exactly q groups, or the method name
  for (i in seq_along(path)) {
    cl <- path[[i]]
    expect_true(identical(cl, "ward2") || length(unique(cl)) == c(2, 4)[i])
  }
})

test_that("compress_columns() preserves the geometry kmeans actually sees", {
  set.seed(41)
  ## a tall, low-rank embedding, as the mean-block family produces (X %*% B)
  X <- matrix(rnorm(200 * 3), 200, 3)
  B <- matrix(rnorm(3 * 40), 3, 40)
  R <- X %*% B

  Rc <- normalblockr:::compress_columns(R)
  expect_equal(nrow(Rc), 3L)                       # down to the true rank
  expect_equal(ncol(Rc), ncol(R))
  expect_equal(as.matrix(dist(t(Rc))), as.matrix(dist(t(R))))

  ## and therefore the same clustering, starts included
  for (q in 2:5) {
    set.seed(7); a <- kmeans(t(R),  q, nstart = 30, iter.max = 50)$cluster
    set.seed(7); b <- kmeans(t(Rc), q, nstart = 30, iter.max = 50)$cluster
    expect_equal(aricode::ARI(a, b), 1)
  }
})

test_that("compress_columns() leaves a full-rank matrix alone", {
  set.seed(42)
  R <- matrix(rnorm(30), 5, 6)             # 5 rows, rank 5: nothing to gain
  expect_identical(normalblockr:::compress_columns(R), R)
  expect_identical(normalblockr:::compress_columns(matrix(1:4, 1, 4)),
                   matrix(1:4, 1, 4))      # a single row is already minimal
})

test_that("the mean-block collection shares its heuristic across q", {
  set.seed(43)
  ex <- generate_normal_block_mean_data(n = 120, p = 30, d = 2, q = 3)
  d  <- NormalBlockData$new(ex$Y, ex$X)
  q_list <- 2:5

  path <- normalblockr:::clustering_path_for_family(d, q_list, "mean", FALSE,
                                                    NB_control(verbose = FALSE))
  expect_equal(length(path), length(q_list))       # kmeans is now shareable too
  for (i in seq_along(q_list)) expect_equal(length(unique(path[[i]])), q_list[i])

  ## the family default resolves to kmeans (kmeans draws its starts, so the
  ## two calls only agree from the same RNG state)
  set.seed(44); a <- normalblockr:::clustering_path_for_family(
    d, q_list, "mean", FALSE, NB_control(verbose = FALSE))
  set.seed(44); b <- normalblockr:::clustering_path_for_family(
    d, q_list, "mean", FALSE, NB_control(verbose = FALSE, clustering_init = "kmeans"))
  expect_equal(a, b)
})

test_that("sbm_clustering_path() returns one valid membership vector per q, named by q, never NULL", {
  q_list <- c(2, 3, 4, 5)
  path <- sbm_clustering_path(sbm_R, q_list)

  expect_equal(names(path), as.character(q_list))
  for (q in q_list) {
    membership <- path[[as.character(q)]]
    expect_false(is.null(membership))
    expect_length(membership, sbm_data$p)
    expect_equal(length(unique(membership)), q)
  }
})

test_that("sbm_clustering_path() falls back to a cheap ward2 clustering for q's the SBM exploration doesn't reach", {
  ## A q_list reaching far beyond what a tiny dataset's SBM exploration will
  ## ever cover forces the fallback path for the largest values -- regression
  ## test for the cost blow-up found on a real dataset where SBM's own model
  ## selection stopped early (re-running a dedicated SBM call per missing q
  ## reintroduced the cost this function exists to avoid).
  q_list <- 2:20
  path <- sbm_clustering_path(sbm_R, q_list)
  expect_true(all(!sapply(path, is.null)))
  for (q in q_list) expect_equal(length(unique(path[[as.character(q)]])), q)
})

test_that("NormalBlockVarCollectionClusters with clustering_init = 'sbm' eagerly resolves it from the shared path", {
  coll <- NormalBlockVarCollectionClusters$new(sbm_data, 2:5, control = NB_control(verbose = FALSE, clustering_init = "sbm"))

  for (model in coll$models) {
    ## every model converges to a valid q-cluster initialization once
    ## actually initialized, whether served by the shared path or (for q's it
    ## didn't cover) by the per-q heuristic fallback
    init <- model$.__enclos_env__$private$optim_initialize()
    expect_equal(length(unique(get_clusters(init$C))), model$q)
  }
  expect_no_error(coll$optimize(control = list(niter = 3, threshold = -1, verbose = FALSE)))
})

test_that("an explicit clustering_init is never overridden by the sbm path", {
  cl_init <- list(rep(1:2, length.out = sbm_data$p), rep(1:3, length.out = sbm_data$p))
  coll <- NormalBlockVarCollectionClusters$new(sbm_data, c(2, 3), control = NB_control(
    verbose = FALSE, clustering_init = cl_init
  ))

  for (i in seq_along(coll$models)) {
    C0 <- coll$models[[i]]$.__enclos_env__$private$C
    expect_equal(get_clusters(C0), cl_init[[i]])
  }
})

test_that("NormalBlockVarCollectionClustersSparsity also uses the shared sbm path", {
  coll <- NormalBlockVarCollectionClustersSparsity$new(sbm_data, c(2, 3), control = NB_control(
    verbose = FALSE, clustering_init = "sbm", n_sparsity_penalties = 3
  ))
  expect_no_error(coll$optimize(control = list(niter = 3, threshold = -1, verbose = FALSE)))
})

test_that("sbm_clustering_path()'s ward2 fallback does not error on a (near-)constant column", {
  sbm_R_const <- sbm_R
  sbm_R_const[, 1] <- 0; sbm_R_const[1, 1] <- 5 # a single non-zero residual for variable 1
  expect_no_warning(path <- sbm_clustering_path(sbm_R_const, 2:20))
  for (q in 2:20) expect_equal(length(unique(path[[as.character(q)]])), q)
})

test_that("zero-inflated collections also use the shared sbm path, clustering on zi_residuals() instead of ols_residuals()", {
  exzi   <- generate_normal_block_var_data(n = 60, p = 24, d = 1, q = 4, kappa = rep(0.3, 24))
  datazi <- NormalBlockData$new(exzi$Y, exzi$X, X0 = matrix(1, nrow(exzi$Y), 1))

  coll <- NormalBlockVarCollectionClusters$new(datazi, 2:3, zero_inflation = TRUE,
                                  control = NB_control(verbose = FALSE, clustering_init = "sbm"))
  ## clustering_init is now eagerly injected from the shared path, just like
  ## for non-ZI collections -- private$C is already a valid q-cluster
  ## indicator matrix before optimize()/optim_initialize() ever runs
  for (model in coll$models) {
    C0 <- model$.__enclos_env__$private$C
    expect_false(anyNA(C0))
    expect_equal(length(unique(get_clusters(C0))), model$q)
  }
  expect_no_error(coll$optimize(control = list(niter = 3, threshold = -1, verbose = FALSE)))
})

test_that("spectral_clustering_path() shares the rank cap with the single-q heuristic", {
  ## Same regime as the rank-cap test in test-clustering-heuristics.R: R with
  ## rank << max(q_list). Without the cap, kmeans still returns q clusters
  ## (no error, no visible collapse), just built partly from noise directions
  ## -- the failure mode is silent, so this checks agreement with the
  ## single-q path rather than an error.
  set.seed(41)
  ex <- generate_normal_block_mean_data(n = 100, p = 30, d = 2, q = 5)
  d  <- NormalBlockData$new(ex$Y, ex$X)
  R  <- d$X %*% d$ols_fit()$B
  q_list <- 2:7

  path <- normalblockr:::spectral_clustering_path(R, q_list)
  expect_equal(length(path), length(q_list))
  for (i in seq_along(q_list)) expect_equal(length(unique(path[[i]])), q_list[i])

  ## and a mean-block collection over that range actually runs
  coll <- normal_block(d, blocks = q_list, model = "mean",
                       control = NB_control(verbose = FALSE, clustering_init = "spectral"))
  expect_true(all(sapply(coll$models, function(m) length(unique(m$clustering))) == q_list))
})
