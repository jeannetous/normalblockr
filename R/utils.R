# turns a list clustering of q cluster labels for N elements into a matrix of
# dimensions (p, q) with a one-hot encoding of the clustering
#
# @param clustering a list of labels
as_indicator <- function(clustering) {
  q <- max(clustering)
  p <- length(clustering)
  Z <- matrix(0, p, q)
  Z[cbind(seq.int(p), clustering)] <- 1
  Z
}

# removes machine's 0 to elements equal to  1 in x
check_one_boundary <- function(x, zero = .Machine$double.eps) {
  x[is.nan(x)] <- zero
  x[x >= 1 - zero] <- 1 - zero
  x
}

# adds machine's 0 to elements equal to 0 in x
check_zero_boundary <- function(x, zero = .Machine$double.eps) {
  x[is.nan(x)] <- zero
  x[x < zero]  <- zero
  x
}

# equivalent of check_zero_boundary(check_one_boundary(x)), used to keep
# variational probabilities (tau) away from the {0, 1} boundaries (mirrors
# clip_probabilities() in src/utils_arma.h). Previously spelled out at both
# call sites (NormalBlockVarUnknownClusters.R/ZINormalBlockVarUnknownClusters.R).
clip_probabilities <- function(x, zero = .Machine$double.eps) {
  check_zero_boundary(check_one_boundary(x, zero), zero)
}

# Projects a symmetric matrix onto the PD cone by flooring its eigenvalues.
# Used by split()/merge()'s new_Omega, hand-edited from an existing
# precision matrix with no general guarantee of staying PD.
ensure_pd <- function(M, floor = 1e-6) {
  M <- (M + t(M)) / 2
  eig <- eigen(M, symmetric = TRUE)
  eig$vectors %*% diag(pmax(eig$values, floor), nrow(M)) %*% t(eig$vectors)
}

# Graphical-lasso estimate of a precision matrix from its covariance estimate,
# with a plain inversion as the fallback when the solver can't produce a finite
# answer. Wraps the in-package solver (src/graphical_lasso.h) so that the R
# reference implementations and the C++ (V)EM cores run the very same code --
# tests/testthat/test-cpp-*.R compare the two at 1e-8, which only holds if they
# share this step exactly.
# Convergence threshold the (V)EM's M-step uses for the graphical lasso.
# MUST match nb_omega::kGlassoThreshold (src/omega_estimation.h): the two
# recursions are compared trace-for-trace at 1e-8 in test-cpp-*.R, so any
# drift between them shows up there rather than silently.
NB_GLASSO_THRESHOLD <- 1e-6

glasso_omega <- function(Sigma, rho) {
  glasso_out <- graphical_lasso_fit(Sigma, rho)
  if (anyNA(glasso_out$wi)) {
    warning("GLasso fails, the penalty is probably too small and the system badly ",
            "conditionned \n reciprocal condition number =", rcond(Sigma),
            "\n We send back the original matrix and its inverse (unpenalized).",
            call. = FALSE)
    return(chol2inv(chol(Sigma)))
  }
  Matrix::symmpart(glasso_out$wi)
}

# computes xlogx, setting it to 0 if x = 0
xlogx <- function(x) ifelse(x < .Machine$double.eps, 0, x * log(x))

# computes xlogy, setting it to 0 if x = 0
xlogy <- function(x,y) ifelse(x < .Machine$double.eps, 0, x * log(y))

# computes softmax
softmax <- function(x) {
  b <- max(x)
  exp(x - b) / sum(exp(x - b))
}

# gets cluster labels from probability matrix
get_clusters <- function(tau) {
  apply(tau, 1, which.max)
}

# for a list of edges, give corresponding (node1, node2) list.
edge_to_node <- function(x, n = max(x)) {
  x <- x - 1 ## easier for arithmetic to number edges starting from 0
  n.node <- round((1 + sqrt(1 + 8*n)) / 2) ## n.node * (n.node -1) / 2 = n (if integer)
  j.grid <- cumsum(0:n.node)
  j <- findInterval(x, vec = j.grid)
  i <- x - j.grid[j]
  ## Renumber i and j starting from 1 to stick with R convention
  data.frame(node1 = i + 1, node2 = j + 1)
}

sigmoid <- function(x){
  return(1 / (1 + exp(-x)))
}

# OLS residuals of Y on X, used to seed the clustering heuristics; shared
# across a collection's models instead of recomputing per q.
ols_residuals <- function(data) {
  B <- data$XtXm1 %*% data$XtY
  data$Y - data$X %*% B
}

# Iteratively reweighted least squares fit of B under a zero-inflation mask
# (weights = zeros_bar, dm1 re-estimated between iterates). Shared by
# zi_residuals() and NormalBlockVarBase's zi_diag_normal_inference(). Each
# column of B is solved independently (the mask varies by row and column);
# ginv() rather than solve() since a design level can be all-zero for some
# variable, making XtWX singular. ssq is floored away from 0 to keep dm1
# finite when a variable has very few non-zero observations.
zi_weighted_fit <- function(data) {
  ssq <- function(B) pmax(colSums(data$zeros_bar * (data$Y - data$X %*% B)^2), .Machine$double.eps)

  B   <- data$XtXm1 %*% data$XtY
  dm1 <- data$nY / ssq(B)
  for (i in 1:3) { # a couple of iterates is enough
    DM1 <- matrix(dm1, data$n, data$p, byrow = TRUE) * data$zeros_bar
    for (j in seq_len(data$p)) {
      w <- DM1[, j]
      XtWX <- crossprod(data$X, data$X * w)
      XtWy <- crossprod(data$X, data$Y[, j] * w)
      B[, j] <- MASS::ginv(XtWX) %*% XtWy
    }
    dm1 <- data$nY / ssq(B)
  }
  list(B = B, dm1 = dm1, R = data$zeros_bar * (data$Y - data$X %*% B))
}

# Zero-inflation analogue of ols_residuals(): kappa isn't computed here since
# the residual only depends on the weighted fit of B, not on kappa.
zi_residuals <- function(data) data$zi_ols_fit()$R

# Ward.D2 clustering tree of the p columns of R by pairwise correlation
# distance (1 - cor); shared by the "ward2" heuristic, its fallback for any
# heuristic that collapses to fewer than q clusters, and sbm_clustering_path().
# cor() is NA for a (near-)constant column: treated as uncorrelated (cor = 0)
# rather than letting dist()/hclust() fail on NA input.
ward2_tree <- function(R) {
  cor_R <- suppressWarnings(stats::cor(R))
  cor_R[is.na(cor_R)] <- 0
  stats::hclust(stats::dist(1 - cor_R), method = "ward.D2")
}

# Clusters R into every q in q_list from a SINGLE sbm::estimateSimpleSBM
# exploration over [min(q_list), max(q_list)] rather than one per q (see
# inst/normal_block_models.qmd, "Implementation notes", for the cost/quality
# comparison). The exploration is adaptive and may stop before exploreMax, or
# return an empty block for some q; either falls back to a shared ward2
# clustering for that q. Returns a list of membership vectors, named by q.
sbm_clustering_path <- function(R, q_list) {
  options <- list(verbosity = 0, exploreMin = min(q_list), exploreMax = max(q_list),
                  plot = FALSE, nbCores = 1)
  mySBM <- sbm::estimateSimpleSBM(stats::cov(R), "gaussian", estimOptions = options)
  explored <- mySBM$storedModels$nbBlocks
  fallback_tree <- ward2_tree(R)

  stats::setNames(
    lapply(q_list, function(q) {
      if (q %in% explored) {
        mySBM$setModel(q)
        memberships <- mySBM$memberships
        if (length(unique(memberships)) == q) return(memberships)
      }
      stats::cutree(fallback_tree, q)
    }),
    q_list
  )
}

# Rewrites R (n x p) with fewer rows but the *same* Euclidean distances
# between its columns, which is all a distance-based clustering of those
# columns can see. R = U D V' with U orthonormal gives
# ||R e_j - R e_k|| = ||D V' (e_j - e_k)||, so D V' (rank(R) x p) is an exact
# stand-in. Worth it because the matrix handed to the heuristics is often a
# tall, low-rank embedding: the mean-block family clusters X %*% B, whose rank
# is d (measured 200 x 60 -> 1 x 60 at d = 1), and OLS residuals still drop
# from n to p rows. Correlation-based heuristics cannot use this -- cor() over
# rank(R) rows is a different quantity -- so it is applied inside the kmeans
# method rather than to every heuristic.
compress_columns <- function(R, tol = 1e-10) {
  if (nrow(R) <= 1) return(R)
  s <- svd(R, nu = 0)
  keep <- s$d > tol * max(s$d, 1)
  if (sum(keep) >= nrow(R)) return(R)
  t(s$v[, keep, drop = FALSE] %*% diag(s$d[keep], sum(keep)))
}

kmeans_columns <- function(R, q) {
  stats::kmeans(t(R), q, nstart = 30, iter.max = 50)$cluster
}

# Clusters R into every q in q_list with kmeans. Only the lossless row
# compression is shared -- Lloyd's algorithm itself depends on q throughout.
kmeans_clustering_path <- function(R, q_list) {
  Rc <- compress_columns(R)
  stats::setNames(lapply(q_list, function(q) kmeans_columns(Rc, q)), q_list)
}

# Clusters R into every q in q_list by cutting a SINGLE hierarchical tree,
# rather than rebuilding it once per q. cutree() can return fewer than q groups
# on exactly tied merge heights; that q is left to the model's own heuristic_clustering().
ward2_clustering_path <- function(R, q_list) {
  tree <- ward2_tree(R)
  stats::setNames(lapply(q_list, function(q) stats::cutree(tree, q)), q_list)
}

# cov(R)'s eigendecomposition plus its numerical rank, shared by the spectral
# heuristic's two call sites (private$clustering_methods$spectral in
# R/NormalBlockBase.R, and spectral_clustering_path() below). Eigenvectors
# past the rank are an arbitrary completion of the null space, not derived
# from the data at all, so capping at the rank rather than at whatever q was
# asked for matters: R = X %*% B for the mean-block family is routinely
# rank << q (d small, q explored well above it).
spectral_eig_rank <- function(R) {
  eig <- eigen(stats::cov(R), symmetric = TRUE)
  list(vectors = eig$vectors, rank = max(1L, sum(eig$values > 1e-8 * max(eig$values))))
}

# The eigendecomposition above is q-independent, only the number of leading
# vectors kept and the kmeans that follows are not -- computed once and
# reused across the whole q_list rather than once per q.
spectral_clustering_path <- function(R, q_list) {
  eig <- spectral_eig_rank(R)
  stats::setNames(
    lapply(q_list, function(q) {
      U <- eig$vectors[, seq_len(min(q, eig$rank)), drop = FALSE]
      U <- U / pmax(sqrt(rowSums(U^2)), 1e-10)
      stats::kmeans(U, q, nstart = 30, iter.max = 50)$cluster
    }),
    q_list
  )
}

# Precomputes, for a collection over q_list, whatever part of the requested
# clustering heuristic is shared across q, and returns one clustering per q
# (named by q). Returns NULL when nothing can be shared. "best_of_inits" is
# the model's own business, and an explicit clustering needs no help -- and
# every model then runs its own heuristic_clustering() as before.
#
# `R` is the matrix the family hands to its clustering heuristics: the
# residuals for variance-block models, the fitted mean trajectory for
# mean-block ones. Any q whose shared computation fails to produce exactly q
# groups is dropped back to the heuristic's name, preserving the per-model
# fallback in heuristic_clustering().
clustering_path_for_collection <- function(R, q_list, method) {
  if (!is.character(method) || length(method) != 1) return(NULL)
  path <- switch(method,
                 "sbm"      = sbm_clustering_path(R, q_list),
                 "ward2"    = ward2_clustering_path(R, q_list),
                 "spectral" = spectral_clustering_path(R, q_list),
                 "kmeans"   = kmeans_clustering_path(R, q_list),
                 NULL)
  if (is.null(path)) return(NULL)
  lapply(stats::setNames(seq_along(q_list), q_list), function(i) {
    cl <- path[[i]]
    if (length(unique(cl)) == q_list[i]) cl else method
  })
}

# Entry point for a collection: resolves the family's default heuristic when
# none was named, and hands clustering_path_for_collection() the matrix that
# family's models cluster -- the residuals for variance-block models, the
# fitted mean trajectory X %*% B for mean-block ones (see the call sites in
# NormalBlockVarUnknownClusters/NormalBlockMeanUnknownClusters, which must
# stay in step with this).
clustering_path_for_family <- function(mydata, q_list, family = c("var", "mean"),
                                       zero_inflation = FALSE, control = NB_control()) {
  family <- match.arg(family)
  method <- control$clustering_init
  if (is.null(method)) method <- if (family == "var") "ward2" else "kmeans"
  if (!is.character(method) || length(method) != 1) return(NULL)
  if (!method %in% c("sbm", "ward2", "spectral", "kmeans")) return(NULL)

  R <- if (family == "var") {
    if (zero_inflation) zi_residuals(mydata) else ols_residuals(mydata)
  } else {
    B <- if (zero_inflation) mydata$zi_ols_fit()$B else mydata$ols_fit()$B
    mydata$X %*% B
  }
  clustering_path_for_collection(R, q_list, method)
}

# Whether NB_control(clustering_init = "best_of_inits") was requested (see
# best_of_inits() in NormalBlockVarBase.R). identical() keeps this safe when
# clustering_init is a list or an explicit clustering.
uses_best_of_inits <- function(control) identical(control$clustering_init, "best_of_inits")
