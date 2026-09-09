## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockVarBase ############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Base Class for Variance-Block Models
#'
#' R6 abstract class for the sparse Normal-Block models, where the clustering
#' structures the latent covariance.
#' @examples
#' # An internal abstract base class, never instantiated directly. See
#' # normal_block() for how concrete models (NormalBlockVarKnownClusters,
#' # NormalBlockVarUnknownClusters, and their zero-inflated variants) are
#' # actually created and fitted.
#' @keywords internal
NormalBlockVarBase <- R6::R6Class(
  classname = "NormalBlockVarBase",
  inherit   =  NormalBlockBase,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(

    #' @description Create a new [`NormalBlockVarBase`] object.
    #' @param data object of NormalBlockData class, with responses and design matrix
    #' @param q number of block/cluster
    #' @param sparsity sparsity penalty on the network density
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @param zero_inflation whether the concrete subclass models zero-inflation;
    #' set by the ZI subclasses themselves, not meant to be set by the end user.
    #' When `FALSE`, the (costly) zero-inflation probability fit (`kappa`/`B0`)
    #' is skipped entirely, since it would otherwise never be used downstream.
    #' @return A new [`NormalBlockVarBase`] object
    initialize = function(data, q, sparsity = 0, control = NB_control(), zero_inflation = FALSE) {
      ## family defaults (see NormalBlockMeanBase for the other ones)
      if (is.null(control$clustering_init)) control$clustering_init <- "ward2"
      if (is.null(control$noise_covariance)) control$noise_covariance <- "diagonal"
      stopifnot("noise_covariance = 'full' only applies to mean-block models, whose Sigma is the full p x p residual covariance" =
                  control$noise_covariance %in% c("diagonal", "spherical"))
      super$initialize(data, q, sparsity, control, zero_inflation)
      ## penalty mask
      private$sparsity_ <- sparsity
      weights <- matrix(1, q, q)
      diag(weights) <- 0
      if (!is.null(control$sparsity_weights)) {
        weights <- control$sparsity_weights
      }
      private$weights <- weights
      ## variant (either diagonal or spherical residuals covariance)
      private$res_covariance <- control$noise_covariance
    },

    #' @description Seed this model's starting parameters from another,
    #' already-optimized model with the same q, instead of a fresh heuristic
    #' clustering. Used by [NormalBlockVarCollectionSparsity] to warm-start
    #' each penalty in a sparsity path from the previous one's solution.
    #' @param other a [NormalBlockVarBase] object, already optimized
    #' @return Update the current object in place with `other`'s parameters
    warm_start_from = function(other) {
      stopifnot("warm_start_from() requires both models to have the same q" = self$q == other$q)
      args <- other$model_par
      args$B0 <- NULL
      if (!is.null(other$var_par)) {
        args$C     <- other$var_par$tau
        args$M     <- other$var_par$M
        args$S     <- other$var_par$S
        args$alpha <- colMeans(other$var_par$tau)
      } else if (!is.null(other$posterior_par)) {
        args$gamma <- other$posterior_par$gamma
        args$mu    <- other$posterior_par$mu
      }
      do.call(self$update, args)
      private$warm_started <- TRUE
      invisible(self)
    },

    #' @description Create a clone of the current [`NormalBlockVarBase`] object after splitting cluster `cl`
    #' We split the cluster according to the species variances
    #' @param index index (integer) of the cluster to split
    #' @param in_place should the split applied to the object itself, or should a copy be sent?
    #' default FALSE (send a copy)
    #' @return A new [`NormalBlockVarBase`] object
    split = function(index, in_place = FALSE) {
      ## update private fields related to group parameters
      ## C, Omega, M, S, sparsity_weights

      ## indices of individuals split within the cluster
      cl  <- self$clustering == index
      var <- 1/private$dm1; var_median <- median(var[cl])
      split1 <- (var > var_median) & cl ;  split2 <- (var <= var_median) & cl

      ## Cluster split
      new_C <- cbind(private$C, .Machine$double.eps)
      new_C[split1, index] <- new_C[split1, index] - .Machine$double.eps
      new_C[split2, self$q + 1] <- new_C[split2, index]
      new_C[split2, index] <- .Machine$double.eps
      new_C <- new_C / rowSums(new_C)

      ## New block starts as a copy of its parent's M/S column; split1/split2
      ## index variables (length p) and must not be reused on M/S, which
      ## are indexed by *individuals* (length n).
      new_M <- cbind(private$M, private$M[, index])

      ## Variational variances
      if (is.matrix(private$S)) {
        new_S <- cbind(private$S, private$S[, index])
      } else {
        new_S <- c(private$S, private$S[index])
      }

      ## Sparsity weights
      if (self$q == 1) {
        new_weights <- matrix(c(0,1,1,0), 2, 2)
      } else {
        weights_cl <-  private$weights[index, setdiff(1:self$q, index)]
        weights_cl <-  c(weights_cl, mean(weights_cl))
        new_weights <- cbind(rbind(private$weights, weights_cl, deparse.level = 0),
                             c(weights_cl, 0))
      }

      ## Precision matrix re-derived from new_M/new_S.
      new_Omega <- private$omega_from_M_S(new_M, new_S, new_weights)

      ## Mark as already initialized so optim_initialize() reuses this state
      ## instead of a fresh heuristic clustering.
      if (in_place) {
        self$update(C = new_C, Omega = new_Omega, M = new_M, S = new_S,
                    alpha = colMeans(new_C), warm_started = TRUE)
        self$sparsity_weights <- new_weights
        return(invisible(self))
      } else {
        new_NB <- self$clone()
        new_NB$update(C = new_C, Omega = new_Omega, M = new_M, S = new_S,
                      alpha = colMeans(new_C), warm_started = TRUE)
        new_NB$sparsity_weights <- new_weights
        return(invisible(new_NB))
      }
    },

    #' @description Create a clone of the current [`NormalBlockVarBase`] object after merging clusters `cl1` and `cl2`
    #' @param indices indices (couple of integer) of the clusters to merge
    #' @param in_place should the split applied to the object itself, or should a copy be sent?
    #' default FALSE (send a copy)
    #' @return A new [`NormalBlockVarBase`] object
    merge = function(indices, in_place=FALSE) {

      ## sorting by increasing group label
      indices <- sort(indices)

      new_C <- private$C[, -indices[2], drop = FALSE]
      new_C[, indices[1]] <- private$C[, indices[1]] + private$C[, indices[2]]

      ## Variational means
      new_M <- private$M[, -indices[2], drop = FALSE]
      new_M[, indices[1]] <- .5 * (private$M[, indices[1]] + private$M[, indices[2]])

      ## Variational variances
      if (is.matrix(private$S)) {
        new_S <- private$S[, -indices[2], drop = FALSE]
        new_S[, indices[1]] <- .5 * (private$S[, indices[1]] + private$S[, indices[2]])
      } else {
        new_S <- private$S[-indices[2]]
        new_S[indices[1]] <- .5 * (private$S[indices[1]] + private$S[indices[2]])
      }

      ## Sparsity weights
      new_weights <-  private$weights[-indices[2], -indices[2], drop = FALSE]

      ## Precision matrix re-derived from new_M/new_S.
      new_Omega <- private$omega_from_M_S(new_M, new_S, new_weights)

      ## See split()'s comment: mark as warm-started so optim_initialize()
      ## reuses this state instead of a fresh heuristic Sigmaq/Omega.
      if (in_place) {
        self$update(C = new_C, Omega = new_Omega, M = new_M, S = new_S,
                    alpha = colMeans(new_C), warm_started = TRUE)
        self$sparsity_weights <- new_weights
        return(self)
      } else {
        new_NB <- self$clone()
        new_NB$update(C = new_C, Omega = new_Omega, M = new_M, S = new_S,
                      alpha = colMeans(new_C), warm_started = TRUE)
        new_NB$sparsity_weights <- new_weights
        return(new_NB)
      }
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    ## heuristics tried by best_of_inits(), best first (see
    ## inst/methods_initialization_and_refine.md for the benchmark)
    default_inits     = c("ward2", "kmeans", "spectral"),
    dm1               = NA, # diagonal vector of inverse variance matrix (variables level)
    gamma             = NA, # variance of  posterior distribution of W
    mu                = NA, # mean for posterior distribution of W
    M                 = NA, # variational mean for posterior distribution of W
    S                 = NA, # variational diagonal of variances for posterior distribution of W
    res_covariance    = NA, # shape of the residuals covariance (diagonal or spherical)

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for integrated (V)EM inference --------------

    ## Closed-form M-step estimate of Omega from M/S (Sigma_hat = M'M/n +
    ## diag(S)); used by split()/merge() to reseed Omega from the new
    ## clustering's own M/S/C rather than the parent's stale Omega.
    ## Cluster pairs are ranked by how strongly related they are in the
    ## current fit: merging two nearly independent blocks is rarely
    ## competitive.
    merge_score = function(pairs) {
      map_dbl(pairs, function(ij) abs(private$Omega[ij[1], ij[2]]))
    },

    omega_from_M_S = function(M, S, weights) {
      s_vec <- if (is.matrix(S)) colMeans(S) else S
      Sigma_hat <- crossprod(M) / self$n + diag(s_vec, ncol(M))
      Omega <- if (private$sparsity_ == 0) {
        chol2inv(chol(Sigma_hat))
      } else {
        glasso_omega(Sigma_hat, private$sparsity_ * weights)
      }
      ensure_pd(Omega)
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for heuristic inference----------------------
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ## Per-variable inverse variance from a residual matrix ("diagonal": one
    ## dm1 per variable; "spherical": one shared value). ddiag is floored
    ## away from 0 to avoid dm1 = Inf on a near-constant residual.
    dm1_from_residuals = function(R) {
      ddiag <- pmax(colMeans(R^2), .Machine$double.eps)
      switch(private$res_covariance,
             "diagonal"  = 1 / as.vector(ddiag),
             "spherical" = rep(1 / mean(ddiag), self$p))
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## MLE of ZI Diagonal Normal distribution
    ## (the weighted-least-squares fit itself lives in zi_weighted_fit(),
    ## R/utils.R, shared with the collection-level zi_residuals())
    zi_diag_normal_inference = function(){
      fit <- self$data$zi_ols_fit()
      list(B = fit$B, dm1 = fit$dm1, kappa = private$kappa, R = fit$R)
    },

    heuristic_optimize = function(control){
      parameters <- private$get_heuristic_parameters()
      c(parameters, list(ll_list = NA))
    },

    heuristic_Sigmaq_from_Sigma = function(Sigma){
      Sigma_q <- (t(private$C) %*% Sigma %*% private$C) / outer(colSums(private$C), colSums(private$C))
      ## NA when a cluster is empty (colSums(C) == 0 there -> division by 0)
      if (anyNA(Sigma_q)) {
        diag(Sigma_q)[is.na(diag(Sigma_q))] <- mean(diag(Sigma_q)[!is.na(diag(Sigma_q))])
        Sigma_q[is.na(Sigma_q)] <- 0
      }
      Sigma_q
    },

    ## Registry of clustering heuristics used to turn the OLS/ZI residuals R
    ## (n x p) into an initial clustering of the p variables into self$q
    ## groups (a vector of length p with values in 1:q), selectable via
    ## NB_control(clustering_init = ...). See
    ## inst/methods_initialization_and_refine.md and
    ## inst/clustering_initialization_benchmark for the rationale and the
    ## empirical comparison behind the "ward2" default.
    clustering_methods = list(
      kmeans   = function(R, q) kmeans(t(R), q, nstart = 30, iter.max = 50)$cluster,
      ## ward2_tree() (R/utils.R) also backs sbm_clustering_path()'s own
      ## fallback. Same computation, shared rather than duplicated.
      ward2    = function(R, q) cutree(ward2_tree(R), q),
      sbm      = function(R, q) {
        options <- list(verbosity = 0, exploreMin = q, exploreMax = q, plot = FALSE, nbCores = 1)
        mySBM <- sbm::estimateSimpleSBM(cov(R), "gaussian", estimOptions = options)
        mySBM$setModel(q)
        mySBM$memberships
      },
      spectral = function(R, q) {
        U <- eigen(cov(R), symmetric = TRUE)$vectors[, seq_len(q), drop = FALSE]
        U <- U / pmax(sqrt(rowSums(U^2)), 1e-10)
        kmeans(U, q, nstart = 30, iter.max = 50)$cluster
      }
    ),

    heuristic_clustering = function(R) {
      clustering <- private$clustering_methods[[private$clustering_approx]](R, self$q)
      if (length(unique(clustering)) < self$q) {
        ## ward2's hclust()/cutree() always yields exactly q groups (modulo
        ## exact tied merge heights), making it a robust fallback whenever
        ## the chosen heuristic collapses to fewer than q clusters.
        clustering <- private$clustering_methods$ward2(R, self$q)
      }
      C <- as_indicator(clustering)
      if (min(colSums(C)) < 1) warning("Initialization failed to place elements in each cluster")
      C
    }

  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field B_original regression coefficients (d x p), converted back to
    #' Y's original units (undoing `NormalBlockData(scale = TRUE)`'s
    #' column-wise rescaling, if any). Use `model_par$B` instead for the
    #' coefficients on the internal fitting scale.
    B_original = function() private$rescale_to_original(private$B, power = 1),
    #' @field d0 number of zi variables (dimensions in X0)
    d0 = function() self$data$d0,
    #' @field model_par a list with the matrices of the model parameters: B (covariates), dm1 (species variance), Omega (groups precision matrix)). On the internal fitting scale (`self$data$Y`, possibly column-rescaled by `NormalBlockData(scale = TRUE)`) -- use `$B_original`/`$dm1_original` for the same quantities converted back to Y's original units.
    model_par = function() list(B = private$B, B0 = private$B0,
                                dm1 = private$dm1, Omega = private$Omega),
    #' @field nb_param number of parameters in the model
    nb_param = function() {
      nb_param_D <- ifelse(private$res_covariance == "diagonal", self$p, 1)
      as.integer(self$p * self$d + self$q + self$n_edges + nb_param_D)
    },
    #' @field sparsity_weights (weights associated to each pair of groups)
    sparsity_weights = function(value) {
      if (missing(value)) {
        private$weights
      } else {
        stopifnot("must be a q x q matrix" =
                    all(is.matrix(value), nrow(value) == ncol(value), ncol(value) == self$q))
        private$weights <- value
      }
    },
    #' @field dm1_original inverse residual variance per variable
    #' (1 / Var(Y_j)), converted back to Y's original units. Use
    #' `model_par$dm1` instead for the internal fitting scale. With
    #' `noise_covariance = "spherical"`, `model_par$dm1` is a single value
    #' repeated p times; once
    #' converted back per-variable, the p values returned here generally
    #' differ from one another whenever Y's columns were rescaled by
    #' different factors.
    dm1_original = function() private$rescale_to_original(private$dm1, power = -2)
  )
)
