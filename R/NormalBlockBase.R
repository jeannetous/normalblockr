## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockBase ############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Root Base Class for Normal-Block Models
#'
#' R6 abstract class shared by the variance-block ([NormalBlockVarBase]) and
#' mean-block ([NormalBlockMeanBase]) model families.
#' @examples
#' # An internal abstract base class, never instantiated directly -- see
#' # normal_block() for how concrete models are created and fitted.
#' @keywords internal
NormalBlockBase <- R6::R6Class(
  classname = "NormalBlockBase",
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field data object of NormalBlockData class, with responses and design matrix
    data  = NULL,

    #' @description Create a new [`NormalBlockBase`] object.
    #' @param data object of NormalBlockData class, with responses and design matrix
    #' @param q number of block/cluster
    #' @param sparsity sparsity penalty on the network density
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @param zero_inflation whether the concrete subclass models zero-inflation;
    #' set by the ZI subclasses themselves, not meant to be set by the end user.
    #' When `FALSE`, the (costly) zero-inflation probability fit (`kappa`/`B0`)
    #' is skipped entirely, since it would otherwise never be used downstream.
    #' @return A new [`NormalBlockBase`] object
    initialize = function(data, q, sparsity = 0, control, zero_inflation = FALSE) {
      self$data <- data
      stopifnot("There cannot be more blocks than there are entities to cluster" = q <= ncol(self$data$Y))

      ## pointer to the chosen optimization function
      private$optimizer <- if (control$heuristic) private$heuristic_optimize else private$EM_optimize
      stopifnot("this model family does not implement NB_control(heuristic = TRUE)" =
                  is.function(private$optimizer))
      private$approx <- control$heuristic


      ## control$clustering_init is either the name of a clustering heuristic
      ## or an actual clustering to use directly
      cl0 <- control$clustering_init
      if (is.character(cl0) && length(cl0) == 1) {
        private$clustering_approx <- cl0
        private$C <- matrix(NA, self$data$n, q)
      } else if (!is.null(cl0)) {
        if (!is.vector(cl0) & !is.matrix(cl0)) stop("Labels must be encoded in vector of labels or indicator matrix")
        if (is.vector(cl0)) {
          if (any(cl0 < 1 | cl0 > q))
            stop("Cluster labels must be between 1 and q")
          if (length(cl0) != self$p)
            stop("Cluster labels must match the number of Y's columns")
          if (length(unique(cl0)) != q)
            stop("The number of clusters in the initial clustering must be equal to q.")
          cl0 <- as_indicator(cl0)
        } else {
          if (nrow(cl0) != self$p)
            stop("Cluster-indicating matrix must have as many rows as Y has columns")
          if (ncol(cl0) != q)
            stop("Cluster-indicating matrix must have q columns")
          if (min(colSums(cl0)) < 1)
            stop("There cannot be empty clusters in the initial clustering matrix.")
        }
        private$C <- cl0
      } else {
        private$C <- matrix(NA, self$data$n, q)
      }

      if (zero_inflation) private$fit_zi_component()
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Setters    ------------------------
    #' @description
    #' Update a [`NormalBlockBase`] object
    #'
    #' All possible parameters of the child classes
    #' @param B regression matrix [all]
    #' @param dm1 diagonal vector of inverse variance matrix (variables level) [NBVar]
    #' @param C the matrix of groups memberships (posterior probabilities) [all]
    #' @param Omega inverse variance matrix (cluster-level for Normal Block
    #' models, variable-level for Normal Mean Block models) [all]
    #' @param gamma  variance of posterior distribution of W [NBVar - known]
    #' @param mu mean for posterior distribution of W [NBVar - known]
    #' @param kappa vector of zero-inflation probabilities [ZINBVar, ZINBMean]
    #' @param alpha vector of groups probabilities [NBVar]
    #' @param M variational mean for posterior distribution of W [NBVar - unknown]
    #' @param S variational diagonal of variances for posterior distribution of W [NBVar - unknown]
    #' @param Psi variational expectation of C'Omega C, intermediary term in calculations [NBMean - unknown]
    #' @param Phi variational correction term used in the Psi/ELBO computations [NBMean - unknown]
    #' @param Lambda variational correction term used in the Sigma-hat update [NBMean - unknown]
    #' @param ll_list  list of log-lik (elbo) values
    #' @param warm_started whether `optim_initialize()` should treat the model as
    #' already initialized (reuse B/Omega/dm1/C/alpha/M/S as they stand)
    #' rather than recomputing a fresh heuristic initialization. Set by
    #' [warm_start_from()] and by [split()]/[merge()].
    #' @param clustering_init initial clustering
    #' @return Update the current [`normal`] object
    update = function(B = NA,
                      dm1 = NA,
                      C = NA,
                      Omega = NA,
                      gamma = NA,
                      mu = NA,
                      kappa = NA,
                      alpha = NA,
                      M = NA,
                      S = NA,
                      Psi = NA,
                      Phi = NA,
                      Lambda = NA,
                      ll_list = NA,
                      warm_started = NA,
                      clustering_init = NA) {
      if (!anyNA(B))       private$B       <- B
      if (!anyNA(dm1))     private$dm1     <- dm1
      if (!anyNA(C))       private$C       <- C
      if (!anyNA(Omega))   private$Omega   <- Omega
      if (!anyNA(gamma))   private$gamma   <- gamma
      if (!anyNA(kappa))   private$kappa   <- kappa
      if (!anyNA(mu))      private$mu      <- mu
      if (!anyNA(alpha))   private$alpha   <- alpha
      if (!anyNA(M))       private$M       <- M
      if (!anyNA(S))       private$S       <- S
      if (!anyNA(Psi))     private$Psi     <- Psi
      if (!anyNA(Phi))     private$Phi     <- Phi
      if (!anyNA(Lambda))  private$Lambda     <- Lambda
      if (!anyNA(ll_list)) private$ll_list <- ll_list
      if (!anyNA(warm_started)) private$warm_started <- warm_started
      if (!anyNA(clustering_init)) {
        stopifnot("clustering_init must be a single heuristic name" =
                    is.character(clustering_init) && length(clustering_init) == 1)
        private$clustering_approx <- clustering_init
        private$C <- matrix(NA, self$p, self$q)
        private$warm_started <- FALSE
      }
    },

    #' @description calls optimization (EM or heuristic) and updates relevant fields
    #' @param control a list for controlling the optimization process
    #' @param warn whether to warn when the (V)EM stops at the `niter` cap
    #' without reaching `threshold` (see `private$warn_if_not_converged()`).
    #' Set to `FALSE` for deliberately-truncated trial fits (cheap candidate
    #' scoring in `candidates_split()`/`candidates_merge()`, the sparsity-path
    #' warm-start probe in [NormalBlockCollectionSparsity]) where stopping at
    #' the cap is expected and not a sign of trouble.
    #' @return optimizes the model and updates its parameters
    optimize = function(control = list(niter = 500, threshold = 1e-4,
                                       fixed_point_niter = 5), warn = TRUE) {
      private$niter_max  <- control$niter
      private$threshold  <- control$threshold
      optim_out <- private$optimizer(control)
      do.call(self$update, optim_out)
      if (warn) private$warn_if_not_converged()
      invisible(self)
    },

    #' @description Try several clustering-initialization heuristics and keep
    #' the best-ELBO converged fit (see `NB_control(clustering_init = )` and
    #' `inst/methods_initialization_and_refine.md` for the rationale). Every
    #' candidate is first screened with a short `trial_niter` run (same idea
    #' as `candidates_split()`/`candidates_merge()`), and only the
    #' `max_training` best-screened ones are fully retrained with `control`.
    #' @param inits vector of clustering-heuristic names to try; defaults to
    #' the model family's own preferred order (`private$default_inits`)
    #' @param trial_niter number of (V)EM iterations used to cheaply screen
    #' every candidate in `inits` before fully retraining the best few
    #' @param max_training how many of the screened candidates (best `loglik`
    #' after `trial_niter` iterations) get fully retrained with `control`
    #' @param control `optimize()` control list (`niter`/`threshold`) used
    #' for the final full retraining of the `max_training` best candidates
    #' @return a new, already-optimized model. Does not mutate the current
    #' object; reassign the result (`model <- model$best_of_inits()`).
    best_of_inits = function(inits = private$default_inits,
                             trial_niter = 10, max_training = 2,
                             control = list(niter = 500, threshold = 1e-4, fixed_point_niter = 5)) {
      stopifnot(
        "best_of_inits() requires (V)EM inference (not the heuristic-only mode, see NB_control(heuristic = ))" =
          !private$approx,
        "best_of_inits() only applies when the initial clustering is inferred by a heuristic (see NB_control(clustering_init = ))" =
          !is.na(private$clustering_approx)
      )
      candidates <- map(inits, function(init) {
        cand <- self$clone()
        cand$update(clustering_init = init)
        cand$optimize(private$trial_control(trial_niter), warn = FALSE)
        cand
      })
      ibest <- order(map_dbl(candidates, "loglik"), decreasing = TRUE)[1:min(max_training, length(candidates))]
      best_candidates <- candidates[ibest]
      map(best_candidates, function(cand) cand$optimize(control, warn = FALSE))
      best_candidates[[which.max(map_dbl(best_candidates, "loglik"))]]
    },

    #' @description generate and select a set of candidate models
    #' by splitting the clusters of the current model
    #' @param trial_niter number of (V)EM iterations used to cheaply score each
    #' candidate before the caller fully re-optimizes the best few -- kept
    #' short on purpose.
    candidates_split = function(trial_niter = 5) {
      # do not split groups with less than 2 guys
      candidates <- map((1:self$q)[self$cluster_sizes > 1], self$split)
      # keep candidates with at least 2 guys per cluster and non empty split.
      # Compared against the number of currently *live* (non-empty) clusters
      # rather than self$q: the (V)EM can converge to an empty cluster,
      # in which case self$q overcounts and every split would otherwise look
      # invalid, silently stalling explore_forward() before n_clusters_range[2].
      clustering_sizes <- map(candidates, "clustering") %>% map(table)
      min_sizes  <- clustering_sizes %>% map_dbl(min)
      n_clusters <- clustering_sizes %>% map_dbl(length)
      n_live <- length(unique(self$clustering))
      candidates <- candidates[min_sizes > 1 & n_clusters == n_live + 1]

      for (i in seq_along(candidates))
        candidates[[i]]$optimize(private$trial_control(trial_niter), warn = FALSE)
      candidates
    },

    #' @description generate and select a set of candidate models
    #' by merging the clusters of the current model
    #' @param max_candidates merge candidates are, unlike split's, quadratic
    #' in q (`choose(q, q-2)` pairs): beyond `max_candidates` pairs, only
    #' the most promising ones are actually built and trial-optimized, ranked
    #' by the family's own `private$merge_score()`. Set to `Inf` to always
    #' try every pair.
    #' @param trial_niter see [candidates_split()]
    candidates_merge = function(max_candidates = 30, trial_niter = 2) {
      stopifnot("need at least two clusters to merge them" = self$q > 1)
      pairs <- combn(self$q, 2, simplify = FALSE)
      if (length(pairs) > max_candidates)
        pairs <- pairs[order(private$merge_score(pairs), decreasing = TRUE)[1:max_candidates]]
      candidates <- map(pairs, self$merge)
      for (i in seq_along(candidates))
        candidates[[i]]$optimize(private$trial_control(trial_niter), warn = FALSE)
      candidates
    },

    #' @description Predicts observations Y for new covariates X, in Y's
    #' original units (like `$fitted`, so that predicting on the training X
    #' reproduces it).
    #' @param new_X new set of covariates.
    #' @return A n*p prediction matrix for new observations
    predict = function(new_X){
      private$rescale_to_original(new_X %*% private$B)
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Extractors ------------------------
    #' @description Extract interaction network in the latent space, as a
    #' matrix rather than a plot -- see `$plot_network()` to plot it instead.
    #' @param type edge value in the network. Can be "support" (binary edges), "precision" (coefficient of the precision matrix) or "partial_cor" (partial correlation between species)
    #' @return a square matrix of size `self$q`
    latent_network = function(type = c("partial_cor", "support", "precision")) {
      net <- switch(
        match.arg(type),
        "support"     = 1 * (private$Omega != 0 & !diag(TRUE, ncol(private$Omega))),
        "precision"   = private$Omega,
        "partial_cor" = {
          tmp <- -private$Omega / tcrossprod(sqrt(diag(private$Omega))); diag(tmp) <- 1
          tmp
        }
      )
      ## Sparse encoding: igraph::graph_from_adjacency_matrix() fails on dsyMatrix
      Matrix::Matrix(net, sparse = TRUE)
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Graphical methods------------------
    #' @param show_increment whether to add a second panel with the (log10)
    #' absolute increment between iterations and the convergence `threshold`
    #' (distinguishes true convergence from a flat-looking objective trace).
    #' @description plots the evolution of the objective (log-likelihood or ELBO)
    #' across the (V)EM iterations of the last call to `optimize()`.
    #' @return a [`ggplot2::ggplot`] graph
    plot_loglik = function(show_increment = TRUE) {
      ll  <- private$ll_list
      obj <- self$objective
      if (length(obj) == 0 || all(is.na(obj))) {
        message("No objective trace to plot (heuristic inference does not compute a log-likelihood/ELBO).")
        return(invisible(NULL))
      }

      ## log10(-obj) spreads out near-convergence iterations on the plot,
      ## same trick as the increment panel below; raw value if not all-negative.
      obj_is_neg <- all(obj < 0)
      obj_facet  <- if (obj_is_neg) "log10(-objective)" else "objective"
      obj_value  <- if (obj_is_neg) log10(-obj) else obj

      dplot <- tibble::tibble(iteration = seq_along(obj), value = obj_value, facet = obj_facet)
      last_increment <- abs(diff(ll))[length(obj)]
      converged <- !is.na(private$threshold) && last_increment < private$threshold
      subtitle  <- if (is.na(private$niter_max)) {
        sprintf("%d iterations", length(obj))
      } else if (converged) {
        sprintf("converged after %d iterations (threshold = %.1e)", length(obj), private$threshold)
      } else {
        sprintf("stopped at the %d-iteration cap -- last increment (%.3g) still above threshold (%.1e)",
                private$niter_max, last_increment, private$threshold)
      }

      if (show_increment) {
        dplot <- dplyr::bind_rows(
          dplot,
          tibble::tibble(iteration = seq_along(obj), value = log10(abs(diff(ll))), facet = "log10(|increment|)")
        )
      }
      dplot$facet <- factor(dplot$facet, levels = c(obj_facet, "log10(|increment|)"))

      p <- ggplot2::ggplot(dplot, ggplot2::aes(x = iteration, y = value)) +
        ggplot2::geom_line() + ggplot2::geom_point(size = .8) +
        ggplot2::ggtitle(label = "(V)EM optimization", subtitle = subtitle) +
        ggplot2::xlab("iteration") + ggplot2::ylab(NULL) +
        ggplot2::theme_bw()

      if (show_increment) {
        p <- p + ggplot2::facet_wrap(~ facet, ncol = 1, scales = "free_y", strip.position = "left")
        if (!is.na(private$threshold))
          p <- p + ggplot2::geom_hline(
            data = data.frame(facet = factor("log10(|increment|)", levels = levels(dplot$facet)),
                              yintercept = log10(private$threshold)),
            ggplot2::aes(yintercept = yintercept),
            linetype = "dashed", colour = "red", alpha = .6
          )
      }
      p
    },

    #' @description plot the latent network. To extract the network as a
    #' matrix instead of plotting it, use `$latent_network()`.
    #' @param type edge value in the network. Either "precision" (coefficient of the precision matrix) or "partial_cor" (partial correlation between species).
    #' @param output Output type. Either `igraph` (for the network) or `corrplot` (for the adjacency matrix)
    #' @param edge.color Length 2 color vector. Color for positive/negative edges. Default is `c("#F8766D", "#00BFC4")`. Only relevant for igraph output.
    #' @param node.labels vector of character. The labels of the nodes. The default will use the column names ot the response matrix.
    #' @param remove.isolated if `TRUE`, isolated node are remove before plotting. Only relevant for igraph output.
    #' @param layout an optional igraph layout. Only relevant for igraph output.
    #' @param plot logical. Should the final network be displayed or only sent back to the user. Default is `TRUE`.
    plot_network = function(type            = c("partial_cor", "support"),
                            output          = c("igraph", "corrplot"),
                            edge.color      = c("#F8766D", "#00BFC4"),
                            remove.isolated = FALSE,
                            node.labels     = NULL,
                            layout          = igraph::layout_in_circle,
                            plot = TRUE) {
      if(anyNA(private$Omega)) stop("NA in the precision matrix")

      type   <- match.arg(type)
      output <- match.arg(output)
      net <- self$latent_network(type)

      if (output == "igraph") {
        G <-  igraph::graph_from_adjacency_matrix(net, mode = "undirected", weighted = TRUE, diag = FALSE)

        if (!is.null(node.labels)) {
          igraph::V(G)$label <- node.labels
        } else {
          igraph::V(G)$label <- unlist(lapply(1:ncol(net), f <- function(x) paste0("Cluster_", x)))
        }
        ## Nice nodes
        V.deg <- igraph::degree(G)/sum(igraph::degree(G))
        igraph::V(G)$label.cex <- V.deg / max(V.deg) + .5
        igraph::V(G)$size <- tabulate(self$clustering, nbins = nrow(self$model_par$Omega)) * 100 / self$p
        igraph::V(G)$label.color <- rgb(0, 0, .2, .8)
        igraph::V(G)$frame.color <- "transparent" # NA triggers an igraph plot() warning; same (invisible) rendering
        ## Nice edges
        igraph::E(G)$color <- ifelse(igraph::E(G)$weight > 0, edge.color[1], edge.color[2])
        if (type == "support")
          igraph::E(G)$width <- abs(igraph::E(G)$weight)
        else
          igraph::E(G)$width <- 15*abs(igraph::E(G)$weight)

        if (remove.isolated) {
          G <- igraph::delete.vertices(G, which(igraph::degree(G) == 0))
        }
        if (plot) plot(G, layout = layout)
      }
      if (output == "corrplot") {
        if (plot) {
          if (ncol(net) > 100)
            colnames(net) <- rownames(net) <- rep(" ", ncol(net))
          G <- net
          diag(net) <- 0
          corrplot::corrplot(as.matrix(net), method = "color", is.corr = FALSE, tl.pos = "td", cl.pos = "n", tl.cex = 0.5, type = "upper")
        } else  {
          G <- net
        }
      }
      invisible(G)
    },

    #' @description plots the evolution of the objective during model optimization
    #' (see `plot_loglik()`)
    plot = function(){
      self$plot_loglik()
    },


    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## S3 methods ----------------------------
    #' @description User friendly print method
    #' @param model First line of the print output
    print = function(model = paste("A", self$who_am_I, ".\n")) {
      cat(model)
      cat("===========================================================================\n")
      print(as.data.frame(round(self$criteria, digits = 3), row.names = ""))
      cat("===========================================================================\n")
      cat("* Useful fields\n")
      cat("    $model_par, $posterior_par / $var_par, $clustering \n")
      cat("    $loglik, $BIC, $ICL, $objective, $nb_param, $criteria\n")
      cat("* Useful S3 methods\n")
      cat("    print(), summary(), plot(), coef(), sigma(), fitted(), predict() \n")
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    kappa             = NA, # matrix of zero-inflation probabilities
    B0                = NA, # zero-inflation regression matrix
    ZI_cond_mean      = NA, # fixed contribution of the ZI component to the log-likelihood
    B                 = NA, # regression matrix
    C                 = NA, # the matrix of posterior probabilities (tau) or group affectation
    Omega             = NA, # precision matrix for clusters or variables
    alpha             = NA, # vector of clusters probabilities
    optimizer         = NA, # a link to the function that perform the optimization
    ll_list           = NA, # list of log-likelihoods or ELBOs
    sparsity_         = NA, # scalar controlling the overall sparsity
    weights           = NA, # sparsity weights specific to each pairs of group
    approx            = NA, # use approximation/heuristic approach or not
    clustering_approx = NA, # name of the clustering heuristic, key into clustering_methods
    niter             = NA, # number of EM iterations required by the inference, if applicable
    niter_max         = NA, # niter cap passed to the last optimize() call (for plot_loglik()'s
    # and warn_if_not_converged()'s diagnostic)
    threshold         = NA, # convergence threshold passed to the last optimize() call (idem)
    warm_started      = FALSE, # set by warm_start_from() and by split()/merge(): tells
    # optim_initialize() to reuse the current B/dm1/Omega (and
    # C/M/S/alpha or gamma/mu) instead of (re-)deriving them
    # from the heuristic clustering.
    # Deliberately NOT inferred from "are these fields non-NA",
    # which would also be true for split()/merge() clones (those
    # keep their own, different, already-tested initialization path).

    ## Warns when the (V)EM hit the niter cap without reaching threshold
    ## (heuristic fits have no ll_list/niter to check, so stay silent).
    ## Control list for the deliberately-truncated trial fits of
    ## best_of_inits()/candidates_split()/candidates_merge(). fixed_point_niter
    ## is only read by the mean-block family; the variance-block one ignores it.
    trial_control = function(trial_niter) {
      list(niter = trial_niter, threshold = 1e-4, fixed_point_niter = 5)
    },

    ## Ranks candidate cluster pairs for candidates_merge(), highest first.
    ## Overridden by each family (see NormalBlockVarBase/NormalBlockMeanBase).
    merge_score = function(pairs) rep(0, length(pairs)),

    warn_if_not_converged = function() {
      if (is.na(private$niter) || private$niter < private$niter_max) return(invisible())
      last_increment <- abs(diff(private$ll_list))[private$niter]
      if (last_increment >= private$threshold)
        warning(sprintf(
          "%s: (V)EM stopped at the niter cap (%d) without reaching the convergence threshold (last increment = %.3g, threshold = %.1e). Consider raising `niter` in NB_control(), especially with many blocks -- see plot_loglik() to check.",
          self$who_am_I, private$niter_max, last_increment, private$threshold), call. = FALSE)
      invisible()
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Zero-inflation probabilities, never revisited by the (V)EM: only their
    ## fixed log-likelihood contribution (ZI_cond_mean) reaches the recursion.
    ## They depend on nothing but the data, so the fit itself lives on NormalBlockData
    fit_zi_component = function() {
      zi <- self$data$zi_fit()
      private$B0           <- zi$B0
      private$kappa        <- zi$kappa
      private$ZI_cond_mean <- zi$ZI_cond_mean
    },

    ## Methods for integrated (V)EM inference --------------
    ## Each concrete subclass overrides EM_optimize() to call its
    ## Rcpp/Armadillo core (see src/exports.cpp and
    ## inst/normal_block_models.qmd); optim_initialize() (heuristic
    ## initialization) stays in R and supplies the starting values.
    ## MLE of MV Normal distribution
    multivariate_normal_inference = function() self$data$ols_fit(),

    optim_initialize = function() {},

    ## Converts a per-variable quantity from the internal (rescaled) fitting
    ## scale back to Y's original units: power = 1 for additive quantities
    ## (B, fitted values), power = -2 for inverse-variance quantities (dm1).
    rescale_to_original = function(M, power = 1) {
      factor <- self$data$Y_scale^power
      if (is.matrix(M)) M * matrix(factor, nrow(M), ncol(M), byrow = TRUE) else M * factor
    },

    get_Omega = function(Sigma) {
      if (private$sparsity_ == 0) {
        chol2inv(chol(Sigma))
      } else {
        glasso_omega(Sigma, private$sparsity_ * self$sparsity_weights)
      }
    },


    ## Registry of clustering heuristics used to turn the OLS/ZI residuals R
    ## (n x p) into an initial clustering of the p variables into self$q
    ## groups (a vector of length p with values in 1:q). One single table
    ## instead of one ad hoc private method per algorithm, selectable via
    ## NB_control(clustering_init = ...) ("ward2"/"kmeans"/"sbm"/"spectral").
    ## Benchmarked separately per family (the two cluster genuinely different
    ## quantities: residual covariance vs. mean trajectory). For the
    ## variance-block family, on three real datasets: no single method
    ## dominates everywhere, but combining each method's BIC rank with how
    ## often its deviance path violates the model's theoretical guarantee
    ## (deviance is non-increasing in q) favors ward2 as the most reliable
    ## default (NormalBlockVarBase$initialize()) -- kmeans has a marginally
    ## better raw BIC rank on average but violates that monotonicity far more
    ## often. For the mean-block family, kmeans is the default instead
    ## (NormalBlockMeanBase$initialize()): it consistently lands in a better
    ## ELBO basin than ward2 or spectral there (inst/mean_block_analyses/).
    ## spectral clusters the eigenvectors of cov(R) (top q, each row rescaled
    ## to unit L2 norm, the classic Ng-Jordan-Weiss normalization).
    clustering_methods = list(
      ## compress_columns() is exact here: kmeans sees only the distances
      ## between R's columns, which it preserves (R/utils.R)
      kmeans   = function(R, q) kmeans_columns(compress_columns(R), q),
      ## ward2_tree() (R/utils.R) also backs sbm_clustering_path()'s own
      ## fallback -- same computation, shared rather than duplicated.
      ward2    = function(R, q) cutree(ward2_tree(R), q),
      sbm      = function(R, q) {
        options <- list(verbosity = 0, exploreMin = q, exploreMax = q, plot = FALSE, nbCores = 1)
        mySBM <- sbm::estimateSimpleSBM(cov(R), "gaussian", estimOptions = options)
        mySBM$setModel(q)
        mySBM$memberships
      },
      ## spectral_eig_rank() (R/utils.R) caps the eigenvectors used at
      ## cov(R)'s numerical rank -- see its comment for why: R = X %*% B for
      ## the mean-block family is routinely rank << q, and eigenvectors past
      ## the rank are arbitrary, not data-derived. Shared with
      ## spectral_clustering_path(), the collection-level analogue.
      spectral = function(R, q) {
        eig <- spectral_eig_rank(R)
        U <- eig$vectors[, seq_len(min(q, eig$rank)), drop = FALSE]
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
    #' @field inference_method inference procedure used (heuristic or integrated with EM)
    inference_method = function() ifelse(private$approx, "heuristic", "integrated"),
    #' @field n number of samples
    n = function() self$data$n,
    #' @field p number of responses per sample
    p = function() self$data$p,
    #' @field d number of variables (dimensions in X)
    d = function() self$data$d,
    #' @field q number of blocks
    q = function() as.integer(ncol(private$C)),
    #' @field n_edges number of edges of the network (non null coefficient of the sparse precision matrix Omega)
    n_edges  = function() sum(private$Omega[upper.tri(private$Omega, diag = FALSE)] != 0),
    #' @field objective evolution of the objective function during (V)EM algorithm
    objective = function() private$ll_list[-1],
    #' @field loglik (or its variational lower bound)
    loglik = function() if (private$approx) NA else private$ll_list[[length(private$ll_list)]] + self$sparsity_term,
    #' @field deviance (or its variational lower bound)
    deviance = function() -2 * self$loglik,
    #' @field BIC (or its variational lower bound)
    #' @field entropy Entropy of the conditional distribution when applicable
    entropy    = function() 0,
    BIC = function() self$deviance + log(self$n) * self$nb_param,
    #' @field ICL variational lower bound of the ICL
    ICL        = function() self$BIC + 2 * self$entropy,
    #' @field EBIC variational lower bound of the EBIC
    EBIC   = function() self$BIC + 2 * ifelse(self$n_edges > 0, self$n_edges * log(self$q), 0),
    #' @field criteria a vector with loglik, BIC and number of parameters
    criteria   = function() {
      data.frame(nb_param = self$nb_param, q = self$q, n_edges = self$n_edges, sparsity = self$sparsity,
                 loglik = self$loglik, deviance = self$deviance, BIC = self$BIC, ICL = self$ICL, EBIC = self$EBIC,
                 niter = private$niter )
    },
    #' @field sparsity (overall sparsity parameter)
    sparsity = function(value) {
      if (missing(value)) {
        private$sparsity_
      } else {
        stopifnot("must be a positive scale" = value >= 0)
        private$sparsity_ <- value
      }
    },
    #' @field sparsity_term (sparsity_term term in log-likelihood due to sparsity)
    sparsity_term = function() self$sparsity * sum(abs(self$sparsity_weights * private$Omega)),
    #' @field memberships cluster memberships
    memberships = function() private$C,
    #' @field clustering given as the list of elements contained in each cluster
    clustering = function() {
      cl <- get_clusters(private$C)
      names(cl) <- colnames(self$data$Y)
      cl
    },
    #' @field cluster_sizes given as a vector of cluster sizes
    cluster_sizes = function() tabulate(self$clustering, nbins = self$q),
    #' @field elements_per_cluster given as the list of elements contained in each cluster
    elements_per_cluster = function() {
      if (is.null(names(self$clustering)))
        base::split(1:self$p, self$clustering)
      else
        base::split(names(self$clustering), self$clustering)
    }
  )
)
