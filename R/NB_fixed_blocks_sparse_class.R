## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_unknown #######################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 generic class for normal-block models
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param C group matrix C_jq = 1 if species j belongs to group q
#' @param penalties list of penalties to be tested
#' @param n_penalties number of penalty values to be tested (if penalties not provided, default is 30)
#' @param min_ratio ratio between min and max penalty to be tested  (if penalties not provided, default is 0.05)
#' @param models uderlying NB_fixed_blocks models for each nb of blocks
#' @param verbose telling if information should be printed during optimization
NB_fixed_blocks_sparse <- R6::R6Class(
  classname = "NB_fixed_blocks_sparse",

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field Y the matrix of responses
    Y  = NULL,
    #' @field X the design matrix
    X = NULL,
    #' @field C the group matrix
    C = NULL,
    #' @field penalties list of penalties values for omegaQ sparsity
    penalties = NULL,
    #' @field n_penalties number of penalty values
    n_penalties = NULL,
    #' @field min_ratio ratio between min and max penalty to be tested
    min_ratio = NULL,
    #' @field models list of NB_fixed_Q models corresponding to each nb_block value
    models = NULL,
    #' @field verbose say whether information should be given about the optimization
    verbose = NULL,
    #' @field latest_niter latest niter value used for optimization
    latest_niter = NULL,
    #' @field latest_threshold latest threshold value used for optimization
    latest_threshold = NULL,

    #' @description Create a new [`NB_fixed_blocks_sparse`] object.
    #' @param Y the matrix of responses (called Y in the model).
    #' @param X design matrix (called X in the model).
    #' @param C group matrix.
    #' @param penalties list of penalties on the network density
    #' @param n_penalties number of penalty values
    #' @param min_ratio ratio between min and max penalty to be tested
    #' @return A new [`nb_fixed_blocks_sparse`] object
    initialize = function(Y, X, C, penalties = NULL,  n_penalties = 30,
                          min_ratio = 0.05, verbose=TRUE) {
      if (!is.matrix(Y) || !is.matrix(X)) {
        stop("Y, X and C must be matrices.")
      }
      if (nrow(Y) != nrow(X)) {
        stop("Y and X must have the same number of rows")
      }
      self$Y <- Y
      self$X <- X
      self$C <- C
      if(!is.null(penalties)){
        self$n_penalties <- length(penalties)
        self$penalties <- penalties[order(penalties)]
      }else{
        self$n_penalties <- n_penalties
        self$min_ratio   <- min_ratio
        init_model <- normal_block(Y, X, C, verbose = FALSE)
        init_model$optimize(5)
        sigmaQ    <- solve(init_model$model_par$omegaQ)
        max_pen   <- max(abs(sigmaQ[upper.tri(sigmaQ, diag = FALSE)]))
        penalties <- 10^seq(log10(max_pen), log10(max_pen* self$min_ratio), len = self$n_penalties)
        self$penalties <- penalties[order(penalties)]
        }
      self$verbose <- verbose
    },


    #' @description optimizes an NB_fixed_Q object for each value of Q
    #' @param niter number of iterations in model optimization
    #' @param threshold loglikelihood threshold under which optimization stops
    optimize = function(niter = 100, threshold = 1e-4) {
      self$latest_niter     <- niter
      self$latest_threshold <- threshold

      self$models <- furrr::future_map(seq_along(self$models), function(m) {
        model <- self$models[[m]]
        if(self$verbose) cat("\t penalty =", self$models[[m]]$sparsity, "          \r")
        flush.console()
        model$optimize(niter, threshold)
        if(m < self$n_penalties){
          self$models[[m + 1]]$update(B      = model$model_par$B,
                                      dm1    = model$model_par$dm1,
                                      omegaQ = model$model_par$omegaQ,
                                      gamma  = model$posterior_par$gamma,
                                      mu     = model$posterior_par$mu)
        }
        model
      }, .options = furrr_options(seed=TRUE))
    },

    #' @description returns the NB_fixed_blocks_sparse model corresponding to given penalty
    #' @param penalty penalty asked by user
    #' @return A NB_fixed_blocks_sparse object with given value penalty
    get_model = function(penalty) {
      if(!(penalty %in% self$penalties)) {
        penalty <-  self$penalties[[which.min(abs(self$penalties - penalty))]]
        cat(paste0("No model with this penalty in the collection. Returning model with closest penalty : ", penalty,  " Collection penalty values can be found via $penalties"))
      }
      penalty_rank <- which(sort(self$penalties) == penalty)
      return(self$models[[penalty_rank]])
    },

    #' @description Extract best model in the collection
    #' @param crit a character for the criterion used to performed the selection.
    #' Either "BIC", "AIC" or "loglik" (-loglik so that criterion to be minimized)
    #' "loglik" is the default criterion
    #' @return a [`NB_fixed_Q`] object
    get_best_model = function(crit = c("loglik", "BIC", "AIC", "ICL", "StARS"),
                              stability = 0.9) {
      crit <- match.arg(crit)
      if (crit == "StARS") {
        if (is.null(private$stab_path)) self$stability_selection()
        max_stab <- max(self$criteria$stability)
        if(max_stab < stability){
          cat(paste0("No model reaches the required stability ", stability, ", returning model with highest stability: ", max_stab))
          stability <- max_stab
        }
        id_stars <- self$criteria %>%
                    dplyr::select(sparsity, stability) %>% dplyr::rename(Stability = stability) %>%
                    dplyr::filter(Stability >= stability) %>%
                    dplyr::pull(sparsity) %>% min() %>% match(self$penalties)
        model <- self$models[[id_stars]]$clone()
      }else{
        stopifnot(!anyNA(self$criteria[[crit]]))
        id <- 1
        if (length(self$criteria[[crit]]) > 1) {
          id <- which.min(self$criteria[[crit]])
        }
        model <- self$models[[id]]$clone()
      }
      model
    },

    #' @description Display various outputs (goodness-of-fit criteria, robustness, diagnostic) associated with a collection of network fits (a [`Networkfamily`])
    #' @param criteria vector of characters. The criteria to plot in `c("loglik", "BIC", "AIC", "ICL")`. Defaults to all of them.
    #' @param log.x logical: should the x-axis be represented in log-scale? Default is `TRUE`.
    #' @return a [`ggplot`] graph
    plot = function(criteria = c("loglik", "BIC", "AIC", "ICL"), log.x = TRUE) {
      vlines <- sapply(intersect(criteria, c("BIC")) , function(crit) self$get_best_model(crit)$sparsity)
      stopifnot(!anyNA(self$criteria[criteria]))

      dplot <- self$criteria %>%
        dplyr::select(dplyr::all_of(c("sparsity", criteria))) %>%
        tidyr::gather(key = "criterion", value = "value", -sparsity) %>%
        dplyr::group_by(criterion)
      if("loglik" %in% criteria){dplot[dplot$criterion == "loglik",]$value <- - dplot[dplot$criterion == "loglik",]$value}
      p <- ggplot2::ggplot(dplot, ggplot2::aes(x = sparsity, y = value, group = criterion, colour = criterion)) +
           ggplot2::geom_line() + ggplot2::geom_point() +
           ggplot2::ggtitle(label    = "Model selection criteria",
                  subtitle = "Lower is better" ) +
            ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = vlines, linetype = "dashed", alpha = 0.25)
      if (log.x) p <- p + ggplot2::coord_trans(x = "log10")
      p
    }
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    stab_path = NULL # a field to store the stability path,
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field n number of samples
    n = function() nrow(self$Y),
    #' @field p number of responses per sample
    p = function() ncol(self$Y),
    #' @field d number of variables (dimensions in X)
    d = function() ncol(self$X),
    #' @field Q number of blocks
    Q = function() ncol(self$C),
    #' @field criteria a data frame with the values of some criteria ((approximated) log-likelihood, BIC, AIC) for the collection of models
    criteria = function(){
      crit <- purrr::map(self$models, "criteria") %>% purrr::reduce(rbind)
      crit$stability <- self$stability
      crit},
    stability_path = function() private$stab_path,
    #' @field stability mean edge stability along the penalty path
    stability = function() {
      if (!is.null(private$stab_path)) {
        stability <- self$stability_path %>%
          dplyr::select(Penalty, Prob) %>%
          dplyr::group_by(Penalty) %>%
          dplyr::summarize(Stability = 1 - mean(4 * Prob * (1 - Prob))) %>%
          dplyr::arrange(desc(Penalty)) %>%
          dplyr::pull(Stability)
      } else {
        stability <- rep(NA, length(self$penalties))
      }
      stability
    }
  )
)


#' R6 generic class for normal-block models
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param C group matrix C_jq = 1 if species j belongs to group q
#' @param penalties list of penalties to be tested
#' @param n_penalties number of penalty values to be tested (if penalties not provided, default is 30)
#' @param min_ratio ratio between min and max penalty to be tested  (if penalties not provided, default is 0.05)
#' @param models uderlying NB_fixed_Q models for each nb of blocks
#' @param verbose telling if information should be printed during optimization
NB_fixed_blocks_diagonal_sparse <- R6::R6Class(
  classname = "NB_fixed_blocks_diagonal_sparse",
  inherit = NB_fixed_blocks_sparse,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(

    #' @description Create a new [`NB_fixed_blocks_diagonal_sparse`] object.
    #' @param Y the matrix of responses (called Y in the model).
    #' @param X design matrix (called X in the model).
    #' @param C group matrix.
    #' @param penalties list of penalties on the network density
    #' @param n_penalties number of penalty values
    #' @param min_ratio ratio between min and max penalty to be tested
    #' @return A new [`NB_fixed_blocks_diagonal_sparse`] object
    initialize = function(Y, X, C, penalties = NULL,  n_penalties = 30,
                          min_ratio = 0.05, verbose=TRUE) {
      super$initialize(Y, X, C, penalties, n_penalties,
                       min_ratio, verbose)
      # instantiates an NB_fixed_blocks_diagonal model for each Q in nb_blocks
      self$models <- map(self$penalties[order(self$penalties)],
                         function(penalty) {
                           model <- NB_fixed_blocks_diagonal$new(self$Y, self$X,
                                                                 self$C,
                                                                 sparsity = penalty)
                         })
    },

    ## Stability -------------------------
    #' @description Compute the stability path by stability selection
    #' @param subsamples a list of vectors describing the subsamples. The number of vectors (or list length) determines the number of subsamples used in the stability selection. Automatically set to 20 subsamples with size `10*sqrt(n)` if `n >= 144` and `0.8*n` otherwise following Liu et al. (2010) recommendations.
    stability_selection = function(subsamples = NULL) {

      ## select default subsamples according to Liu et al. (2010) recommendations.
      if (is.null(subsamples)) {
        subsample.size <- round(ifelse(self$n >= 144, 10*sqrt(self$n), 0.8*private$n))
        subsamples <- replicate(20, sample.int(self$n, subsample.size), simplify = FALSE)
      }

      ## got for stability selection
      cat("\nStability Selection for NB_fixed_blocks_sparse: ")
      cat("\nsubsampling: ")

      stabs_out <- future.apply::future_lapply(subsamples, function(subsample) {
        cat("+")
        inception_ <- self$get_model(self$penalties[[1]])$clone()
        inception_$update(
          mu     = inception_$posterior_par$mu[subsample, ]
        )


        data <- list(
          Y  = self$Y  [subsample, , drop = FALSE],
          X  = self$X [subsample, , drop = FALSE])

        myNB <- NB_fixed_blocks_diagonal_sparse$new(data$Y, data$X, self$C, self$penalties)
        myNB$optimize(niter = self$latest_niter, threshold = self$latest_threshold)

        nets <- do.call(cbind, lapply(myNB$models, function(model) {
          as.matrix(model$latent_network("support"))[upper.tri(diag(self$Q))]
        }))
        nets
      }, future.seed = TRUE, future.scheduling = structure(TRUE, ordering = "random"))

      prob <- Reduce("+", stabs_out, accumulate = FALSE) / length(subsamples)
      ## formatting/tyding
      node_set <- lapply(1:self$Q, f <- function(g){paste0("group_", g)})
      colnames(prob) <- self$penalties
      private$stab_path <- prob %>%
        as.data.frame() %>%
        dplyr::mutate(Edge = 1:dplyr::n()) %>%
        tidyr::gather(key = "Penalty", value = "Prob", -Edge) %>%
        dplyr::mutate(Penalty = as.numeric(Penalty),
                      Node1   = node_set[edge_to_node(Edge)$node1],
                      Node2   = node_set[edge_to_node(Edge)$node2],
                      Edge    = paste0(Node1, "|", Node2))
      invisible(subsamples)
    }
  )
)

#' R6 generic class for normal-block models
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @param C group matrix C_jq = 1 if species j belongs to group q
#' @param penalties list of penalties to be tested
#' @param n_penalties number of penalty values to be tested (if penalties not provided, default is 30)
#' @param min_ratio ratio between min and max penalty to be tested  (if penalties not provided, default is 0.05)
#' @param models uderlying NB_fixed_Q models for each nb of blocks
#' @param verbose telling if information should be printed during optimization
NB_fixed_blocks_spherical_sparse <- R6::R6Class(
  classname = "NB_fixed_blocks_spherical_sparse",
  inherit = NB_fixed_blocks_sparse,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(

    #' @description Create a new [`NB_fixed_blocks_sparse`] object.
    #' @param Y the matrix of responses (called Y in the model).
    #' @param X design matrix (called X in the model).
    #' @param C group matrix.
    #' @param penalties list of penalties on the network density
    #' @param n_penalties number of penalty values
    #' @param min_ratio ratio between min and max penalty to be tested
    #' @return A new [`nb_fixed_blocks_sparse`] object
    initialize = function(Y, X, C, penalties = NULL,  n_penalties = 30,
                          min_ratio = 0.05, verbose=TRUE) {
      super$initialize(Y, X, C, penalties, n_penalties,
                       min_ratio, verbose)
      # instantiates an NB_fixed_blocks_spherical model for each Q in nb_blocks
      self$models <- map(self$penalties[order(self$penalties)],
                         function(penalty) {
                           model <- NB_fixed_blocks_spherical$new(self$Y, self$X,
                                                                 self$C, penalty)
                         })
    },


    ## Stability -------------------------
    #' @description Compute the stability path by stability selection
    #' @param subsamples a list of vectors describing the subsamples. The number of vectors (or list length) determines the number of subsamples used in the stability selection. Automatically set to 20 subsamples with size `10*sqrt(n)` if `n >= 144` and `0.8*n` otherwise following Liu et al. (2010) recommendations.
    stability_selection = function(subsamples = NULL) {

      ## select default subsamples according to Liu et al. (2010) recommendations.
      if (is.null(subsamples)) {
        subsample.size <- round(ifelse(self$n >= 144, 10*sqrt(self$n), 0.8*private$n))
        subsamples <- replicate(20, sample.int(self$n, subsample.size), simplify = FALSE)
      }

      ## got for stability selection
      cat("\nStability Selection for NB_fixed_blocks_sparse: ")
      cat("\nsubsampling: ")

      stabs_out <- future.apply::future_lapply(subsamples, function(subsample) {
        cat("+")
        inception_ <- self$get_model(self$penalties[[1]])$clone()
        inception_$update(
          mu     = inception_$posterior_par$mu[subsample, ]
        )


        data <- list(
          Y  = self$Y  [subsample, , drop = FALSE],
          X  = self$X [subsample, , drop = FALSE])

        myNB <- NB_fixed_blocks_spherical_sparse$new(data$Y, data$X, self$C, self$penalties)
        myNB$optimize(niter = self$latest_niter, threshold = self$latest_threshold)
        nets <- do.call(cbind, lapply(myNB$models, function(model) {
          as.matrix(model$latent_network("support"))[upper.tri(diag(self$Q))]
        }))
        nets
      }, future.seed = TRUE, future.scheduling = structure(TRUE, ordering = "random"))

      prob <- Reduce("+", stabs_out, accumulate = FALSE) / length(subsamples)
      ## formatting/tyding
      node_set <- lapply(1:self$Q, f <- function(g){paste0("group_", g)})
      colnames(prob) <- self$penalties
      private$stab_path <- prob %>%
        as.data.frame() %>%
        dplyr::mutate(Edge = 1:dplyr::n()) %>%
        tidyr::gather(key = "Penalty", value = "Prob", -Edge) %>%
        dplyr::mutate(Penalty = as.numeric(Penalty),
                      Node1   = node_set[edge_to_node(Edge)$node1],
                      Node2   = node_set[edge_to_node(Edge)$node2],
                      Edge    = paste0(Node1, "|", Node2))
      invisible(subsamples)
    }
  )
)
