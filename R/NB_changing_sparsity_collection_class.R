## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_changing_sparsity #########################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 class for normal-block model with unknown Q (number of groups)
#' @param data contains the matrix of responses (Y) and the design matrix (X).
#' @param zero_inflation whether the models should be zero-inflated or not
#' @param blocks group matrix or number of blocks.
#' @param control structured list for specific parameters (including initial clustering proposal)
NB_changing_sparsity <- R6::R6Class(
  classname = "NB_changing_sparsity",
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field data object of normal_data class, with responses and design matrix
    data  = NULL,
    #' @field blocks group matrix or number of blocks.
    blocks  = NULL,
    #' @field zero_inflation whether the object is a collection of zero-inflated models or not
    zero_inflation = NULL,
    #' @field models list of NB_fixed_Q models corresponding to each nb_block value
    models = NULL,
    #' @field verbose say whether information should be given about the optimization
    verbose = NULL,
    #' @field latest_niter latest niter value used for optimization
    latest_niter = NULL,
    #' @field latest_threshold latest threshold value used for optimization
    latest_threshold = NULL,

    #' @description Create a new [`NB_changing_sparsity`] object.
    #' @param data object of normal_data class, with responses and design matrix
    #' @param blocks group matrix or number of blocks.
    #' @param zero_inflation boolean to specify whether data is zero-inflated
    #' @param control structured list of parameters to handle sparsity control
    #' @return A new [`NB_changing_sparsity`] object
    initialize = function(data, blocks, zero_inflation = FALSE,
                          control = NB_control()) {
      stopifnot("blocks must be either a clustering matrix or a fixed number of blocks" =
                  is.matrix(blocks) | length(blocks) == 1)
      if (!is.null(control$penalties)) {
        if (min(control$penalties) <= 0) stop("All penalties must be strictly positive")
      }
      if(!is.matrix(blocks)){
        if(blocks == 1){
          stop("No penalty can be applied with one cluster as there is no network.")
        }
      }

      private$penalties <- c(0)
      self$data   <- data
      self$blocks <- blocks
      self$zero_inflation   <- zero_inflation
      self$verbose <- control$verbose

      private$penalties        <- control$penalties
      private$n_penalties      <- control$n_penalties
      private$min_ratio        <- control$min_ratio
      private$approx           <- control$heuristic

      if(!is.null(control$penalties)){
        private$n_penalties <- length(control$penalties)
      }else{
          init_model <- get_model(data, blocks,
                                  zero_inflation = self$zero_inflation,
                                  control = control)
          init_model$optimize(control = list(niter = 5, threshold = 1e-4))
          SigmaQ    <- solve(init_model$model_par$OmegaQ)
          diag_pen  <- max(diag(init_model$penalty_weights)) > 0
          weights   <- init_model$penalty_weights ; weights[weights == 0] <- 1
          max_pen   <- max(abs((SigmaQ / weights)[upper.tri(SigmaQ, diag = diag_pen)]))
          private$penalties <- 10^seq(log10(max_pen), log10(max_pen * private$min_ratio), len = private$n_penalties)
      }
      private$penalties <- sort(private$penalties, decreasing = TRUE)
      self$models <- map(private$penalties,
                         function(lambda) {
                           model <- get_model(data, blocks,
                                              sparsity = lambda,
                                              zero_inflation = zero_inflation,
                                              control = control)
                         })
    },

    #' @description optimizes a model for each penalty
    #' @param control optimization parameters (niter and threshold)
    optimize = function(control = list(niter = 100, threshold = 1e-4)) {
      self$latest_niter     <- control$niter
      self$latest_threshold <- control$threshold

      self$models <- lapply(seq_along(self$models), function(m) {
        model <- self$models[[m]]
        if(self$verbose) cat("\t penalty =", self$models[[m]]$penalty, "          \r")
        flush.console()
        model$optimize(control)
        model
      })
    },

    #' @description returns the NB_fixed_block model corresponding to given penalty
    #' @param penalty penalty asked by user
    #' @return A NB_fixed_blocks_sparse object with given value penalty
    get_model = function(penalty) {
      if(!(penalty %in% private$penalties)) {
        penalty <-  private$penalties[[which.min(abs(private$penalties - penalty))]]
        cat(paste0("No model with this penalty in the collection. Returning model with closest penalty : ", penalty,  " Collection penalty values can be found via $penalties_list \n"))
      }
      penalty_rank <- which(sort(private$penalties, decreasing = TRUE) == penalty)
      return(self$models[[penalty_rank]])
    },

    #' @description Extract best model in the collection
    #' @param crit a character for the criterion used to performed the selection.
    #' @param stability if criterion = "StARS" gives level of stability required.
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
          dplyr::select(penalty, stability) %>% dplyr::rename(Stability = stability) %>%
          dplyr::filter(Stability >= stability) %>%
          dplyr::pull(penalty) %>% min() %>% match(private$penalties)
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
    #' @param criteria vector of characters. The criteria to plot in `c("deviance", BIC", "AIC", "ICL")`. Defaults to all of them.
    #' @param log.x logical: should the x-axis be represented in log-scale? Default is `TRUE`.
    #' @importFrom tidyr gather
    #' @return a [`ggplot`] graph
    plot = function(criteria = c("deviance", "BIC", "AIC", "ICL"), log.x = TRUE) {
      vlines <- sapply(intersect(criteria, c("BIC")) , function(crit) self$get_best_model(crit)$penalty)
      stopifnot(!anyNA(self$criteria[criteria]))

      dplot <- self$criteria %>%
        dplyr::select(dplyr::all_of(c("penalty", criteria))) %>%
        tidyr::gather(key = "criterion", value = "value", -penalty) %>%
        dplyr::group_by(criterion)
      p <- ggplot2::ggplot(dplot, ggplot2::aes(x = penalty, y = value, group = criterion, colour = criterion)) +
        ggplot2::geom_line() + ggplot2::geom_point() +
        ggplot2::ggtitle(label    = "Model selection criteria",
                         subtitle = "Lower is better" ) +
        ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = vlines, linetype = "dashed", alpha = 0.25)
      if (log.x) p <- p + ggplot2::coord_trans(x = "log10")
      p
    },

    ## Stability -------------------------
    #' @description Compute the stability path by stability selection
    #' @param subsamples a list of vectors describing the subsamples. The number of vectors (or list length) determines the number of subsamples used in the stability selection. Automatically set to 20 subsamples with size `10*sqrt(n)` if `n >= 144` and `0.8*n` otherwise following Liu et al. (2010) recommendations.
    #' @param n_subsamples number of subsamples to create if the subsamples are not given
    stability_selection = function(subsamples = NULL, n_subsamples = 10) {

      ## select default subsamples according to Liu et al. (2010) recommendations.
      if (is.null(subsamples)) {
        subsample.size <- round(ifelse(self$n >= 144, 10*sqrt(self$n), 0.8*self$n))
        subsamples <- replicate(n_subsamples, sample.int(self$n, subsample.size), simplify = FALSE)
      }

      ## got for stability selection
      if(self$verbose){
        cat("\nStability Selection for NB_fixed_blocks_sparse: ")
        cat("\nsubsampling: ")
      }

      stabs_out <- lapply(subsamples, function(subsample) {
        data <- normal_data$new(Y  = self$data$Y  [subsample, , drop = FALSE],
                                X  = self$data$X  [subsample, , drop = FALSE])

        if(is.matrix(self$blocks)){
          blocks          <- self$blocks
          clustering_init <- NULL
        }else{
          blocks          <- self$Q
          clustering_init <- self$get_model(min(private$penalties))$var_par$tau
        }
        myNB <- get_model(data, blocks = blocks,
                          sparsity = T,
                          zero_inflation = self$zero_inflation,
                          control = NB_control(sparsity_weights = self$penalty_weights,
                                               penalties = private$penalties,
                                               fixed_tau = TRUE,
                                               clustering_init = clustering_init,
                                               noise_covariance = self$get_res_covariance,
                                               verbose = FALSE))
        if(is.null(self$latest_niter)){
          stop("The model must be optimized before running a stability selection.")
        }

        myNB$optimize(control = list(niter = self$latest_niter, threshold = self$latest_threshold))

        nets <- do.call(cbind, lapply(myNB$models, function(model) {
          as.matrix(model$latent_network("support"))[upper.tri(diag(self$Q))]
        }))
        nets
      })

      prob <- Reduce("+", stabs_out, accumulate = FALSE) / length(subsamples)
      ## formatting/tyding
      node_set <- lapply(1:self$Q, f <- function(g){paste0("group_", g)})
      colnames(prob) <- private$penalties
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
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    stab_path            = NULL, # a field to store the stability path,
    penalties            = NA, # penalty values in the collection for sparsity (lambda)
    n_penalties          = NA, # number of distinct penalty values in the collection
    min_ratio            = NA, # ratio between min and max penalties when their values are not explicitly given as input
    approx               = NA # use approximation/heuristic approach or not
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field n number of samples
    n = function() self$data$n,
    #' @field p number of responses per sample
    p = function() self$data$p,
    #' @field d number of variables (dimensions in X)
    d = function() self$data$d,
    #' @field Q number of blocks
    Q = function(value) ifelse(is.matrix(self$blocks), ncol(self$blocks), self$blocks),
    #' @field inference_method inference procedure used (heuristic or integrated with EM)
    inference_method = function(value) ifelse(private$approx, "heuristical", "integrated"),
    #' @field penalties_list list of lambda values for sparsity penalties
    penalties_list = function() private$penalties,
    #' @field penalties_details list of information about model's penalties
    penalties_details = function() list("n_penalties" = private$n_penalties, "min_ratio" = private$min_ratio,
                                      "min_penalty" = min(private$penalties), "max_penalty" = max(private$penalties)),
    #' @field penalty_weights (weights associated to each pair of groups)
    penalty_weights = function(value) self$models[[1]]$penalty_weights,
    #' @field get_res_covariance whether the residual covariance is diagonal or spherical
    get_res_covariance = function(value) self$models[[1]]$get_res_covariance,
    #' @field criteria a data frame with the values of some criteria ((approximated) log-likelihood, BIC, AIC) for the collection of models
    criteria = function(){
      crit <- purrr::map(self$models, "criteria") %>% purrr::reduce(rbind)
      crit$stability <- self$stability
      crit},
    #' @field stability_path measure of edges stability based on StARS method
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
        stability <- rep(NA, length(private$penalties))
      }
      stability
    },
    #' @field who_am_I a method to print what model is being fitted
    who_am_I  = function(value){paste0("Collection of ", ifelse(self$zero_inflation, " zero-inflated ", ""),
                                       self$noise_cov, "normal-block models with ",
                                       ifelse(is.matrix(self$blocks), "fixed blocks", "fixed Q"), ", with different sparsity penalties.")}
  )
)
