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
    #' @field models list of NB_fixed_Q models corresponding to each nb_block value
    models = NA,
    #' @field control store the list of user-defined model settings and optimization parameters
    control = NA,
    #' @field data object of NBData class, with responses and design matrix
    data   = NA,

    #' @description Create a new [`NB_changing_sparsity`] object.
    #' @param mydata object of NBData class, with responses and design matrix
    #' @param zero_inflation boolean to specify whether data is zero-inflated
    #' @param control structured list of parameters to handle sparsity control
    #' @return A new [`NB_changing_sparsity`] object
    initialize = function(mydata, blocks, zero_inflation = FALSE,
                          control = NB_control()) {

      ## Store user-defined fields
      self$data    <- mydata
      self$control <- control
      self$control$zero_inflation <- zero_inflation
      private$blocks_ <- blocks

      ## Check block format
      stopifnot("blocks must be either a clustering matrix or a fixed number of blocks" =
                  is.matrix(blocks) | length(blocks) == 1)

      ## extract a collection of sparsifying penalties
      if (!is.null(control$sparsity_penalties)){
        stopifnot("All penalties must be strictly positive" =
                    (min(control$sparsity_penalties) > 0))
        sparsity <- control$sparsity_penalties
      } else {
        init_model <- get_model(mydata, blocks, 0, zero_inflation, control)
        init_model$optimize(control = list(niter=5, threshold=1e-4, verbose=FALSE))
        SigmaQ   <- solve(init_model$model_par$OmegaQ)
        diag_pen <- max(diag(init_model$sparsity_weights)) > 0
        weights  <- init_model$sparsity_weights
        weights  <- abs((SigmaQ / weights)[upper.tri(SigmaQ, diag = diag_pen)])
        if (length(weights) > 0) {
          max_pen  <- max(weights)
          sparsity <- 10^seq(log10(max_pen), log10(max_pen * control$min_ratio), len = control$n_sparsity_penalties)
        } else {
          sparsity <- rep(0, len = control$n_sparsity_penalties)
        }
      }

      private$sparsity_ <- sort(sparsity, decreasing = TRUE)
      ## Instantiation of the models in the collection
      self$models <- map(private$sparsity_, function(lambda)
        model <- get_model(mydata, blocks, lambda, zero_inflation, control)
      )
    },

    #' @description optimizes a model for each penalty
    #' @param control optimization parameters (niter and threshold)
    optimize = function(control = list(niter=100, threshold=1e-4, verbose=TRUE)) {
      self$models <- lapply(self$models, function(model) {
        if (control$verbose) cat("\t penalty =", model$sparsity, "          \r")
        flush.console()
        model$optimize(control)
        model
      })
    },

    #' @description returns the NB_fixed_block model corresponding to given penalty
    #' @param sparsity sparsity penalty asked by user
    #' @return A NB_fixed_blocks_sparse object with given value penalty
    get_model = function(sparsity) {
      if (!(sparsity %in% private$sparsity_)) {
        sparsity <-  private$sparsity_[[which.min(abs(private$sparsity_ - sparsity))]]
        cat(paste0("No model with this penalty in the collection. Returning model with closest penalty: ",
                   sparsity,  " Collection penalty values can be found via $sparsity \n"))
      }
      self$models[[which(private$sparsity_ == sparsity)]]
    },

    #' @description Extract best model in the collection
    #' @param crit a character for the criterion used to performed the selection.
    #' @param stability if criterion = "StARS" gives level of stability required.
    #' Either "BIC", "EBIC", "ICL" or "StARS". Default is BIC
    #' @return a [`NB_fixed_Q`] object
    get_best_model = function(crit = c("BIC", "EBIC", "ICL", "StARS"),
                              stability = 0.9) {
      stopifnot("Log-likelihood based criteria do not apply to the heuristic method" = (self$models[[1]]$inference_method == "integrated" | crit == "StARS" ))
      crit <- match.arg(crit)
      if (crit == "StARS") {
        if (is.null(private$stab_path)) self$stability_selection()
        max_stab <- max(self$criteria$stability)
        if (max_stab < stability) {
          cat(paste0("No model reaches the required stability ", stability, ", returning model with highest stability: ", max_stab))
          stability <- max_stab
        }
        id_stars <- self$criteria %>%
          dplyr::select(sparsity, stability) %>% dplyr::rename(Stability = stability) %>%
          dplyr::filter(Stability >= stability) %>%
          dplyr::pull(sparsity) %>% min() %>% match(private$sparsity_)
        model <- self$models[[id_stars]]$clone()
      } else {
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
    #' @param criteria vector of characters. The criteria to plot in `c("deviance", BIC", "EBIC", "ICL")`. Defaults to all of them.
    #' @param log.x logical: should the x-axis be represented in log-scale? Default is `TRUE`.
    #' @importFrom tidyr gather
    #' @return a [`ggplot`] graph
    plot = function(criteria = c("deviance", "BIC", "EBIC", "ICL"), log.x = TRUE) {
      vlines <- sapply(intersect(criteria, c("BIC")) , function(crit) self$get_best_model(crit)$sparsity)
      stopifnot(!is.null(self$criteria[criteria]))

      dplot <- self$criteria %>%
        dplyr::select(dplyr::all_of(c("sparsity", criteria))) %>%
        tidyr::gather(key = "criterion", value = "value", -sparsity) %>%
        dplyr::group_by(criterion)
      p <- ggplot2::ggplot(dplot, ggplot2::aes(x = sparsity, y = value, group = criterion, colour = criterion)) +
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
        subsample.size <- round(ifelse(self$data$n >= 144, 10*sqrt(self$data$n), 0.8*self$data$n))
        subsamples <- replicate(n_subsamples, sample.int(self$data$n, subsample.size), simplify = FALSE)
      }

      ## retrieve all the appropriate control parameters for calling stabselection
      control_stabs <- self$control
      control_stabs$sparsity <- private$sparsity_
      control_stabs$verbose <- FALSE
      control_stabs$fixed_tau <- TRUE
      blocks  <- private$blocks_
      if (is.matrix(blocks)){
        control_stabs$clustering_init <- blocks
      } else {
        control_stabs$clustering_init <- self$get_model(min(private$sparsity_))$var_par$tau
      }

      ## got for stability selection
      if(self$control$verbose)
        cat("\nStability Selection for NB_fixed_blocks_sparse: \nsubsampling: ")

      stabs_out <- lapply(subsamples, function(subsample) {
        mydata <- NBData$new(Y = self$data$Y[subsample, , drop = FALSE],
                                  X = self$data$X[subsample, , drop = FALSE])
        myNB <- NB_changing_sparsity$new(
          mydata, blocks, self$control$zero_inflation, control_stabs)
        myNB$optimize(control_stabs)

        upper_tri <- upper.tri(diag(self$Q))
        nets <- do.call(cbind, lapply(myNB$models, function(model) {
          model$latent_network("support")[upper_tri]
        }))
        nets
      })

      prob <- Reduce("+", stabs_out, accumulate = FALSE) / length(subsamples)
      ## formatting/tyding
      node_set <- lapply(1:self$Q, f <- function(g){paste0("group_", g)})
      colnames(prob) <- private$sparsity_
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
    stab_path    = NULL, # a field to store the stability path,
    sparsity_    = NA,   # penalty values in the collection for sparsity (lambda)
    blocks_      = NA    # blocks (either a scalar or an indicator matrix)
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field Q number of blocks
    Q = function(value) ifelse(is.matrix(private$blocks_), ncol(private$blocks_), private$blocks_),
    #' @field blocks group matrix or number of blocks.
    blocks = function(value) private$blocks_,
    #' @field sparsity list of sparsity penalties
    sparsity = function() private$sparsity_,
    #' @field sparsity_details list of information about model's penalties
    sparsity_details = function()
      list("n_penalties" = length(self$sparsity), "min_ratio" = min(self$sparsity)/max(self$sparsity),
           "min_sparsity" = min(self$sparsity), "max_sparsity" = max(self$sparsity)),
    #' @field criteria a data frame with the values of some criteria ((approximated) log-likelihood, BIC) for the collection of models
    criteria = function() map_df(self$models, "criteria") %>% mutate(stability = self$stability),
    #' @field stability_path measure of edges stability based on StARS method
    stability_path = function() private$stab_path,
    #' @field stability mean edge stability along the sparsity penalties path
    stability = function() {
      if (!is.null(private$stab_path)) {
        stability <- self$stability_path %>%
          dplyr::select(Penalty, Prob) %>%
          dplyr::group_by(Penalty) %>%
          dplyr::summarize(Stability = 1 - mean(4 * Prob * (1 - Prob))) %>%
          dplyr::arrange(desc(Penalty)) %>%
          dplyr::pull(Stability)
      } else {
        stability <- rep(NA, length(private$sparsity_))
      }
      stability
    },
    #' @field who_am_I a method to print what model is being fitted
    who_am_I  = function(value){
      paste0("Collection of ",
             ifelse(self$control$zero_inflation, " zero-inflated ", ""),
                    self$control$noise_covariance, " normal-block models with ",
             ifelse(is.matrix(private$blocks_), "fixed blocks", "fixed Q"),
        ", with different sparsity penalties.")}
  )
)
