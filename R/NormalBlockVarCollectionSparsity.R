## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockVarCollectionSparsity #########################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Collection of Normal-Block Models over a Sparsity Path
#'
#' R6 class for a collection of normal-block models with a fixed clustering
#' (blocks) and different sparsity levels.
#' @examples
#' ex <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' models <- normal_block(data, blocks = ex$parameters$C, sparsity = TRUE,
#'                        control = NB_control(verbose = FALSE, n_sparsity_penalties = 5))
#' models$get_best_model("BIC")$sparsity
#' @export
NormalBlockVarCollectionSparsity <- R6::R6Class(
  classname = "NormalBlockVarCollectionSparsity",
  inherit   = NormalBlockCollectionSparsity,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`NormalBlockVarCollectionSparsity`] object.
    #' @param mydata object of NormalBlockData class, with responses and design matrix
    #' @param blocks either a clustering matrix (known, fixed clustering) or a single integer (number of blocks to infer)
    #' @param zero_inflation boolean to specify whether data is zero-inflated
    #' @param control structured list of parameters to handle sparsity control
    #' @return A new [`NormalBlockVarCollectionSparsity`] object
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
        init_model$optimize(control = list(niter=5, threshold=1e-4, verbose=FALSE), warn = FALSE)
        Sigmaq   <- chol2inv(chol(init_model$model_par$Omega))
        diag_pen <- max(diag(init_model$sparsity_weights)) > 0
        weights  <- init_model$sparsity_weights
        weights  <- abs((Sigmaq / weights)[upper.tri(Sigmaq, diag = diag_pen)])
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
      private$progress_field <- "sparsity"
      private$progress_label <- "penalty"
    },

    #' @description Extract best model in the collection
    #' @param crit a character for the criterion used to performed the selection.
    #' @param stability if criterion = "StARS" gives level of stability required.
    #' Either "BIC", "EBIC", "ICL" or "StARS". Default is BIC
    #' @return a [`NormalBlockVarUnknownClusters`] object
    get_best_model = function(crit = c("BIC", "EBIC", "ICL", "StARS"),
                              stability = 0.9) {
      crit <- match.arg(crit)
      if (crit == "StARS") {
        if (is.null(private$stab_path)) self$stability_selection()
        max_stab <- max(self$criteria$stability)
        if (max_stab < stability) {
          message("No model reaches the required stability ", stability, ", returning model with highest stability: ", max_stab)
          stability <- max_stab
        }
        id_stars <- self$criteria %>%
          dplyr::select(sparsity, stability) %>% dplyr::rename(Stability = stability) %>%
          dplyr::filter(Stability >= stability) %>%
          dplyr::pull(sparsity) %>% min() %>% match(private$sparsity_)
        model <- self$models[[id_stars]]$clone()
      } else {
        id <- private$best_id(crit)
        model <- self$models[[id]]$clone()
      }
      model
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
        cat("\nStability Selection for NormalBlockVarKnownClusters (sparse): \nsubsampling: ")

      stabs_out <- lapply(subsamples, function(subsample) {
        mydata <- NormalBlockData$new(Y = self$data$Y[subsample, , drop = FALSE],
                                  X = self$data$X[subsample, , drop = FALSE])
        myNB <- NormalBlockVarCollectionSparsity$new(
          mydata, blocks, self$control$zero_inflation, control_stabs)
        myNB$optimize(control_stabs)

        upper_tri <- upper.tri(diag(self$q))
        nets <- do.call(cbind, lapply(myNB$models, function(model) {
          model$latent_network("support")[upper_tri]
        }))
        nets
      })

      prob <- Reduce("+", stabs_out, accumulate = FALSE) / length(subsamples)
      ## formatting/tyding
      node_set <- lapply(1:self$q, f <- function(g){paste0("group_", g)})
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
    stab_path = NULL # stability path, filled in by stability_selection()
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field sparsity_details list of information about model's penalties
    sparsity_details = function()
      list("n_penalties" = length(self$sparsity), "min_ratio" = min(self$sparsity)/max(self$sparsity),
           "min_sparsity" = min(self$sparsity), "max_sparsity" = max(self$sparsity)),
    #' @field criteria a data frame with the values of some criteria ((approximated) log-likelihood, BIC) for the collection of models
    criteria = function() super$criteria %>% dplyr::mutate(stability = self$stability),
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
    who_am_I  = function(){
      nc <- if (is.null(self$control$noise_covariance)) "diagonal" else self$control$noise_covariance
      zi <- if (isTRUE(self$control$zero_inflation)) "zero-inflated " else ""
      paste0("Collection of ", zi, nc, " normal-block-var models with ",
             ifelse(is.matrix(private$blocks_), "fixed blocks", "fixed q"),
        ", with different sparsity penalties.")}
  )
)
