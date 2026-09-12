## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockCollectionSparsity ############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Base Class for a Collection of Models over a Sparsity Path
#'
#' Shared scaffolding for [NormalBlockVarCollectionSparsity] and
#' [NormalBlockMeanCollectionSparsity]: the warm-started path traversal, the
#' penalty lookup and the criteria plot. Concrete subclasses derive the
#' penalty grid in their `initialize()` and provide their own
#' `get_best_model()` (only the variance-block family offers StARS).
#' @examples
#' # An internal abstract base class, never instantiated directly -- see
#' # normal_block() for how collections are created and fitted.
#' @keywords internal
NormalBlockCollectionSparsity <- R6::R6Class(
  classname = "NormalBlockCollectionSparsity",
  inherit   = NormalBlockCollection,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field data object of NormalBlockData class, with responses and design matrix
    data = NA,

    #' @description optimizes every model in the sparsity path, warm-starting
    #' each one (after the first) from the previous, adjacent penalty's
    #' converged parameters (see the family's base class's `warm_start_from()`)
    #' instead of re-deriving everything from the heuristic clustering, the
    #' way the generic [NormalBlockCollection] `optimize()` would. `blocks`
    #' (hence q) is fixed across the whole path, only the sparsity penalty
    #' changes, so the warm start is always between models of matching shape.
    #' @param control optimization parameters (niter and threshold)
    optimize = function(control = list(niter = 500, threshold = 1e-4, verbose = TRUE)) {
      previous <- NULL
      self$models <- lapply(self$models, function(model) {
        if (control$verbose)
          cat("\t", private$progress_label, "=", model[[private$progress_field]], "          \r")
        flush.console()
        if (!is.null(previous)) model$warm_start_from(previous)
        model$optimize(control)
        previous <<- model
        model
      })
      invisible(self)
    },

    #' @description returns the NormalBlockVarKnownClusters model corresponding to given penalty
    #' @param sparsity sparsity penalty asked by user
    #' @return A NormalBlockVarKnownClusters (sparse) object with given value penalty
    get_model = function(sparsity) {
      if (!(sparsity %in% private$sparsity_)) {
        sparsity <-  private$sparsity_[[which.min(abs(private$sparsity_ - sparsity))]]
        message("No model with this penalty in the collection. Returning model with closest penalty: ",
                sparsity,  " Collection penalty values can be found via $sparsity")
      }
      self$models[[which(private$sparsity_ == sparsity)]]
    },

    #' @description Display various outputs (goodness-of-fit criteria, robustness, diagnostic) associated with a collection of network fits (a collection)
    #' @param criteria vector of characters. The criteria to plot in `c("deviance", BIC", "EBIC", "ICL")`. Defaults to all of them.
    #' @param log.x logical: should the x-axis be represented in log-scale? Default is `TRUE`.
    #' @return a [`ggplot2::ggplot`] graph
    plot = function(criteria = c("deviance", "BIC", "EBIC", "ICL"), log.x = TRUE) {
      stopifnot(!is.null(self$criteria[criteria]))
      p <- private$plot_criteria_path("sparsity", criteria)
      if (log.x) p <- p + ggplot2::coord_transform(x = "log10")
      p
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    sparsity_ = NA, # penalty values in the collection (lambda)
    blocks_   = NA  # blocks (either a scalar or an indicator matrix)
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field q number of blocks
    q = function() ifelse(is.matrix(private$blocks_), ncol(private$blocks_), private$blocks_),
    #' @field blocks group matrix or number of blocks
    blocks = function() private$blocks_,
    #' @field sparsity list of sparsity penalties
    sparsity = function() private$sparsity_
  )
)
