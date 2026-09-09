## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockVarCollectionClusters #################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Collection of Normal-Block Models over a Range of Cluster Counts
#'
#' R6 class for a collection of normal-block models with different number of
#' clusters (q) and a fixed sparsity level.
#' @examples
#' ex <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' models <- normal_block(data, blocks = 2:4, control = NB_control(verbose = FALSE))
#' models$get_best_model("ICL")$q
#' @export
NormalBlockVarCollectionClusters <- R6::R6Class(
  classname = "NormalBlockVarCollectionClusters",
  inherit   = NormalBlockCollectionClusters,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`NormalBlockVarCollectionClusters`] object.
    #' @param mydata object of NormalBlockData class, with responses and design matrix
    #' @param q_list list of q values (number of groups) in the collection
    #' @param zero_inflation whether the models in the collection should be zero-inflated or not
    #' @param sparsity sparsity penalty on the network density
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return A new [`NormalBlockVarCollectionClusters`] object
    initialize = function(mydata, q_list, zero_inflation = FALSE,
                          sparsity = 0, control = NB_control()) {
      stopifnot("each nb_blocks value can only be present once in nb_blocks" =
                  length(q_list) == length(unique(q_list)))
      stopifnot("There cannot be more blocks than there are entities to cluster." =
                  max(q_list) <= ncol(mydata$Y))

      self$control <- control
      self$control$zero_inflation <- zero_inflation

      ## Whatever part of the requested clustering heuristic doesn't depend on
      ## q is computed once here instead of once per model: a single wide SBM
      ## exploration covers every intermediate block count, one hierarchical
      ## tree is cut at each q, one eigendecomposition serves every spectral
      ## cut. See clustering_path_for_family() in R/utils.R; NULL when nothing is
      ## shareable ("kmeans", "best_of_inits", an explicit clustering), and
      ## every model then runs its own heuristic_clustering() as before.
      clustering_path <- clustering_path_for_family(mydata, q_list, "var", zero_inflation, control)

      # instantiates an NormalBlockVarUnknownClusters model for each q in nb_blocks
      this_control <- control
      self$models <- map(rank(q_list),
          function(r) {
            this_control$clustering_init <-
              if (!is.null(clustering_path)) clustering_path[[as.character(q_list[r])]]
              ## a list means one explicit clustering per q; anything else
              ## (a heuristic name, or NULL) applies identically to every q
              else if (is.list(control$clustering_init)) control$clustering_init[[r]]
              else control$clustering_init
            model <- get_model(mydata,
                               q_list[r],
                               sparsity = sparsity,
                               zero_inflation = zero_inflation,
                               control = this_control)
        })
      private$progress_field <- "q"
      private$progress_label <- "number of blocks"
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field who_am_I a method to print what model is being fitted
    ## self$control$noise_covariance stays NULL unless the caller set it
    ## explicitly: each model resolves its own "diagonal" default locally,
    ## never writing it back to the collection's control.
    who_am_I  = function(){
      nc <- if (is.null(self$control$noise_covariance)) "diagonal" else self$control$noise_covariance
      paste0(nc, " normal-block-var model with unknown q")
    }
  )
)
