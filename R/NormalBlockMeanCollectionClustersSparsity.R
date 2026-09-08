## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockMeanCollectionClustersSparsity ##############
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Collection of Mean-Block Models over Cluster Counts and Sparsity Levels
#'
#' R6 class for a collection of mean-block models ([NormalBlockMeanBase]) over
#' both a range of cluster counts (q) and a sparsity path, i.e. one
#' [NormalBlockMeanCollectionSparsity] per q. Mirrors
#' [NormalBlockVarCollectionClustersSparsity], minus its SBM-path shortcut for
#' the initial clustering (ill-suited to this family, see
#' [NormalBlockMeanCollectionClusters]).
#' @examples
#' ex <- generate_normal_block_mean_data(n = 60, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' models <- normal_block(data, blocks = 2:4, sparsity = TRUE, model = "mean",
#'                        control = NB_control(n_sparsity_penalties = 4))
#' models$plot("BIC")
#' models$get_best_model("BIC")
#' @export
NormalBlockMeanCollectionClustersSparsity <- R6::R6Class(
  classname = "NormalBlockMeanCollectionClustersSparsity",
  inherit   = NormalBlockCollectionClustersSparsity,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`NormalBlockMeanCollectionClustersSparsity`] object.
    #' @param mydata object of NormalBlockData class, with responses and design matrix
    #' @param q_list list of q values (number of groups) in the collection
    #' @param control structured list of parameters to handle sparsity control
    #' @return A new [`NormalBlockMeanCollectionClustersSparsity`] object
    initialize = function(mydata, q_list, control = NB_control()) {
      self$control <- control
      control_ <- control

      self$models <- map(seq_along(q_list), function(rank) {
        ## a list means one explicit clustering per q; anything else applies identically to every q
        if (is.list(control$clustering_init))
          control_$clustering_init <- control$clustering_init[[rank]]
        NormalBlockMeanCollectionSparsity$new(mydata, q_list[rank], control_)
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
    who_am_I = function()
      "collection of normal-block-mean models with different values of q and different penalties"
  )
)
