## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockMeanCollectionClusters ############
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Collection of Mean-Block Models over a Range of Cluster Counts
#'
#' R6 class for a collection of mean-block models ([NormalBlockMeanBase])
#' with different numbers of clusters (q). Inherits its scaffolding
#' (`print()`/`summary()`/`plot()`/`optimize()`, the `criteria` table) from
#' [NormalBlockCollection]. Unlike [NormalBlockVarCollectionClusters], there
#' is no SBM-path shortcut here: the shared clustering-heuristic registry's
#' cov()/correlation-based methods are ill-suited to the mean-block family
#' (see [NormalBlockMeanBase]'s own default), so each q is fit
#' independently.
#' @examples
#' ex <- generate_normal_block_mean_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' models <- normal_block(data, blocks = 2:5, model = "mean")
#' models$plot(c("BIC", "ICL"))
#' models$get_best_model()
#' @export
NormalBlockMeanCollectionClusters <- R6::R6Class(
  classname = "NormalBlockMeanCollectionClusters",
  inherit   = NormalBlockCollectionClusters,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`NormalBlockMeanCollectionClusters`] object.
    #' @param mydata object of NormalBlockData class, with responses and design matrix
    #' @param q_list list of q values (number of groups) in the collection
    #' @param zero_inflation whether Y carries structural zeros; every model in
    #' the collection is then zero-inflated (Sigma diagonal or spherical only,
    #' see [NormalBlockMeanBase])
    #' @param sparsity sparsity penalty on the network density
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return A new [`NormalBlockMeanCollectionClusters`] object
    initialize = function(mydata, q_list, zero_inflation = FALSE, sparsity = 0,
                          control = NB_control()) {
      stopifnot("each nb_blocks value can only be present once in nb_blocks" =
                  length(q_list) == length(unique(q_list)))
      stopifnot("There cannot be more blocks than there are entities to cluster." =
                  max(q_list) <= ncol(mydata$Y))

      self$control <- control

      ## Whatever part of the requested heuristic doesn't depend on q, computed
      ## once instead of once per model. For this family's default (kmeans)
      ## that is the lossless row compression of the mean trajectory, whose
      ## rank is only d. See clustering_path_for_family() in R/utils.R.
      clustering_path <- clustering_path_for_family(mydata, q_list, "mean",
                                                    zero_inflation, control)

      this_control <- control
      self$models <- map(rank(q_list),
          function(r) {
            ## a list means one explicit clustering per q; anything else
            ## applies identically to every q
            this_control$clustering_init <-
              if (!is.null(clustering_path)) clustering_path[[as.character(q_list[r])]]
              else if (is.list(control$clustering_init)) control$clustering_init[[r]]
              else control$clustering_init
            get_model(mydata, q_list[r], sparsity = sparsity,
                     zero_inflation = zero_inflation,
                     model = "mean", control = this_control)
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
    who_am_I = function() "normal-block-mean model with unknown q"
  )
)
