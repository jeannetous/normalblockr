## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_unknown_Q #################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 class for normal-block model with unknown Q (number of groups)
#' @param data contains the matrix of responses (Y) and the design matrix (X).
#' @param zero_inflation whether the models should be zero-inflated or not
#' @param control structured list for specific parameters (including initial clustering proposal)
NB_unknown_Q <- R6::R6Class(
  classname = "NB_unknown_Q",

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field models list of NB_fixed_Q models corresponding to each nb_block value
    models = NULL,
    #' @field control store the list of user-defined model settings and optimization parameters
    control = NA,

    #' @description Create a new [`NB_unknown_Q`] object.
    #' @param mydata object of normal_data class, with responses and design matrix
    #' @param Q_list list of Q values (number of groups) in the collection
    #' @param zero_inflation whether the models in the collection should be zero-inflated or not
    #' @param penalty penalty on the network density
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return A new [`NB_unknown_Q`] object
    initialize = function(mydata, Q_list, zero_inflation = FALSE,
                          penalty = 0, control = NB_control()) {
      stopifnot("each nb_blocks value can only be present once in nb_blocks" =
                  length(Q_list) == length(unique(Q_list)))
      stopifnot("There cannot be more blocks than there are entities to cluster." =
                  max(Q_list) <= ncol(mydata$Y))

      self$control <- control
      self$control$zero_inflation <- zero_inflation

      # instantiates an NB_fixed_Q model for each Q in nb_blocks
      this_control <- control
      self$models <- map(rank(Q_list),
          function(r) {
            this_control$clustering_init <- control$clustering_init[[r]]
            model <- get_model(mydata,
                               Q_list[r],
                               sparsity = penalty,
                               zero_inflation = zero_inflation,
                               control = this_control)
        })
    },

    #' @description optimizes an NB_fixed_Q object for each value of Q
    #' @param control optimization parameters (niter and threshold)
    optimize = function(control = list(niter = 100, threshold = 1e-4, verbose=TRUE)) {
      self$models <- map(self$models, function(model) {
        if(control$verbose) cat("\tnumber of blocks =", model$Q, "          \r")
        flush.console()
        model$optimize(control)
        model
      })
    },

    #' @description returns the NB_fixed_Q model corresponding to given Q
    #' @param Q number of blocks asked by user
    #' @return A NB_fixed_Q object with given value Q
    get_model = function(Q) {
      stopifnot("No such model in the collection. Acceptable values can be found via $Q" = Q %in% self$Q_list)
      self$models[[which(self$Q_list == Q)]]
    },

    #' @description Extract best model in the collection
    #' @param crit a character for the criterion used to performed the selection.
    #' Either "ICL" or "BIC". "ICL" is the default criterion
    #' @return a [`NB_fixed_Q`] object
    get_best_model = function(crit = c("ICL", "BIC")) {
      stopifnot("Log-likelihood based criteria do not apply to the heuristic method" = self$models[[1]]$inference_method == "integrated")
      crit <- match.arg(crit)
      stopifnot(!anyNA(self$criteria[[crit]]))
      id <- 1
      if (length(self$criteria[[crit]]) > 1) {
        id <- which.min(self$criteria[[crit]])
      }
      model <- self$models[[id]]$clone()
      model
    },

    #' @description Display various outputs (goodness-of-fit criteria, robustness, diagnostic) associated with a collection of network fits (a [`Networkfamily`])
    #' @param criteria vector of characters. The criteria to plot in `c("deviance", "BIC", "ICL")`. Defaults to all of them.
    #' @return a [`ggplot2::ggplot`] graph
    plot = function(criteria = c("deviance", "ICL", "BIC", "EBIC")) {
      vlines <- sapply(intersect(criteria, c("ICL")) , function(crit) self$get_best_model(crit)$Q)
      stopifnot(!anyNA(self$criteria[criteria]))

      dplot <- self$criteria %>%
        dplyr::select(dplyr::all_of(c("Q", criteria))) %>%
        tidyr::gather(key = "criterion", value = "value", -Q) %>%
        dplyr::group_by(criterion)
      p <- ggplot2::ggplot(dplot, ggplot2::aes(x = Q, y = value, group = criterion, colour = criterion)) +
        ggplot2::geom_line() + ggplot2::geom_point() +
        ggplot2::ggtitle(label    = "Model selection criteria",
                         subtitle = "Lower is better" ) +
        ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = vlines, linetype = "dashed", alpha = 0.25)
      p
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field Q_list number of blocks
    Q_list = function(value) map_dbl(self$models, "Q"),
    #' @field criteria a data frame with the values of some criteria ((approximated) log-likelihood, BIC, ICL) for the collection of models
    criteria = function() purrr::map_df(self$models, "criteria"),
    #' @field who_am_I a method to print what model is being fitted
    who_am_I  = function(value){paste0(self$control$noise_covariance, " normal-block model with unknown Q")}
  )
)
