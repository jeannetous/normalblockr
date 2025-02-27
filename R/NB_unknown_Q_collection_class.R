## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_unknown_Q #################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 class for normal-block model with unknown Q (number of groups)
#' @param data contains the matrix of responses (Y) and the design matrix (X).
#' @param zero_inflation whether the models should be zero-inflated or not
#' @param verbose telling if information should be printed during optimization
#' @param control structured list for specific parameters (including initial clustering proposal)
NB_unknown_Q <- R6::R6Class(
  classname = "NB_unknown_Q",

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field data object of normal_data class, with responses and design matrix
    data  = NULL,
    #' @field zero_inflation whether the object is a collection of zero-inflated models or not
    zero_inflation = NULL,
    #' @field models list of NB_fixed_Q models corresponding to each nb_block value
    models = NULL,
    #' @field verbose say whether information should be given about the optimization
    verbose = NULL,

    #' @description Create a new [`NB_unknown_Q`] object.
    #' @param data object of normal_data class, with responses and design matrix
    #' @param Q_list list of Q values (number of groups) in the collection
    #' @param zero_inflation whether the models in the collection should be zero-inflated or not
    #' @param penalty penalty on the network density
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return A new [`NB_unknown_Q`] object
    initialize = function(data, Q_list, zero_inflation = FALSE,
                          penalty = 0, control = NB_control()) {
      stopifnot("each nb_blocks value can only be present once in nb_blocks" =
                  length(Q_list) == length(unique(Q_list)))
      stopifnot("each nb_blocks value can only be present once in nb_blocks" =
                  length(Q_list) == length(unique(Q_list)))
      stopifnot("There cannot be more blocks than there are entities to cluster." =
                  max(Q_list) <= ncol(data$Y))
      self$data      <- data
      self$verbose   <- control$verbose
      private$Q_list <- Q_list

      # instantiates an NB_fixed_Q model for each Q in nb_blocks
      self$models <- map(order(Q_list),
          function(Q_rank) {
            this_control <- control
            if(length(this_control$clustering_init) > 1){
              this_control$clustering_init <- control$clustering_init[[Q_rank]]
            }
            model <- get_model(data,
                               Q_list[[Q_rank]],
                               sparsity = penalty,
                               zero_inflation = zero_inflation,
                               control = this_control)
        })
    },

    #' @description optimizes an NB_fixed_Q object for each value of Q
    #' @param niter number of iterations in model optimization
    #' @param threshold loglikelihood threshold under which optimization stops
    optimize = function(control = list(niter = 100, threshold = 1e-4)) {
      self$models <- map(self$models, function(model) {
        if(self$verbose) cat("\tnumber of blocks =", model$Q, "          \r")
        flush.console()
        model$optimize(control)
        model
      })
    },

    #' @description returns the NB_fixed_Q model corresponding to given Q
    #' @param Q number of blocks asked by user
    #' @return A NB_fixed_Q object with given value Q
    get_model = function(Q) {
      if(!(Q %in% self$nb_blocks)) {
        stop("No such model in the collection. Acceptable parameter values can be found via $nb_blocks")
      }
      Q_rank <- which(sort(self$nb_blocks) == Q)
      self$models[[Q_rank]]
    },

    #' @description Extract best model in the collection
    #' @param crit a character for the criterion used to performed the selection.
    #' Either "ICL", "BIC" or AIC". "ICL" is the default criterion
    #' @return a [`NB_fixed_Q`] object
    get_best_model = function(crit = c("ICL", "BIC", "AIC")) {
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
    #' @param criteria vector of characters. The criteria to plot in `c("deviance", "BIC", "AIC", "ICL")`. Defaults to all of them.
    #' @return a [`ggplot2::ggplot`] graph
    plot = function(criteria = c("deviance", "BIC", "AIC", "ICL")) {
      vlines <- sapply(intersect(criteria, c("BIC")) , function(crit) self$get_best_model(crit)$Q)
      stopifnot(!anyNA(self$criteria[criteria]))

      dplot <- self$criteria %>%
        dplyr::select(dplyr::all_of(c("Q", criteria))) %>%
        tidyr::gather(key = "criterion", value = "value", -Q) %>%
        dplyr::group_by(criterion)
      if("loglik" %in% criteria){dplot[dplot$criterion == "loglik",]$value <- - dplot[dplot$criterion == "loglik",]$value}
      p <- ggplot2::ggplot(dplot, ggplot2::aes(x = Q, y = value, group = criterion, colour = criterion)) +
        ggplot2::geom_line() + ggplot2::geom_point() +
        ggplot2::ggtitle(label    = "Model selection criteria",
                         subtitle = "Lower is better" ) +
        ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = vlines, linetype = "dashed", alpha = 0.25)
      p
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    Q_list            = NA # list of Q values (number of groups) in the collection
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field Q number of blocks
    Q = function(value) private$Q_list,
    #' @field n number of samples
    n = function() self$data$n,
    #' @field p number of responses per sample
    p = function() self$data$p,
    #' @field d number of variables (dimensions in X)
    d = function() self$data$d,
    #' @field criteria a data frame with the values of some criteria ((approximated) log-likelihood, BIC, AIC) for the collection of models
    criteria = function() purrr::map(self$models, "criteria") %>% purrr::reduce(rbind),
    #' @field who_am_I a method to print what model is being fitted
    who_am_I  = function(value){paste0("sparse ", self$noise_cov, " normal-block model with unknown Q")}
  )
)
