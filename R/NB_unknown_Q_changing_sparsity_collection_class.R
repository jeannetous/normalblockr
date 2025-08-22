## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_unknown_Q_changing_sparsity ###############
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 class for a collection of normal-block models with different number of clusters (Q) and different sparsity levels.
#' @param data contains the matrix of responses (Y) and the design matrix (X).
#' @param zero_inflation whether the models should be zero-inflated or not
#' @param control structured list for specific parameters (including initial clustering proposal)
NB_unknown_Q_changing_sparsity <- R6::R6Class(
  classname = "NB_unknown_Q_changing_sparsity",
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field models list of NB_fixed_Q models corresponding to each nb_block value
    models = NA,
    #' @field control store the list of user-defined model settings and optimization parameters
    control = NA,

    #' @description Create a new [`NB_unknown_Q_changing_sparsity`] object.
    #' @param mydata object of NB_data class, with responses and design matrix
    #' @param Q_list list of Q values (number of groups) in the collection
    #' @param zero_inflation boolean to specify whether data is zero-inflated
    #' @param control structured list of parameters to handle sparsity control
    #' @return A new [`NB_unknown_Q_changing_sparsity`] object
    initialize = function(mydata, Q_list, zero_inflation = FALSE,
                          control = NB_control()) {

      ## Store user-defined fields
      self$control <- control
      self$control$zero_inflation <- zero_inflation
      control_ <- control
      self$models <- purrr::map(seq_along(Q_list), function(rank) {
        if (!is.null(control$clustering_init))
          control_$clustering_init <- control$clustering_init[[rank]]
        model <- NB_changing_sparsity$new(mydata, Q_list[rank], zero_inflation, control_)
      })
    },

    #' @description optimizes an NB_changing_sparsity object for each penalty value
    #' @param control optimization parameters (niter and threshold)
    optimize = function(control = list(niter=100, threshold=1e-4, verbose=TRUE)) {
      self$models <- map(self$models, function(model) {
        if(control$verbose) cat("\tnumber of blocks =", model$Q, "          \r")
        flush.console()
        model$optimize(control)
        model
      })
    },

    #' @description returns a collection of NB_unknown models corresponding to given Q
    #' or one single model if penalty is also given
    #' @param Q number of blocks asked by user.
    #' @param sparsity sparsity penalty penalty asked by user
    #' @return either a NB_changing_sparsity or a NB_fixed_Q object
    get_model = function(Q, sparsity = NA) {
      stopifnot("No such model in the collection. Acceptable values can be found via $Q" = Q %in% self$Q_list)
      model <- self$models[[which(self$Q_list == Q)]]
      if (!is.na(sparsity)) {
        if (!(sparsity %in% model$sparsity)) {
          sparsity <-  model$sparsity[[which.min(abs(model$sparsity - sparsity))]]
          cat(paste0("No model with this penalty in the collection. Returning model with closest penalty: ",
                     sparsity,  " Collection penalty values can be found via $sparsity \n"))
        }
        model <- model$models[[which(model$sparsity == sparsity)]]
      }
      model
    },

    #' @description Extract best model in the collection
    #' @param crit a character for the criterion used to performed the selection.
    #' Either "BIC", "EBIC" or "ICL. "ICL" is the default criterion
    #' @return a [`NB_unknown`] object
    get_best_model = function(crit = c("ICL", "BIC", "EBIC")) {
      crit <- match.arg(crit)
      stopifnot(!anyNA(self$criteria[[crit]]))
      id <- 1
      if (length(self$criteria[[crit]]) > 1) {
        id       <- which.min(self$criteria[[crit]])
        best_pen <- self$criteria$sparsity[[id]]
        best_Q   <- self$criteria$Q[[id]]
      }
      model <- self$get_model(best_Q, best_pen)$clone()
      model
    },

    #' @description Display various outputs (goodness-of-fit criteria, robustness, diagnostic) associated with a collection of network fits (a [`Networkfamily`])
    #' @param criterion The criteria to plot in `c("deviance", BIC", "EBIC", "ICL")`. Defaults deviance.
    #' @param n_intervals number of intervals into which the penalties range should be splitted
    #' @importFrom tidyr gather
    #' @return a [`ggplot`] heatmap
    plot = function(criterion = c("deviance", "ICL", "BIC", "EBIC"),
                    n_intervals = NULL) {
      criterion   <- match.arg(criterion)
      if(is.null(n_intervals)) n_intervals <- round(0.1 * length(unique(self$criteria$sparsity )))
      df <- self$criteria %>% mutate(pen_binned = cut(sparsity, breaks = n_intervals)) %>%
        group_by(pen_binned, Q) %>% summarize(avg_crit = mean(.data[[criterion]]), .groups = "drop")
      p  <- ggplot2::ggplot(df, aes(x = pen_binned, y = Q, fill = avg_crit)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::theme_minimal() +
        ggplot2::ggtitle(label    = criterion,
                         subtitle = "Lower is better" ) +
        ggplot2::labs(x = "Penalties (Binned)", y = "Q", fill = paste0("Average ", criterion),
                      title = criterion) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      p
    }

  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field Q_list  number of blocks
    Q_list = function(value) map_dbl(self$models, "Q"),
    #' @field criteria a data frame with the values of some criteria ((approximated) log-likelihood, BIC, ICL) for the collection of models
    criteria = function() map_df(self$models, "criteria") ,
    #' @field sparsity list of penalties used for each Q
    sparsity = function(){
      self$criteria %>%
        dplyr::group_by(Q) %>%
        dplyr::summarize(sparsity = paste(round(sparsity, 2), collapse = ", "))
    },
    #' @field who_am_I a method to print what model is being fitted
    who_am_I  = function(value){
      paste("Collection of ",
            ifelse(self$control$zero_inflation, " zero-inflated", ""),
            self$control$noise_covariance,
            "normal-block models with different values of Q and different penalties.")
    }

  )
)
