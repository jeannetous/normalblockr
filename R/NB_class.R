## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB #######################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 generic class for normal-block models
#' @param Y the matrix of responses (called Y in the model).
#' @param X design matrix (called X in the model).
#' @export
NB <- R6::R6Class(
  classname = "NB",
  inherit = MVEM,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field penalty penalty on the network density
    penalty = NULL,
    #' @field sparsity_weights distribution of sparsity on the elements of the network
    sparsity_weights = NULL,
    #' @field Q number of blocks
    Q = NULL,

    #' @description Create a new [`NB`] object.
    #' @param Y the matrix of responses (called Y in the model).
    #' @param X design matrix (called X in the model).
    #' @param Q number of groups
    #' @param penalty penalty on the network density
    #' @param control structured list of parameters, including sparsity_weights
    #' @return A new [`nb_fixed`] object
    initialize = function(Y, X, Q, penalty = 0, control = NB_param()) {
      super$initialize(Y, X)
      self$Q <- Q
      private$omegaQ <- diag(1, Q, Q)
      self$penalty <- penalty
      if (penalty > 0) {
        sparsity_weights  <- control$sparsity_weights
        if(is.null(sparsity_weights)){
          sparsity_weights <- matrix(1, self$Q, self$Q)
          diag(sparsity_weights) <- 0
        }
        self$sparsity_weights  <- sparsity_weights
      }
    },

    #' @description
    #' Update a [`NB`] object
    #' @param B regression matrix
    #' @param dm1 diagonal vector of species inverse variance matrix
    #' @param omegaQ groups inverse variance matrix
    #' @param ll_list log-likelihood during optimization
    #' @return Update the current [`NB`] object
    update = function(B = NA, dm1 = NA, omegaQ = NA, ll_list = NA) {
      super$update(B=B, dm1=dm1, ll_list=ll_list)
      if (!anyNA(omegaQ)) private$omegaQ  <- omegaQ
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Extractors ------------------------
    #' @description Extract interaction network in the latent space
    #' @param type edge value in the network. Can be "support" (binary edges), "precision" (coefficient of the precision matrix) or "partial_cor" (partial correlation between species)
    #' @importFrom Matrix Matrix
    #' @return a square matrix of size `NB_fixed_blocks_class$Q`
    latent_network = function(type = c("partial_cor", "support", "precision")) {
      net <- switch(
        match.arg(type),
        "support"     = 1 * (private$omegaQ != 0 & !diag(TRUE, ncol(private$omegaQ))),
        "precision"   = private$omegaQ,
        "partial_cor" = {
          tmp <- -private$omegaQ / tcrossprod(sqrt(diag(private$omegaQ))); diag(tmp) <- 1
          tmp
        }
      )
      ## Enforce sparse Matrix encoding to avoid downstream problems with igraph::graph_from_adjacency_matrix
      ## as it fails when given dsyMatrix objects
      Matrix(net, sparse = TRUE)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    omegaQ    = NA,   # groups variance matrix
    ## funciton to optimize OmegaQ (in the sparse and non sparse case)
    get_omegaQ = function(sigmaQ) {
      if (self$penalty == 0) {
        omegaQ <- solve(sigmaQ)
      } else {
        glasso_out <- glassoFast::glassoFast(sigmaQ, rho = self$penalty * self$sparsity_weights)
        if (anyNA(glasso_out$wi)) stop("GLasso fails")
        omegaQ <- Matrix::symmpart(glasso_out$wi)
      }
      omegaQ
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field nb_param number of parameters in the model
    nb_param = function() as.integer(super$nb_param + self$Q + self$n_edges),
    #' @field n_edges number of edges of the network (non null coefficient of the sparse precision matrix OmegaQ)
    n_edges  = function() {sum(private$omegaQ[upper.tri(private$omegaQ, diag = FALSE)] != 0)},
    #' @field model_par a list with the matrices of the model parameters: B (covariates), dm1 (species variance), omegaQ (groups precision matrix))
    model_par = function() list(B = private$B, dm1 = private$dm1, omegaQ = private$omegaQ),
    #' @field penalty_term (penalty term in log-likelihood due to sparsity)
    penalty_term = function() self$penalty * sum(abs(self$sparsity_weights * private$omegaQ)),
    # #' @field EBIC variational lower bound of the EBIC
    # EBIC      = function() {self$BIC + 2 * ifelse(self$n_edges > 0, self$n_edges * log(.5 * self$Q*(self$Q - 1)/self$n_edges), 0)},
    #' @field criteria a vector with loglik, BIC and number of parameters
    criteria   = function() {
      res <- super$criteria
      ## res$EBIC <- self$EBIC
      res$Q <- self$Q
      res$n_edges <- self$n_edges
      res$penalty <- self$penalty
      res
    }
  )
)


#' NB_fixed_Q_sparse_param
#' @param sparsity_weights weights with which penalty should be applied in case
#' sparsity is required, non-0 values on the diagonal mean diagonal shall be
#' penalized too (default is non-penalized diagonal)
#' Generates control parameters for NB models
#' @export
NB_param <- function(sparsity_weights = NULL){
  structure(list(
    sparsity_weights = sparsity_weights
  ))
}
