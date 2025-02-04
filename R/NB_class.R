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
    initialize = function(Y, X, Q, penalty = 0, control = normal_block_control()) {
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
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Graphical methods------------------
    #' @description plot the latent network.
    #' @param type edge value in the network. Either "precision" (coefficient of the precision matrix) or "partial_cor" (partial correlation between species).
    #' @param output Output type. Either `igraph` (for the network) or `corrplot` (for the adjacency matrix)
    #' @param edge.color Length 2 color vector. Color for positive/negative edges. Default is `c("#F8766D", "#00BFC4")`. Only relevant for igraph output.
    #' @param node.labels vector of character. The labels of the nodes. The default will use the column names ot the response matrix.
    #' @param remove.isolated if `TRUE`, isolated node are remove before plotting. Only relevant for igraph output.
    #' @param layout an optional igraph layout. Only relevant for igraph output.
    #' @param plot logical. Should the final network be displayed or only sent back to the user. Default is `TRUE`.
    plot_network = function(type            = c("partial_cor", "support"),
                            output          = c("igraph", "corrplot"),
                            edge.color      = c("#F8766D", "#00BFC4"),
                            remove.isolated = FALSE,
                            node.labels     = NULL,
                            layout          = igraph::layout_in_circle,
                            plot = TRUE){
      if(anyNA(private$omegaQ)) stop("NA in the precision matrix")

      type   <- match.arg(type)
      output <- match.arg(output)

      net <- self$latent_network(type)

      if (output == "igraph") {
        G <-  igraph::graph_from_adjacency_matrix(net, mode = "undirected", weighted = TRUE, diag = FALSE)

        if (!is.null(node.labels)) {
          igraph::V(G)$label <- node.labels
        } else {
          igraph::V(G)$label <- unlist(lapply(1:ncol(net), f <- function(x) paste0("Cluster_", x)))
        }
        ## Nice nodes
        V.deg <- igraph::degree(G)/sum(igraph::degree(G))
        igraph::V(G)$label.cex <- V.deg / max(V.deg) + .5
        igraph::V(G)$size <- V.deg * 100
        igraph::V(G)$label.color <- rgb(0, 0, .2, .8)
        igraph::V(G)$frame.color <- NA
        ## Nice edges
        igraph::E(G)$color <- ifelse(igraph::E(G)$weight > 0, edge.color[1], edge.color[2])
        if (type == "support")
          igraph::E(G)$width <- abs(igraph::E(G)$weight)
        else
          igraph::E(G)$width <- 15*abs(igraph::E(G)$weight)

        if (remove.isolated) {
          G <- delete.vertices(G, which(degree(G) == 0))
        }
        if (plot) plot(G, layout = layout)
      }
      if (output == "corrplot") {
        if (plot) {
          if (ncol(net) > 100)
            colnames(net) <- rownames(net) <- rep(" ", ncol(net))
          G <- net
          diag(net) <- 0
          corrplot(as.matrix(net), method = "color", is.corr = FALSE, tl.pos = "td", cl.pos = "n", tl.cex = 0.5, type = "upper")
        } else  {
          G <- net
        }
      }
      invisible(G)
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
        if (anyNA(glasso_out$wi)) {
          warning(
            "GLasso fails, the penalty is probably too small and the system badly conditionned \n reciprocal condition number =",
            rcond(sigmaQ), "\n We send back the original matrix and its inverse (unpenalized)."
          )
          omegaQ <- solve(sigmaQ)
        } else {
          omegaQ <- Matrix::symmpart(glasso_out$wi)
        }
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
    #' @field EBIC variational lower bound of the EBIC
    EBIC      = function() {self$BIC + 2 * ifelse(self$n_edges > 0, self$n_edges * log(.5 * self$Q*(self$Q - 1)/self$n_edges), 0)},
    #' @field criteria a vector with loglik, BIC and number of parameters
    criteria   = function() {
      res   <- super$criteria
      res$Q <- self$Q
      res$n_edges <- self$n_edges
      res$penalty <- self$penalty
      res$EBIC    <- self$EBIC
      res
    }
  )
)

