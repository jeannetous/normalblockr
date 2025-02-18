## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_fixed_sparsity ############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 class for a generic normal model
#' @param data contains the matrix of responses (Y) and the design matrix (X).
#' @param penalty to apply on variance matrix when calling GLASSO
#' @param Q number of clusters
#' @param control structured list of more specific parameters, to generate with normal_control
NB_fixed_sparsity <- R6::R6Class(
  classname = "NB_fixed_sparsity",
  inherit   = normal_fixed_sparsity,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(

    #' @field Q number of blocks
    Q = NULL,

    #' @description Create a new [`NB_fixed_sparsity`] object.
    #' @param data object of normal_data class, with responses and design matrix
    #' @param penalty penalty on the network density
    #' @param control structured list of more specific parameters, to generate with normal_control
    #' @return A new [`NB_fixed_sparsity`] object
    initialize = function(data, Q, penalty = 0,
                          control = normal_control()) {
      super$initialize(data,  penalty, control)
      self$Q <- Q
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
    #' Update a [`NB_fixed_sparsity`] object
    #' @param B regression matrix
    #' @param dm1 diagonal vector of species inverse variance matrix
    #' @param OmegaQ groups inverse variance matrix
    #' @param ll_list log-likelihood during optimization
    #' @return Update the current [`NB`] object
    update = function(B = NA, OmegaQ = NA, dm1 = NA,  ll_list = NA) {
      super$update(B = B, ll_list = ll_list)
      if (!anyNA(dm1))        private$dm1     <- dm1
      if (!anyNA(OmegaQ))     private$OmegaQ  <- OmegaQ
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
        "support"     = 1 * (private$OmegaQ != 0 & !diag(TRUE, ncol(private$OmegaQ))),
        "precision"   = private$OmegaQ,
        "partial_cor" = {
          tmp <- -private$OmegaQ / tcrossprod(sqrt(diag(private$OmegaQ))); diag(tmp) <- 1
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
      if(anyNA(private$OmegaQ)) stop("NA in the precision matrix")

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
    C          = NA, # the matrix of species groups
    dm1        = NA, # diagonal of variables' inverse variance matrix
    OmegaQ     = NA, # precision matrix for clusters

    heuristic_SigmaQ_from_Sigma = function(Sigma){
      Sigma_Q <- (t(private$C) %*% Sigma %*% private$C) / outer(colSums(private$C), colSums(private$C))
      if(anyNA(Sigma_Q)){
        diag(Sigma_Q)[is.na(diag(Sigma_Q))] <- mean(diag(Sigma_Q)[!is.na(diag(Sigma_Q))])
        Sigma_Q[is.na(Sigma_Q)] <- 0
      }
      Sigma_Q
    },

    heuristic_loglik = function(B, OmegaQ){
      Sigma_tilde <- private$C %*% solve(OmegaQ) %*% t(private$C) + diag(1e-6, self$p)
      log_det_Sigma_tilde <- as.numeric(determinant(Sigma_tilde, logarithm = TRUE)$modulus)
      R <- self$data$Y - self$data$X %*% B
      J <- -.5 * self$p * log(2 * pi)
      J <- J + .5 * self$n  * log_det_Sigma_tilde
      J <- J - .5 * sum(diag((R %*% solve(Sigma_tilde) %*% t(R))))
      J
    },

    get_clustering = function(Sigma = NA, R = NA){
      if(self$clustering_method ==  "cluster_sigma"){
        C <- private$cluster_sigma(Sigma)
      }else{
        C <- private$cluster_residuals(R)
      }
    },

    cluster_sigma = function(Sigma){
      sink('/dev/null')
      mySBM <- Sigma %>%
        sbm::estimateSimpleSBM("gaussian",
                               estimOption=list(verbosity=0, exploreMin=self$Q, verbosity=0, plot=FALSE, nbCores=1)
        )
      sink()
      mySBM$setModel(self$Q)
      return(as_indicator(mySBM$memberships))
    },

    cluster_residuals = function(R){
      return(as_indicator(kmeans(t(R), self$Q, nstart = 30, iter.max = 50)$cluster))
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field nb_param number of parameters in the model
    nb_param = function() as.integer(super$nb_param + self$Q + self$n_edges),
    #' @field n_edges number of edges of the network (non null coefficient of the sparse precision matrix OmegaQ)
    n_edges  = function() sum(private$OmegaQ[upper.tri(private$OmegaQ, diag = FALSE)] != 0),
    #' @field model_par a list with the matrices of the model parameters: B (covariates), dm1 (species variance), OmegaQ (groups precision matrix))
    model_par = function(){
      par <- super$model_par
      par$dm1 <- private$dm1 ; par$OmegaQ <- private$OmegaQ
      par
    },
    #' @field penalty_term (penalty term in log-likelihood due to sparsity)
    penalty_term = function() self$penalty * sum(abs(self$sparsity_weights * private$OmegaQ)),
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
    },
    #' @field clustering given as the list of elements contained in each cluster
    clustering = function()  get_clusters(private$C),
    #' @field elements_per_cluster given as the list of elements contained in each cluster
    elements_per_cluster = function() split(names(self$clustering), self$clustering)
  )
)
