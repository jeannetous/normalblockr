## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB ############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 class for a generic sparse Normal Block model
#' @param data contains the matrix of responses (Y) and the design matrix (X).
#' @param Q number of clusters
#' @param control structured list of more specific parameters, to generate with NB_control
NB <- R6::R6Class(
  classname = "NB",
  inherit   = normal_models,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(

    #' @description Create a new [`NB`] object.
    #' @param data object of normal_data class, with responses and design matrix
    #' @param penalty penalty on the network density
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return A new [`NB`] object
    initialize = function(data, Q, penalty = 0, control = NB_control()) {
      super$initialize(data, control)
      private$C <- matrix(NA, self$data$n, Q)

      ## variant (either diaongal or spherical residuals covariance)
      private$res_covariance <- control$noise_covariance

      ## pointer to the chosen optimization function
      private$optimizer <- ifelse(control$heuristic,
                                  private$heuristic_optimize,
                                  private$EM_optimize)
      ## pointer to the chosen clustering function for heuristic approach
      private$approx <- control$heuristic
      private$clustering_approx <-
        switch(control$clustering_approx,
               "residuals"  = private$heuristic_cluster_residuals,
               "covariance" = private$heuristic_cluster_sigma
        )
      ## penalty mask
      private$lambda <- penalty
      weights <- matrix(1, self$Q, self$Q)
      diag(weights) <- 0
      if (!is.null(control$sparsity_weights)) {
        weights <- control$sparsity_weights
      }
      private$weights <- weights

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
                            plot = TRUE) {
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
    lambda            = NA, # scalar controlling the overall sparsity
    weights           = NA, # sparsity weights specific to each pairs of group
    res_covariance    = NA, # shape of the residuals covariance (diagonal or)
    approx            = NA, # use approximation/heuristic approach or not
    clustering_approx = NA, # clustering function in the heuristic approach

    get_OmegaQ = function(Sigma) {
      if (self$penalty == 0) {
        Omega <- solve(Sigma)
      } else {
        glasso_out <- glassoFast::glassoFast(Sigma, rho = self$penalty * self$penalty_weights)
        if (anyNA(glasso_out$wi)) {
          warning(
            "GLasso fails, the penalty is probably too small and the system badly conditionned \n reciprocal condition number =",
            rcond(Sigma), "\n We send back the original matrix and its inverse (unpenalized)."
          )
          Omega <- solve(Sigma)
        } else {
          Omega <- Matrix::symmpart(glasso_out$wi)
        }
      }
      Omega
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for heuristic inference----------------------
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    heuristic_optimize = function(control){
      parameters <- private$get_heuristic_parameters()
      c(parameters, list(ll_list = NA))
    },

    heuristic_SigmaQ_from_Sigma = function(Sigma){
      Sigma_Q <- (t(private$C) %*% Sigma %*% private$C) / outer(colSums(private$C), colSums(private$C))
### TODO: why is there any NA?
      if (anyNA(Sigma_Q)) {
        diag(Sigma_Q)[is.na(diag(Sigma_Q))] <- mean(diag(Sigma_Q)[!is.na(diag(Sigma_Q))])
        Sigma_Q[is.na(Sigma_Q)] <- 0
      }
      Sigma_Q
    },

    heuristic_cluster_sigma = function(R){
      sink('/dev/null')
      mySBM <- cov(R) %>%
        sbm::estimateSimpleSBM("gaussian",
                               estimOption=list(verbosity=0, exploreMin=self$Q, verbosity=0, plot=FALSE, nbCores=1)
        )
      sink()
      mySBM$setModel(self$Q)
      mySBM$memberships |> as_indicator()
    },

    heuristic_cluster_residuals = function(R){
      kmeans(t(R), self$Q, nstart = 30, iter.max = 50)$cluster |>
        as_indicator()
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field inference_method inference procedure used (heuristic or integrated with EM)
    inference_method = function(value) ifelse(private$approx, "heuristic", "integrated"),
    #' @field Q number of blocks
    Q = function(value) as.integer(ncol(private$C)),
    #' @field nb_param number of parameters in the model
    nb_param = function(value) {
      nb_param_D <- ifelse(private$res_covariance == "diagonal", self$p, 1)
      as.integer(super$nb_param + self$Q + self$n_edges + nb_param_D)
      }, # adding OmegaQ and dm1
    #' @field n_edges number of edges of the network (non null coefficient of the sparse precision matrix OmegaQ)
    n_edges  = function(value) sum(private$OmegaQ[upper.tri(private$OmegaQ, diag = FALSE)] != 0),
    #' @field model_par a list with the matrices of the model parameters: B (covariates), dm1 (species variance), OmegaQ (groups precision matrix))
    model_par = function(value) c(super$model_par, list(OmegaQ = private$OmegaQ)),
    #' @field penalty (overall sparsity parameter)
    penalty = function(value) private$lambda,
    #' @field penalty_weights (weights associated to each pair of groups)
    penalty_weights = function(value) private$weights,
    #' @field penalty_term (penalty term in log-likelihood due to sparsity)
    penalty_term = function(value) self$penalty * sum(abs(self$penalty_weights * private$OmegaQ)),
    #' @field loglik (or its variational lower bound)
    loglik = function(value) if (private$approx) NA else super$loglik + self$penalty_term,
    #' @field EBIC variational lower bound of the EBIC
    EBIC      = function(value) self$BIC + 2 * ifelse(self$n_edges > 0, self$n_edges * log(self$Q), 0),
    #' @field criteria a vector with loglik, BIC and number of parameters
    criteria   = function(value) c(Q = self$Q, n_edges = self$n_edges, penalty = self$penalty, super$criteria, EBIC = self$EBIC),
    #' @field get_res_covariance whether the residual covariance is diagonal or spherical
    get_res_covariance = function(value) private$res_covariance,
    #' @field clustering given as the list of elements contained in each cluster
    clustering = function(value) get_clusters(private$C),
    #' @field elements_per_cluster given as the list of elements contained in each cluster
    elements_per_cluster = function(value) split(names(self$clustering), self$clustering)
  )
)
