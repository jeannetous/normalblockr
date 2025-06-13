## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB ############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 abstract class for a generic sparse Normal Block model
NB <- R6::R6Class(
  classname = "NB",
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field data object of NBData class, with responses and design matrix
    data  = NULL,

    #' @description Create a new [`NB`] object.
    #' @param data object of NBData class, with responses and design matrix
    #' @param Q number of block/cluster
    #' @param sparsity sparsity penalty on the network density
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return A new [`NB`] object
    initialize = function(data, Q, sparsity = 0, control = NB_control()) {
      self$data <- data

      stopifnot("There cannot be more blocks than there are entities to cluster" = Q <= ncol(self$data$Y))

      ## variant (either diagonal or spherical residuals covariance)
      private$res_covariance <- control$noise_covariance

      ## pointer to the chosen optimization function
      private$optimizer <- ifelse(control$heuristic,
                                  private$heuristic_optimize,
                                  private$EM_optimize)
      ## pointer to the chosen clustering function for heuristic approach
      private$approx <- control$heuristic
      private$clustering_approx <-
        switch(control$clustering_approx,
               "kmeans"    = private$heuristic_cluster_residuals_kmeans,
               "ward2"     = private$heuristic_cluster_sigma_ward2,
               "sbm"       = private$heuristic_cluster_sigma_sbm
        )

      ## penalty mask
      private$sparsity_ <- sparsity
      weights <- matrix(1, Q, Q)
      diag(weights) <- 0
      if (!is.null(control$sparsity_weights)) {
        weights <- control$sparsity_weights
      }
      private$weights <- weights

      cl0 <- control$clustering_init
      if (!is.null(cl0)) {
        if (!is.vector(cl0) & !is.matrix(cl0)) stop("Labels must be encoded in vector of labels or indicator matrix")
        if (is.vector(cl0)) {
          if (any(cl0 < 1 | cl0 > Q))
            stop("Cluster labels must be between 1 and Q")
          if (length(cl0) != self$p)
            stop("Cluster labels must match the number of Y's columns")
          if (length(unique(cl0)) != Q)
            stop("The number of clusters in the initial clustering must be equal to Q.")
          cl0 <- as_indicator(cl0)
        } else {
          if (nrow(cl0) != self$p)
            stop("Cluster-indicating matrix must have as many rows as Y has columns")
          if (ncol(cl0) != Q)
            stop("Cluster-indicating matrix must have Q columns")
          if ((min(colSums(cl0)) < 1) & !Q)
            stop("The number of clusters in the initial clustering must be equal to Q.")
        }
        private$C <- cl0
      } else {
        private$C <- matrix(NA, self$data$n, Q)
      }

      ## ZI parameters that will remain fixed
      if(control$zero_inflation_type == "column"){
        private$kappa <- colMeans(data$zeros)
        private$ZI_cond_mean <-
          sum(xlogy(data$zeros, rep(private$kappa, each = data$n) )) +
          sum(xlogy(data$zeros_bar, 1 - rep(private$kappa, each = data$n) ))
      }else{
        B0_list <- lapply(1:self$data$p,
                          f <- function(j){
                            df <- data.frame("zeros" = data$zeros[,j], self$data$X0)
                            model <- glm(zeros ~ 0 + ., family=binomial(link = "logit"), data=df)
                            return(model$coefficients)})
        private$B0 <- t(sapply(B0_list, unlist))
        private$kappa <- apply(self$data$X0 %*% private$B0, MARGIN = c(1,2), FUN = sigmoid)
        private$ZI_cond_mean <-
          sum(xlogy(data$zeros, private$kappa)) +
          sum(xlogy(data$zeros_bar, 1 - private$kappa))
      }
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Setters    ------------------------
    #' @description
    #' Update a [`NB`] object
    #'
    #' All possible parameters of the child classes
    #' @param B regression matrix
    #' @param dm1 diagonal vector of inverse variance matrix (variables level)
    #' @param C the matrix of groups memberships (posterior probabilities)
    #' @param OmegaQ groups inverse variance matrix
    #' @param gamma  variance of posterior distribution of W
    #' @param mu mean for posterior distribution of W
    #' @param kappa vector of zero-inflation probabilities
    #' @param alpha vector of groups probabilities
    #' @param M variational mean for posterior distribution of W
    #' @param S variational diagonal of variances for posterior distribution of W
    #' @param ll_list  list of log-lik (elbo) values
    #' @return Update the current [`normal`] object
    update = function(B = NA,
                      dm1 = NA,
                      C = NA,
                      OmegaQ = NA,
                      gamma = NA,
                      mu = NA,
                      kappa = NA,
                      alpha = NA,
                      M = NA,
                      S = NA,
                      ll_list = NA) {
      if (!anyNA(B))       private$B       <- B
      if (!anyNA(dm1))     private$dm1     <- dm1
      if (!anyNA(C))       private$C       <- C
      if (!anyNA(OmegaQ))  private$OmegaQ  <- OmegaQ
      if (!anyNA(gamma))   private$gamma   <- gamma
      if (!anyNA(kappa))   private$kappa   <- kappa
      if (!anyNA(mu))      private$mu      <- mu
      if (!anyNA(alpha))   private$alpha   <- alpha
      if (!anyNA(M))       private$M       <- M
      if (!anyNA(S))       private$S       <- S
      if (!anyNA(ll_list)) private$ll_list <- ll_list
    },

    #' @description calls optimization (EM or heuristic) and updates relevant fields
    #' @param control a list for controlling the optimization proces
    #' @return optimizes the model and updates its parameters
    optimize = function(control = list(niter = 100, threshold = 1e-4)) {
      optim_out <- private$optimizer(control)
      do.call(self$update, optim_out)
    },

    #' @description Create a clone of the current [`NB`] object after splitting cluster `cl`
    #' We split the cluster according to the species variances
    #' @param index index (integer) of the cluster to split
    #' @param in_place should the split applied to the object itself, or should a copy be sent?
    #' default FALSE (send a copy)
    #' @return A new [`NB`] object
    split = function(index, in_place = FALSE) {
      ## update private fields related to group parameters
      ## C, OmegaQ, M, S, sparsity_weights

      ## indices of individuals split within the cluster
      cl  <- self$clustering == index
      var <- 1/private$dm1; var_median <- median(var[cl])
      split1 <- (var > var_median) & cl ;  split2 <- (var <= var_median) & cl

      ## Cluster split
      new_C <- cbind(private$C, .Machine$double.eps)
      new_C[split1, index] <- new_C[split1, index] - .Machine$double.eps
      new_C[split2, self$Q + 1] <- new_C[split2, index]
      new_C[split2, index] <- .Machine$double.eps
      new_C <- new_C / rowSums(new_C)

      ## Variational means
      new_M <- cbind(private$M, 0)
      new_M[split2, self$Q + 1] <- new_M[split2, index]
      new_M[split2, index] <- 0

      ## Variational variances
      if (is.matrix(private$S)) {
        new_S <- cbind(private$S, 0.1)
        new_S[split2, self$Q + 1] <- new_C[split2, index]
        new_S[split2, index] <- 0.1
      } else {
        new_S <- c(private$S, mean(private$S))
      }

      ## Precision matrix
      new_OmegaQ <- cbind(rbind(private$OmegaQ,  0), 0)
      new_OmegaQ[index, index] <- private$OmegaQ[index, index]/2
      new_OmegaQ[self$Q + 1, self$Q + 1] <- private$OmegaQ[index, index]/2

      ## Sparsity weights
      if (self$Q == 1) {
        new_weights <- matrix(c(0,1,1,0), 2, 2)
      } else {
        weights_cl <-  private$weights[index, setdiff(1:self$Q, index)]
        weights_cl <-  c(weights_cl, mean(weights_cl))
        new_weights <- cbind(rbind(private$weights, weights_cl, deparse.level = 0),
                             c(weights_cl, 0))
      }

      if (in_place) {
        self$update(C = new_C, OmegaQ = new_OmegaQ, M = new_M, S = new_S)
        self$sparsity_weights <- new_weights
        return(invisible(self))
      } else {
        new_NB <- self$clone()
        new_NB$update(C = new_C, OmegaQ = new_OmegaQ, M = new_M, S = new_S)
        new_NB$sparsity_weights <- new_weights
        return(invisible(new_NB))
      }
    },

    #' @description generate and select a set of candidate models
    #' by splitting the clusters of the current model
    candidates_split = function() {
      # do not split groups with less than 2 guys
      candidates <- map((1:self$Q)[self$cluster_sizes > 1], self$split)
      # keep candidates with at least 2 guys per cluster and non empty split
      clustering_sizes <- map(candidates, "clustering") %>% map(table)
      min_sizes  <- clustering_sizes %>% map_dbl(min)
      n_clusters <- clustering_sizes %>% map_dbl(length)
      candidates <- candidates[min_sizes > 1 & n_clusters == self$Q + 1]

      for (i in seq_along(candidates))
        candidates[[i]]$optimize(list(niter = 5, threshold = 1e-4))
      candidates
    },

    #' @description generate and select a set of candidate models
    #' by merging the clusters of the current model
    candidates_merge = function() {
      stopifnot("need at least two clusters to merge them" = self$Q > 1)
      candidates <- map(combn(self$Q, 2, simplify = FALSE), self$merge)
      for (i in seq_along(candidates))
        candidates[[i]]$optimize(list(niter = 5, threshold = 1e-4))
      candidates
    },

    #' @description Create a clone of the current [`NB`] object after merging clusters `cl1` and `cl2`
    #' @param indices indices (couple of integer) of the clusters to merge
    #' @param in_place should the split applied to the object itself, or should a copy be sent?
    #' default FALSE (send a copy)
    #' @return A new [`NB`] object
    merge = function(indices, in_place=FALSE) {

      ## sorting by increasing group label
      indices <- sort(indices)

      ## Cluster merge
      new_C <- private$C[, -indices[2]]
      new_C[, indices[1]] <- private$C[, indices[1]] + private$C[, indices[2]]

      ## Variational means
      new_M <- private$M[, -indices[2]]
      new_M[, indices[1]] <- .5 * (private$M[, indices[1]] + private$M[, indices[2]])

      ## Variational variances
      if (is.matrix(private$S)) {
        new_S <- private$S[, -indices[2]]
        new_S[, indices[1]] <- .5 * (private$S[, indices[1]] + private$S[, indices[2]])
      } else {
        new_S <- private$S[-indices[2]]
        new_S[indices[1]] <- .5 * (private$S[indices[1]] + private$S[indices[2]])
      }

      ## Precision matrix
      new_OmegaQ <- private$OmegaQ[-indices[2], -indices[2]]
      new_OmegaQ[indices[1], indices[1]] <-
        .5 * (private$OmegaQ[indices[1], indices[1]] + private$OmegaQ[indices[2], indices[2]])

      ## Sparsity weights
      new_weights <-  private$weights[-indices[2], -indices[2]]

      if (in_place) {
        self$update(C = new_C, OmegaQ = new_OmegaQ, M = new_M, S = new_S)
        self$sparsity_weights <- new_weights
        return(self)
      } else {
        new_NB <- self$clone()
        new_NB$update(C = new_C, OmegaQ = new_OmegaQ, M = new_M, S = new_S)
        new_NB$sparsity_weights <- new_weights
        return(new_NB)
      }
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Graphical methods------------------
    #' @param type char for line type (see plot.default)
    #' @param log char for logarithmic axes (see plot.default)
    #' @param neg boolean plot negative log-likelihood (useful when log="y")
    #' @description plots log-likelihood values during model optimization
    plot_loglik = function(type = "b", log = "xy", neg = TRUE) {
      neg <- ifelse(neg, -1, 1)
      plot(seq_along(self$objective), neg * self$objective, type = type, log = log)
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
        igraph::V(G)$size <- table(self$clustering) * 100 / self$p
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
    B                 = NA, # regression matrix
    dm1               = NA, # diagonal vector of inverse variance matrix (variables level)
    C                 = NA, # the matrix of posterior probabilities (tau) or group affectation
    OmegaQ            = NA, # precision matrix for clusters
    kappa             = NA, # vector of zero-inflation probabilities
    B0                = NA, # vector of zero-inflation regression matrix
    alpha             = NA, # vector of groups probabilities
    gamma             = NA, # variance of  posterior distribution of W
    mu                = NA, # mean for posterior distribution of W
    M                 = NA, # variational mean for posterior distribution of W
    S                 = NA, # variational diagonal of variances for posterior distribution of W
    optimizer         = NA, # a link to the function that perform the optimization
    ll_list           = NA, # list of log-likelihoods or ELBOs
    sparsity_         = NA, # scalar controlling the overall sparsity
    weights           = NA, # sparsity weights specific to each pairs of group
    res_covariance    = NA, # shape of the residuals covariance (diagonal or spherical)
    approx            = NA, # use approximation/heuristic approach or not
    clustering_approx = NA, # clustering function in the heuristic approach
    ZI_cond_mean      = NA, # conditional mean of the ZI component (fixed)

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for integrated EM inference------------------
    EM_optimize = function(control) {
      parameters <- private$EM_initialize()
      ll_list    <- do.call(private$compute_loglik, parameters)
      for (h in 1:control$niter) {
        parameters <- do.call(private$EM_step, parameters)
        ll_list    <- c(ll_list, do.call(private$compute_loglik, parameters))
        if (abs(ll_list[h + 1] - ll_list[h]) < control$threshold)
          break
      }
      c(parameters, list(ll_list = ll_list))
    },
    EM_step = function() {},
    EM_initialize = function() {},
    compute_loglik  = function() {},

    get_OmegaQ = function(Sigma) {
      if (private$sparsity_ == 0) {
        Omega <- solve(Sigma)
      } else {
        glasso_out <- glassoFast::glassoFast(Sigma, rho = private$sparsity_ * self$sparsity_weights)
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

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## MLE of MV Normal distribution
    multivariate_normal_inference = function(){
      B     <- self$data$XtXm1 %*% self$data$XtY
      R     <- self$data$Y - self$data$X %*% B
      Sigma <- cov(R)
      list(B = B, R = R, Sigma = Sigma)
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## MLE of ZI Diagonal Normal distribution
    zi_diag_normal_inference = function(){
      B     <- self$data$XtXm1 %*% self$data$XtY
      dm1   <- self$data$nY / colSums(self$data$zeros_bar * (self$data$Y - self$data$X %*% B)^2)
      for (i in 1:3) { # a couple of iterates is enough
        B     <- private$zi_diag_normal_optim_B(B, dm1)
        dm1   <- self$data$nY / colSums(self$data$zeros_bar * (self$data$Y - self$data$X %*% B)^2)
      }
      R <- self$data$zeros_bar * (self$data$Y - self$data$X %*% B)
      list(B = B, dm1 = dm1, kappa = private$kappa, R = R)
    },

    zi_diag_normal_obj_grad_B = function(B_vec, DM1) {
      R <- self$data$Y - self$data$X %*% matrix(B_vec, nrow = self$d, ncol = self$p)
      grad <- crossprod(self$data$X, DM1 * R)
      obj <- -.5 * sum(DM1 * R^2)
      res <- list("objective" = -obj, "gradient"  = -grad)
      res
    },

    zi_diag_normal_optim_B = function(B0, dm1) {
      DM1 <- matrix(dm1, self$data$n, self$data$p, byrow = TRUE) * self$data$zeros_bar
      res <- nloptr::nloptr(
        x0 = as.vector(B0),
        eval_f = private$zi_diag_normal_obj_grad_B,
        opts = list(
          algorithm = "NLOPT_LD_LBFGS",
          maxeval = 100
        ),
        DM1 = DM1
      )
      newB <- matrix(res$solution, nrow = self$d, ncol = self$p)
      newB
    },

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

    heuristic_clustering = function(R) {
      clustering <- private$clustering_approx(R)
      if (length(unique(clustering)) < self$Q) {
        clustering <- cutree(ClustOfVar::hclustvar(R), self$Q)
      }
      C <- as_indicator(clustering)
      if (min(colSums(C)) < 1) warning("Initialization failed to place elements in each cluster")
      C
    },

    heuristic_cluster_sigma_ward2 = function(R){
      hc <- hclust(dist(1 - cor(R)), method = "ward.D2")
      cutree(hc, self$Q)
    },

    heuristic_cluster_sigma_sbm = function(R){
      options <- list(verbosity=0, exploreMin=self$Q, verbosity=0, plot=FALSE, nbCores=1)
      mySBM <- sbm::estimateSimpleSBM(cov(R), "gaussian", estimOptions = options)
      mySBM$setModel(self$Q)
      mySBM$memberships
    },

    heuristic_cluster_residuals_kmeans = function(R) {
      kmeans(t(R), self$Q, nstart = 30, iter.max = 50)$cluster
    }

  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field inference_method inference procedure used (heuristic or integrated with EM)
    inference_method = function(value) ifelse(private$approx, "heuristic", "integrated"),
    #' @field n number of samples
    n = function() self$data$n,
    #' @field p number of responses per sample
    p = function() self$data$p,
    #' @field d number of variables (dimensions in X)
    d = function() self$data$d,
    #' @field Q number of blocks
    Q = function(value) as.integer(ncol(private$C)),
    #' @field n_edges number of edges of the network (non null coefficient of the sparse precision matrix OmegaQ)
    n_edges  = function(value) sum(private$OmegaQ[upper.tri(private$OmegaQ, diag = FALSE)] != 0),
    #' @field model_par a list with the matrices of the model parameters: B (covariates), dm1 (species variance), OmegaQ (groups precision matrix))
    model_par = function(value) list(B = private$B, B0 = private$B0,
                                     dm1 = private$dm1, OmegaQ = private$OmegaQ),
    #' @field nb_param number of parameters in the model
    nb_param = function(value) {
      nb_param_D <- ifelse(private$res_covariance == "diagonal", self$p, 1)
      as.integer(self$p * self$d + self$Q + self$n_edges + nb_param_D)
    },
    #' @field objective evolution of the objective function during (V)EM algorithm
    objective = function() private$ll_list[-1],
    #' @field loglik (or its variational lower bound)
    loglik = function(value) if (private$approx) NA else private$ll_list[[length(private$ll_list)]] + self$sparsity_term,
    #' @field deviance (or its variational lower bound)
    deviance = function() -2 * self$loglik,
    #' @field BIC (or its variational lower bound)
    #' @field entropy Entropy of the conditional distribution when applicable
    entropy    = function() 0,
    BIC = function() self$deviance + log(self$n) * self$nb_param,
    #' @field ICL variational lower bound of the ICL
    ICL        = function() self$BIC + 2 * self$entropy,
    #' @field EBIC variational lower bound of the EBIC
    EBIC   = function(value) self$BIC + 2 * ifelse(self$n_edges > 0, self$n_edges * log(self$Q), 0),
    #' @field criteria a vector with loglik, BIC and number of parameters
    criteria   = function() {
      data.frame(nb_param = self$nb_param, Q = self$Q, n_edges = self$n_edges, sparsity = self$sparsity,
                 loglik = self$loglik, deviance = self$deviance, BIC = self$BIC, ICL = self$ICL, EBIC = self$EBIC)
    },
    #' @field sparsity (overall sparsity parameter)
    sparsity = function(value) {
      if (missing(value)) {
        private$sparsity_
      } else {
        stopifnot("must be a positive scale" = value >= 0)
        private$sparsity_ <- value
      }
    },
    #' @field sparsity_weights (weights associated to each pair of groups)
    sparsity_weights = function(value) {
      if (missing(value)) {
        private$weights
      } else {
        stopifnot("must be a Q x Q matrix" =
                    all(is.matrix(value), nrow(value) == ncol(value), ncol(value) == self$Q))
        private$weights <- value
      }
    },
    #' @field sparsity_term (sparsity_term term in log-likelihood due to sparsity)
    sparsity_term = function(value) self$sparsity * sum(abs(self$sparsity_weights * private$OmegaQ)),
    #' @field get_res_covariance whether the residual covariance is diagonal or spherical
    get_res_covariance = function(value) private$res_covariance,
    #' @field memberships cluster memberships
    memberships = function(value) private$C,
    #' @field clustering given as the list of elements contained in each cluster
    clustering = function(value) get_clusters(private$C),
    #' @field cluster_sizes given as a vector of cluster sizes
    cluster_sizes = function(value) table(self$clustering),
    #' @field elements_per_cluster given as the list of elements contained in each cluster
    elements_per_cluster = function(value) {
      if (is.null(names(self$clustering)))
        base::split(1:self$p, self$clustering)
      else
        base::split(names(self$clustering), self$clustering)
    }
  )
)
