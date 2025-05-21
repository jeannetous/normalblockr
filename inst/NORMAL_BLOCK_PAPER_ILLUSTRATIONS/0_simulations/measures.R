library(combinat)

#'Function to compute ARI between two clusters
#'@param true_clusters true clustering to compare inferred clustering to
#'@param inferred_clusters inferred clustering to compare the true clustering to
ARI <- function(true_clusters, inferred_clusters){
  pdfCluster::adj.rand.index(true_clusters, inferred_clusters)}

#' Function to compute the RMSE between two vectors / matrices
#' @param x vector / matrix to compare to x
#' @param y vector / matrix to compare to y
rmse <- function(x, y){
  return(Metrics::rmse(x, y))
}

#' Function to compute the RMSE between true and inferred Omega
#' @param omega_true true value of omega
#' @param omega_hat estimated value of omega
#' @param permutation permutation to apply to omega_hat before computing the rmse
rmse_omega <- function(omega_true, omega_hat, permutation = NULL) {
  if (is.null(permutation)) {
    res <- rmse(omega_true, omega_hat)
  } else if (anyNA(permutation)) {
    res <- NA
  } else {
    res <- rmse(omega_true, omega_hat[permutation,permutation])
  }
  res
}

#' Function to compute recall fallout and precision in network inference
#' @param omega_true true value of omega
#' @param omega_hat estimated value of omega
#' @param permutation permutation to apply to omega_hat before computing the rmse
roc_metrics <- function(omega_true, omega_estimate, permutation = NULL){
  if (!is.null(permutation)) omega_estimate <- omega_estimate[permutation, permutation]

  diag(omega_true) <- 0 ; diag(omega_estimate) <- 0

  true.nzero <- which(omega_true != 0)
  true.zero  <- which(omega_true == 0)

  nzero <- which(omega_estimate != 0)
  zero  <- which(omega_estimate == 0)

  TP <- sum(nzero %in% true.nzero)
  TN <- sum(zero %in%  true.zero) - nrow(omega_true) # removing diagonal values that do not count
  FP <- sum(nzero %in% true.zero)
  FN <- sum(zero %in%  true.nzero)

  recall    <- TP/(TP + FN) ## also recall and sensitivity
  fallout   <- FP/(FP + TN) ## also 1 - specificity
  precision <- TP/(TP + FP) ## also PPR
  recall[TP + FN == 0] <- NA
  fallout[TN + FP == 0] <- NA
  precision[TP + FP == 0] <- NA

  res <-  round(c(fallout,recall,precision),3)
  res[is.nan(res)] <- 0
  names(res) <- c("fallout","recall", "precision")

  return(res)
}

#' Function to compute the Area Under the Curve from recall and fallout lists
#' @param recall (True Positive) / (Number of Positive)
#' @param fallout (True Negative) / (Number of Negative)
auc <- function(recall, fallout){
  return(sum(diff(fallout) * (recall[-1] + recall[-length(recall)]) / 2))
}

#' Function to compute the Area Under the Curve from true omega value and Normal-block model
#' @param omega_true true value of omega
#' @param NB_sparse sparse Normal-Block model (collection of models with different penalties)
#' @param permutation permutation to apply to omega_hat before comparing it with the true omega
get_auc <- function(omega_true, NB_sparse, permutation = NULL){
  fallout <- c() ; recall <- c()
  for(pen in NB_sparse$penalties){
    omega_estimate <- NB_sparse$get_model(pen)$model_par$omega
    res <- roc_metrics(omega_true, omega_estimate, permutation)
   if(!is.na(res[["fallout"]]) && !is.na(res[["recall"]])){
     fallout <- c(fallout, res[["fallout"]]) ; recall <- c(recall, res[["recall"]])
   }
  }
  recall <- rev(recall) ; fallout <- rev(fallout)
  # One value of fallout may correspond to different recall values depending on the penalty
  fallout_unique <- unique(fallout) ; recall_unique <- c()
  for(x in fallout_unique){
    recall_unique <- c(recall_unique, max(recall[which(fallout == x)]))
  }
  return(auc(recall_unique, fallout_unique))
}

#' Function to compute different measures to assess the quality of the inference by a given Normal-Block model
#' @param NB_object Normal-Block object to assess
#' @param data NBData object
#' @param param list of true parameters to compare the results to
#' @param model_selection criterion for model selection if applicable
#' @param observed_blocks boolean, whether the clustering is observed or not
#' @param heuristic boolean, whether a heuristic method was used for the inference
get_measures <- function(NB_object, param, data, model_selection = "None",
                         observed_blocks = FALSE, heuristic = FALSE) {
  # Select best sparsity level according to the chosen criterion
  if(heuristic){
    model <- NB_object
  }else{model <- NB_object$get_best_model(model_selection)}

  # Boolean indicating wether blocks are known or estimated
  observed_blocks <- (inherits(model, "NB_fixed_blocks") | observed_blocks)

  # Get best permutation of Omega according to rmse when possible
  omega_hat <- model$model_par$OmegaQ
  if (observed_blocks) {
    best_perm <- NULL
  } else if (model$Q < 7) {
    perms <- permn(model$Q)
    ibest <- perms %>%
      purrr::map_dbl(\(P) rmse(param$Omega, omega_hat[P,P])) %>%
      which.min()
    best_perm <- perms[[ibest]]
  } else {
    best_perm <- NA
  }

  ## get metrics
  res <- as.data.frame(c(list(
    criterion = model_selection,
    fixed_blocks = observed_blocks,
    AUC = get_auc(param$Omega, NB_object, best_perm),
    ARI = ARI(get_clusters(param$C), model$clustering),
    rmse_B = rmse(param$B, model$model_par$B),
    rmse_D = rmse(diag(param$D), 1/model$model_par$dm1),
    rmse_kappa = ifelse(is.null(model$model_par$kappa), NA, rmse(param$kappa, model$model_par$kappa)),
    rmse_fit = rmse(model$fitted, data$Y),
    rmse_omega = rmse_omega(param$Omega, omega_hat, best_perm)
  ), roc_metrics(param$Omega, omega_hat, best_perm)))
  res
}
