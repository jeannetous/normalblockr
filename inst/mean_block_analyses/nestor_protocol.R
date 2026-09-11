## ---------------------------------------------------------------------------
## Nestor Ngalala Manguitini's M2 simulation protocol, run on the package's
## mean-block implementation.
##
## The point is NOT to reproduce his numbers. His prototype
## (inst/normalblockmean/code_nestor) and this package's implementation are
## different code; where they disagree the package is the one that is tested,
## profiled and validated against its own C++/R equivalence suite. The point is
## to put our estimator through the same experiment, so that its behaviour can
## be read against the report.
##
## What is taken from his code, verbatim: the data-generating process
## (01_generation_donnees.R -- X Gaussian, B random, every cluster non-empty,
## Sigma from an Erdos-Renyi precision matrix rescaled to a target SNR, and the
## exp/log round trip). What is ours: the estimator, the initialisations, the
## criteria.
##
## Deliberate departures from his settings, each for a reason:
##
##  * `scale = FALSE` on NormalBlockData. He does not standardise, and the
##    RMSEs below compare to true parameters on their original scale, which
##    column-scaling would invalidate.
##  * Initialisations are ours (kmeans/ward2/spectral/sbm), not his five
##    (kmeans/clustofvar/hclust/gmm/pam). ClustOfVar was benchmarked worst of
##    all and dropped from this package, taking its dependency with it; hclust
##    is ward2 here. Only "kmeans" is common to both lists by name.
##  * nu = d*q + p(p+1)/2 + (q-1) matches his `calculer_nu()` exactly, as do
##    BIC = -2*ELBO + nu*log(n) and ICL = BIC + 2*entropy -- but our nb_param
##    counts the *non-zero* off-diagonal terms of Omega, so under a penalty it
##    is smaller than his fixed p(p+1)/2. Reported below.
##
## Run from the package root with devtools::load_all() already done, or just
## source it: it loads the package itself.
## ---------------------------------------------------------------------------

suppressMessages({
  library(igraph); library(MASS); library(aricode); library(dplyr)
})
if (!"normalblockr" %in% loadedNamespaces()) suppressMessages(devtools::load_all(quiet = TRUE))
source("inst/normalblockmean/code_nestor/01_generation_donnees.R")

OUT <- "inst/mean_block_analyses/nestor_protocol_results.rds"

## his grid, verbatim (parametres.R section 5)
n_list   <- c(80, 150, 220, 300)
p_list   <- 60
q_list   <- c(3, 5)
d_list   <- 2
snr_list <- c(0.3, 1.5)
n_simu   <- 10
lambda   <- 0.1                # his lambda_glasso, diagonal unpenalised
q_select <- 2:6                # his q_selection_list
inits    <- c("kmeans", "ward2", "spectral", "sbm")

rmse <- function(a, b) sqrt(mean((a - b)^2))

## Every metric he computes, plus a note on which are permutation-invariant.
## Cluster labels are only identified up to a permutation, so RMSE(B) and
## RMSE(alpha) compare quantities whose columns may be in any order: they are
## reported because the report reports them, but RMSE(BC') is the one that
## means anything.
metrics_for <- function(fit, truth, n_warn) {
  Sigma_hat <- solve(fit$model_par$Omega)
  data.frame(
    ARI        = aricode::ARI(fit$clustering, truth$cl),
    RMSE_BCt   = rmse(truth$B %*% t(truth$C), fit$model_par$B %*% t(fit$memberships)),
    RMSE_Sigma = rmse(truth$Sigma, Sigma_hat),
    RMSE_B     = rmse(truth$B, fit$model_par$B),                       # label-dependent
    RMSE_alpha = rmse(colMeans(truth$C), colMeans(fit$memberships)),   # label-dependent
    ELBO       = fit$loglik,
    BIC        = fit$BIC,
    ICL        = fit$ICL,
    nb_param   = fit$nb_param,
    niter      = length(fit$objective),
    glasso_warn = n_warn
  )
}

quiet_fit <- function(expr) {
  n <- 0
  val <- withCallingHandlers(expr, warning = function(w) { n <<- n + 1; invokeRestart("muffleWarning") })
  list(fit = val, n_warn = n)
}

grid <- expand.grid(n = n_list, p = p_list, q = q_list, d = d_list,
                    snr = snr_list, simu = seq_len(n_simu),
                    KEEP.OUT.ATTRS = FALSE)
cat(sprintf("%d configurations x %d initialisations, plus q selection over %s\n",
            nrow(grid), length(inits), paste(range(q_select), collapse = ":")))

set.seed(42)                                   # his graine_simulations
rows_unknown <- list(); rows_known <- list(); rows_qsel <- list()
t0 <- Sys.time()

for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  dd <- generer_donnees(n = g$n, p = g$p, d = g$d, q = g$q, snr = g$snr)
  truth <- list(B = dd$B_vrai, C = dd$C, Sigma = dd$sigma_vrai,
                cl = apply(dd$C, 1, which.max))
  ## he fits log(Y), which is exactly the Gaussian Z he generated
  data <- NormalBlockData$new(dd$Z, dd$X, scale = FALSE)

  ## --- unknown clustering, one fit per initialisation ---------------------
  for (init in inits) {
    r <- quiet_fit(normal_block(data, blocks = g$q, sparsity = lambda, model = "mean",
                                control = NB_control(verbose = FALSE, clustering_init = init)))
    rows_unknown[[length(rows_unknown) + 1]] <-
      cbind(g, init = init, metrics_for(r$fit, truth, r$n_warn),
            snr_effectif = dd$snr_effectif)
  }

  ## --- known clustering (his "clustering observe") -------------------------
  r <- quiet_fit(normal_block(data, blocks = dd$C, sparsity = lambda, model = "mean",
                              control = NB_control(verbose = FALSE)))
  rows_known[[length(rows_known) + 1]] <- cbind(
    g,
    RMSE_B     = rmse(truth$B, r$fit$model_par$B),
    RMSE_Sigma = rmse(truth$Sigma, solve(r$fit$model_par$Omega)),
    loglik     = r$fit$loglik,
    glasso_warn = r$n_warn)

  ## --- automatic selection of q by BIC and ICL -----------------------------
  r <- quiet_fit(normal_block(data, blocks = q_select, sparsity = lambda, model = "mean",
                              control = NB_control(verbose = FALSE)))
  rows_qsel[[length(rows_qsel) + 1]] <- cbind(
    g,
    q_BIC = r$fit$get_best_model("BIC")$q,
    q_ICL = r$fit$get_best_model("ICL")$q,
    glasso_warn = r$n_warn)

  if (i %% 20 == 0)
    cat(sprintf("  %3d/%d  (%.1f min)\n", i, nrow(grid),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

unknown <- bind_rows(rows_unknown)
known   <- bind_rows(rows_known)
qsel    <- bind_rows(rows_qsel)
saveRDS(list(unknown = unknown, known = known, qsel = qsel,
             settings = list(grid = grid, lambda = lambda, inits = inits,
                             q_select = q_select)), OUT)

cat(sprintf("\ndone in %.1f min -> %s\n\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), OUT))

## ---------------------------------------------------------------------------
## Summary
## ---------------------------------------------------------------------------
cat("== unknown clustering: ARI by initialisation and SNR ==\n")
print(unknown |>
        group_by(snr, init) |>
        summarise(ARI = mean(ARI), RMSE_BCt = mean(RMSE_BCt), .groups = "drop") |>
        arrange(snr, desc(ARI)) |> as.data.frame(), digits = 3)

cat("\n== unknown clustering: does more data help? (kmeans, the default) ==\n")
print(unknown |> filter(init == "kmeans") |>
        group_by(snr, q, n) |>
        summarise(ARI = mean(ARI), RMSE_BCt = mean(RMSE_BCt),
                  RMSE_Sigma = mean(RMSE_Sigma), .groups = "drop") |>
        as.data.frame(), digits = 3)

cat("\n== known clustering: RMSE against n ==\n")
print(known |> group_by(snr, q, n) |>
        summarise(RMSE_B = mean(RMSE_B), RMSE_Sigma = mean(RMSE_Sigma), .groups = "drop") |>
        as.data.frame(), digits = 3)

cat("\n== selection of q (true q in the row label) ==\n")
print(qsel |> group_by(snr, q) |>
        summarise(pct_BIC_exact = 100 * mean(q_BIC == q),
                  pct_ICL_exact = 100 * mean(q_ICL == q),
                  median_q_BIC = median(q_BIC), median_q_ICL = median(q_ICL),
                  .groups = "drop") |> as.data.frame(), digits = 3)

cat(sprintf("\nglasso fallbacks: %d of %d penalised fits\n",
            sum(unknown$glasso_warn) + sum(known$glasso_warn) + sum(qsel$glasso_warn),
            nrow(unknown) + nrow(known) + nrow(qsel)))
