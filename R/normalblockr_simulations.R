library(reshape2)
library(tidyverse)
library(viridis)
# variance matrices for simulations
Sigmas <- list()
for(Q in 1:20){
  if(Q == 5 | Q==8){
    Sigma <- matrix(unlist(toeplitz(lapply(0:(Q - 1), f<-function(x){0.75 ** x}))), nrow=Q, ncol=Q)
    Sigmas <- c(Sigmas, list(Sigma))
  }else{
    Sigma <- matrix(0, Q, Q)
    diag(Sigma) <- 1
    Sigma[row(Sigma) != col(Sigma)] = 0.1
    Sigmas <- c(Sigmas, list(Sigma))
  }
}


# code for one simulation
one_simulation <- function(n, p, d, Q, minX = 0, maxX = 1, minB = 0, maxB = 1,
                          minD = 0.5, maxD = 3, zero_inflation = FALSE,
                          minKappa = 0, maxKappa = 0.2){
  SigmaQ = Sigmas[[Q]]
  B <- matrix(rep(1, d * p), nrow = d)
  for (dim in 1:d) B[dim, ] <- runif(p, min = 0, max = 1)
  X <- matrix(rep(1, d * n), nrow = d)
  for (dim in 1:d) X[dim, ] <- runif(n, min = minX[[dim]], max = maxX[[dim]])
  C <- matrix(rep(0, p * Q), nrow = p)
  groups <- sample(1 : Q, size <- p, replace <- TRUE)
  for(dim in 1:p) C[dim, groups[[dim]]] <- 1
  D <- diag(runif(p, min <- minD, max <- maxD))
  dm1 <- 1 / diag(D)
  W <- t(MASS::mvrnorm(n, mu = matrix(rep(0, Q), Q, 1), Sigma = SigmaQ))
  epsilon <- t(MASS::mvrnorm(n, mu = matrix(rep(0, p), p, 1), Sigma = D))
  Y <- t(B) %*% X + C %*% W + epsilon
  if(zero_inflation){
    kappa <- runif(p, min = minKappa, max = maxKappa)
    Y[matrix(rbinom(n * p, size = 1, prob = rep(kappa, n)), nrow = p, ncol = n) == 1] <- 0
  }
  Y    <- t(Y)
  X    <- t(X)
  return(list(B = B, SigmaQ = SigmaQ, dm1 = dm1, Y = Y, C = C, X = X,
              kappa = kappa))
}


multiple_optimizations <- function(n_simus, n, p, d, Q, minX, maxX, minB, maxB,
                                   minD, maxD, zero_inflation, minKappa,
                                   maxKappa, fixed_blocks = FALSE,
                                   niter = 100, threshold = 1e-4){
  calc_time  = c()
  if(! fixed_blocks){clustering = list("true_blocks"=c(), "model_blocks"=c())}
  rmse_B     = c()
  rmse_dm1   = c()
  if(zero_inflation){rmse_kappa = c()}

  for(h in 1:n_simus){
    simu <- one_simulation(n, p, d, Q, minX = 0, maxX = 1, minB = 0, maxB = 1,
                           minD = 0.5, maxD = 3, zero_inflation = FALSE,
                           minKappa = 0, maxKappa = 0.2)
    Y = simu$Y
    X = simu$X
    B = simu$B
    C = simu$C
    if(fixed_blocks) {
      blocks = C
    }else{blocks = NULL}
    t0 = Sys.time()
    res_optim <- normal_block(Y, X, nb_blocks = Q, blocks = blocks,
                             zero_inflation = zero_inflation, niter = niter,
                             threshold = threshold)
    calc_time = c(calc_time, Sys.time() - t0)
    rmse_B = c(rmse_B, Metrics::rmse(res_optim$model_par$B, B))
    rmse_dm1 = c(rmse_dm1, Metrics::rmse(res_optim$model_par$dm1, simu$dm1))
    if(! fixed_blocks){
      clustering$true_blocks  = c(clustering$true_blocks, get_clusters(C))
      clustering$model_blocks = c(clustering$model_blocks, res_optim$clustering) }
    if(zero_inflation){rmse_kappa = c(rmse_kappa,
                                      Metrics::rmse(res_optim$model_par$kappa, simu$kappa))}

  }
  res <- list("calc_time" = calc_time, "rmse_B" = rmse_B, "rmse_dm1" = rmse_dm1)
  if(!fixed_blocks){finalRes[["clustering"]] = clustering}
  if(zero_inflation){finalRes[["rmse_kappa"]] = rmse_kappa}
  res
}


# grid of simulations over different sets of parameters
meta_grid <- function(ns, ps, Qs, d, minXs, maxXs, minDs, maxDs, minBs, maxBs,
                      fixed_blocks, zero_inflation, minKappas=c(0),
                      maxKappas=c(0), n_simus, niter, threshold,
                      verbose = FALSE){
  calc_time <- data.frame("n"=c(), "p"=c(), "Q" = c(), "minX"=c(), "maxX"=c(),
                         "minKappa"=c(), "maxKappa"=c(),"minD"=c(), "maxD"=c(),
                         "minB"=c(), "maxB"=c(),  "time"=c())
  rmse_B    <- data.frame("n"=c(), "p"=c(), "Q" = c(), "minX"=c(), "maxX"=c(),
                         "minKappa"=c(), "maxKappa"=c(),"minD"=c(),
                         "maxD"=c(), "minB"=c(), "maxB"=c(), "rmse"=c())
  rmse_dm1  <- data.frame("n"=c(), "p"=c(), "Q" = c(), "minX"=c(), "maxX"=c(),
                         "minKappa"=c(), "maxKappa"=c(), "minD"=c(),
                         "maxD"=c(), "minB"=c(), "maxB"=c(), "rmse"=c())

  if(!fixed_blocks){
    clustering <- data.frame("n"=c(), "p"=c(), "Q" = c(), "minD"=c(),
                             "maxD"=c(), "minB"=c(), "maxB"=c(), "ARI"=c())
  }
  if(zero_inflation){
    rmse_Kappa <- data.frame("n"=c(), "p"=c(), "Q" = c(), "minX"=c(), "maxX"=c(),
                             "minKappa"=c(), "maxKappa"=c(), "minD"=c(),
                             "maxD"=c(), "minB"=c(), "maxB"=c(),"rmse"=c())
  }

  for(n in ns){
    if(verbose){cat("---------","\n", "n : ", n,"\n", "---------","\n")}
    for(p in ps){
      if(n >= p){
        if(verbose){cat("---------","\n", "p : ", p,"\n", "---------","\n")}
        for(Qind in seq_along(Qs)){
          Q <- Qs[[Qind]]
          if(p >= 2 * Q){
            if(verbose){cat("---------","\n", "Q : ", Q,"\n", "---------","\n")}
            for(minXind in seq_along(minXs)){
              minX = minXs[[minXind]] ; maxX = maxXs[[minXind]]
              if(verbose){cat("---------","\n", minX, " < X < ", maxX,"\n", "---------","\n")}
              for(minDind in seq_along(minDs)){
                minD = minDs[[minDind]] ; maxD = maxDs[[minDind]]
                for(minBind in seq_along(minBs)){
                  minB = minBs[[minBind]] ; maxB = maxBs[[minBind]]
                  for(minKappaind in seq_along(minKappas)){
                      minKappa = minKappas[[minKappaind]] ; maxKappa = maxKappas[[minKappaind]]

                      res <- multiple_optimizations(n_simus, n, p, d, Q, minX,
                                                    maxX, minB, maxB, minD,
                                                    maxD, zero_inflation,
                                                    minKappa, maxKappa,
                                                    fixed_blocks, niter,
                                                    threshold)


                    this_calc_time = data.frame("Q" = rep(Q, n_simus),
                                                "n"=rep(n, n_simus),
                                                "p"=rep(p, n_simus),
                                                "minX"=rep(minX, n_simus),
                                                "maxX"=rep(maxX, n_simus),
                                                "minKappa"=rep(minKappa, n_simus),
                                                "maxKappa"=rep(maxKappa, n_simus),
                                                "minD"=rep(minD, n_simus),
                                                "maxD"=rep(maxD, n_simus),
                                                "minB"=rep(minB, n_simus),
                                                "maxB"=rep(maxB, n_simus),
                                                "time"=res$calc_time)
                      calc_time = rbind(calc_time, this_calc_time)

                      if(!fixed_blocks){
                        this_clustering = data.frame("Q" = rep(Q, n_simus),
                                                    "n"=rep(n, n_simus),
                                                    "p"=rep(p, n_simus),
                                                    "minX"=rep(minX, n_simus),
                                                    "maxX"=rep(maxX, n_simus),
                                                    "minKappa"=rep(minKappa, n_simus),
                                                    "maxKappa"=rep(maxKappa, n_simus),
                                                    "minD"=rep(minD, n_simus),
                                                    "maxD"=rep(maxD, n_simus),
                                                    "minB"=rep(minB, n_simus),
                                                    "maxB"=rep(maxB, n_simus),
                                                    "ARI"=c(sapply(1:n_simus, function(i){matchingGroupScores(res$clustering$true_blocks[[i]], res$clustering$model_blocks[[i]])})))
                        clustering = rbind(clustering, this_clustering)
                      }


                      this_rmse_B = data.frame("Q" = rep(Q, n_simus),
                                               "n"=rep(n, n_simus),
                                               "p"=rep(p, n_simus),
                                               "minX"=rep(minX, n_simus),
                                               "maxX"=rep(maxX, n_simus),
                                               "minKappa"=rep(minKappa, n_simus),
                                               "maxKappa"=rep(maxKappa, n_simus),
                                               "minD"=rep(minD, n_simus),
                                               "maxD"=rep(maxD, n_simus),
                                               "minB"=rep(minB, n_simus),
                                               "maxB"=rep(maxB, n_simus),
                                               "rmse"= res$rmse_B)
                      rmse_B = rbind(rmse_B, this_rmse_B)
                      this_rmse_dm1 = data.frame("Q" = rep(Q, n_simus),
                                                 "n"=rep(n, n_simus),
                                                 "p"=rep(p, n_simus),
                                                 "minX"=rep(minX, n_simus),
                                                 "maxX"=rep(maxX, n_simus),
                                                 "minKappa"=rep(minKappa, n_simus),
                                                 "maxKappa"=rep(maxKappa, n_simus),
                                                 "minD"=rep(minD, n_simus),
                                                 "maxD"=rep(maxD, n_simus),
                                                 "minB"=rep(minB, n_simus),
                                                 "maxB"=rep(maxB, n_simus),
                                                 "rmse"= res$rmse_dm1)
                      rmse_dm1 = rbind(rmse_dm1, this_rmse_dm1)
                      if(zero_inflation){
                        this_rmse_kappa= data.frame("Q" = rep(Q, n_simus),
                                                    "n"=rep(n, n_simus),
                                                    "p"=rep(p, n_simus),
                                                    "minX"=rep(minX, n_simus),
                                                    "maxX"=rep(maxX, n_simus),
                                                    "minKappa"=rep(minKappa, n_simus),
                                                    "maxKappa"=rep(maxKappa, n_simus),
                                                    "minD"=rep(minD, n_simus),
                                                    "maxD"=rep(maxD, n_simus),
                                                    "minB"=rep(minB, n_simus),
                                                    "maxB"=rep(maxB, n_simus),
                                                    "rmse"= res$rmse_kappa)
                        rmse_kappa = rbind(rmse_kappa, this_rmse_kappa)
                      }
                    }
                }
              }
            }
          }
        }
      }
    }
  }
  res = list("calc_time"=calc_time, "rmse_B"=rmse_B, "rmse_dm1"=rmse_dm1)
  if(!fixed_blocks){res[["clustering"]] = clustering}
  if(zero_inflation){res[["rmse_kappa"]] = rmse_kappa}
  res
}


# heatmaps for each pair of parameters
pair_heatmap <- function(param1, param2, df){
  indicator <- colnames(df)[[length(colnames(df))]]
  agg_df <- df %>%
    group_by(!!sym(param1), !!sym(param2)) %>%
    summarize(!!indicator := mean(!!sym(indicator)))

  # Step 2: Create the formula dynamically
  formula <- as.formula(paste(param1, "~", param2))

  # Reshape the data into a matrix
  agg_matrix <- acast(agg_df, formula, value.var = indicator, fill = NA_real_)
  agg_long <- as.data.frame(as.table(agg_matrix))
  colnames(agg_long) <- c(param1, param2, indicator)

  ggplot(agg_long, aes_string(x = param2, y = param1, fill = indicator)) +
    geom_tile() +
    scale_fill_viridis(option = "viridis") + # Using viridis color scale
    theme_minimal()+
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(), # To remove major grid lines
      panel.grid.minor = element_blank()  # To remove minor grid lines
    )
}

# saving meta_grid results
save_meta_grid <- function(sub_dir, rank_simu, res, ns, ps, Qs, minXs, maxXs,
                          minDs, maxDs, minBs, maxBs,  fixed_blocks,
                          zero_inflation, minKappas, maxKappas,
                          n_simus, parameters = c("n", "Q", "p", "minD")){
  main_dir <- getwd()
  save_dir <- paste0(sub_dir, "/simulation_", as.character(rank_simu))
  if (file.exists(save_dir)){
    print("Given rank simu already used - no saving done.")
    return()
  } else {
    dir.create(file.path(main_dir, save_dir))
    setwd(file.path(main_dir, save_dir))
  }

  # Saving results data frames as csv
  write.csv(res$calc_time, paste("simu_", as.character(rank_simu), "_calc_time.csv"))
  if(!fixed_blocks){write.csv(res$clustering, paste("simu_", as.character(rank_simu), "_clustering.csv"))}
  write.csv(res$rmse_B, paste("simu_", as.character(rank_simu), "_rmse_b.csv"))
  write.csv(res$rmse_dm1, paste("simu_", as.character(rank_simu), "_rmse_dm1.csv"))
  if(zero_inflation){write.csv(res$rmse_kappa, paste("simu_", as.character(rank_simu), "_rmse_kappa.csv"))}

  # Saving ARI and errors heatmaps
  for(param_ind1 in (1:(length(parameters) - 1))){
    param1 = parameters[[param_ind1]]
    for(param_ind2 in ((param_ind1 + 1):length(parameters))){
      param2 = parameters[[param_ind2]]

      if(!fixed_blocks){
        df <- res$clustering
        H_clustering <- pair_heatmap(param1, param2, df)
        ggsave(paste0("metagrid_simu_", param1, "_", param2, "_simu_", as.character(rank_simu), "_ari.jpg"),
               H_clustering, width=21, height=14, units="cm")
      }

      df <- res$rmse_B
      H_rmse_B <- pair_heatmap(param1, param2, df)
      ggsave(paste0("metagrid_simu_", param1, "_", param2, "_simu_", as.character(rank_simu), "_rmse_B.jpg"),
             H_rmse_B, width=21, height=14, units="cm")

      df <- res$rmse_dm1
      H_rmse_dm1 <- pair_heatmap(param1, param2, df)
      ggsave(paste0("metagrid_simu_", param1, "_", param2, "_simu_", as.character(rank_simu), "_rmse_dm1.jpg"),
             H_rmse_dm1, width=21, height=14, units="cm")

      if(zero_inflation){
        df <- res$rmse_kappa
        Hrmse_dm1 <- pair_heatmap(param1, param2, df)
        ggsave(paste0("metagrid_simu_", param1, "_", param2, "_simu_", as.character(rank_simu), "_rmse_kappa.jpg"),
               Hrmse_dm1, width=21, height=14, units="cm")
      }

    }
  }

  setwd(main_dir)

  # Saving simulation parameters
  param_filename = paste0(save_dir, "/simu_", as.character(rank_simu), "_parameters.txt")
  file_connexion <- file(param_filename, open = "w")
  cat(paste0("n = ", paste(ns, collapse=" ; "), "\n"), file = file_connexion, append = TRUE)
  cat(paste0("p = ", paste(ps, collapse=" ; "), "\n"), file = file_connexion, append = TRUE)
  cat(paste0("Q = ", paste(Qs, collapse=" ; "), "\n"), file = file_connexion, append = TRUE)
  cat(paste0("minD = ", paste(minDs, collapse=" ; "), "\n"), file = file_connexion, append = TRUE)
  cat(paste0("maxD = ", paste(maxDs, collapse=" ; "), "\n"), file = file_connexion, append = TRUE)
  cat(paste0("minB = ", paste(minBs, collapse=" ; "), "\n"), file = file_connexion, append = TRUE)
  cat(paste0("maxB = ", paste(maxBs, collapse=" ; "), "\n"), file = file_connexion, append = TRUE)
  cat(paste0("minX = ", paste(minXs, collapse=" ; "), "\n"), file = file_connexion, append = TRUE)
  cat(paste0("maxX = ", paste(maxXs, collapse=" ; "), "\n"), file = file_connexion, append = TRUE)

  cat(" ", file = file_connexion, append = TRUE)
  cat("Sigmas \n", file = file_connexion, append = TRUE)
  for(Q in Qs){
    Sigma = Sigmas[Q]
    cat(paste0("Sigma", as.character(Q), " = ", as.character(Sigma), "\n"), file = file_connexion, append = TRUE)
  }
  close(file_connexion)
}



ns             = c(200, 300, 500)
ps             = c(50, 100, 200 )
Qs             = c(3, 4, 6)
d              = 1
minXs          = c(0, 0, 0)
maxXs          = c(1, 10, 100)
minDs          = c(0.5, 1, 2)
maxDs          = c(1, 2, 4)
minBs          = c(2)
maxBs          = c(4)
fixed_blocks   = TRUE
zero_inflation = FALSE
minKappas      = c(0)
maxKappas      = c(0)
n_simus        = 2
niter          = 75
threshold      = 1e-5
verbose        = TRUE


meta_grid_res <- meta_grid(ns, ps, Qs, d, minXs, maxXs, minDs, maxDs, minBs,
                           maxBs, fixed_blocks, zero_inflation,
                           minKappas, maxKappas, n_simus, niter, threshold,
                           verbose)

sub_dir = "simulations_results/simulations_fixed_blocks"
parameters = c("n", "Q", "p", "minD", "maxX")
rank_simu = 1
save_meta_grid(sub_dir, rank_simu, meta_grid_res, ns, ps, Qs, minXs, maxXs, minDs,
               maxDs, minBs, maxBs,  fixed_blocks, zero_inflation, minKappas,
               maxKappas, n_simus, parameters)
