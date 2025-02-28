###############################################################################
###############################################################################
## Use pre-save testdata (seed are hard to handle in testhat)
testdata <- readRDS("testdata/testdata_normal_zi.RDS")
Y <- testdata$Y ; X <- testdata$X
C <- testdata$parameters$C ; Q <- ncol(C)

###############################################################################

test_that("normal block with diagonal residual covariance and unknown Q", {
  clustering_init_3 <- get_clusters(C)
  clustering_init_2 <- clustering_init_3 ; clustering_init_2[which(clustering_init_2 == 3) ] = 2
  clustering_init_1 <- clustering_init_2 ; clustering_init_1[which(clustering_init_1 == 2) ] = 1
  clustering_init   <- list(clustering_init_1, clustering_init_2, clustering_init_3)

  data  <- normal_data$new(Y, X)
  model <- NB_unknown_Q$new(data, c(1, 2, 3),zero_inflation = FALSE,
                            control = NB_control(clustering_init = clustering_init))
  model$optimize()
  model_BIC <- model$get_best_model("BIC")
  expect_lt(model_BIC$BIC, 7455)
})

test_that("normal block with unknown Q, heuristic", {
  data  <- normal_data$new(Y, X)
  model <- NB_unknown_Q$new(data, c(2, 3, 4), penalty = 0.05,
                            control = NB_control(heuristic = TRUE))
  model$optimize()
  model_3 <- model$get_model(3)
  expect_lt(Metrics::rmse(model_3$fitted, Y), 2.9)
})
