###############################################################################
###############################################################################

testdata <- readRDS("testdata/testdata_normal.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C ; Q <- ncol(C)
data <- NB_data$new(Y, X)

test_that("Robustness of starting clustering with NB", {
  ## Diagonal model
  model_kmeans <- NB_fixed_Q$new(data, Q, control = NB_control(clustering_approx = "kmeans"))
  model_kmeans$optimize()

  model_ward2 <- NB_fixed_Q$new(data, Q, control = NB_control(clustering_approx = "ward2"))
  model_ward2$optimize()

  model_sbm <- NB_fixed_Q$new(data, Q, control = NB_control(clustering_approx = "sbm"))
  model_sbm$optimize()

  expect_gt(model_kmeans$loglik, -2650)
  expect_gt(model_ward2$loglik, -2650)
  expect_gt(model_sbm$loglik, -2650)

  expect_gte(model_kmeans$loglik, model_sbm$loglik)
  expect_gte(model_ward2$loglik, model_sbm$loglik)

  expect_gt(aricode::ARI(model_kmeans$clustering, model_sbm$clustering), 0.75)
  expect_gt(aricode::ARI(model_kmeans$clustering, model_ward2$clustering), 0.75)
})


testdata <- readRDS("testdata/testdata_normal_zi.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C ; Q <- ncol(C)
data <- NB_data$new(Y, X)

test_that("Robustness of starting clustering with ZINB", {
  ## Diagonal model
  model_ward2 <- ZINB_fixed_Q$new(data, Q, control = NB_control(clustering_approx = "ward2"))
  model_ward2$optimize()

  model_kmeans <- ZINB_fixed_Q$new(data, Q, control = NB_control(clustering_approx = "kmeans"))
  model_kmeans$optimize()

  model_sbm <- ZINB_fixed_Q$new(data, Q, control = NB_control(clustering_approx = "sbm"))
  model_sbm$optimize()

  expect_gt(model_kmeans$loglik, -2750)
  expect_gt(model_ward2$loglik, -2750)
  expect_gt(model_sbm$loglik, -2750)

  expect_gte(model_kmeans$loglik, model_sbm$loglik)
  expect_gte(model_ward2$loglik, model_sbm$loglik)

  expect_gte(aricode::ARI(model_kmeans$clustering, model_sbm$clustering)  , 0.5)
  expect_gte(aricode::ARI(model_kmeans$clustering, model_ward2$clustering), 0.99)
})
