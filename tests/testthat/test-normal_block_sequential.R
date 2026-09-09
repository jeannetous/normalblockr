###############################################################################
## normal_block_sequential() fits a mean-block model, then a variance-block one
## on its residuals. The positive control below carries two genuinely distinct
## structures (unrelated partitions in the mean and in the covariance); the
## estimator must recover both. See
## inst/mean_block_analyses/sequential_mean_then_variance.R.

test_that("normal_block_sequential() recovers two distinct structures", {
  skip_if_not_installed("aricode")
  set.seed(42)
  n <- 300; p <- 100; d <- 3; q_mean <- 5; q_var <- 4
  X  <- cbind(1, matrix(rnorm(n * (d - 1)), n, d - 1))
  Cm <- normalblockr:::as_indicator(sample(rep(1:q_mean, length.out = p)))
  Cv <- normalblockr:::as_indicator(sample(rep(1:q_var,  length.out = p)))
  Y  <- X %*% matrix(rnorm(d * q_mean, sd = 2), d, q_mean) %*% t(Cm) +
        matrix(rnorm(n * q_var, sd = 1.5), n, q_var) %*% t(Cv) +
        matrix(rnorm(n * p, sd = 0.5), n, p)

  fit <- normal_block_sequential(NormalBlockData$new(Y, X),
                                 blocks_mean = seq(3, 15, by = 2), blocks_var = 1:8)

  expect_s3_class(fit, "normal_block_sequential")
  expect_equal(fit$mean$q, q_mean)
  expect_equal(fit$var$q,  q_var)
  expect_gt(aricode::ARI(fit$mean$clustering, normalblockr:::get_clusters(Cm)), 0.95)
  expect_gt(aricode::ARI(fit$var$clustering,  normalblockr:::get_clusters(Cv)), 0.95)
  ## the two structures are unrelated by construction, and must stay so
  expect_lt(aricode::ARI(fit$mean$clustering, fit$var$clustering), 0.2)
  expect_equal(dim(fit$residuals), dim(Y))
})

test_that("normal_block_sequential() accepts fixed clusterings and prints", {
  testdata <- readRDS("testdata/testdata_normal_mean_block.RDS")
  data <- NormalBlockData$new(testdata$Y, testdata$X)
  C    <- testdata$parameters$C

  fit <- normal_block_sequential(data, blocks_mean = C, blocks_var = 2)
  expect_s3_class(fit$mean, "NormalBlockMeanKnownClusters")
  expect_equal(fit$var$q, 2)
  expect_output(print(fit), "sequential mean-then-variance")
  expect_invisible(print(fit))
})

test_that("an explicit zeros mask overrides the default Y == 0", {
  ex <- generate_normal_block_mean_data(n = 40, p = 10, d = 1, q = 2)
  mask <- matrix(0, 40, 10); mask[1:5, 1:2] <- 1

  default  <- NormalBlockData$new(ex$Y, ex$X)
  explicit <- NormalBlockData$new(ex$Y, ex$X, zeros = mask)

  expect_equal(sum(default$zeros), 0)
  expect_equal(explicit$zeros, mask)
  expect_equal(explicit$zeros_bar, 1 - mask)
  expect_equal(explicit$npY, 40 * 10 - 10)
  expect_equal(explicit$nY, colSums(1 - mask))

  expect_error(NormalBlockData$new(ex$Y, ex$X, zeros = mask[, -1]),
               "same dimensions")
})

test_that("the zero-inflated chain runs and shares one mask across both stages", {
  set.seed(7)
  ex <- generate_normal_block_mean_data(n = 80, p = 20, d = 1, q = 3)
  Y  <- ex$Y; Y[runif(length(Y)) < 0.2] <- 0
  data <- NormalBlockData$new(Y, ex$X)

  fit <- normal_block_sequential(data, blocks_mean = 3, blocks_var = 2,
                                 zero_inflation = TRUE)

  expect_s3_class(fit, "normal_block_sequential")
  expect_true(inherits(fit$mean, "ZINormalBlockMeanUnknownClusters"))
  expect_true(inherits(fit$var, "ZINormalBlockVarUnknownClusters"))
  ## the residuals carry no zero of their own, so only the inherited mask can
  ## have told the second stage where the structural zeros are
  expect_equal(fit$var$data$zeros, data$zeros)
  expect_true(all(fit$residuals[data$zeros == 1] == 0))
})
