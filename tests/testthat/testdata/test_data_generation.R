set.seed(3)

testdata <- generate_normal_block_data(
  n = 80,
  p = 20,
  d = 1,
  Q = 3
)

saveRDS(testdata, file = "tests/testthat/testdata/testdata_normal.RDS")

testdata <- generate_normal_block_data(
  n = 80,
  p = 20,
  d = 1,
  Q = 3,
  kappa = runif(20, min = .2, max = .6)
)

saveRDS(testdata, file = "tests/testthat/testdata/testdata_normal_zi.RDS")
