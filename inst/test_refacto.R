testdata <- readRDS("tests/testthat/testdata/testdata_normal.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C ; Q <- ncol(C)

data  <- normalblockr:::normal_data$new(Y, X)


model_fixed_block <- normalblockr:::NB_fixed_blocks_diagonal$new(data, C)
model_fixed_block$optimize()
model_fixed_block <- normalblockr:::NB_fixed_blocks_diagonal$new(data, C, penalty = 0.05)
model_fixed_block$optimize()
model_fixed_block <- normalblockr:::NB_fixed_blocks_diagonal$new(data, C, penalty = 0.05,
                                                                 control = NB_control(heuristic = TRUE))
model_fixed_block$optimize()

model_unkwn_block <- normalblockr:::NB_fixed_Q_diagonal$new(data, Q)
model_unkwn_block$optimize()
model_unkwn_block <- normalblockr:::NB_fixed_Q_diagonal$new(data, Q, penalty = 0.05)
model_unkwn_block$optimize()
model_unkwn_block <- normalblockr:::NB_fixed_Q_diagonal$new(data, Q, penalty = 0.05,
                                                            control = NB_control(heuristic = TRUE))
model_unkwn_block$optimize()


