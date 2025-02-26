testdata <- readRDS("tests/testthat/testdata/testdata_normal.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C ; Q <- ncol(C)

data  <- normalblockr:::normal_data$new(Y, X)

## =======================================================
## FIXED BLOCKS

## Diagonal model
model_fixed_block <- normalblockr:::NB_fixed_blocks$new(data, C)
model_fixed_block$optimize()
model_fixed_block <- normalblockr:::NB_fixed_blocks$new(data, C, penalty = 0.05)
model_fixed_block$optimize()
model_fixed_block <- normalblockr:::NB_fixed_blocks$new(data, C, penalty = 0.05,
                                                                 control = NB_control(heuristic = TRUE))
model_fixed_block$optimize()

## Spherical model
ctrl <- NB_control(noise_covariance = "spherical")
model_fixed_block <- normalblockr:::NB_fixed_blocks$new(data, C, control = ctrl)
model_fixed_block$optimize()
model_fixed_block <- normalblockr:::NB_fixed_blocks$new(data, C, penalty = 0.05, control = ctrl)
model_fixed_block$optimize()
model_fixed_block <- normalblockr:::NB_fixed_blocks$new(data, C, penalty = 0.05,
                                                        control = NB_control(heuristic = TRUE, noise_covariance = "spherical"))
model_fixed_block$optimize()

## =======================================================
## UNKNOWN BLOCKS

## Diagonal model
model_unkwn_block <- normalblockr:::NB_fixed_Q$new(data, Q)
model_unkwn_block$optimize()
model_unkwn_block <- normalblockr:::NB_fixed_Q$new(data, Q, penalty = 0.05)
model_unkwn_block$optimize()
model_unkwn_block <- normalblockr:::NB_fixed_Q$new(data, Q, penalty = 0.05,
                                                            control = NB_control(heuristic = TRUE))
model_unkwn_block$optimize()

## Spherical model
ctrl <- NB_control(noise_covariance = "spherical")
model_unkwn_block <- normalblockr:::NB_fixed_Q$new(data, Q, control = ctrl)
model_unkwn_block$optimize()
model_unkwn_block <- normalblockr:::NB_fixed_Q$new(data, Q, penalty = 0.05, control = ctrl)
model_unkwn_block$optimize()
model_unkwn_block <- normalblockr:::NB_fixed_Q$new(data, Q, penalty = 0.05,
                                                   control = NB_control(heuristic = TRUE,
                                                                        noise_covariance = "spherical"))
model_unkwn_block$optimize()

