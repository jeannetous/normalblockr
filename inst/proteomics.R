library(tidyverse)
library(future)

proteomics <- read_csv(file = "itraq/77_cancer_proteomes_CPTAC_itraq.csv")
metadata <- dplyr::select(proteomics, c(RefSeq_accession_number, gene_symbol, gene_name))
proteom <- dplyr::select(proteomics, -c(RefSeq_accession_number, gene_symbol, gene_name))

## 1st filter with NA
to_keep <- rowSums(is.na(proteom)) == 0
proteom <- filter(proteom, to_keep)
metadata <- filter(metadata, to_keep)

## 2nd filter with variances
variances <- apply(proteom, 1, var)
proteom  <- filter(proteom , variances > quantile(variances, 0.80))
metadata <- filter(metadata, variances > quantile(variances, 0.80))
Y <- as.matrix(t(proteom))
X <- matrix(1, nrow(Y), 1)


plan(multisession, workers = 15)

## Adjust NORMAL Block for varying number of groups
# n = 87, so need to regularize to go beyond 80 groups
# nb_blocks <- c(2:9,seq(10,85,by=5))
# res <- normal_block(Y, X, nb_blocks)
nb_blocks <- c(2:9,seq(10,150,by=10))
res <- normal_block(Y, X, nb_blocks, sparsity = 0.1)

plan("sequential")

# take 10 seconds on my computer

par(mfrow = c(1,1))
plot(nb_blocks, -2*res$criteria$loglik, type = 'l', log = "y")
lines(nb_blocks, -2 *res$criteria$pen_loglik, col="green")
lines(nb_blocks, res$criteria$BIC, col="red")
lines(nb_blocks, res$criteria$EBIC, col="blue")
legend("topright", legend=c("-2 loglik", "-2 pen_loglik", "BIC", "EBIC"),
       col=c("black", "green", "red", "blue"), lty=1, cex=0.8)

K_star <- nb_blocks[which.min(res$criteria$BIC)]
best_model <- res$models[[which.min(res$criteria$BIC)]]

corrplot::corrplot(cov2cor(solve(best_model$model_par$omegaQ)), order = "hclust")

gene_groups <- split(metadata, best_model$clustering)
