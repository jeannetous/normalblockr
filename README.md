
# normalblockr

<!-- badges: start -->
  [![R-CMD-check](https://github.com/jeannetous/normalblockr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jeannetous/normalblockr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Normal-Block is a graphical model designed for the multivariate analysis of continuous data.
It clusters variables and builds on Graphical-Lasso to infer a network of statistical dependencies 
between clusters.

## Installation

You can install the development version of normalblockr from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("jeannetous/normalblockr")
```

## Usage 

All fitting is done using the function normal_block. An object of the class
normal_data must be created (using observations Y and covariates X) and used as 
an input to normal_block. normal_block parameters allow to choose between sparse
/ not sparse models, fixed or unknown blocks... 


The package comes with a small artificial data set simulated under the model.

```r
testdata <- readRDS("testdata/testdata_normal.RDS")
Y        <- testdata$Y ; X <- testdata$X
C        <- testdata$parameters$C ; Q <- ncol(C)
data     <- normal_data$new(Y, X)
model    <- normal_block(data, blocks = C)
```

### Sparse model

```r
model    <- normal_block(data, blocks = C, sparsity = TRUE)
```

### With unobserved clustering

```r
model    <- normal_block(data, blocks = 2:5, sparsity = TRUE)
model$plot_criteria()
model_BIC <- model$get_best_model("BIC")
```


