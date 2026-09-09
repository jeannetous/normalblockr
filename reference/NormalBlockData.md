# Data Container for Normal-Block Models

R6 class holding the responses and design matrix used to fit a
normal-block model.

## Public fields

- `Y`:

  the matrix of responses (rescaled column-wise if \`scale = TRUE\`)

- `Y_scale`:

  the per-column standard deviation Y was divided by (all 1's if \`scale
  = FALSE\`)

- `X`:

  the matrix of covariates

- `X0`:

  the matrix of zero-inflation covariates, if applicable

- `formula`:

  describes the relationship between Y and X, and X0 if applicable,
  useful if not all of X's or X0's covariates should be used, should be
  formatted ~ X1 + X2... \| Z1 + Z2... with the Normal formula before
  the \| and the ZI formula after the \|

- `n`:

  sample size

- `d`:

  number of covariates

- `d0`:

  number of zero-inflation covariates, if applicable

- `p`:

  number of variables

- `XtX`:

  useful for inference in some cases

- `XtXm1`:

  inverse of XtX, useful for inference

- `XtY`:

  useful for inference

- `npY`:

  total number of non zeros in Y

- `nY`:

  total number of non zeros for each column/variable in Y

- `zeros`:

  where are the zero in Y

- `zeros_bar`:

  where are the non-zeros in Y

## Methods

### Public methods

- [`NormalBlockData$new()`](#method-NormalBlockData-initialize)

- [`NormalBlockData$ols_fit()`](#method-NormalBlockData-ols_fit)

- [`NormalBlockData$zi_ols_fit()`](#method-NormalBlockData-zi_ols_fit)

- [`NormalBlockData$zi_fit()`](#method-NormalBlockData-zi_fit)

- [`NormalBlockData$clone()`](#method-NormalBlockData-clone)

------------------------------------------------------------------------

### `NormalBlockData$new()`

Create a new \[\`NormalBlockData\`\] object.

#### Usage

    NormalBlockData$new(
      Y,
      X,
      X0 = NULL,
      formula = NULL,
      scale = TRUE,
      zeros = NULL
    )

#### Arguments

- `Y`:

  the matrix of responses (called Y in the model).

- `X`:

  design matrix (called X in the model).

- `X0`:

  zero-inflation design matrix, if applicable.

- `formula`:

  describes the relationship between Y and X, useful if not all of X's
  covariates should be used.

- `scale`:

  whether to rescale each column of Y by its own standard deviation (no
  centering). Default TRUE, see the class-level documentation for the
  rationale and its limits.

- `zeros`:

  an optional 0/1 matrix of structural zeros, overriding the default \`Y
  == 0\`.

------------------------------------------------------------------------

### `NormalBlockData$ols_fit()`

Ordinary-least-squares fit of Y on X, with its residuals and their
covariance. Computed once and memoized.

#### Usage

    NormalBlockData$ols_fit()

#### Returns

a list with \`B\` (d x p), \`R\` (n x p residuals) and \`Sigma\` (p x p
residual covariance)

------------------------------------------------------------------------

### `NormalBlockData$zi_ols_fit()`

Masked counterpart of \`ols_fit()\`: a per-variable weighted
least-squares fit of B under the zero-inflation mask, with its inverse
residual variances and residuals (see \`zi_weighted_fit()\`). Computed
once and memoized.

#### Usage

    NormalBlockData$zi_ols_fit()

#### Returns

a list with \`B\` (d x p), \`dm1\` (p) and \`R\` (n x p masked
residuals)

------------------------------------------------------------------------

### `NormalBlockData$zi_fit()`

Zero-inflation component: \`p\` independent logistic regressions of each
variable's zero pattern on \`X0\`, and the fixed contribution they make
to the log-likelihood. The (V)EM never revisits these, so they are a
property of the data rather than of a model.

#### Usage

    NormalBlockData$zi_fit()

#### Returns

a list with \`B0\` (d0 x p), \`kappa\` (n x p zero-inflation
probabilities) and \`ZI_cond_mean\` (a scalar)

------------------------------------------------------------------------

### `NormalBlockData$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockData$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
ex <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)
c(n = data$n, p = data$p, d = data$d)
#>  n  p  d 
#> 50 20  1 
```
