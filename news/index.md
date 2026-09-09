# Changelog

## normalblockr 0.3.0

### New model family: clustering in the mean

The package now fits two complementary families. The original
Normal-Block model clusters variables by how they *covary*; the new
mean-block family clusters them by how they *respond* to the covariates,
constraining the mean as `mu_i = C B' X_i` with one regression profile
per cluster.

- `normal_block(..., model = "mean")` fits it, for a known clustering
  (`NormalBlockMeanKnownClusters`) or an unknown one inferred by
  variational EM (`NormalBlockMeanUnknownClusters`), over a range of
  cluster counts (`NormalBlockMeanCollectionClusters`), a penalty path
  (`NormalBlockMeanCollectionSparsity`) or both
  (`NormalBlockMeanCollectionClustersSparsity`).
- Zero-inflated counterparts (`ZINormalBlockMeanKnownClusters`,
  `ZINormalBlockMeanUnknownClusters`), and collections of those over a
  range of cluster counts.
- `noise_covariance` gains a `"full"` shape for this family, whose Sigma
  is the p x p residual covariance rather than a diagonal noise term.
  `"diagonal"` is the default: a full Sigma spends p(p+1)/2 parameters
  that drown the mean structure the criteria are weighing, and selects
  the number of clusters markedly worse as p approaches n. Asking for
  `sparsity > 0` implies `"full"`, since a diagonal precision matrix has
  nothing for the graphical lasso to penalize.
- [`normal_block_sequential()`](../reference/normal_block_sequential.md)
  chains the two families: a mean-block fit, then a variance-block fit
  on its residuals. A heuristic two-stage estimator, not a joint model.
- `NormalBlockData$new()` gains a `zeros` argument, to carry an explicit
  zero-inflation mask when the matrix handed to a model is no longer the
  one carrying the zeros (the residuals of a first stage, say).

### The graphical lasso is now in-package

This also fixes the segfault CRAN reported on
`r-devel-linux-x86_64-fedora-gcc` when re-building
`breast-cancer-proteomics.Rmd`
(`address 0x580, cause 'memory not mapped'`), which our own CI had been
hitting intermittently on Linux.

- `glassoFast` is no longer a dependency. The C++ (V)EM used to call
  back into R once per M-step to reach it – thousands of round trips
  from inside a single `.Call` for one sparsity path, and a recurring
  source of intermittent crashes. The solver now lives in
  `src/graphical_lasso.h`.
- Two bugs in `glassoFast` are fixed on the way. When the empirical
  covariance has no off-diagonal mass the problem separates and the
  solution is diagonal; it returned `1 / max(rho_ii, eps)` there,
  dropping the variance term, so with an unpenalized diagonal it
  returned about 9e15 instead of `1 / S_ii` – and every `q = 1` problem
  takes that branch. Separately, its inner loop could not terminate on
  non-finite input.
- Each M-step warm-starts from the previous one, at a tightened
  threshold: faster *and* closer to the optimum than the previous cold
  starts.

### Bug fixes

- Zero-inflation with more than one zero-inflation covariate (`X0` with
  two or more columns) failed with “non-conformable arguments”: `B0` was
  built transposed. A second bug behind it made `kappa` non-finite when
  the response had no zeros at all. Both were unreachable before, so no
  working fit changes.
- Fixed two crashes on the penalized path, one of them a long-standing
  memory-safety bug (an R object cached in a C++ `static`).

### Other changes

- A collection over cluster counts no longer repeats work that does not
  depend on the cluster count: the OLS fit, the zero-inflation logistic
  regressions, and the shareable part of each clustering heuristic (one
  hierarchical tree cut at each q, one eigendecomposition, one lossless
  row compression) are computed once. Fitted values are unchanged.
- Printed model descriptions now name the residual-covariance shape for
  every family, and the zero-inflated variance-block models name their
  family.

## normalblockr 0.2.1

CRAN release: 2026-09-03

First CRAN submission

- S3 methods [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) and
  [`BIC()`](https://rdrr.io/r/stats/AIC.html) for fitted models (any
  `NormalBlockVarBase` subclass), and
  [`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)/[`logLik()`](https://rdrr.io/r/stats/logLik.html)/[`BIC()`](https://rdrr.io/r/stats/AIC.html)
  for collections of models; accessing `$loglik` on a collection now
  raises an informative error instead of silently returning `NULL`.
- Shortened/title-cased man page titles, added missing `@examples`,
  cross-referenced `$plot_network()`/`$latent_network()` in each other’s
  documentation, and replaced a few inefficient matrix operations
  ([`solve()`](https://rdrr.io/r/base/solve.html) on symmetric
  positive-definite matrices, `M %*% t(C)`) with
  `chol2inv(chol())`/[`tcrossprod()`](https://rdrr.io/r/base/crossprod.html).

## normalblockr 0.2.0

### New features

- Zero-inflation extension
  (`ZINormalBlockVarKnownClusters`/`ZINormalBlockVarUnknownClusters`)
  for data with an excess of exact zeros.
- Sparsity path on the cluster-level precision matrix (graphical lasso),
  with warm-starting across penalties
  (`NormalBlockVarCollectionSparsity`, `sparsity = TRUE` in
  [`normal_block()`](../reference/normal_block.md)).
- Accelerated variational EM (SQUAREM-style extrapolation) for both
  known- and unknown-clustering models.
- Several clustering-initialization heuristics (`ward2`, `kmeans`,
  `spectral`, `sbm`, selectable via `NB_control(clustering_init = )`),
  and `best_of_inits()` to try several and keep the best-ELBO fit.
- `refine()` on `NormalBlockVarCollectionClusters`: post-hoc split/merge
  search seeded from neighboring cluster counts, to escape mediocre
  local optima left by independent per-q cold starts.
- New real datasets: `brca_rppa` (breast cancer proteomics), `onema`
  (French stream fish biomass, zero-inflated), `university` (WebKB text
  data), each with a dedicated vignette.

### Other changes

- Renamed the model classes (`NormalBlock*` to `NormalBlockVar*`) for
  consistency; `NormalBlockData` is unchanged.
- `NormalBlockData` rescales columns of `Y` by default; fitted values
  and regression coefficients are reported back on the original scale.
- Removed the `ClustOfVar` dependency (the `kmeansvar` clustering
  heuristic was dropped after benchmarking showed it was both the
  worst-ranked and least reliable of the available heuristics).
- Package cleanup for CRAN submission: license, documentation, and
  package structure.

## normalblockr 0.1.0

- Initial implementation of the Normal-Block model: a Gaussian graphical
  model with a latent clustering structure, for known or unknown
  clusterings, fit by variational EM.
