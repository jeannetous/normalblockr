# Mean-block analyses

Exploratory analyses and simulation studies backing the design choices of the
mean-block family (`normal_block(..., model = "mean")`). Kept out of the built
package (`.Rbuildignore`); the scripts are meant to be sourced from the package
root with `devtools::load_all()` already done.

Not to be confused with `inst/CSDA_analyses/`, which holds the analyses of the
published article on the *variance*-block model.

| script | what it establishes |
|---|---|
| `deviance_monotonicity.R` | Deviance must not increase with `q` (nested models), so every violation measures an optimization failure. On `brca_rppa`, a full Sigma violates it 4/14 times and `refine()` repairs only 2, while diagonal and spherical violate it 0/14. Going from `q = 1` to `q = 15` also buys the full Sigma only 0.65% of deviance, against 38% for the diagonal one: its ~13,000 covariance parameters absorb the structure the clustering would otherwise explain. |
| `covariance_shape_selection_study.R` | Simulation study behind the diagonal default, 12 replicates, both generative regimes. At `n/p = 1.3`, BIC picks the true `q` 10/12 times with a diagonal Sigma against 6/12 with a full one -- even when the data *were* generated with a full Sigma. The ARI at fixed `q` is the same either way: what a full Sigma costs is the choice of `q`. The honest counter-example is `n/p = 3.3` with a genuinely full Sigma, where "full" does select better (9/12 vs 6/12). |
| `sequential_mean_then_variance.R` | Validates `normal_block_sequential()`. A positive control carrying two genuinely distinct structures recovers both exactly (`q` and ARI = 1.000 on each). On `brca_rppa` the two partitions are unrelated (ARI 0.018), so a model forcing them equal would be badly misspecified; but the covariance clustering is largely unchanged by removing the mean structure first (ARI 0.712 against fitting it directly), which suggests the two structures are separable rather than competing. |
| `nestor_protocol.R` | Nestor Ngalala Manguitini's M2 simulation protocol (`inst/normalblockmean/code_nestor`), re-run on this package's estimator: his data-generating process verbatim, our estimator, his grid (n in 80/150/220/300, p = 60, q in 3/5, SNR in 0.3/1.5, 10 replicates, lambda = 0.1). 160 configurations, ~43 min. Results in `nestor_protocol_results.rds`. Four findings: (i) `kmeans` is the best of our four initialisations at both SNRs (mean ARI 0.90 at SNR 0.3, 0.97 at SNR 1.5), with `sbm` a close second and `ward2`/`spectral` clearly behind -- which is what makes `kmeans` this family's default; (ii) RMSE on `BC'` falls monotonically with `n` in every cell, so the estimator behaves; (iii) BIC and ICL put the median selected `q` exactly on the truth in all four (SNR, q) cells, though only 45-80% of individual replicates land exactly right, worse at SNR 0.3 as expected; (iv) RMSE on Sigma does *not* improve with `n`, which is his protocol rather than our estimator -- lambda is fixed at 0.1 while Sigma is rescaled per dataset to hit the target SNR, so the penalty bias dominates. Note his `RMSE(B)` and `RMSE(alpha)` compare quantities identified only up to a permutation of the cluster labels; `RMSE(BC')` is the one that means anything. |
| `nestor_protocol_qsel_unpenalised.R` | Follow-up to the above, separating criterion from estimator. Nestor penalises Omega at lambda = 0.1 but his `calculer_nu()` counts the *dense* p(p+1)/2 regardless, so his BIC is monotone in q by construction while ours -- which counts the non-zero off-diagonal terms, as lasso degrees-of-freedom theory prescribes -- is not: 977, 884, 887, 890, 893 at q = 2..6. Re-running his grid with a full but **unpenalised** Sigma, where our `nb_param` is identical to his nu at every q (1835, 1838, 1841, 1844, 1847), reverses the gap: exact q selection by BIC 80.0/62.5/87.5/87.5% against his 72.5/60.0/90.0/87.5%, with over-selection down from 35/27.5/25/20% to 7.5/10/10/7.5%. The earlier deficit was the criterion, not the estimator. Caveat: one lambda, one data-generating process; the non-monotonicity is expected whenever a larger q leaves Omega sparser, but its size has not been mapped. |
| `diagnostic_mean_block_brca.R` | Initialization, convergence and `refine()` on `brca_rppa`. |
| `diagnostic_mean_block_onema.R` | The same on `onema`, to check the findings are not specific to one dataset. |
| `diagnostic_mean_block_scaling.R` | How cost and quality scale with `n`, `p` and `q`, with a profiling pass. |

## Open questions

- **Joint mean + covariance clustering.** The sequential estimator above is a
  heuristic, not a joint model. Free partitions look more plausible than linked
  ones (the two selected `q` differ a lot: 70 and 20 on `brca_rppa`), and the
  separability measured there is an argument for identifiability. What a joint
  fit would buy over two stages -- statistical efficiency, a single criterion --
  is still open, as is the cost of selecting over a `q_mean` x `q_var` grid.
- **When does a full Sigma become viable again?** It wins at a comfortable
  `n/p`; the transition point is not mapped.
