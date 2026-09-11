# normalblockr: Gaussian Graphical Models with Latent Clustering Structure

Implements the Normal-Block model, a Gaussian graphical model with a
latent clustering structure for the multivariate analysis of continuous
data. The model clusters variables and, building on the graphical lasso,
infers a network of statistical dependencies between clusters rather
than between individual variables, for known or unknown clusterings,
with an optional zero-inflation extension for data with an excess of
exact zeros. A complementary family clusters variables by their
regression response to covariates rather than by their covariance,
sharing one profile per cluster. See Tous & Chiquet (2026)
[doi:10.1016/j.csda.2026.108347](https://doi.org/10.1016/j.csda.2026.108347)
for the model itself and its variational expectation-maximization
estimation procedure.

## See also

Useful links:

- <https://github.com/jchiquet/normalblockr>

- <https://jchiquet.github.io/normalblockr/>

- Report bugs at <https://github.com/jchiquet/normalblockr/issues>

## Author

**Maintainer**: Julien Chiquet <julien.chiquet@inrae.fr>
([ORCID](https://orcid.org/0000-0002-3629-3429))

Authors:

- Julien Chiquet <julien.chiquet@inrae.fr>
  ([ORCID](https://orcid.org/0000-0002-3629-3429))

- Jeanne Tous <jeanne.tous@inrae.fr>

Other contributors:

- Nestor Ngalala Manguitini \[contributor\]
