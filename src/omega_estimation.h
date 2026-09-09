#ifndef NORMALBLOCKR_OMEGA_ESTIMATION_H
#define NORMALBLOCKR_OMEGA_ESTIMATION_H

#include <RcppArmadillo.h>
#include "graphical_lasso.h"

// Equivalent of the shared private method `NormalBlockVarBase$get_Omega` (R/NormalBlockVarBase.R),
// used identically by the M-step of both NormalBlockVarKnownClusters and NormalBlockVarUnknownClusters:
// estimate the precision matrix of the blocks from its covariance estimate
// Sigma_hat (q x q). When `sparsity <= 0`, this is a plain matrix inversion;
// otherwise it runs the in-package graphical lasso (src/graphical_lasso.h).
namespace nb_omega {

// MUST match NB_GLASSO_THRESHOLD (R/utils.R), which the R reference recursion uses: the two are compared trace-for-trace at 1e-8 in test-cpp-normal-block-mean.R
constexpr double kGlassoThreshold = 1e-6;

// Armadillo equivalent of ensure_pd() (R/utils.R): the graphical lasso can
// return a precision matrix that is not quite positive definite (an EM
// iterate's Sigma can be badly conditioned, especially for the p x p Sigma of
// the mean-block family), which would then make log_det_sympd() throw
inline arma::mat ensure_pd(const arma::mat& M, double floor_value = 1e-6) {
  arma::mat sym = arma::symmatu(M), R;
  if (arma::chol(R, sym)) return sym;
  arma::vec eigval;
  arma::mat eigvec;
  if (!arma::eig_sym(eigval, eigvec, sym)) return sym;
  eigval = arma::clamp(eigval, floor_value, eigval.max());
  return eigvec * arma::diagmat(eigval) * eigvec.t();
}

// `warm`, when given, carries the previous M-step's (W, X) so the graphical
// lasso can resume from it. Consecutive M-steps solve nearly the same problem,
// and the solver's stopping rule measures per-sweep progress rather than
// distance to the optimum, so a warm start both converges in fewer sweeps
// and, at the tightened threshold this affords, lands closer to the optimum
// than a cold start at the looser one.
inline arma::mat estimate(const arma::mat& Sigma_hat, double sparsity,
                          const arma::mat& sparsity_weights,
                          nb_glasso::State* warm = nullptr,
                          double thr = 1e-4) {
  if (sparsity <= 0.0) {
    return arma::inv_sympd(Sigma_hat);
  }

  const arma::mat penalty = sparsity * sparsity_weights;
  nb_glasso::Result glasso_out = nb_glasso::solve(Sigma_hat, penalty, thr, 10000, warm);

  // A warm start is only ever a starting point, and a bad one is not free: on
  // an ill-conditioned Sigma it can send the coordinate descent off to
  // infinity where a cold start on the very same problem converges in a
  // handful of sweeps (observed on the mean-block p x p Sigma, and it is why
  // this retry exists rather than a direct fall-through to the unpenalized
  // inverse.
  if (glasso_out.X.has_nan() && warm != nullptr && warm->usable_for(Sigma_hat.n_rows)) {
    glasso_out = nb_glasso::solve(Sigma_hat, penalty, thr, 10000, nullptr);
  }

  if (warm != nullptr) {
    if (glasso_out.X.is_finite() && glasso_out.W.is_finite()) warm->store(glasso_out);
    else warm->reset();
  }

  if (glasso_out.X.has_nan()) {
    Rcpp::warning("GLasso fails, the penalty is probably too small and the system badly "
                  "conditionned (reciprocal condition number = %f). "
                  "Sending back the original matrix and its inverse (unpenalized).",
                  arma::rcond(Sigma_hat));
    return arma::inv_sympd(Sigma_hat);
  }
  return ensure_pd(glasso_out.X); // symmpart(), plus a guard against a non-PD glasso output
}

} // namespace nb_omega

#endif // NORMALBLOCKR_OMEGA_ESTIMATION_H
