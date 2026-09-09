#ifndef NORMALBLOCKR_NORMAL_BLOCK_VAR_KNOWN_CLUSTERS_H
#define NORMALBLOCKR_NORMAL_BLOCK_VAR_KNOWN_CLUSTERS_H

#include <RcppArmadillo.h>
#include "normal_block_var_base.h"
#include "normal_block_data.h"
#include "utils_arma.h"

// Normal-block model with known clusters (fixed indicator matrix C), templated
// on the residual-noise policy (see noise_models.h). Equivalent of the R6
// class NormalBlockVarKnownClusters (R/NormalBlockVarKnownClusters.R); see
// inst/normal_block_models.qmd §1/§2 for the EM recursion.
template <typename NoisePolicy>
class NormalBlockVarKnownClusters : public NormalBlockVarBase {
  arma::mat C_;      // p x q, fixed cluster-indicator matrix
  arma::mat Gamma_;  // q x q, posterior covariance of W | Y
  arma::mat Mu_;     // n x q, posterior mean of W | Y

  void E_step() override {
    // Gamma = (Omega + C^T diag(dm1) C)^{-1}, Mu = R diag(dm1) C Gamma
    arma::mat dm1C = C_;
    dm1C.each_col() %= dm1_;
    Gamma_ = arma::inv_sympd(Omega_ + arma::diagmat(C_.t() * dm1_));
    arma::mat R = data_.Y - XB();
    Mu_ = R * dm1C * Gamma_;
  }

  void M_step() override {
    arma::mat YmMuCT = data_.Y - Mu_ * C_.t();
    set_B(data_.XtXm1 * (data_.X.t() * YmMuCT));
    arma::mat resid = YmMuCT - XB();
    arma::vec ddiag = arma::vectorise(arma::mean(arma::square(resid), 0));
    ddiag += C_ * Gamma_.diag();
    dm1_ = NoisePolicy::update_dm1(ddiag);
    arma::mat Sigma_hat = Mu_.t() * Mu_ / data_.n + Gamma_;
    Omega_ = estimate_omega(Sigma_hat);
  }

public:
  NormalBlockVarKnownClusters(const NormalBlockData& data, const arma::mat& C,
                              const arma::mat& B0, const arma::vec& dm1_0, const arma::mat& Omegaq0,
                              double sparsity, const arma::mat& sparsity_weights) :
    NormalBlockVarBase(data, C.n_cols, B0, dm1_0, Omegaq0, sparsity, sparsity_weights),
    C_(C), Gamma_(arma::eye(C.n_cols, C.n_cols)), Mu_(arma::zeros(data.n, C.n_cols)) {}

  // General (non-profiled) marginal log-likelihood of Y, valid at any
  // (B_, dm1_, Omega_), not just an M-step optimum, via the matrix
  // determinant lemma and Woodbury identity on the q x q posterior
  // precision Gamma^{-1}. See inst/normal_block_models.qmd ("Criterion",
  // §1/§2) for the derivation and the sign-bug fix this replaced.
  double objective() const override {
    arma::mat R = data_.Y - XB();
    arma::mat Gamma_inv = Omega_ + arma::diagmat(C_.t() * dm1_);
    auto Gamma_fresh = nb_utils::inv_and_log_det_sympd(Gamma_inv);
    arma::mat dm1C = C_;
    dm1C.each_col() %= dm1_;
    arma::mat Mu_fresh = R * dm1C * Gamma_fresh.inv;

    double log_det_Omega = arma::log_det_sympd(Omega_);
    double sum_log_dm1    = arma::sum(arma::log(dm1_));
    double SSQ_w = arma::dot(dm1_, arma::vectorise(arma::sum(arma::square(R), 0)));
    double quad  = SSQ_w - arma::trace(Gamma_inv * (Mu_fresh.t() * Mu_fresh));

    double J = -0.5 * data_.n * data_.p * std::log(2.0 * arma::datum::pi);
    J += 0.5 * data_.n * sum_log_dm1;
    J += 0.5 * data_.n * log_det_Omega;
    J -= 0.5 * data_.n * Gamma_fresh.log_det; // -log|Gamma^{-1}| = +log|Gamma|
    J -= 0.5 * quad;
    return J;
  }

  const arma::mat& Gamma() const { return Gamma_; }
  const arma::mat& Mu() const { return Mu_; }

  std::unique_ptr<NormalBlockEMBase> clone() const override {
    return std::make_unique<NormalBlockVarKnownClusters<NoisePolicy>>(*this);
  }
  void restore_from(const NormalBlockEMBase& other) override {
    const auto& o = static_cast<const NormalBlockVarKnownClusters<NoisePolicy>&>(other);
    copy_tracked_state_from(o);
    Gamma_ = o.Gamma_;
    Mu_ = o.Mu_;
  }
  // Excludes sparsity_ > 0: the graphical lasso is an approximate iterative
  // solver,
  // so even plain EM's M-step isn't a guaranteed ascent step there, and
  // SQUAREM's larger jumps amplify that instability for a smaller speedup
  // than the unpenalized case. See inst/normal_block_models.qmd ("SQUAREM").
  bool supports_acceleration() const override { return sparsity_ <= 0.0; }
};

#endif // NORMALBLOCKR_NORMAL_BLOCK_VAR_KNOWN_CLUSTERS_H
