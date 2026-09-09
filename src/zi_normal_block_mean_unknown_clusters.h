#ifndef NORMALBLOCKR_ZI_NORMAL_BLOCK_MEAN_UNKNOWN_CLUSTERS_H
#define NORMALBLOCKR_ZI_NORMAL_BLOCK_MEAN_UNKNOWN_CLUSTERS_H

#include <RcppArmadillo.h>
#include "zi_normal_block_mean_base.h"
#include "utils_arma.h"

// Zero-inflated mean-block model with unknown clusters, inferred by
// variational EM. Equivalent of the R6 class ZINormalBlockMeanUnknownClusters
// (R/ZINormalBlockMeanUnknownClusters.R).
class ZINormalBlockMeanUnknownClusters : public ZINormalBlockMeanBase {
  arma::mat tau_;   // p x q, variational membership probabilities
  arma::vec alpha_; // q, prior cluster probabilities
  bool fixed_tau_;  // if true, tau is left at its initial value (stability selection)

  void M_step() override {
    update_B(tau_);
    refresh_stats();
    update_omega(s_from(tau_));
    alpha_ = arma::vectorise(arma::mean(tau_, 0));
  }

  // With a diagonal Sigma the ELBO is linear in each row of tau (up to the
  // entropy) and the rows are independent, so this softmax is the exact
  // maximizer in one shot.
  void E_step() override {
    if (fixed_tau_) return;
    arma::vec w = Omega_.diag();
    arma::mat eta = P_ - 0.5 * A_;
    eta.each_col() %= w;
    eta.each_row() += arma::log(arma::clamp(alpha_, 1e-300, arma::datum::inf)).t();
    tau_ = nb_utils::softmax_rows(eta);
    tau_ = nb_utils::check_one_boundary(nb_utils::check_zero_boundary(tau_));
    tau_.each_col() /= arma::sum(tau_, 1);
  }

public:
  ZINormalBlockMeanUnknownClusters(const ZINormalBlockData& data,
                                   const arma::mat& B0, const arma::mat& Omega0,
                                   const arma::mat& tau0, bool accelerate, bool fixed_tau,
                                   const std::string& cov_structure) :
    ZINormalBlockMeanBase(data, tau0.n_cols, B0, Omega0, accelerate, cov_structure),
    tau_(tau0), fixed_tau_(fixed_tau) {
    alpha_ = arma::vectorise(arma::mean(tau_, 0));
  }

  double objective() const override {
    refresh_stats();
    return gaussian_objective(tau_)
      + arma::accu(tau_ * arma::log(alpha_)) - arma::accu(tau_ % arma::log(tau_));
  }

  const arma::mat& tau() const { return tau_; }
  const arma::vec& alpha() const { return alpha_; }

  std::unique_ptr<NormalBlockEMBase> clone() const override {
    return std::make_unique<ZINormalBlockMeanUnknownClusters>(*this);
  }
  void restore_from(const NormalBlockEMBase& other) override {
    const auto& o = static_cast<const ZINormalBlockMeanUnknownClusters&>(other);
    copy_tracked_state_from(o);
    tau_ = o.tau_; alpha_ = o.alpha_;
  }
};

#endif // NORMALBLOCKR_ZI_NORMAL_BLOCK_MEAN_UNKNOWN_CLUSTERS_H
