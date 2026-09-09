#ifndef NORMALBLOCKR_ZI_NORMAL_BLOCK_VAR_KNOWN_CLUSTERS_H
#define NORMALBLOCKR_ZI_NORMAL_BLOCK_VAR_KNOWN_CLUSTERS_H

#include <RcppArmadillo.h>
#include "normal_block_var_base.h"
#include "zi_normal_block_data.h"
#include "zi_closed_form_solvers.h"
#include "utils_arma.h"

// Zero-inflated normal-block-var model with known clusters (fixed indicator
// matrix C), templated on the residual-noise policy (see zi_noise_models.h).
// Equivalent of the R6 class ZINormalBlockVarKnownClusters
// (R/ZINormalBlockVarKnownClusters.R). Unlike the non-ZI counterpart, the
// posterior covariance of W | Y is row-dependent (the ZI mask varies row by
// row), so Gamma_ is a cube (q x q x n) rather than a single matrix, and B's
// normal equations are solved one column at a time (zi_closed_form_solvers.h).
template <typename NoisePolicy>
class ZINormalBlockVarKnownClusters : public NormalBlockVarBase {
  const ZINormalBlockData& zi_data_;
  arma::mat C_;       // p x q, fixed cluster-indicator matrix
  arma::cube Gamma_;  // q x q x n, posterior covariance of W | Y (one matrix per row)
  arma::mat Mu_;      // n x q, posterior mean of W | Y
  arma::mat R_;       // n x p, Y - X*B computed at the start of the (V)EM step,
                       // reused unchanged by both E_step() and M_step()
  arma::mat dm1_mat_; // n x p, repmat(dm1_, n) % zeros_bar, computed once in
                      // E_step() and reused in M_step() (dm1_ does not change
                      // in between -- it is only updated at the end of M_step())
  arma::mat sum_slices() const {
    arma::mat out(q_, q_, arma::fill::zeros);
    for (arma::uword i = 0; i < Gamma_.n_slices; ++i) out += Gamma_.slice(i);
    return out;
  }

  void E_step() override {
    R_ = zi_data_.Y - XB();
    dm1_mat_ = arma::repmat(dm1_.t(), zi_data_.n, 1) % zi_data_.zeros_bar;
    arma::mat dm1C = dm1_mat_ * C_;             // n x q
    arma::mat Rdm1C = (R_ % dm1_mat_) * C_;     // n x q

    // Same shape as solve_M_ridge()'s loop (zi_closed_form_solvers.h): one
    // q x q system per row, Omega fixed, only its diagonal moving. The full
    // inverse is genuinely needed here: M_step() reads both each slice's
    // diagonal and their sum, but the two temporaries the naive form
    // allocated per row are not.
    const arma::vec omega_diag = Omega_.diag();
    arma::mat A = Omega_;
    for (arma::uword i = 0; i < zi_data_.n; ++i) {
      A.diag() = omega_diag + dm1C.row(i).t();
      Gamma_.slice(i) = arma::inv_sympd(A);
      Mu_.row(i) = Rdm1C.row(i) * Gamma_.slice(i);
    }
  }

  void M_step() override {
    arma::mat MuCT = Mu_ * C_.t();
    set_B(nb_optim::solve_wls(dm1_mat_, zi_data_.Y, zi_data_.X, MuCT));
    R_ = zi_data_.Y - XB();

    arma::mat RmMuCT = R_ - MuCT;
    arma::mat CgC(zi_data_.n, zi_data_.p);
    for (arma::uword i = 0; i < zi_data_.n; ++i) CgC.row(i) = (C_ * Gamma_.slice(i).diag()).t();
    arma::mat A = arma::square(RmMuCT) + CgC;

    arma::vec weighted_ssq = arma::vectorise(arma::sum(zi_data_.zeros_bar % A, 0));
    dm1_ = NoisePolicy::update_dm1(weighted_ssq, zi_data_.nY, zi_data_.npY);

    arma::mat Sigma_hat = (Mu_.t() * Mu_ + sum_slices()) / zi_data_.n;
    Omega_ = estimate_omega(Sigma_hat);
  }

public:
  ZINormalBlockVarKnownClusters(const ZINormalBlockData& data, const arma::mat& C,
                              const arma::mat& B0, const arma::vec& dm1_0, const arma::mat& Omega0,
                              double sparsity, const arma::mat& sparsity_weights) :
    NormalBlockVarBase(data, C.n_cols, B0, dm1_0, Omega0, sparsity, sparsity_weights),
    zi_data_(data), C_(C), Gamma_(C.n_cols, C.n_cols, data.n), Mu_(arma::zeros(data.n, C.n_cols)) {
    for (arma::uword i = 0; i < data.n; ++i) Gamma_.slice(i) = arma::eye(C.n_cols, C.n_cols);
  }

  // General (non-profiled) marginal log-likelihood given the fixed ZI mask,
  // valid at any (B_, dm1_, Omega_); same Woodbury/determinant-lemma trick
  // as NormalBlockVarKnownClusters::objective(), applied per row since the
  // ZI mask makes each row's marginal (and posterior precision Gamma^{(i)})
  // distinct. See inst/normal_block_models.qmd §6/§7.
  double objective() const override {
    arma::mat R = zi_data_.Y - XB();
    arma::mat dm1_mat = arma::repmat(dm1_.t(), zi_data_.n, 1) % zi_data_.zeros_bar;
    arma::mat dm1C = dm1_mat * C_;          // n x q
    arma::mat Rdm1C = (R % dm1_mat) * C_;   // n x q

    double sum_log_Gamma_inv = 0.0;
    double quad_mu = 0.0;
    for (arma::uword i = 0; i < zi_data_.n; ++i) {
      arma::mat Gamma_inv_i = Omega_ + arma::diagmat(dm1C.row(i).t());
      auto Gamma_fresh_i = nb_utils::inv_and_log_det_sympd(Gamma_inv_i);
      sum_log_Gamma_inv += Gamma_fresh_i.log_det;
      arma::vec mu_i = (Rdm1C.row(i) * Gamma_fresh_i.inv).t();
      quad_mu += arma::as_scalar(mu_i.t() * Gamma_inv_i * mu_i);
    }

    double log_det_Omega = arma::log_det_sympd(Omega_);
    double weighted_sum_log_dm1 = arma::accu(zi_data_.nY % arma::log(dm1_));
    double SSQ_w = arma::accu(dm1_mat % arma::square(R));

    double J = -0.5 * zi_data_.npY * std::log(2.0 * arma::datum::pi);
    J += 0.5 * weighted_sum_log_dm1;
    J += 0.5 * zi_data_.n * log_det_Omega;
    J -= 0.5 * sum_log_Gamma_inv; // -log|Gamma^{-1}| = +log|Gamma|, summed over rows
    J -= 0.5 * (SSQ_w - quad_mu);
    J += zi_data_.zi_cond_mean;
    return J;
  }

  const arma::cube& Gamma() const { return Gamma_; }
  const arma::mat& Mu() const { return Mu_; }

  std::unique_ptr<NormalBlockEMBase> clone() const override {
    return std::make_unique<ZINormalBlockVarKnownClusters<NoisePolicy>>(*this);
  }
  void restore_from(const NormalBlockEMBase& other) override {
    const auto& o = static_cast<const ZINormalBlockVarKnownClusters<NoisePolicy>&>(other);
    copy_tracked_state_from(o);
    Gamma_ = o.Gamma_;
    Mu_ = o.Mu_;
    R_ = o.R_;
    dm1_mat_ = o.dm1_mat_;
  }
  // Same general/always-valid objective() as NormalBlockVarKnownClusters, so
  // the same SQUAREM gate applies; see its supports_acceleration().
  bool supports_acceleration() const override { return sparsity_ <= 0.0; }
};

#endif // NORMALBLOCKR_ZI_NORMAL_BLOCK_VAR_KNOWN_CLUSTERS_H
