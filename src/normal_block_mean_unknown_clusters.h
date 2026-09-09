#ifndef NORMALBLOCKR_NORMAL_BLOCK_MEAN_UNKNOWN_CLUSTERS_H
#define NORMALBLOCKR_NORMAL_BLOCK_MEAN_UNKNOWN_CLUSTERS_H

#include <RcppArmadillo.h>
#include "normal_block_mean_base.h"
#include "normal_block_data.h"
#include "utils_arma.h"

// Mean-block model with unknown clusters, inferred through a variational EM
// (mean-field on the cluster membership). Equivalent of the R6 class
// NormalBlockMeanUnknownClusters (R/NormalBlockMeanUnknownClusters.R); see
// the intern's report (inst/normalblockmean) for the derivations, eq.
// (2.32)-(2.40).
class NormalBlockMeanUnknownClusters : public NormalBlockMeanBase {
  arma::mat tau_;      // p x q, variational membership probabilities
  arma::vec alpha_;    // q, prior cluster probabilities
  arma::mat Psi_;      // q x q, E_Q[C' Omega C], set by the M-step
  arma::vec lambda_;   // p, diagonal of the M-step's variational correction
  int fixed_point_niter_;
  bool fixed_tau_;     // if true, tau is left at its initial value (stability selection)
  mutable arma::mat Phi_; // Psi - tau' Omega tau, refreshed by objective()

  // Psi = tau' Omega tau + Diag(tau' w) - tau' Diag(w) tau, w = diag(Omega),
  // plus the same jitter the R implementation adds (eq. 2.32).
  arma::mat psi_from(const arma::mat& Omega, const arma::mat& tau) const {
    arma::vec w = Omega.diag();
    arma::mat Psi = arma::diagmat(tau.t() * w);
    // with a diagonal Omega the two remaining terms cancel exactly
    if (!omega_is_diagonal()) {
      arma::mat tw = tau;
      tw.each_col() %= w;
      Psi += (tau.t() * Omega) * tau - tau.t() * tw;
    }
    Psi.diag() += 1e-8;
    return Psi;
  }

  arma::mat tau_omega_tau(const arma::mat& tau) const {
    if (!omega_is_diagonal()) return (tau.t() * Omega_) * tau;
    arma::mat tw = tau;
    tw.each_col() %= Omega_.diag();
    return tau.t() * tw;
  }

  // Diagonal of the correction term of eq. (2.35), kept as a vector (the R
  // version materializes two dense p x p diagonal matrices instead).
  arma::vec lambda_from(const arma::mat& B, const arma::mat& tau) const {
    arma::mat M = (B.t() * data_.XtX) * B;
    return (tau * M.diag() - arma::sum((tau * M) % tau, 1)) / data_.n;
  }

  void M_step() override {
    Psi_ = psi_from(Omega_, tau_);
    B_   = ((((data_.XtXm1 * data_.XtY) * Omega_) * tau_) * arma::inv(Psi_));
    lambda_ = lambda_from(B_, tau_);

    arma::mat R_bar = data_.Y - (data_.X * B_) * tau_.t();
    Omega_ = omega_from_residuals(R_bar, lambda_);
    alpha_ = arma::vectorise(arma::mean(tau_, 0));
  }

  // Sequential (Gauss-Seidel) sweep over the rows of tau: each row's softmax
  // maximizes the ELBO exactly (its quadratic terms cancel), so a sweep can't
  // decrease it -- updating all rows at once can cycle.
  void E_step() override {
    if (fixed_tau_) return;
    arma::mat M = (B_.t() * data_.XtX) * B_;
    arma::rowvec diag_M = M.diag().t();
    arma::vec w = Omega_.diag();
    arma::mat ZtXB = data_.Y.t() * (data_.X * B_);
    arma::rowvec log_alpha = arma::log(arma::clamp(alpha_, 1e-300, arma::datum::inf)).t();

    for (int sweep = 0; sweep < fixed_point_niter_; ++sweep) {
      arma::mat G = ZtXB - tau_ * M; // kept in sync with tau_, row by row
      for (arma::uword j = 0; j < tau_.n_rows; ++j) {
        arma::rowvec omega_row_G;
        if (omega_is_diagonal()) omega_row_G = w(j) * G.row(j); else omega_row_G = Omega_.row(j) * G;
        arma::rowvec expo = log_alpha + omega_row_G
          + w(j) * (tau_.row(j) * M) - 0.5 * w(j) * diag_M;
        arma::rowvec t = arma::exp(expo - expo.max());
        t /= arma::accu(t);
        tau_.row(j) = t;
        G.row(j) = ZtXB.row(j) - t * M;
      }
      tau_ = nb_utils::check_one_boundary(nb_utils::check_zero_boundary(tau_));
      tau_.each_col() /= arma::sum(tau_, 1);
    }
  }

public:
  NormalBlockMeanUnknownClusters(const NormalBlockData& data,
                                 const arma::mat& B0, const arma::mat& Omega0,
                                 const arma::mat& tau0,
                                 double sparsity, const arma::mat& sparsity_weights,
                                 int fixed_point_niter, bool accelerate, bool fixed_tau,
                                 const std::string& cov_structure) :
    NormalBlockMeanBase(data, tau0.n_cols, B0, Omega0, sparsity, sparsity_weights, accelerate,
                        cov_structure),
    tau_(tau0), fixed_point_niter_(fixed_point_niter), fixed_tau_(fixed_tau) {
    alpha_  = arma::vectorise(arma::mean(tau_, 0));
    Psi_    = psi_from(Omega_, tau_);
    lambda_ = lambda_from(B_, tau_);
  }

  // Penalized ELBO of eq. (2.33). Phi is recomputed here rather than reused from the
  // M-step: it depends on the Omega/tau that step has just updated.
  double objective() const override {
    Phi_ = psi_from(Omega_, tau_) - tau_omega_tau(tau_);
    arma::mat M = (B_.t() * data_.XtX) * B_;
    arma::mat R_bar = data_.Y - (data_.X * B_) * tau_.t();

    double J = -0.5 * data_.n * data_.p * std::log(2.0 * arma::datum::pi)
           + 0.5 * data_.n * log_det_omega()
           + arma::accu(tau_ * arma::log(alpha_)) - arma::accu(tau_ % arma::log(tau_))
           - 0.5 * (arma::accu(R_bar % rmult_omega(R_bar)) + arma::trace(Phi_ * M));
    // penalized objective, as in the variance-block core: the R side adds the
    // penalty back through the `sparsity_term` active binding, so $loglik
    // stays the plain ELBO. Unlike there, no trace correction is needed: this
    // ELBO is written out in full and never assumes Omega = Sigma^-1.
    if (sparsity_ > 0.0) J -= sparsity_ * arma::accu(arma::abs(sparsity_weights_ % Omega_));
    return J;
  }

  const arma::mat& tau() const { return tau_; }
  const arma::vec& alpha() const { return alpha_; }
  const arma::mat& Psi() const { return Psi_; }
  const arma::mat& Phi() const { return Phi_; }
  arma::mat Lambda() const { return arma::diagmat(lambda_); }

  std::unique_ptr<NormalBlockEMBase> clone() const override {
    return std::make_unique<NormalBlockMeanUnknownClusters>(*this);
  }
  void restore_from(const NormalBlockEMBase& other) override {
    const auto& o = static_cast<const NormalBlockMeanUnknownClusters&>(other);
    copy_tracked_state_from(o);
    tau_ = o.tau_; alpha_ = o.alpha_; Psi_ = o.Psi_; lambda_ = o.lambda_; Phi_ = o.Phi_;
  }
};

#endif // NORMALBLOCKR_NORMAL_BLOCK_MEAN_UNKNOWN_CLUSTERS_H
