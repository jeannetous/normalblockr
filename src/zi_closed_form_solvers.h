#ifndef NORMALBLOCKR_ZI_CLOSED_FORM_SOLVERS_H
#define NORMALBLOCKR_ZI_CLOSED_FORM_SOLVERS_H

#include <RcppArmadillo.h>

// The B- and M-update subproblems of the zero-inflated normal-block models
// have no closed form in general nonlinear optimization. However, both
// objectives are exactly quadratic in the parameter being optimized (the
// zero-inflation mask only reweights residuals, it does not introduce any
// nonlinearity), so each one reduces to a (weighted) linear system: no
// iterative optimizer is needed at all.
namespace nb_optim {

// Solves argmax_Theta -.5 * sum(W % (Y - Design*Theta - Offset)^2), i.e. the
// B-update subproblem shared by every zero-inflated variant. W varies both
// per row (zero-inflation mask) and per column (dm1), so the normal
// equations do not share a common cross-product across columns: each column
// of Theta is solved independently via its own weighted least squares.
inline arma::mat solve_wls(const arma::mat& W, const arma::mat& Y,
                            const arma::mat& Design, const arma::mat& Offset) {
  arma::uword d = Design.n_cols, p = Y.n_cols;
  arma::mat Theta(d, p);
  arma::mat Z = Y - Offset;
  // Buffers hoisted out of the loop, and X'WX solved as the symmetric
  // positive definite system it is (W >= 0, Design full rank): Armadillo's
  // general solver re-checks conditioning on every call, and at these sizes
  // that overhead is the cost, not the arithmetic. Same reasoning as
  // solve_M_ridge() below.
  arma::mat WX(Design.n_rows, d), XtWX(d, d);
  arma::vec XtWz(d), theta(d);
  for (arma::uword j = 0; j < p; ++j) {
    WX = Design;
    WX.each_col() %= W.col(j);
    XtWX = Design.t() * WX;
    XtWz = Design.t() * (W.col(j) % Z.col(j));
    if (!arma::solve(theta, XtWX, XtWz, arma::solve_opts::likely_sympd + arma::solve_opts::fast))
      theta = arma::solve(XtWX, XtWz); // degenerate column: take the safe path
    Theta.col(j) = theta;
  }
  return Theta;
}

// Solves argmax_M -.5 * ( sum(DM1 % ((M^2)*C^T - 2*R%(M*C^T))) + sum((M*Omega)%M) ),
// i.e. the variational-mean (M) subproblem of the unknown-cluster
// zero-inflated model. The objective is row-separable (mirrors the
// known-cluster E-step's per-row Gamma_i/Mu_i): each row solves an
// independent q x q ridge system M_i = (Omega + diag(DM1C_i))^{-1} (DM1*R*C)_i.
inline arma::mat solve_M_ridge(const arma::mat& DM1, const arma::mat& R,
                                  const arma::mat& C, const arma::mat& Omega) {
  arma::uword n = DM1.n_rows, q = C.n_cols;
  arma::mat DM1C = DM1 * C;          // n x q
  arma::mat DM1RC = (DM1 % R) * C;   // n x q
  arma::mat M(n, q);
  // One q x q system per row, all sharing Omega and differing only in the
  // diagonal added to it. There is no factorization to share across rows, but
  // no reason to pay for a fresh allocation and a general solver either: the
  // system is symmetric positive definite (Omega is, and DM1C is
  // non-negative).
  const arma::vec omega_diag = Omega.diag();
  arma::mat A(q, q); arma::vec b(q), x(q);
  for (arma::uword i = 0; i < n; ++i) {
    A = Omega;
    A.diag() = omega_diag + DM1C.row(i).t();
    b = DM1RC.row(i).t();
    if (!arma::solve(x, A, b, arma::solve_opts::likely_sympd + arma::solve_opts::fast))
      x = arma::solve(A, b); // not observed in practice; stay correct if it happens
    M.row(i) = x.t();
  }
  return M;
}

} // namespace nb_optim

#endif // NORMALBLOCKR_ZI_CLOSED_FORM_SOLVERS_H
