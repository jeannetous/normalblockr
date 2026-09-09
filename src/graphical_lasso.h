#ifndef NORMALBLOCKR_GRAPHICAL_LASSO_H
#define NORMALBLOCKR_GRAPHICAL_LASSO_H

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <cstddef>

// In-package graphical lasso: minimize over a positive definite Theta
//
//   -log det(Theta) + tr(S Theta) + || L o Theta ||_1
//
// by the block coordinate descent of Friedman, Hastie & Tibshirani (2008),
// in the bookkeeping of Sustik & Calderhead (2012), "GLASSOFAST: An efficient
// GLASSO implementation" (TR-12-29, UT Austin).
//
// This is a direct port of the `glassofast` Fortran subroutine the glassoFast
// package ships, which the package used to call back into R for, once per
// M-step, from inside the C++ (V)EM loop (thousands of R round trips for a
// single sparsity path). Two things motivated bringing it in-house:
//
//  * that callback was a live memory-safety hazard -- one bug there was
//    already found and fixed (a `static Rcpp::Function` holding an R object
//    outside R's protection discipline), and CI kept crashing with heap
//    corruption while rebuilding the sparsity-path vignette;
//  * having the solver in-house lets it be warm-started (nb_glasso::State),
//    which the R interface made awkward.
//
// The port is faithful except on two points, both deliberate:
//
//  1. When S carries no off-diagonal mass the problem separates exactly and
//     Theta is diagonal. glassoFast returns 1 / max(L_ii, eps) there, dropping
//     S_ii entirely: with an unpenalized diagonal (our default) it hands
//     back ~9.09e15 instead of 1 / S_ii. That is a bug, not a convention, and
//     it is reachable: any 1 x 1 problem (q = 1) takes this branch. We return
//     the correct 1 / (S_ii + L_ii).
//  2. The inner coordinate descent is guaranteed to terminate. In the Fortran
//     it is an unbounded `do` loop whose only exit is `dlx < thrLasso`, which
//     is never true once a NaN reaches `dlx`, hence, a non-finite input hangs the
//     process rather than failing. We reject non-finite input and
//     non-positive S_ii + L_ii up front (together the only ways to reach that
//     state), break on a non-finite `dlx`, and keep a far-off backstop cap, so
//     a caller gets a result it can test. `converged` reports which happened.
// Fortran guarantees non-aliasing arrays and vectorizes the rank-1 updates
// below on that basis; C++ has to be told, or the compiler keeps them scalar.
// The pointers it is applied to are always distinct allocations (a standalone
// arma::vec and the columns of a matrix).
#if defined(__GNUC__) || defined(__clang__)
  #define NB_RESTRICT __restrict__
#elif defined(_MSC_VER)
  #define NB_RESTRICT __restrict
#else
  #define NB_RESTRICT
#endif

namespace nb_glasso {

// The Fortran's EPS parameter, kept to the digit for comparability.
constexpr double kEps = 1.1e-16;

// Last-resort backstop on the inner coordinate descent, not a working limit:
// termination is guaranteed structurally instead (non-finite input and a
// non-positive S_ii + L_ii are both rejected up front, and a non-finite dlx
// breaks the loop). It is set far above what a well-posed problem needs:
// weak penalties genuinely take thousands of passes, up to ~15k measured over
// a 432-case sweep, and an earlier 10k cap silently degraded two of them.
constexpr int kMaxInner = 500000;

struct Result {
  arma::mat W;            // covariance estimate (glassoFast's `w`)
  arma::mat X;            // precision estimate  (glassoFast's `wi`)
  int niter = 0;          // outer sweeps performed
  bool converged = true;  // false if a sweep cap was hit
};

// Persistent (W, X) carried between calls to warm-start the next solve. The
// pair is only reused when it still has the right size; anything else (a
// changed q, a first call) silently falls back to a cold start.
//
// Used by the (V)EM between M-steps (nb_omega::estimate), which is where it
// pays: the outer loop stops on `dw <= shr`, how much a whole sweep moved W
// rather than how far W still is from the optimum, so resuming both takes
// fewer sweeps and, at the tightened threshold that affords, lands closer to
// the optimum than a cold start at the looser one.
//
// A warm start is never load-bearing, though. On an ill-conditioned Sigma a
// bad one can send the coordinate descent below to infinity where a cold
// start on the same problem converges in a few sweeps, so nb_omega::estimate()
// retries cold on a non-finite result rather than treating it as failure.
struct State {
  arma::mat W;
  arma::mat X;
  bool filled = false;

  bool usable_for(arma::uword n) const {
    return filled && W.n_rows == n && W.n_cols == n && X.n_rows == n && X.n_cols == n
           && W.is_finite() && X.is_finite();
  }
  void store(const Result& r) { W = r.W; X = r.X; filled = true; }
  void reset() { filled = false; }
};

// `warm` is used only when it is `usable_for(S.n_rows)`; pass nullptr for a
// cold start. Mirrors glassoFast's defaults (thr = 1e-4, max_iter = 10000).
inline Result solve(const arma::mat& S, const arma::mat& L,
                    double thr = 1e-4, int max_iter = 10000,
                    const State* warm = nullptr) {
  const arma::uword n = S.n_rows;

  Result out;
  out.W.zeros(n, n);
  out.X.zeros(n, n);
  if (n == 0) return out;

  if (!S.is_finite() || !L.is_finite()) {
    out.W.fill(arma::datum::nan);
    out.X.fill(arma::datum::nan);
    out.converged = false;
    return out;
  }

  // S_ii + L_ii is what the soft-threshold below divides by, and the only way
  // a finite input could manufacture an infinity there; a non-positive entry
  // is a zero-variance coordinate, i.e. degenerate input. Checked before the
  // separable branch so that both paths answer the same way rather than that
  // one quietly returning 1 / eps.
  const arma::vec diag_sum = S.diag() + L.diag();
  if (!arma::all(diag_sum > 0.0)) {
    out.W.fill(arma::datum::nan);
    out.X.fill(arma::datum::nan);
    out.converged = false;
    return out;
  }

  arma::mat& W = out.W;
  arma::mat& X = out.X;

  // Total off-diagonal absolute mass of S; sets both convergence thresholds.
  const double off_mass = arma::accu(arma::abs(S)) - arma::accu(arma::abs(S.diag()));

  if (off_mass <= 0.0) {
    // Separable: no coupling to estimate, so the exact solution is diagonal.
    // (See note 1 above -- this is where we part with glassoFast.)
    W.diag() = diag_sum;
    X.diag() = 1.0 / diag_sum;
    return out;
  }

  const double shr = thr * off_mass / static_cast<double>(n - 1);
  const double thr_lasso = std::max(shr / static_cast<double>(n), 2.0 * kEps);

  if (warm != nullptr && warm->usable_for(n)) {
    // The recursion carries X as the negated normalized regression
    // coefficients of each column, not as the precision matrix; a warm start
    // has to be pushed back into that representation first.
    W = warm->W;
    X = warm->X;
    for (arma::uword i = 0; i < n; ++i) {
      const double xii = X(i, i);
      X.col(i) /= -xii;
      X(i, i) = 0.0;
    }
    if (!X.is_finite()) { // a singular warm start (some X_ii == 0)
      W = S;
      X.zeros();
    }
  } else {
    W = S;
    X.zeros();
  }

  arma::vec Wd(n);
  for (arma::uword i = 0; i < n; ++i) {
    Wd(i) = S(i, i) + L(i, i);
    W(i, i) = Wd(i);
  }

  arma::vec WXj(n);
  int iter = 0;
  bool outer_converged = false;

  // The rest of this function goes through raw column pointers rather than
  // Armadillo element access. This is the hot loop -- a rank-1 update of WXj
  // per accepted coordinate, and the weakest penalties take thousands of
  // passes over it -- and the package does not set ARMA_NO_DEBUG, so every
  // `X(i, j)` here would otherwise carry a bounds check. Measured 7-10x on
  // n = 50-150. Bounds checking stays on everywhere else.
  const double* const NB_RESTRICT Wd_p = Wd.memptr();
  double* const NB_RESTRICT WXj_p = WXj.memptr();

  for (iter = 1; iter <= max_iter; ++iter) {
    double dw = 0.0;

    for (arma::uword j = 0; j < n; ++j) {
      double* const NB_RESTRICT Xj = X.colptr(j);
      const double* const NB_RESTRICT Sj = S.colptr(j);
      const double* const NB_RESTRICT Lj = L.colptr(j);

      // WXj = W * X.col(j), skipping the zeros X is expected to be full of
      std::fill(WXj_p, WXj_p + n, 0.0);
      for (arma::uword i = 0; i < n; ++i) {
        const double xij = Xj[i];
        if (xij != 0.0) {
          const double* const NB_RESTRICT Wi = W.colptr(i);
          for (arma::uword k = 0; k < n; ++k) WXj_p[k] += Wi[k] * xij;
        }
      }

      int inner = 0;
      for (;;) {
        double dlx = 0.0;
        for (arma::uword i = 0; i < n; ++i) {
          if (i == j) continue;
          const double a = Sj[i] - WXj_p[i] + Wd_p[i] * Xj[i];
          const double b = std::fabs(a) - Lj[i];
          // soft-threshold; matches Fortran's sign(b, a), which returns +|b|
          // when a is zero (including -0.0)
          const double c = (b > 0.0) ? ((a >= 0.0 ? b : -b) / Wd_p[i]) : 0.0;
          const double delta = c - Xj[i];
          if (delta != 0.0) {
            Xj[i] = c;
            const double* const NB_RESTRICT Wi = W.colptr(i);
            for (arma::uword k = 0; k < n; ++k) WXj_p[k] += Wi[k] * delta;
            const double ad = std::fabs(delta);
            if (ad > dlx) dlx = ad;
          }
        }
        if (dlx < thr_lasso) break;
        if (!std::isfinite(dlx) || ++inner >= kMaxInner) { out.converged = false; break; }
      }

      WXj_p[j] = Wd_p[j];
      double* const NB_RESTRICT Wj = W.colptr(j);
      double acc = 0.0;
      for (arma::uword k = 0; k < n; ++k) acc += std::fabs(WXj_p[k] - Wj[k]);
      if (acc > dw) dw = acc;
      for (arma::uword k = 0; k < n; ++k) Wj[k] = WXj_p[k];
      for (arma::uword k = 0; k < n; ++k) W.colptr(k)[j] = WXj_p[k]; // W(j, :)
    }

    if (dw <= shr) { outer_converged = true; break; }
  }

  out.niter = std::min(iter, max_iter);
  if (!outer_converged) out.converged = false;

  // Back out the precision matrix from the regression coefficients. X(i,i) is
  // still 0 here, so it drops out of the dot product on its own.
  for (arma::uword i = 0; i < n; ++i) {
    const double tmp = 1.0 / (Wd(i) - arma::dot(X.col(i), W.col(i)));
    X.col(i) *= -tmp;
    X(i, i) = tmp;
  }

  const arma::mat Xt = X.t(); // averaging the two triangles leaves the diagonal as is
  X = 0.5 * (X + Xt);

  return out;
}

} // namespace nb_glasso

#undef NB_RESTRICT

#endif // NORMALBLOCKR_GRAPHICAL_LASSO_H
