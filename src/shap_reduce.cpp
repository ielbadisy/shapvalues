// Parallel marginal-contribution reduction for sampling (permutation) SHAP.
//
// The R side issues a single batched prediction call over all coalition rows
// and hands the results here as a matrix. This file only does the arithmetic
// of turning those predictions into per-observation, per-feature attributions,
// which is the part worth parallelising.

#include <Rcpp.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

// Layout of `preds` (N rows x O outputs), with N = n_obs * nsim * (p + 1):
//   block index         b   = i * nsim + m           (i = obs, m = permutation)
//   row of coalition k  is  b * (p + 1) + k          (k = 0 .. p)
//   perm_order(b, k-1)  is  the 0-based feature added at step k
//
// phi(i, f, o) = mean over m of sum over steps k where perm_order(b,k-1)==f of
//                ( preds(b*(p+1)+k, o) - preds(b*(p+1)+k-1, o) )
//
// Output is a flat vector of length n_obs * p * O in column-major (i fastest,
// then f, then o) so R can wrap it with dim = c(n_obs, p, O).

struct ShapReducer : public Worker {
  const RMatrix<double> preds;
  const RMatrix<int> perm_order;   // (n_obs*nsim) x p
  const std::size_t n_obs;
  const std::size_t nsim;
  const std::size_t p;
  const std::size_t n_out;
  RVector<double> phi;             // length n_obs * p * n_out

  ShapReducer(const NumericMatrix preds, const IntegerMatrix perm_order,
              std::size_t n_obs, std::size_t nsim, std::size_t p,
              std::size_t n_out, NumericVector phi)
    : preds(preds), perm_order(perm_order), n_obs(n_obs), nsim(nsim),
      p(p), n_out(n_out), phi(phi) {}

  inline std::size_t phi_idx(std::size_t i, std::size_t f, std::size_t o) const {
    return i + n_obs * (f + p * o);
  }

  void operator()(std::size_t begin, std::size_t end) {
    const std::size_t L = p + 1;
    for (std::size_t i = begin; i < end; ++i) {
      for (std::size_t m = 0; m < nsim; ++m) {
        const std::size_t b = i * nsim + m;
        const std::size_t row0 = b * L;
        for (std::size_t k = 1; k <= p; ++k) {
          const int f = perm_order(b, k - 1);
          if (f < 0 || static_cast<std::size_t>(f) >= p) continue;
          const std::size_t r_hi = row0 + k;
          const std::size_t r_lo = row0 + k - 1;
          for (std::size_t o = 0; o < n_out; ++o) {
            phi[phi_idx(i, static_cast<std::size_t>(f), o)] +=
              preds(r_hi, o) - preds(r_lo, o);
          }
        }
      }
      // average over permutations
      for (std::size_t f = 0; f < p; ++f) {
        for (std::size_t o = 0; o < n_out; ++o) {
          phi[phi_idx(i, f, o)] /= static_cast<double>(nsim);
        }
      }
    }
  }
};

// [[Rcpp::export]]
NumericVector shap_reduce_cpp(NumericMatrix preds, IntegerMatrix perm_order,
                              int n_obs, int nsim, int p, int n_out,
                              int grain_size = 1) {
  NumericVector phi(static_cast<R_xlen_t>(n_obs) * p * n_out);
  ShapReducer worker(preds, perm_order,
                     static_cast<std::size_t>(n_obs),
                     static_cast<std::size_t>(nsim),
                     static_cast<std::size_t>(p),
                     static_cast<std::size_t>(n_out), phi);
  parallelFor(0, static_cast<std::size_t>(n_obs), worker,
              static_cast<std::size_t>(grain_size < 1 ? 1 : grain_size));
  phi.attr("dim") = IntegerVector::create(n_obs, p, n_out);
  return phi;
}

// Column means of the all-background rows (coalition step 0 of every block):
// used for the SHAP baseline E[f(x)]. Returns length n_out.
// [[Rcpp::export]]
NumericVector shap_baseline_cpp(NumericMatrix preds, int n_obs, int nsim, int p) {
  const std::size_t L = static_cast<std::size_t>(p) + 1;
  const std::size_t nblk = static_cast<std::size_t>(n_obs) * nsim;
  const std::size_t n_out = preds.ncol();
  NumericVector out(n_out);
  for (std::size_t b = 0; b < nblk; ++b) {
    for (std::size_t o = 0; o < n_out; ++o) out[o] += preds(b * L, o);
  }
  for (std::size_t o = 0; o < n_out; ++o) out[o] /= static_cast<double>(nblk);
  return out;
}
