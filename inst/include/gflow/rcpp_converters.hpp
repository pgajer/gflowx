#pragma once

#include <Rcpp.h>
#include <vector>
#include <algorithm>

namespace gflow {
namespace rcpp {

// R matrix (n x d) -> std::vector<std::vector<double>> size n, each length d
inline std::vector<std::vector<double>>
matrix_to_vecvec(const Rcpp::NumericMatrix& X) {
  const int n = X.nrow();
  const int d = X.ncol();
  std::vector<std::vector<double>> out;
  out.resize(n);
  for (int i = 0; i < n; ++i) {
    out[i].resize(d);
    for (int j = 0; j < d; ++j)
      out[i][j] = X(i, j);
  }
  return out;
}

// 3-level trajectory (steps x n x d) -> list of n x d numeric matrices
inline Rcpp::List
traj_to_list_of_mats(const std::vector<std::vector<std::vector<double>>>& X_traj) {
  const std::size_t T = X_traj.size();
  Rcpp::List out(T);
  for (std::size_t t = 0; t < T; ++t) {
    const auto& Xt = X_traj[t];
    if (Xt.empty()) {
      out[t] = Rcpp::NumericMatrix(0, 0);
      continue;
    }
    const int n = static_cast<int>(Xt.size());
    const int d = static_cast<int>(Xt[0].size());
    Rcpp::NumericMatrix M(n, d);
    for (int i = 0; i < n; ++i) {
      const int dd = std::min(d, static_cast<int>(Xt[i].size()));
      for (int j = 0; j < dd; ++j) M(i, j) = Xt[i][j];
    }
    out[t] = M;
  }
  return out;
}

// Wrap results struct into an R list
template <class ResultsT>
inline Rcpp::List results_to_R_list(const ResultsT& res) {
  Rcpp::NumericVector med(res.median_kdistances.begin(), res.median_kdistances.end());
  Rcpp::List traj = traj_to_list_of_mats(res.X_traj);
  return Rcpp::List::create(
    Rcpp::Named("X_traj") = traj,
    Rcpp::Named("median_kdistances") = med
  );
}

} // namespace rcpp
} // namespace gflow
