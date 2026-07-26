#include <Rcpp.h>
#include <vector>
#include <memory>

#include "mean_shift_smoother.hpp"
#include "rcpp_converters.hpp"

using gflow::rcpp::matrix_to_vecvec;
using gflow::rcpp::traj_to_list_of_mats;
using gflow::rcpp::results_to_R_list;

using Rcpp::IntegerVector;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::List;
using Rcpp::Named;

//' Fully Adaptive Mean Shift with Gradient Field Averaging (C++ engine)
//'
//' Low-level Rcpp wrapper. Prefer calling the R wrapper
//' \code{adaptive_mean_shift_gfa()} for a friendlier interface.
//'
//' @param X numeric matrix (n x d).
//' @param k integer, k-NN for gradient estimation (excl. self).
//' @param density_k integer, k-NN for density estimation (excl. self).
//' @param n_steps integer, number of iterations.
//' @param initial_step_size positive numeric, initial per-point step size.
//' @param ikernel Type of kernel function to use (default: 1 - Epanechnikov).
//'               Available kernels:
//'               - 0-Constant,
//'               - 1-Epanechnikov,
//'               - 2-Triangular,
//'               - 3-TrExponential,
//'               - 4-Laplace,
//'               - 5-Normal,
//'               - 6-Biweight,
//'               - 7-Tricube,
//'               - 8-Cosine
//'               - 9-Hyperbolic
//' @param dist_normalization_factor numeric > 0.
//' @param average_direction_only logical; if TRUE, average directions only.
//' @param momentum numeric, gradient momentum (default 0.9).
//' @param increase_factor numeric, step up factor (default 1.2).
//' @param decrease_factor numeric, step down factor (default 0.5).
//'
//' @return list with \code{X_traj} (list of n x d matrices over iterations)
//'   and \code{median_kdistances} (numeric).
// [[Rcpp::export]]
Rcpp::List rcpp_adaptive_mean_shift_gfa(
    const NumericMatrix& X,
    int k,
    int density_k,
    int n_steps,
    double initial_step_size,
    int ikernel = 1,
    double dist_normalization_factor = 1.01,
    bool average_direction_only = false,
    double momentum = 0.9,
    double increase_factor = 1.2,
    double decrease_factor = 0.5
) {
  if (X.nrow() == 0 || X.ncol() == 0)
    Rcpp::stop("X must be a non-empty numeric matrix.");
  if (k <= 0 || density_k <= 0)
    Rcpp::stop("k and density_k must be positive.");
  if (n_steps <= 0)
    Rcpp::stop("n_steps must be positive.");
  if (initial_step_size <= 0.0)
    Rcpp::stop("initial_step_size must be positive.");

  auto X_cpp = matrix_to_vecvec(X);

  std::unique_ptr<mean_shift_smoother_results_t> res =
    adaptive_mean_shift_data_smoother_with_grad_field_averaging(
      X_cpp, k, density_k, n_steps, initial_step_size,
      ikernel, dist_normalization_factor, average_direction_only,
      momentum, increase_factor, decrease_factor
    );

  if (!res)
    Rcpp::stop("C++ returned null results pointer.");

  return results_to_R_list(*res);
}

//' kNN-Adaptive Mean Shift with Gradient Field Averaging (C++ engine)
//'
//' Low-level Rcpp wrapper. Prefer calling the R wrapper
//' \code{knn_adaptive_mean_shift_gfa()} for a friendlier interface.
//'
//' @param X numeric matrix (n x d).
//' @param k integer, k-NN for gradient estimation (excl. self).
//' @param density_k integer, k-NN for density estimation (excl. self).
//' @param n_steps integer, number of iterations.
//' @param step_size positive numeric step size (constant over iterations).
//' @param ikernel Type of kernel function to use (default: 1 - Epanechnikov).
//'               Available kernels:
//'               - 0-Constant,
//'               - 1-Epanechnikov,
//'               - 2-Triangular,
//'               - 3-TrExponential,
//'               - 4-Laplace,
//'               - 5-Normal,
//'               - 6-Biweight,
//'               - 7-Tricube,
//'               - 8-Cosine
//'               - 9-Hyperbolic
//' @param dist_normalization_factor numeric > 0.
//' @param average_direction_only logical; if TRUE, average directions only.
//'
//' @return list with \code{X_traj} (list of n x d matrices over iterations)
//'   and \code{median_kdistances} (numeric).
// [[Rcpp::export]]
Rcpp::List rcpp_knn_adaptive_mean_shift_gfa(
    const NumericMatrix& X,
    int k,
    int density_k,
    int n_steps,
    double step_size,
    int ikernel = 1,
    double dist_normalization_factor = 1.01,
    bool average_direction_only = false
) {
  if (X.nrow() == 0 || X.ncol() == 0)
    Rcpp::stop("X must be a non-empty numeric matrix.");
  if (k <= 0 || density_k <= 0)
    Rcpp::stop("k and density_k must be positive.");
  if (n_steps <= 0)
    Rcpp::stop("n_steps must be positive.");
  if (step_size <= 0.0)
    Rcpp::stop("step_size must be positive.");

  auto X_cpp = matrix_to_vecvec(X);

  std::unique_ptr<mean_shift_smoother_results_t> res =
    knn_adaptive_mean_shift_data_smoother_with_grad_field_averaging(
      X_cpp, k, density_k, n_steps, step_size,
      ikernel, dist_normalization_factor, average_direction_only
    );

  if (!res)
    Rcpp::stop("C++ returned null results pointer.");

  return results_to_R_list(*res);
}
