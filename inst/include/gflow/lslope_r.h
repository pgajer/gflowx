#ifndef LSLOPE_R_H
#define LSLOPE_R_H

#include <R.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

	/**
	 * @brief SEXP interface for gradient-restricted local slope (instrumented)
	 *
	 * Computes asymmetric association measures along gradient direction with
	 * full diagnostic output.
	 *
	 * @param s_adj_list R list of integer vectors (0-based adjacency)
	 * @param s_weight_list R list of numeric vectors (edge weights)
	 * @param s_y R numeric vector of directing function values
	 * @param s_z R numeric vector of response function values
	 * @param s_type R character: "slope", "normalized", or "sign"
	 * @param s_y_diff_type R character: "difference" or "logratio"
	 * @param s_z_diff_type R character: "difference" or "logratio"
	 * @param s_epsilon R numeric: pseudocount (0 = adaptive)
	 * @param s_sigmoid_alpha R numeric: sigmoid scale (0 = auto-calibrate)
	 * @param s_ascending R logical: use ascending (TRUE) or descending (FALSE) gradient
	 *
	 * @return R list with coefficients and diagnostics
	 */
	SEXP S_lslope_gradient_instrumented(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y,
		SEXP s_z,
		SEXP s_type,
		SEXP s_y_diff_type,
		SEXP s_z_diff_type,
		SEXP s_epsilon,
		SEXP s_sigmoid_alpha,
		SEXP s_sigmoid_type,
		SEXP s_ascending
		);

	/**
	 * @brief SEXP interface for gradient-restricted local slope (production)
	 *
	 * Streamlined version returning only coefficient vector.
	 */
	SEXP S_lslope_gradient(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y,
		SEXP s_z,
		SEXP s_type,
		SEXP s_y_diff_type,
		SEXP s_z_diff_type,
		SEXP s_epsilon,
		SEXP s_sigmoid_alpha,
		SEXP s_sigmoid_type,
		SEXP s_ascending
		);

	/**
	 * @brief SEXP interface for neighborhood local regression coefficient
	 *
	 * Computes local regression coefficient using all neighborhood edges.
	 *
	 * @param s_adj_list R list of integer vectors (0-based adjacency)
	 * @param s_weight_list R list of numeric vectors (edge weights)
	 * @param s_y R numeric vector of directing function values
	 * @param s_z R numeric vector of response function values
	 * @param s_weight_type R character: "unit" or "derivative"
	 * @param s_y_diff_type R character: "difference" or "logratio"
	 * @param s_z_diff_type R character: "difference" or "logratio"
	 * @param s_epsilon R numeric: pseudocount (0 = adaptive)
	 *
	 * @return R list with coefficients and diagnostics
	 */
	SEXP S_lslope_neighborhood(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y,
		SEXP s_z,
		SEXP s_weight_type,
		SEXP s_y_diff_type,
		SEXP s_z_diff_type,
		SEXP s_epsilon
		);

	/**
	 * @brief Compute local slope between vector y and matrix Z
	 *
	 * @param s_adj_list R list of integer vectors (0-based adjacency)
	 * @param s_weight_list R list of numeric vectors (edge weights)
	 * @param s_y R numeric vector of directing function values
	 * @param s_Z R numeric matrix of response function values
	 * @param s_type R character: "slope", "normalized", or "sign"
	 * @param s_y_diff_type R character: "difference" or "logratio"
	 * @param s_z_diff_type R character: "difference" or "logratio"
	 * @param s_epsilon R numeric: pseudocount (0 = adaptive)
	 * @param s_sigmoid_alpha R numeric: sigmoid scale (0 = auto-calibrate)
	 * @param s_ascending R logical: use ascending gradient
	 * @param s_n_threads R integer: number of OpenMP threads
	 *
	 * @return R list with coefficients matrix and metadata
	 */
	SEXP S_lslope_vector_matrix(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y,
		SEXP s_Z,
		SEXP s_type,
		SEXP s_y_diff_type,
		SEXP s_z_diff_type,
		SEXP s_epsilon,
		SEXP s_sigmoid_alpha,
		SEXP s_ascending,
		SEXP s_n_threads
		);

	/**
	 * @brief SEXP interface for memory-efficient lslope with posterior
	 *
	 * Computes local slope (lslope) statistics with posterior uncertainty
	 * propagation from spectral smoothing. Memory-efficient implementation
	 * that processes features sequentially.
	 *
	 * @param s_adj_list Adjacency list (0-based indices)
	 * @param s_weight_list Edge weight list
	 * @param s_y_hat Smoothed response values (length n)
	 * @param s_Z_abundances Original feature matrix (n x p)
	 * @param s_V Eigenvector matrix (n x m)
	 * @param s_eigenvalues Eigenvalue vector (length m)
	 * @param s_filter_type Filter type string
	 * @param s_eta_default Default eta value
	 * @param s_lslope_type Integer: 0=slope, 1=normalized, 2=sign
	 * @param s_ascending Logical for gradient direction
	 * @param s_y_diff_type Edge difference type for y ("difference" or "logratio")
	 * @param s_z_diff_type Edge difference type for z ("difference" or "logratio")
	 * @param s_per_column_gcv Logical for per-column GCV selection
	 * @param s_n_samples Number of posterior samples
	 * @param s_credible_level Credible interval level
	 * @param s_seed Random seed
	 * @param s_n_cores Number of OpenMP threads
	 * @param s_verbose Logical for progress output
	 *
	 * @return R list with components:
	 *   - mean: Matrix (p x n) of posterior mean lslope
	 *   - sd: Matrix (p x n) of posterior SD
	 *   - lower: Matrix (p x n) of lower credible bounds
	 *   - upper: Matrix (p x n) of upper credible bounds
	 *   - eta.used: Vector (p) of smoothing parameters
	 *   - effective.df: Vector (p) of effective degrees of freedom
	 */
	SEXP S_lslope_with_posterior_internal(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y_hat,
		SEXP s_Z_abundances,
		SEXP s_V,
		SEXP s_eigenvalues,
		SEXP s_filter_type,
		SEXP s_eta_default,
		SEXP s_lslope_type,
		SEXP s_ascending,
		SEXP s_y_diff_type,
		SEXP s_z_diff_type,
		SEXP s_per_column_gcv,
		SEXP s_n_samples,
		SEXP s_credible_level,
		SEXP s_seed,
		SEXP s_n_cores,
		SEXP s_verbose
		);

#ifdef __cplusplus
}
#endif

#endif // LSLOPE_R_H
