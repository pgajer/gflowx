/**
 * @file lcor_posterior_r.h
 * @brief Declaration of S_lcor_with_posterior_internal for R registration
 */

#ifndef LCOR_POSTERIOR_R_H
#define LCOR_POSTERIOR_R_H

#include <R.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Memory-efficient lcor computation with posterior uncertainty propagation
 *
 * Computes local correlation statistics with full posterior uncertainty
 * quantification without materializing all posterior samples in memory. For
 * each feature, posterior samples of the smoothed values are generated, lcor
 * is computed for each sample, summary statistics are computed, and samples
 * are discarded.
 *
 * @param s_adj_list Adjacency list (0-based indices)
 * @param s_weight_list Edge weight list
 * @param s_y_hat Smoothed response vector (length n)
 * @param s_Z Original (unsmoothed) feature matrix (n x p)
 * @param s_V Eigenvector matrix from fitted model (n x m)
 * @param s_eigenvalues Raw eigenvalues (length m)
 * @param s_eta_fixed Fixed eta value (used when per_column_gcv = FALSE)
 * @param s_lcor_type "derivative", "unit", or "sign"
 * @param s_filter_type Filter type string ("heat_kernel", "tikhonov", etc.)
 * @param s_per_column_gcv Logical: select eta per column via GCV
 * @param s_n_gcv_candidates Number of GCV candidates
 * @param s_n_posterior_samples Number of posterior samples per feature
 * @param s_credible_level Credible level (e.g., 0.95)
 * @param s_seed Base random seed
 * @param s_n_cores Number of OpenMP threads
 * @param s_verbose Logical: print progress
 *
 * @return R list with:
 *   - mean: Matrix (p x n) of posterior mean lcor values
 *   - sd: Matrix (p x n) of posterior standard deviations
 *   - lower: Matrix (p x n) of lower credible bounds
 *   - upper: Matrix (p x n) of upper credible bounds
 *   - eta.used: Vector (length p) of eta values used
 *   - effective.df: Vector (length p) of effective degrees of freedom
 *   - n.samples: Number of posterior samples used
 *   - credible.level: Credible level used
 */
SEXP S_lcor_with_posterior_internal(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y_hat,
    SEXP s_Z,
    SEXP s_V,
    SEXP s_eigenvalues,
    SEXP s_eta_fixed,
    SEXP s_lcor_type,
    SEXP s_filter_type,
    SEXP s_per_column_gcv,
    SEXP s_n_gcv_candidates,
    SEXP s_n_posterior_samples,
    SEXP s_credible_level,
    SEXP s_seed,
    SEXP s_n_cores,
    SEXP s_verbose
);

#ifdef __cplusplus
}
#endif

#endif /* LCOR_POSTERIOR_R_H */
