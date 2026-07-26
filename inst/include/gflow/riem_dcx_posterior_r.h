/**
 * @file riem_dcx_posterior_r.h
 * @brief Declaration of S_compute_posterior_summary for R registration
 */

#ifndef RIEM_DCX_POSTERIOR_R_H
#define RIEM_DCX_POSTERIOR_R_H

#include <R.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Compute posterior summary for spectral-filtered regression
 *
 * @param s_V Eigenvector matrix (n x m)
 * @param s_eigenvalues Raw eigenvalues (length m)
 * @param s_filtered_eigenvalues Filter weights (length m)
 * @param s_y Original response vector (length n)
 * @param s_y_hat Fitted values (length n)
 * @param s_eta Smoothing parameter (scalar)
 * @param s_credible_level Coverage probability (scalar)
 * @param s_n_samples Number of Monte Carlo samples (integer)
 * @param s_seed Random seed (integer)
 * @param s_return_samples Whether to return samples (logical)
 *
 * @return R list with posterior summary components:
 *   - lower: Vector of lower credible bounds
 *   - upper: Vector of upper credible bounds
 *   - sd: Vector of posterior standard deviations
 *   - credible.level: Coverage probability
 *   - sigma: Estimated residual standard deviation
 *   - samples: (optional) Matrix of posterior samples
 */
SEXP S_compute_posterior_summary(
    SEXP s_V,
    SEXP s_eigenvalues,
    SEXP s_filtered_eigenvalues,
    SEXP s_y,
    SEXP s_y_hat,
    SEXP s_eta,
    SEXP s_credible_level,
    SEXP s_n_samples,
    SEXP s_seed,
    SEXP s_return_samples
);

#ifdef __cplusplus
}
#endif

#endif /* RIEM_DCX_POSTERIOR_R_H */
