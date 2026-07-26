/**
 * @file riem_dcx_posterior_r.cpp
 * @brief SEXP interface for posterior summary computation
 *
 * Provides R-callable interface for Bayesian posterior inference on
 * spectral-filtered regression estimates. This standalone implementation
 * enables posterior sampling from refit.rdgraph.regression() without
 * requiring a full riem_dcx_t object.
 */

#include <R.h>
#include <Rinternals.h>

#include <Eigen/Core>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

// Type aliases consistent with gflow
using vec_t = Eigen::VectorXd;

/**
 * @brief Posterior summary structure
 *
 * Contains vertex-wise credible interval bounds and posterior standard deviations
 * from Bayesian inference on the smoothed response surface.
 */
struct posterior_summary_result_t {
    vec_t lower;              // Lower credible bounds (length n)
    vec_t upper;              // Upper credible bounds (length n)
    vec_t posterior_sd;       // Posterior standard deviations (length n)
    double credible_level;    // Coverage probability
    double sigma_hat;         // Estimated residual SD

    // Optional: full posterior samples (n x B matrix)
    Eigen::MatrixXd samples;
    bool has_samples;
};

/**
 * @brief Compute posterior summary for spectral-filtered regression
 *
 * Generates Monte Carlo samples from the posterior distribution of the
 * smoothed response and computes vertex-wise credible intervals.
 *
 * The posterior distribution follows from a Bayesian interpretation of
 * spectral smoothing: coefficients in the eigenbasis have posterior
 * variance inversely proportional to (1 + eta * lambda_j), where lambda_j
 * is the j-th eigenvalue and eta is the smoothing parameter.
 *
 * @param V Eigenvector matrix (n x m) from spectral decomposition
 * @param eigenvalues Raw eigenvalues (length m)
 * @param filtered_eigenvalues Filter weights f(lambda; eta) (length m)
 * @param y Original response vector (length n)
 * @param y_hat Fitted values (length n)
 * @param eta Smoothing parameter
 * @param credible_level Coverage probability (e.g., 0.95)
 * @param n_samples Number of Monte Carlo samples
 * @param seed Random seed for reproducibility
 * @param return_samples If true, store full sample matrix
 *
 * @return posterior_summary_result_t with credible bounds and optional samples
 */
static posterior_summary_result_t compute_posterior_summary_impl(
    const Eigen::MatrixXd& V,
    const vec_t& eigenvalues,
    const vec_t& filtered_eigenvalues,
    const vec_t& y,
    const vec_t& y_hat,
    double eta,
    double credible_level,
    int n_samples,
    unsigned int seed,
    bool return_samples
) {
    const int n = V.rows();
    const int m = V.cols();

    // ================================================================
    // STEP 0: VALIDATE INPUTS
    // ================================================================

    if (y.size() != n) {
        throw std::invalid_argument(
            "Response vector y length (" + std::to_string(y.size()) +
            ") does not match number of vertices (" + std::to_string(n) + ")"
        );
    }

    if (y_hat.size() != n) {
        throw std::invalid_argument(
            "Fitted values y_hat length (" + std::to_string(y_hat.size()) +
            ") does not match number of vertices (" + std::to_string(n) + ")"
        );
    }

    if (eigenvalues.size() != m) {
        throw std::invalid_argument(
            "Eigenvalue vector length (" + std::to_string(eigenvalues.size()) +
            ") does not match number of eigenvectors (" + std::to_string(m) + ")"
        );
    }

    if (filtered_eigenvalues.size() != m) {
        throw std::invalid_argument(
            "Filtered eigenvalue vector length (" + std::to_string(filtered_eigenvalues.size()) +
            ") does not match number of eigenvectors (" + std::to_string(m) + ")"
        );
    }

    // ================================================================
    // STEP 1: ESTIMATE RESIDUAL VARIANCE
    // ================================================================

    vec_t residuals = y - y_hat;
    double eff_df = filtered_eigenvalues.sum();  // Effective degrees of freedom

    // Guard against degenerate case
    if (n - eff_df < 1.0) {
        throw std::runtime_error(
            "Cannot estimate residual variance: effective degrees of freedom (" +
            std::to_string(eff_df) + ") too close to sample size (" +
            std::to_string(n) + ")"
        );
    }

    double sigma_hat = std::sqrt(residuals.squaredNorm() / (n - eff_df));

    // ================================================================
    // STEP 2: COMPUTE POSTERIOR PARAMETERS IN SPECTRAL DOMAIN
    // ================================================================

    // Posterior mean: alpha_mean = filtered_eigenvalues .* (V^T y)
    vec_t Vt_y = V.transpose() * y;
    vec_t alpha_mean = filtered_eigenvalues.array() * Vt_y.array();

    // Posterior standard deviation for each coefficient
    // SD(alpha_j | y) = sigma_hat / sqrt(1 + eta * lambda_j)
    vec_t alpha_sd(m);
    for (int j = 0; j < m; ++j) {
        double posterior_precision = 1.0 + eta * eigenvalues[j];
        alpha_sd[j] = sigma_hat / std::sqrt(posterior_precision);
    }

    // ================================================================
    // STEP 3: GENERATE POSTERIOR SAMPLES
    // ================================================================

    // Setup random number generator
    std::mt19937 rng(seed);
    std::normal_distribution<double> std_normal(0.0, 1.0);

    // Storage for samples at each vertex
    Eigen::MatrixXd y_samples(n, n_samples);

    for (int b = 0; b < n_samples; ++b) {
        // Draw spectral coefficients from posterior
        vec_t alpha_sample(m);
        for (int j = 0; j < m; ++j) {
            alpha_sample[j] = alpha_mean[j] + alpha_sd[j] * std_normal(rng);
        }

        // Transform to vertex domain: y_sample = V * alpha_sample
        y_samples.col(b) = V * alpha_sample;
    }

    // ================================================================
    // STEP 4: COMPUTE QUANTILES AT EACH VERTEX
    // ================================================================

    vec_t lower(n);
    vec_t upper(n);
    vec_t posterior_sd(n);

    // Quantile positions
    double alpha_lower = (1.0 - credible_level) / 2.0;
    double alpha_upper = (1.0 + credible_level) / 2.0;

    for (int v = 0; v < n; ++v) {
        // Extract samples for this vertex
        std::vector<double> vertex_samples(n_samples);
        for (int b = 0; b < n_samples; ++b) {
            vertex_samples[b] = y_samples(v, b);
        }

        // Sort for quantile computation
        std::sort(vertex_samples.begin(), vertex_samples.end());

        // Compute empirical quantiles
        int idx_lower = static_cast<int>(std::floor(alpha_lower * (n_samples - 1)));
        int idx_upper = static_cast<int>(std::floor(alpha_upper * (n_samples - 1)));

        // Ensure indices are valid
        idx_lower = std::max(0, std::min(idx_lower, n_samples - 1));
        idx_upper = std::max(0, std::min(idx_upper, n_samples - 1));

        lower[v] = vertex_samples[idx_lower];
        upper[v] = vertex_samples[idx_upper];

        // Compute posterior standard deviation
        double mean_v = 0.0;
        for (double val : vertex_samples) {
            mean_v += val;
        }
        mean_v /= n_samples;

        double var_v = 0.0;
        for (double val : vertex_samples) {
            var_v += (val - mean_v) * (val - mean_v);
        }
        posterior_sd[v] = std::sqrt(var_v / (n_samples - 1));
    }

    // ================================================================
    // STEP 5: PACKAGE RESULTS
    // ================================================================

    posterior_summary_result_t summary;
    summary.lower = lower;
    summary.upper = upper;
    summary.posterior_sd = posterior_sd;
    summary.credible_level = credible_level;
    summary.sigma_hat = sigma_hat;

    // Conditionally store samples
    if (return_samples) {
        summary.samples = std::move(y_samples);
        summary.has_samples = true;
    } else {
        summary.has_samples = false;
    }

    return summary;
}

/**
 * @brief Convert posterior_summary_result_t to SEXP R list
 */
static SEXP create_posterior_sexp(const posterior_summary_result_t& summary) {
    const int n = summary.lower.size();

    // Determine list size based on whether samples are included
    int n_components = summary.has_samples ? 6 : 5;
    SEXP s_posterior = PROTECT(Rf_allocVector(VECSXP, n_components));
    SEXP s_names = PROTECT(Rf_allocVector(STRSXP, n_components));

    // Component 1: Lower credible bounds
    SEXP s_lower = PROTECT(Rf_allocVector(REALSXP, n));
    double* p_lower = REAL(s_lower);
    for (int i = 0; i < n; ++i) {
        p_lower[i] = summary.lower[i];
    }
    SET_VECTOR_ELT(s_posterior, 0, s_lower);
    SET_STRING_ELT(s_names, 0, Rf_mkChar("lower"));

    // Component 2: Upper credible bounds
    SEXP s_upper = PROTECT(Rf_allocVector(REALSXP, n));
    double* p_upper = REAL(s_upper);
    for (int i = 0; i < n; ++i) {
        p_upper[i] = summary.upper[i];
    }
    SET_VECTOR_ELT(s_posterior, 1, s_upper);
    SET_STRING_ELT(s_names, 1, Rf_mkChar("upper"));

    // Component 3: Posterior standard deviations
    SEXP s_sd = PROTECT(Rf_allocVector(REALSXP, n));
    double* p_sd = REAL(s_sd);
    for (int i = 0; i < n; ++i) {
        p_sd[i] = summary.posterior_sd[i];
    }
    SET_VECTOR_ELT(s_posterior, 2, s_sd);
    SET_STRING_ELT(s_names, 2, Rf_mkChar("sd"));

    // Component 4: Credible level
    SEXP s_level = PROTECT(Rf_allocVector(REALSXP, 1));
    REAL(s_level)[0] = summary.credible_level;
    SET_VECTOR_ELT(s_posterior, 3, s_level);
    SET_STRING_ELT(s_names, 3, Rf_mkChar("credible.level"));

    // Component 5: Estimated sigma
    SEXP s_sigma = PROTECT(Rf_allocVector(REALSXP, 1));
    REAL(s_sigma)[0] = summary.sigma_hat;
    SET_VECTOR_ELT(s_posterior, 4, s_sigma);
    SET_STRING_ELT(s_names, 4, Rf_mkChar("sigma"));

    // Component 6: Posterior samples (optional)
    if (summary.has_samples) {
        const int n_samples = summary.samples.cols();

        SEXP s_samples = PROTECT(Rf_allocMatrix(REALSXP, n, n_samples));
        double* p_samples = REAL(s_samples);

        // Copy column-major (R and Eigen both use column-major)
        for (int j = 0; j < n_samples; ++j) {
            for (int i = 0; i < n; ++i) {
                p_samples[i + j * n] = summary.samples(i, j);
            }
        }

        SET_VECTOR_ELT(s_posterior, 5, s_samples);
        SET_STRING_ELT(s_names, 5, Rf_mkChar("samples"));

        UNPROTECT(1);  // s_samples
    }

    Rf_setAttrib(s_posterior, R_NamesSymbol, s_names);

    // UNPROTECT: s_posterior, s_names, s_lower, s_upper, s_sd, s_level, s_sigma
    UNPROTECT(7);

    return s_posterior;
}


// ================================================================
// SEXP INTERFACE FUNCTION
// ================================================================

/**
 * @brief R-callable interface for posterior summary computation
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
 * @return R list with posterior summary components
 */
extern "C" SEXP S_compute_posterior_summary(
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
) {
    // ================================================================
    // EXTRACT DIMENSIONS
    // ================================================================

    SEXP s_dim = Rf_getAttrib(s_V, R_DimSymbol);
    if (Rf_isNull(s_dim) || Rf_length(s_dim) != 2) {
        Rf_error("V must be a matrix");
    }
    const int n = INTEGER(s_dim)[0];
    const int m = INTEGER(s_dim)[1];

    // ================================================================
    // CONVERT R OBJECTS TO EIGEN
    // ================================================================

    // Eigenvector matrix V
    Eigen::MatrixXd V(n, m);
    const double* p_V = REAL(s_V);
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < n; ++i) {
            V(i, j) = p_V[i + j * n];  // Column-major
        }
    }

    // Eigenvalues
    if (Rf_length(s_eigenvalues) != m) {
        Rf_error("eigenvalues length (%d) must match ncol(V) (%d)",
                 Rf_length(s_eigenvalues), m);
    }
    vec_t eigenvalues(m);
    const double* p_eigenvalues = REAL(s_eigenvalues);
    for (int j = 0; j < m; ++j) {
        eigenvalues[j] = p_eigenvalues[j];
    }

    // Filtered eigenvalues
    if (Rf_length(s_filtered_eigenvalues) != m) {
        Rf_error("filtered_eigenvalues length (%d) must match ncol(V) (%d)",
                 Rf_length(s_filtered_eigenvalues), m);
    }
    vec_t filtered_eigenvalues(m);
    const double* p_filtered = REAL(s_filtered_eigenvalues);
    for (int j = 0; j < m; ++j) {
        filtered_eigenvalues[j] = p_filtered[j];
    }

    // Response vector y
    if (Rf_length(s_y) != n) {
        Rf_error("y length (%d) must match nrow(V) (%d)",
                 Rf_length(s_y), n);
    }
    vec_t y(n);
    const double* p_y = REAL(s_y);
    for (int i = 0; i < n; ++i) {
        y[i] = p_y[i];
    }

    // Fitted values y_hat
    if (Rf_length(s_y_hat) != n) {
        Rf_error("y_hat length (%d) must match nrow(V) (%d)",
                 Rf_length(s_y_hat), n);
    }
    vec_t y_hat(n);
    const double* p_y_hat = REAL(s_y_hat);
    for (int i = 0; i < n; ++i) {
        y_hat[i] = p_y_hat[i];
    }

    // Scalar parameters
    double eta = Rf_asReal(s_eta);
    double credible_level = Rf_asReal(s_credible_level);
    int n_samples = Rf_asInteger(s_n_samples);
    unsigned int seed = static_cast<unsigned int>(Rf_asInteger(s_seed));
    bool return_samples = Rf_asLogical(s_return_samples);

    // ================================================================
    // COMPUTE POSTERIOR SUMMARY
    // ================================================================

    posterior_summary_result_t summary;

    try {
        summary = compute_posterior_summary_impl(
            V, eigenvalues, filtered_eigenvalues,
            y, y_hat, eta,
            credible_level, n_samples, seed, return_samples
        );
    } catch (const std::exception& e) {
        Rf_error("Posterior computation failed: %s", e.what());
    }

    // ================================================================
    // CONVERT TO SEXP AND RETURN
    // ================================================================

    return create_posterior_sexp(summary);
}
