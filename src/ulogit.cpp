#include <R.h>
#include <Rinternals.h>

#include <vector>
#include <array>
#include <cmath>      // for fabs()
#include <algorithm>  // for std::find,

#include <Eigen/Dense>

#include "ulogit.hpp"
#include "error_utils.h"

extern "C" {
    SEXP S_ulogit(SEXP x_sexp,
                  SEXP y_sexp,
                  SEXP w_sexp,
                  SEXP max_iterations_sexp,
                  SEXP ridge_lambda_sexp,
                  SEXP max_beta_sexp,
                  SEXP tolerance_sexp,
                  SEXP verbose_sexp);

    SEXP S_eigen_ulogit(SEXP x_sexp,
                        SEXP y_sexp,
                        SEXP w_sexp,
                        SEXP fit_quadratic_sexp,
                        SEXP with_errors_sexp,
                        SEXP max_iterations_sexp,
                        SEXP ridge_lambda_sexp,
                        SEXP tolerance_sexp);
}

/**
 * @brief Fits a weighted logistic regression model on a window of one-dimensional data
 *
 * @details This function implements a locally weighted logistic regression for binary classification
 * using Newton-Raphson optimization. The model fits a logistic curve to a window of data points,
 * with each point weighted by both the provided weights and its position in the window.
 *
 * The function performs the following steps:
 * 1. Centers the x values for numerical stability
 * 2. Fits a logistic model using Newton-Raphson iteration
 * 3. Computes predictions and Leave-One-Out Cross-Validation (LOOCV) errors
 *
 * The logistic model has the form:
 * P(y=1|x) = 1 / (1 + exp(-(β₀ + β₁(x - x̄))))
 * where x̄ is the weighted mean of x values in the window
 *
 * @param x Pointer to array of x values (predictor variable) in the window
 * @param y Pointer to array of binary y values (0 or 1) corresponding to x
 * @param w Vector of observation weights for the window
 * @param max_iterations Maximum number of iterations for Newton-Raphson optimization (default: 100)
 * @param ridge_lambda  Ridge regularization parameter for stability (default: 0.1)
 * @param max_beta      Maximum allowed absolute value for coefficient estimates (default: 100.0)
 * @param tolerance     Convergence tolerance for optimization (default: 1e-8)
 * @param verbose       Flag to enable detailed output during optimization (default: false)
 *
 * @return ulogit_t struct containing:
 *         - predictions: fitted probabilities
 *         - errors: leave-one-out cross-validation errors
 *         - weights: observation weights used
 *         - x_min_index: index of smallest x value
 *         - x_max_index: index of largest x value
 *         If verbose is true, also includes:
 *         - iteration_count: number of iterations until convergence
 *         - converged: whether optimization converged
 *
 *
 * @note
 * - Input y values must be binary (0 or 1)
 * - The window size must be at least 2 points
 * - The function uses weighted maximum likelihood estimation
 * - LOOCV errors are computed using log loss: -log(p) for y=1, -log(1-p) for y=0
 *
 * @Rf_warning
 * - The function assumes x values are sorted in ascending order
 * - Numerical instability may occur with extremely imbalanced weights
 * - The Newton-Raphson method may not converge for pathological data
 *
 * @see ulogit_t for the return type structure
 * @see wmabilo() for the main algorithm using this function
 *
 * @example
 * ```cpp
 * const double x[] = {1.0, 2.0, 3.0, 4.0, 5.0};
 * const double y[] = {0.0, 0.0, 1.0, 1.0, 1.0};
 * std::vector<double> w = {0.2, 0.5, 1.0, 0.5, 0.2};
 * int x_min_idx = 0;
 * int x_max_idx = 4;
 *
 * auto result = ulogit(x, y, w, x_min_idx, x_max_idx);
 * // Access predictions with result.predictions
 * // Access LOOCV errors with result.errors
 * ```
 */
ulogit_t ulogit(const double* x,
                const double* y,
                const std::vector<double>& w,
                int max_iterations,
                double ridge_lambda,
                double max_beta,
                double tolerance,
                bool verbose) {

    int window_size = w.size();
    ulogit_t result;
    result.w = w;
    result.predictions.resize(window_size);
    result.errors.resize(window_size);

    if (verbose) {
        result.iteration_count = 0;
        result.converged = false;
    }

    // Weight validation
    double total_weight = 0.0;
    for (const auto& weight : w) {
        total_weight += weight;
    }
    if (total_weight <= 0) REPORT_ERROR("total_weight: %.2f   Sum of weights must be positive", total_weight);

    // Check for effective complete separation considering weights
    double weighted_sum_y = 0.0;
    double total_nonzero_weight = 0.0;
    for (int i = 0; i < window_size; ++i) {
        if (w[i] > tolerance) {  // Only consider points with non-negligible weights
            weighted_sum_y += w[i] * y[i];
            total_nonzero_weight += w[i];
        }
    }

    double prior_prob = weighted_sum_y / total_nonzero_weight;

    // If all significant weights are associated with same y value
    if (total_nonzero_weight > 0) {  // Ensure we have some meaningful data
        double weighted_mean_y = weighted_sum_y / total_nonzero_weight;
        if (weighted_mean_y > 1.0 - tolerance || weighted_mean_y < tolerance) {
            // We have effective complete separation
            double pred_value = (weighted_mean_y > 0.5) ? 1.0 : 0.0;
            std::fill(result.predictions.begin(), result.predictions.end(), pred_value);
            std::fill(result.errors.begin(), result.errors.end(), -std::log(1.0 - tolerance));
            return result;
        }
    }

    // Define Newton-Raphson optimization as a lambda function with ridge regularization
    auto newton_raphson = [&](const std::vector<double>& weights, double x_mean, bool& converged) -> std::array<double, 2> {
        double weights_sum = 0.0;
        double y_mean = 0.0;
        for (int i = 0; i < window_size; ++i) {
            y_mean += weights[i] * y[i];
            weights_sum += weights[i];
        }
        y_mean /= weights_sum;

        // Apply logit function to get initial b0
        y_mean = std::clamp(y_mean, tolerance, 1.0 - tolerance); // If we get here, we need to protect the logit transformation in beta initialization
        std::array<double, 2> beta = {std::log(y_mean / (1.0 - y_mean)), 0.0};  // {intercept, slope}

        converged = false;
        int iter = 0;
        while (iter < max_iterations && !converged) {
            double g0 = 0.0, g1 = 0.0;  // Gradients
            double h00 = 0.0, h01 = 0.0, h11 = 0.0;  // Hessian

            for (int i = 0; i < window_size; ++i) {
                double x_centered = x[i] - x_mean;
                double xbeta = beta[0] + beta[1] * x_centered;

                // Prevent overflow in exp
                xbeta = std::clamp(xbeta, -max_beta, max_beta);
                double p = 1.0 / (1.0 + std::exp(-xbeta));

                double w_i = weights[i];

                // Gradient components
                double Rf_error = p - y[i];
                g0 += w_i * Rf_error;
                g1 += w_i * Rf_error * x_centered;

                // Hessian components
                double p_1mp = p * (1.0 - p);
                h00 += w_i * p_1mp;
                h01 += w_i * p_1mp * x_centered;
                h11 += w_i * p_1mp * x_centered * x_centered;
            }

            // Add ridge regularization terms
            g0 += ridge_lambda * beta[0];
            g1 += ridge_lambda * beta[1];
            h00 += ridge_lambda;
            h11 += ridge_lambda;

            // Solve 2x2 system using direct inverse
            double det = h00 * h11 - h01 * h01;
            double step = 1.0; //std::min(1.0, 1.0 / (1.0 + std::exp(-std::abs(det))));

            if (std::abs(det) < tolerance) {
                beta[0] -= g0 * step / (h00 + ridge_lambda);
                beta[1] -= g1 * step / (h11 + ridge_lambda);
            } else {
                double d_beta0 = (h11 * g0 - h01 * g1) / det;
                double d_beta1 = (-h01 * g0 + h00 * g1) / det;
                beta[0] -= step * d_beta0;
                beta[1] -= step * d_beta1;
            }

            // Clamp coefficients to prevent explosion
            beta[0] = std::clamp(beta[0], -max_beta, max_beta);
            beta[1] = std::clamp(beta[1], -max_beta, max_beta);

            // Check convergence using provided tolerance
            if (std::abs(g0) < tolerance && std::abs(g1) < tolerance) {
                converged = true;
                if (verbose) {
                    result.iteration_count = iter + 1;
                    result.converged = true;
                    Rprintf("Newton-Raphson converged at iter: %d\n", iter + 1);
                }
            }

            iter++;
        }
        return beta;
    };

    // Fit full model
    double x_mean = 0.0;
    for (int i = 0; i < window_size; ++i) {
        x_mean += w[i] * x[i];
    }
    x_mean /= total_weight;

    bool converged;
    std::array<double, 2> beta = newton_raphson(w, x_mean, converged);

    // Compute predictions for full model
    for (int i = 0; i < window_size; ++i) {
        double x_centered = x[i] - x_mean;
        double xbeta = beta[0] + beta[1] * x_centered;
        xbeta = std::clamp(xbeta, -100.0, 100.0);  // Prevent overflow
        result.predictions[i] = 1.0 / (1.0 + std::exp(-xbeta));
    }

    // Compute LOOCV errors, handling potential complete separation in subsets
    for (int i = 0; i < window_size; ++i) {
        // Create leave-one-out weights
        std::vector<double> loo_weights = w;
        loo_weights[i] = 0.0;

        double weighted_sum_y = 0.0;
        double total_weight = 0.0;
        for (int j = 0; j < window_size; ++j) {
            if (j != i && loo_weights[j] > tolerance) {
                weighted_sum_y += loo_weights[j] * y[j];
                total_weight += loo_weights[j];
            }
        }

        double loo_pred;
        if (total_weight > 0) {
            double weighted_mean_y = weighted_sum_y / total_weight;
            if (weighted_mean_y > 1.0 - tolerance || weighted_mean_y < tolerance) {
                // In case of complete separation in LOO sample
                loo_pred = (weighted_mean_y > 1.0 - tolerance) ? 1.0 - tolerance : tolerance;
            } else {
                // Compute new weighted mean without point i
                double loo_total_weight = 0.0;
                double loo_x_mean = 0.0;
                for (int j = 0; j < window_size; ++j) {
                    if (j != i) {
                        loo_x_mean += loo_weights[j] * x[j];
                        loo_total_weight += loo_weights[j];
                    }
                }
                loo_x_mean /= loo_total_weight;

                // Fit model without point i
                bool loo_converged;
                std::array<double, 2> beta_loo = newton_raphson(loo_weights, loo_x_mean, loo_converged);

                // Compute prediction for held-out point
                double x_centered = x[i] - loo_x_mean;
                double xbeta = beta_loo[0] + beta_loo[1] * x_centered;
                xbeta = std::clamp(xbeta, -100.0, 100.0);  // Prevent overflow
                loo_pred = 1.0 / (1.0 + std::exp(-xbeta));
            }
        } else {
            loo_pred = prior_prob;
        }

        // Compute cross-entropy Rf_error
        //result.errors[i] = -y[i] * std::log(std::max(loo_pred, tolerance)) - (1-y[i]) * std::log(std::max(1.0 - loo_pred, tolerance));
        //result.errors[i] = std::abs(y[i] - loo_pred); // absolute deviation Rf_error
        result.errors[i] = (y[i] - loo_pred) * (y[i] - loo_pred);
    }

    return result;
}

inline bool r_logical_to_bool_true_only(SEXP s) {
    if (TYPEOF(s) != LGLSXP || XLENGTH(s) < 1)
        Rf_error("Expected logical(1)");
    return (LOGICAL(s)[0] == 1) == 1;  // TRUE->true; FALSE/NA->false
}

/**
 * @brief R interface wrapper for univariate logistic regression function
 *
 * @details SEXP wrapper function that provides an R interface to the C++ ulogit function.
 * Handles conversion between R and C++ data types, memory management, and returns results in an R-compatible format.
 * The function expects sorted x values and binary y values (0 or 1).
 *
 * @param x_sexp SEXP (numeric vector) Input predictor values
 *               Must be sorted in ascending order
 * @param y_sexp SEXP (numeric vector) Binary response values (0 or 1)
 *               Must be same length as x_sexp
 * @param w_sexp SEXP (numeric vector) Window weights
 *               Length must equal x_max_index - x_min_index + 1
 * @param max_iterations_sexp SEXP containing integer for maximum iterations
 * @param ridge_lambda_sexp  SEXP containing numeric ridge regularization parameter
 * @param max_beta_sexp     SEXP containing numeric maximum coefficient value
 * @param tolerance_sexp     SEXP containing numeric convergence tolerance
 * @param verbose_sexp      SEXP containing logical verbose flag
 *
 * @return SEXP list containing:
 *         - predictions: numeric vector of fitted probabilities
 *         - errors: numeric vector of leave-one-out cross-validation errors
 *         - weights: numeric vector of weights used in fitting
 *
 * @note
 * - All indices are converted between R's 1-based and C++'s 0-based indexing
 * - The function uses PROTECT/UNPROTECT for R's garbage collection
 * - Return value is a named list for easy access in R
 *
 * @Rf_warning
 * - No input validation is performed in the C++ code
 * - R code should validate inputs before calling this function
 * - Memory management relies on R's garbage collection system
 *
 * @see ulogit() for the underlying C++ implementation
 * @see ulogit_t for the C++ return type structure
 *
 * Usage in R:
 * ```r
 * result <- .Call("S_ulogit",
 *                 as.double(x),
 *                 as.double(y),
 *                 as.double(w),
 *                 as.double(tolerance))
 * ```
 */
SEXP S_ulogit(SEXP x_sexp,
              SEXP y_sexp,
              SEXP w_sexp,
              SEXP max_iterations_sexp,
              SEXP ridge_lambda_sexp,
              SEXP max_beta_sexp,
              SEXP tolerance_sexp,
              SEXP verbose_sexp) {
    // Basic type/length checks (optional but nice)
    if (TYPEOF(x_sexp) != REALSXP || TYPEOF(y_sexp) != REALSXP || TYPEOF(w_sexp) != REALSXP)
        Rf_error("x, y, w must be numeric");
    R_xlen_t n = XLENGTH(w_sexp);
    if (XLENGTH(x_sexp) < n || XLENGTH(y_sexp) < n)
        Rf_error("x and y must have length at least length(w)");

    const double* x = REAL(x_sexp);
    const double* y = REAL(y_sexp);
    const double* w_r = REAL(w_sexp);

    int max_iterations = INTEGER(max_iterations_sexp)[0];
    double ridge_lambda = REAL(ridge_lambda_sexp)[0];
    double max_beta = REAL(max_beta_sexp)[0];
    double tolerance = REAL(tolerance_sexp)[0];

    // SAFE: treat only TRUE as true; FALSE/NA => false
    const int* vptr = LOGICAL(verbose_sexp);
    bool verbose = (vptr && vptr[0] == 1);

    // Copy weights
    std::vector<double> w;
    w.assign(w_r, w_r + n);

    // Compute
    ulogit_t result = ulogit(x, y, w,
                             max_iterations,
                             ridge_lambda,
                             max_beta,
                             tolerance,
                             verbose);

    // Build return list
    const int N_COMPONENTS = 3;
    SEXP out = PROTECT(Rf_allocVector(VECSXP, N_COMPONENTS));
    {
        // names
        SEXP names = PROTECT(Rf_allocVector(STRSXP, N_COMPONENTS));
        SET_STRING_ELT(names, 0, Rf_mkChar("predictions"));
        SET_STRING_ELT(names, 1, Rf_mkChar("errors"));
        SET_STRING_ELT(names, 2, Rf_mkChar("weights"));
        Rf_setAttrib(out, R_NamesSymbol, names);
        UNPROTECT(1);
    }

    // predictions
    {
        SEXP predictions = PROTECT(Rf_allocVector(REALSXP, n));
        for (R_xlen_t i = 0; i < n; ++i)
            REAL(predictions)[i] = result.predictions[i];
        SET_VECTOR_ELT(out, 0, predictions);
        UNPROTECT(1);
    }

    // errors
    {
        SEXP errors = PROTECT(Rf_allocVector(REALSXP, n));
        for (R_xlen_t i = 0; i < n; ++i)
            REAL(errors)[i] = result.errors[i];
        SET_VECTOR_ELT(out, 1, errors);
        UNPROTECT(1);
    }


    // weights
    {
        SEXP weights = PROTECT(Rf_allocVector(REALSXP, n));
        for (R_xlen_t i = 0; i < n; ++i)
            REAL(weights)[i] = result.w[i];
        SET_VECTOR_ELT(out, 2, weights);
        UNPROTECT(1);
    }

    UNPROTECT(1);
    return out;
}

/**
 * @brief Fits a univariate logistic regression model and returns predictions
 *
 * @details This function fits a logistic regression model using weighted maximum likelihood
 * estimation with ridge regularization. It handles complete separation cases and includes
 * safeguards against numerical instability. The optimization is performed using
 * Newton-Raphson iteration.
 *
 * The model fitted is: logit(p) = β₀ + β₁(x - x̄), where x̄ is the weighted mean of x.
 *
 * @param x Pointer to predictor values array
 * @param y Pointer to binary response values array (should contain only 0s and 1s)
 * @param w Vector of observation weights (must be non-negative)
 * @param max_iterations Maximum number of Newton-Raphson iterations (default: 100)
 * @param ridge_lambda Ridge regularization parameter (default: 0.002)
 * @param max_beta Maximum absolute value for coefficient estimates (default: 100.0)
 * @param tolerance Convergence tolerance for Newton-Raphson iteration (default: 1e-8)
 *
 * @return Vector of predicted probabilities, one for each input observation
 *
 * @throws std::invalid_argument If the sum of weights is not positive
 *
 * @note The function handles complete separation by returning appropriate constant predictions
 * @note Ridge regularization is applied to both intercept and slope coefficients
 */
std::vector<double> ulogit_predict(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    int max_iterations,
    double ridge_lambda,
    double max_beta,
    double tolerance) {

    int window_size = w.size();
    std::vector<double> predictions(window_size);

    // Weight validation
    double total_weight = 0.0;
    for (const auto& weight : w) {
        total_weight += weight;
    }
    if (total_weight <= 0) {
        Rf_error("Sum of weights must be positive");
    }

    // Check for effective complete separation considering weights
    double weighted_sum_y = 0.0;
    double total_nonzero_weight = 0.0;
    for (int i = 0; i < window_size; ++i) {
        if (w[i] > tolerance) {  // Only consider points with non-negligible weights
            weighted_sum_y += w[i] * y[i];
            total_nonzero_weight += w[i];
        }
    }

    // Handle complete separation case
    if (total_nonzero_weight > 0) {
        double weighted_mean_y = weighted_sum_y / total_nonzero_weight;
        if (weighted_mean_y > 1.0 - tolerance || weighted_mean_y < tolerance) {
            double pred_value = (weighted_mean_y > 0.5) ? 1.0 - tolerance : tolerance;
            std::fill(predictions.begin(), predictions.end(), pred_value);
            return predictions;
        }
    }

    // Define Newton-Raphson optimization
    auto newton_raphson = [&](const std::vector<double>& weights, double x_mean, bool& converged) -> std::array<double, 2> {
        double weights_sum = 0.0;
        double y_mean = 0.0;
        for (int i = 0; i < window_size; ++i) {
            y_mean += weights[i] * y[i];
            weights_sum += weights[i];
        }
        y_mean /= weights_sum;

        // Initialize coefficients
        y_mean = std::clamp(y_mean, tolerance, 1.0 - tolerance);
        std::array<double, 2> beta = {std::log(y_mean / (1.0 - y_mean)), 0.0};

        converged = false;
        int iter = 0;
        while (iter < max_iterations && !converged) {
            double g0 = 0.0, g1 = 0.0;  // Gradients
            double h00 = 0.0, h01 = 0.0, h11 = 0.0;  // Hessian

            for (int i = 0; i < window_size; ++i) {
                double x_centered = x[i] - x_mean;
                double xbeta = beta[0] + beta[1] * x_centered;
                xbeta = std::clamp(xbeta, -max_beta, max_beta);
                double p = 1.0 / (1.0 + std::exp(-xbeta));

                double w_i = weights[i];
                double Rf_error = p - y[i];
                double p_1mp = p * (1.0 - p);

                // Update gradient and Hessian
                g0 += w_i * Rf_error;
                g1 += w_i * Rf_error * x_centered;
                h00 += w_i * p_1mp;
                h01 += w_i * p_1mp * x_centered;
                h11 += w_i * p_1mp * x_centered * x_centered;
            }

            // Add ridge regularization
            g0 += ridge_lambda * beta[0];
            g1 += ridge_lambda * beta[1];
            h00 += ridge_lambda;
            h11 += ridge_lambda;

            // Solve system
            double det = h00 * h11 - h01 * h01;
            if (std::abs(det) < tolerance) {
                beta[0] -= g0 / (h00 + ridge_lambda);
                beta[1] -= g1 / (h11 + ridge_lambda);
            } else {
                beta[0] -= (h11 * g0 - h01 * g1) / det;
                beta[1] -= (-h01 * g0 + h00 * g1) / det;
            }

            beta[0] = std::clamp(beta[0], -max_beta, max_beta);
            beta[1] = std::clamp(beta[1], -max_beta, max_beta);

            converged = std::abs(g0) < tolerance && std::abs(g1) < tolerance;
            iter++;
        }
        return beta;
    };

    // Compute weighted mean of x
    double x_mean = 0.0;
    for (int i = 0; i < window_size; ++i) {
        x_mean += w[i] * x[i];
    }
    x_mean /= total_weight;

    // Fit model
    bool converged;
    std::array<double, 2> beta = newton_raphson(w, x_mean, converged);

    // Compute predictions
    for (int i = 0; i < window_size; ++i) {
        double x_centered = x[i] - x_mean;
        double xbeta = beta[0] + beta[1] * x_centered;
        xbeta = std::clamp(xbeta, -max_beta, max_beta);
        predictions[i] = 1.0 / (1.0 + std::exp(-xbeta));
    }

    return predictions;
}

/**
 * @brief Fits univariate logistic regression using Eigen with LOOCV Rf_error estimation
 *
 * @details Implements weighted logistic regression using iteratively reweighted least
 * squares (IRLS) with the Newton-Raphson algorithm. The implementation includes several
 * numerical stability safeguards:
 *
 * Algorithm Details:
 * 1. Initialization:
 *    - Initial beta values set using log-odds of weighted mean response
 *    - Working weights initialized using initial probabilities
 *
 * 2. IRLS Implementation:
 *    - Uses working response $$\eta + (y - \mu)/var$$ where:
 *      - η is the linear predictor
 *      - μ is the fitted probability
 *      - var is μ(1-μ)
 *    - Working weights are w * var where w are the input weights
 *    - Includes step halving for better convergence
 *    - Convergence checked using deviance change
 *
 * 3. Numerical Stability Features:
 *    - Bounds probabilities away from 0/1 using tolerance parameter
 *    - Adds ridge penalty to diagonal of XᵀWX matrix
 *    - Uses stable matrix decompositions (LDLT, QR, or SVD as fallbacks)
 *    - Includes step halving in Newton-Raphson updates
 *
 * 4. Leave-One-Out Cross-Validation (LOOCV):
 *    The function computes three types of LOOCV prediction errors:
 *
 *    a) Deviance errors:
 *       $$-[y_i\log(\hat{p}_{(-i)}) + (1-y_i)\log(1-\hat{p}_{(-i)})]$$
 *
 *    b) Brier score errors:
 *       $$(\hat{p}_{(-i)} - y_i)^2$$
 *
 *    c) Absolute errors:
 *       $$|y_i - \hat{p}_{(-i)}|$$
 *
 *    where $$\hat{p}_{(-i)}$$ is approximated using:
 *    $$\hat{p}_{(-i)} \approx \hat{p}_i - \frac{h_i(y_i - \hat{p}_i)}{1 - h_i}$$
 *
 *    For high leverage points (h_i > 3p/n), a more stable Williams' form is used:
 *    $$\hat{p}_{(-i)} = \hat{p}_i - \frac{h_i r_i}{1-h_i}\sqrt{\frac{1-h_i^*}{1-h_i}}$$
 *    where r_i is the Pearson residual and h_i* = h_i(1-h_i)
 *
 * @param x Pointer to predictor values
 * @param y Pointer to binary response values (0/1)
 * @param w Vector of observation weights (must sum to positive value)
 * @param fit_quadratic If true, fits quadratic term in addition to linear
 * @param max_iterations Maximum number of Newton-Raphson iterations (default: 25)
 * @param ridge_lambda Ridge regularization parameter for stability (default: 1e-6)
 * @param tolerance Numerical tolerance for:
 *                 - Probability bounds (clamped to [tolerance, 1-tolerance])
 *                 - Convergence checking in deviance
 *                 - Minimum value for working weights
 * @param with_errors If true, compute LOOCV prediction errors
 *
 * @return eigen_ulogit_t structure containing:
 *         - predictions: Fitted probabilities
 *         - loocv_brier_errors: Brier score LOOCV errors
 *         - beta: Model coefficients (intercept, linear [, quadratic])
 *         - converged: Whether algorithm converged
 *         - iterations: Number of iterations used
 *
 * @throws std::invalid_argument If:
 *         - Weights sum to zero or negative value
 *         - x or y pointers are null
 *         - Weight vector size doesn't match data
 *
 * @note
 * 1. Handles complete separation through ridge regularization
 * 2. Uses multiple matrix decomposition methods for stability:
 *    - LDLT decomposition (primary)
 *    - QR with column pivoting (first fallback)
 *    - SVD (second fallback)
 * 3. If all decompositions fail for Rf_error estimation, uses fitted
 *    probabilities as approximation for LOOCV estimates
 * 4. Particularly suitable for local logistic regression fitting where
 *    numerical stability is crucial
 */
eigen_ulogit_t eigen_ulogit_fit(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    bool fit_quadratic,
    int max_iterations,
    double ridge_lambda,
    double tolerance,
    bool with_errors) {

    int n = w.size();
    int p = fit_quadratic ? 3 : 2;

    // Initialize return structure
    eigen_ulogit_t result;
    result.beta = Eigen::VectorXd::Zero(p);
    result.converged = false;
    result.iterations = 0;
    result.predictions.resize(n);
    result.fit_quadratic = fit_quadratic;

    // Initialize beta using weighted mean of y
    double sum_w = 0.0, sum_wy = 0.0;
    for (int i = 0; i < n; ++i) {
        sum_w += w[i];
        sum_wy += w[i] * y[i];
    }
    double y_bar = sum_wy / sum_w;

    // Set initial values using log-odds
    if (y_bar > 0 && y_bar < 1) {
        result.beta(0) = std::log(y_bar / (1.0 - y_bar));
    }

    double prev_deviance = std::numeric_limits<double>::max();

    auto compute_probability_bound = [](int n_pts, double C = 0.1) -> double {
        return std::max(1e-15, C / static_cast<double>(n_pts));
    };

    double prob_bound = compute_probability_bound(n);
    double log_prob_bound = std::log(prob_bound);

    for (int iter = 0; iter < max_iterations; ++iter) {
        // First pass: compute probabilities and working variables
        std::vector<double> mu(n);
        std::vector<double> working_y(n);
        std::vector<double> working_w(n);
        double current_deviance = 0.0;

        // Compute current probabilities and deviance
        for (int i = 0; i < n; ++i) {
            // Compute linear predictor
            double eta = result.beta(0) + result.beta(1) * x[i];
            if (fit_quadratic) {
                eta += result.beta(2) * x[i] * x[i];
            }

            // Compute probability with bounds
            if (eta > -log_prob_bound) {
                mu[i] = 1.0 - prob_bound;
            } else if (eta < log_prob_bound) {
                mu[i] = prob_bound;
            } else {
                mu[i] = 1.0 / (1.0 + std::exp(-eta));
            }

            // Compute working response
            double var = mu[i] * (1.0 - mu[i]);
            if (var < 1e-10) var = 1e-10;
            working_y[i] = eta + (y[i] - mu[i]) / var;
            working_w[i] = w[i] * var;

            // Compute deviance
            if (y[i] > 0) {
                current_deviance -= 2 * w[i] * y[i] * std::log(mu[i]);
            }
            if (y[i] < 1) {
                current_deviance -= 2 * w[i] * (1 - y[i]) * std::log(1 - mu[i]);
            }
        }

        // Set up weighted least squares matrices
        Eigen::MatrixXd X(n, p);
        Eigen::VectorXd z(n);
        Eigen::VectorXd weights(n);

        for (int i = 0; i < n; ++i) {
            X(i, 0) = 1.0;
            X(i, 1) = x[i];
            if (fit_quadratic) {
                X(i, 2) = x[i] * x[i];
            }
            z(i) = working_y[i];
            weights(i) = std::sqrt(working_w[i]);
        }

        // Weight the matrices
        Eigen::MatrixXd X_w = weights.asDiagonal() * X;
        Eigen::VectorXd z_w = weights.cwiseProduct(z);

        // Add tiny ridge penalty for stability
        Eigen::MatrixXd XtX = X_w.transpose() * X_w;

        // Solve weighted least squares with proper Rf_error checking
        Eigen::LDLT<Eigen::MatrixXd> ldlt(XtX);
        Eigen::VectorXd new_beta;
        bool solve_succeeded = false;

        if (ldlt.info() == Eigen::Success) {
            new_beta = ldlt.solve(X_w.transpose() * z_w);
            solve_succeeded = true;
        } else {
            // Try QR with pivoting
            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(XtX);
            qr.setThreshold(1e-10);
            if (qr.rank() == XtX.cols()) {
                new_beta = qr.solve(X_w.transpose() * z_w);
                solve_succeeded = true;
            } else {
                // Try SVD as last resort
                Eigen::JacobiSVD<Eigen::MatrixXd> svd(XtX, Eigen::ComputeThinU | Eigen::ComputeThinV);
                double threshold = 1e-10 * svd.singularValues()(0);
                if (svd.singularValues()(svd.singularValues().size()-1) > threshold) {
                    new_beta = svd.solve(X_w.transpose() * z_w);
                    solve_succeeded = true;
                }
            }
        }

        if (!solve_succeeded) {
            // Use previous beta values
            new_beta = result.beta;
            // Add Rf_warning
            result.warnings.push_back("Numerical issues in IWLS step - using previous estimates");

            // Add extra ridge penalty and try again
            XtX.diagonal() += Eigen::VectorXd::Constant(p, 1e-4);  // Increased ridge penalty
            Eigen::LDLT<Eigen::MatrixXd> ldlt_retry(XtX);
            if (ldlt_retry.info() == Eigen::Success) {
                new_beta = ldlt_retry.solve(X_w.transpose() * z_w);
                result.warnings.push_back("Succeeded with additional ridge penalty");
            }
        }

        // Check for extreme fitted probabilities
        bool has_extreme_probs = false;
        for (int i = 0; i < n; ++i) {
            if (mu[i] < 1e-10 || mu[i] > 1 - 1e-10) {
                has_extreme_probs = true;
                break;
            }
        }
        if (has_extreme_probs) {
            result.warnings.push_back("fitted probabilities numerically 0 or 1 occurred");
        }

        // Step halving if needed
        double step = 1.0;
        Eigen::VectorXd beta_old = result.beta;

        while (step > 1e-10) {
            Eigen::VectorXd trial_beta = beta_old + step * (new_beta - beta_old);

            // Compute trial deviance
            double trial_deviance = 0.0;
            for (int i = 0; i < n; ++i) {
                double eta = trial_beta(0) + trial_beta(1) * x[i];
                if (fit_quadratic) {
                    eta += trial_beta(2) * x[i] * x[i];
                }

                double p = 1.0 / (1.0 + std::exp(-eta));
                if (y[i] > 0) trial_deviance -= 2 * w[i] * y[i] * std::log(p);
                if (y[i] < 1) trial_deviance -= 2 * w[i] * (1 - y[i]) * std::log(1 - p);
            }

            if (trial_deviance <= current_deviance * (1 + 1e-4)) {
                result.beta = trial_beta;
                current_deviance = trial_deviance;
                break;
            }

            step *= 0.5;
        }

        double deviance_change = std::abs(current_deviance - prev_deviance);

        // Check convergence using deviance change
        if (iter > 0 && deviance_change < 1e-8) {
            result.converged = true;
            break;
        }

        prev_deviance = current_deviance;
        result.iterations = iter + 1;
    }

    // Compute final predictions
    for (int i = 0; i < n; ++i) {
        double eta = result.beta(0) + result.beta(1) * x[i];
        if (fit_quadratic) {
            eta += result.beta(2) * x[i] * x[i];
        }
        result.predictions[i] = 1.0 / (1.0 + std::exp(-eta));
    }

    if (with_errors) {
        // Compute design matrix
        Eigen::MatrixXd X(n, p);
        Eigen::MatrixXd W_sqrt = Eigen::MatrixXd::Zero(n, n);

        for (int i = 0; i < n; ++i) {
            X(i, 0) = 1.0;
            X(i, 1) = x[i];
            if (fit_quadratic) {
                X(i, 2) = x[i] * x[i];
            }

            // Square root of weight matrix diagonal elements
            W_sqrt(i, i) = std::sqrt(w[i] * result.predictions[i] * (1.0 - result.predictions[i]));
        }

        // Compute hat matrix with added stability
        Eigen::MatrixXd WX = W_sqrt * X;
        Eigen::MatrixXd XtWX = X.transpose() * W_sqrt.transpose() * W_sqrt * X;

        // Add ridge penalty for numerical stability
        XtWX.diagonal() += Eigen::VectorXd::Constant(p, ridge_lambda);

        bool decomposition_succeeded = false;
        Eigen::MatrixXd H;

        // Try LDLT first
        Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
        if (ldlt.info() == Eigen::Success) {
            H = WX * ldlt.solve(X.transpose() * W_sqrt.transpose());
            decomposition_succeeded = true;
        } else {
            // Try QR with pivoting
            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(XtWX);
            qr.setThreshold(1e-10);
            if (qr.rank() == XtWX.cols()) {
                H = WX * qr.solve(X.transpose() * W_sqrt.transpose());
                decomposition_succeeded = true;
            } else {
                // Try SVD as last resort
                Eigen::JacobiSVD<Eigen::MatrixXd> svd(XtWX, Eigen::ComputeThinU | Eigen::ComputeThinV);
                double threshold = 1e-10 * svd.singularValues()(0);
                if (svd.singularValues()(svd.singularValues().size()-1) > threshold) {
                    H = WX * svd.matrixV() *
                        (svd.singularValues().array().abs() > threshold).select(
                            svd.singularValues().array().inverse(), 0).matrix().asDiagonal() *
                        svd.matrixU().transpose() * X.transpose() * W_sqrt.transpose();
                    decomposition_succeeded = true;
                }
            }
        }

        result.loocv_brier_errors.resize(n);
        result.loocv_deviance_errors.resize(n);

        for (int i = 0; i < n; ++i) {
            // Compute linear predictor
            double eta = result.beta(0) + result.beta(1) * x[i];
            if (fit_quadratic) {
                eta += result.beta(2) * x[i] * x[i];
            }

            double y_i = y[i];
            double p_i = result.predictions[i];
            double p_leave_one_out;

            if (eta > -log_prob_bound) {  // p_i ≈ 1
                if (y_i == 1) {
                    p_leave_one_out = 1 - prob_bound;  // minimal correction
                } else {
                    p_leave_one_out = std::max(0.5, 1 - 2*prob_bound);  // larger correction
                }
            } else if (eta < log_prob_bound) {  // p_i ≈ 0
                if (y_i == 0) {
                    p_leave_one_out = prob_bound;  // minimal correction
                } else {
                    p_leave_one_out = std::min(0.5, 2*prob_bound);  // larger correction
                }
            } else { // non-extreme probabilities
                if (decomposition_succeeded) {
                    double h_i = H(i, i);

                    // Calculate p_{(-i)}
                    double denom = std::max(1e-10, 1.0 - h_i);
                    double resid = y_i - p_i;

                    // Basic term
                    double base_correction = (h_i * resid / denom);

                    // Second-order term
                    double second_order = (h_i * h_i * resid * resid * (1 - 2*p_i)) /
                        (2 * denom * denom);

                    // Bias correction term (simplified version)
                    double w_i = p_i * (1 - p_i);
                    double B_i = (h_i * h_i * h_i * resid * resid * resid) /
                        (6 * w_i * denom * denom * denom);

                    p_leave_one_out = p_i - base_correction + second_order + B_i;
                } else {
                    p_leave_one_out = p_i;  // Use fitted value as approximation
                }
            }

            // Ensure valid probabilities
            p_leave_one_out = std::clamp(p_leave_one_out, tolerance, 1.0 - tolerance);

            // Compute Brier score
            result.loocv_brier_errors[i] = std::pow(p_leave_one_out - y_i, 2);

            // deviance Rf_error
            result.loocv_deviance_errors[i] = -(y[i] * std::log(p_leave_one_out) + (1 - y[i]) * std::log(1 - p_leave_one_out));
        }
    }

    return result;
}



/**
 * @brief Fits linear logistic regression and returns predictions
 *
 * @details Wrapper around eigen_ulogit_fit() for linear model case.
 * Returns only predictions, discarding other fit information.
 *
 * @param x Pointer to predictor values
 * @param y Pointer to binary response values (0/1)
 * @param w Vector of observation weights
 * @param max_iterations Maximum number of Newton-Raphson iterations
 * @param ridge_lambda Ridge regularization parameter
 * @param tolerance Convergence tolerance for Newton-Raphson algorithm
 *
 * @return Vector of fitted probabilities
 *
 * @throws std::invalid_argument if weights sum to zero or negative
 *
 * @see eigen_ulogit_fit()
 */
std::vector<double> eigen_ulogit_predict(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    int max_iterations = 100,
    double ridge_lambda = 0.002,
    double tolerance = 1e-8) {

    auto fit = eigen_ulogit_fit(x, y, w,
                                false,
                                max_iterations,
                                ridge_lambda,
                                tolerance,
                                false);
    return fit.predictions;
}

/**
 * @brief Fits quadratic logistic regression and returns predictions
 *
 * @details Wrapper around eigen_ulogit_fit() for quadratic model case.
 * Returns only predictions, discarding other fit information.
 * Includes both linear and quadratic terms in the model.
 *
 * @param x Pointer to predictor values
 * @param y Pointer to binary response values (0/1)
 * @param w Vector of observation weights
 * @param max_iterations Maximum number of Newton-Raphson iterations
 * @param ridge_lambda Ridge regularization parameter
 * @param tolerance Convergence tolerance for Newton-Raphson algorithm
 *
 * @return Vector of fitted probabilities
 *
 * @throws std::invalid_argument if weights sum to zero or negative
 *
 * @see eigen_ulogit_fit()
 */
std::vector<double> eigen_ulogit_quadratic_predict(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    int max_iterations = 100,
    double ridge_lambda = 0.002,
    double tolerance = 1e-8) {

    auto fit = eigen_ulogit_fit(x, y, w,
                                true,
                                max_iterations,
                                ridge_lambda,
                                tolerance,
                                false);
    return fit.predictions;
}


/**
 * @brief R interface function for Eigen-based univariate logistic regression with leave-one-out cross-validation
 *
 * @details Converts R objects to C++ types, calls eigen_ulogit_fit(), and returns results as an R list.
 * Supports both linear and quadratic models, with optional leave-one-out cross-validation (LOOCV)
 * Rf_error estimation using influence-based approximations. The LOOCV errors are computed using
 * the Brier score (squared Rf_error) metric and leverage-based approximations for computational
 * efficiency.
 *
 * @param x_sexp Numeric vector of predictor values
 * @param y_sexp Numeric vector of binary response values (0/1)
 * @param w_sexp Numeric vector of observation weights
 * @param fit_quadratic_sexp Logical scalar indicating whether to fit quadratic term
 * @param with_errors_sexp Logical scalar indicating whether to compute LOOCV errors
 * @param max_iterations_sexp Integer scalar of maximum Newton-Raphson iterations
 * @param ridge_lambda_sexp Numeric scalar of ridge regularization parameter
 * @param tolerance_sexp Numeric scalar of convergence tolerance
 *
 * @return An R list with components:
 *   - predictions: Numeric vector of fitted probabilities
 *   - converged: Logical indicating convergence status
 *   - iterations: Integer giving number of iterations used
 *   - beta: Numeric vector of fitted coefficients (length 2 for linear, 3 for quadratic)
 *   - loocv_brier_errors: Numeric vector of leave-one-out Brier scores (when with_errors=TRUE)
 *
 * @note Protects R objects from garbage collection during computation. The LOOCV errors
 * are computed using efficient influence-based approximations that avoid refitting the
 * model for each observation.
 */
SEXP S_eigen_ulogit(SEXP x_sexp,
                    SEXP y_sexp,
                    SEXP w_sexp,
                    SEXP fit_quadratic_sexp,
                    SEXP with_errors_sexp,
                    SEXP max_iterations_sexp,
                    SEXP ridge_lambda_sexp,
                    SEXP tolerance_sexp) {
    // Convert inputs from R to C++
    double* x = REAL(x_sexp);
    double* y = REAL(y_sexp);
    double* w_r = REAL(w_sexp);
    bool fit_quadratic = (LOGICAL(fit_quadratic_sexp)[0] == 1);
    bool with_errors = (LOGICAL(with_errors_sexp)[0] == 1);
    int max_iterations = INTEGER(max_iterations_sexp)[0];
    double ridge_lambda = REAL(ridge_lambda_sexp)[0];
    double tolerance = REAL(tolerance_sexp)[0];

    // Convert R vector to std::vector
    int window_size = Rf_length(w_sexp);
    std::vector<double> w(w_r, w_r + window_size);

    // Call the actual function with all parameters
    eigen_ulogit_t result = eigen_ulogit_fit(x, y, w,
                                           fit_quadratic,
                                           max_iterations,
                                           ridge_lambda,
                                           tolerance,
                                           with_errors);

    // Creating return list - adjust size based on whether errors are included
    const int N_COMPONENTS = with_errors ? 6 : 4;
    SEXP r_result = PROTECT(Rf_allocVector(VECSXP, N_COMPONENTS));
    {
        // Set list names
        SEXP r_names = PROTECT(Rf_allocVector(STRSXP, N_COMPONENTS));
        SET_STRING_ELT(r_names, 0, Rf_mkChar("predictions"));
        SET_STRING_ELT(r_names, 1, Rf_mkChar("converged"));
        SET_STRING_ELT(r_names, 2, Rf_mkChar("iterations"));
        SET_STRING_ELT(r_names, 3, Rf_mkChar("beta"));
        if (with_errors) {
            SET_STRING_ELT(r_names, 4, Rf_mkChar("loocv_brier_errors"));
            SET_STRING_ELT(r_names, 5, Rf_mkChar("loocv_deviance_errors"));
        }
        Rf_setAttrib(r_result, R_NamesSymbol, r_names);
        UNPROTECT(1); // r_names
    }

    // Convert predictions to R vector
    {
        SEXP predictions = PROTECT(Rf_allocVector(REALSXP, window_size));
        for(int i = 0; i < window_size; i++) {
            REAL(predictions)[i] = result.predictions[i];
        }
        SET_VECTOR_ELT(r_result, 0, predictions);
        UNPROTECT(1);
    }

    // Convert convergence status
    {
        SEXP converged = PROTECT(Rf_allocVector(LGLSXP, 1));
        LOGICAL(converged)[0] = result.converged;
        SET_VECTOR_ELT(r_result, 1, converged);
        UNPROTECT(1);
    }

    // Convert iteration count
    {
        SEXP iterations = PROTECT(Rf_allocVector(INTSXP, 1));
        INTEGER(iterations)[0] = result.iterations;
        SET_VECTOR_ELT(r_result, 2, iterations);
        UNPROTECT(1);
    }

    // Convert Eigen::VectorXd beta to R vector
    {
        int beta_size = result.beta.size();
        SEXP beta = PROTECT(Rf_allocVector(REALSXP, beta_size));
        for(int i = 0; i < beta_size; i++) {
            REAL(beta)[i] = result.beta(i);
        }
        SET_VECTOR_ELT(r_result, 3, beta);
        UNPROTECT(1);
    }

    // Convert LOOCV errors if computed
    if (with_errors) {
        SEXP loocv_brier_errors = PROTECT(Rf_allocVector(REALSXP, window_size));
        for(int i = 0; i < window_size; i++) {
            REAL(loocv_brier_errors)[i] = result.loocv_brier_errors[i];
        }
        SET_VECTOR_ELT(r_result, 4, loocv_brier_errors);
        UNPROTECT(1);

        SEXP loocv_deviance_errors = PROTECT(Rf_allocVector(REALSXP, window_size));
        for(int i = 0; i < window_size; i++) {
            REAL(loocv_deviance_errors)[i] = result.loocv_deviance_errors[i];
        }
        SET_VECTOR_ELT(r_result, 5, loocv_deviance_errors);
        UNPROTECT(1);
    }

    UNPROTECT(1);
    return r_result;
}

