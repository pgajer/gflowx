// ==========================================================================
// TIE-BREAKING AND CONDITIONAL EXPECTATION PREPARATION METHODS
// ==========================================================================
//
// These methods prepare conditional expectation estimates for extremality
// analysis by enforcing proper bounds, applying winsorization to handle
// outliers from slow diffusion convergence, and breaking ties to ensure
// unique function values required for Morse-Smale complex construction.

#include "riem_dcx.hpp"
#include <R.h>
#include <Rinternals.h>
#include <random>          // For std::mt19937, std::uniform_real_distribution
#include <unordered_map>   // For std::unordered_map
#include <unordered_set>   // For std::unordered_set
#include <algorithm>       // For std::sort, std::max
#include <cmath>           // For std::abs, std::isfinite

/**
 * @brief Break ties in function values with adaptive noise
 *
 * Adds minimal noise to break ties while preserving global extrema. The noise
 * scale adapts to the data range, ensuring perturbations are negligible relative
 * to the function's variation. For boundary values (global min/max), the method
 * ensures at least one point remains at the exact boundary when preserve_bounds
 * is enabled.
 *
 * @param y Function values (will be modified in place)
 * @param noise_scale Relative noise as fraction of range (default 1e-10)
 * @param min_abs_noise Minimum absolute noise magnitude (default 1e-12)
 * @param preserve_bounds Keep exact min/max (default true)
 * @param seed Random seed for reproducibility (use negative for random)
 * @param verbose Print diagnostic information
 *
 * @details
 * The method handles three cases for tied values:
 * - Global minimum: first instance stays at exact value, others perturb upward
 * - Global maximum: first instance stays at exact value, others perturb downward
 * - Interior values: symmetric perturbation centered at original value
 *
 * Noise magnitude is max(range * noise_scale, min_abs_noise) to ensure both
 * relative scaling and absolute minimum perturbation size.
 */
void riem_dcx_t::break_ties(
    vec_t& y,
    double noise_scale,
    double min_abs_noise,
    bool preserve_bounds,
    int seed,
    verbose_level_t verbose_level
) const {

    const Eigen::Index n = y.size();

    if (n == 0) {
        return;
    }

    // Check for any NaN values
    for (Eigen::Index i = 0; i < n; ++i) {
        if (!std::isfinite(y[i])) {
            Rf_error("break_ties: y contains non-finite values at index %lld",
                     static_cast<long long>(i));
        }
    }

    // Compute range
    const double y_min = y.minCoeff();
    const double y_max = y.maxCoeff();
    const double y_range = y_max - y_min;

    if (y_range == 0.0) {
        Rf_error("break_ties: y is constant, cannot break ties");
    }

    // Determine noise magnitude
    const double noise_magnitude = std::max(y_range * noise_scale, min_abs_noise);

    // Initialize RNG
    std::mt19937 rng;
    if (seed >= 0) {
        rng.seed(static_cast<unsigned int>(seed));
    } else {
        rng.seed(std::random_device{}());
    }
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    // Find tied values
    std::unordered_map<double, std::vector<Eigen::Index>> tied_groups;
    for (Eigen::Index i = 0; i < n; ++i) {
        tied_groups[y[i]].push_back(i);
    }

    // Count tied values
    size_t n_tied_values = 0;
    for (const auto& pair : tied_groups) {
        if (pair.second.size() > 1) {
            n_tied_values++;
        }
    }

    if (vl_at_least(verbose_level, verbose_level_t::TRACE) && n_tied_values == 0) {
        Rprintf("break_ties: No ties to break\n");
        return;
    }

    if (vl_at_least(verbose_level, verbose_level_t::TRACE)) {
        Rprintf("Breaking ties:\n");
        Rprintf("  Tied values: %zu\n", n_tied_values);
        Rprintf("  Noise magnitude: %.3e\n", noise_magnitude);
        Rprintf("  As %% of range: %.4f%%\n",
                100.0 * noise_magnitude / y_range);
    }

    // Process each tied value
    for (auto& pair : tied_groups) {
        const double val = pair.first;
        std::vector<Eigen::Index>& indices = pair.second;
        const size_t n_tied = indices.size();

        if (n_tied == 1) {
            continue;  // Not tied
        }

        const bool is_global_min = (std::abs(val - y_min) < 1e-14);
        const bool is_global_max = (std::abs(val - y_max) < 1e-14);

        // Generate perturbations based on location
        if (preserve_bounds && is_global_min) {
            // Keep first at exact min, perturb others upward
            std::vector<double> perturbations(n_tied);
            perturbations[0] = 0.0;

            for (size_t j = 1; j < n_tied; ++j) {
                perturbations[j] = uniform(rng) * noise_magnitude;
            }

            // Sort perturbations (except first which stays at 0)
            std::sort(perturbations.begin() + 1, perturbations.end());

            for (size_t j = 0; j < n_tied; ++j) {
                y[indices[j]] = val + perturbations[j];
            }

        } else if (preserve_bounds && is_global_max) {
            // Keep first at exact max, perturb others downward
            std::vector<double> perturbations(n_tied);
            perturbations[0] = 0.0;

            for (size_t j = 1; j < n_tied; ++j) {
                perturbations[j] = -uniform(rng) * noise_magnitude;
            }

            // Sort perturbations in decreasing order (except first which stays at 0)
            std::sort(perturbations.begin() + 1, perturbations.end(),
                     std::greater<double>());

            for (size_t j = 0; j < n_tied; ++j) {
                y[indices[j]] = val + perturbations[j];
            }

        } else {
            // Symmetric perturbation for interior values
            std::vector<double> perturbations(n_tied);
            double sum = 0.0;

            for (size_t j = 0; j < n_tied; ++j) {
                perturbations[j] = (uniform(rng) - 0.5) * noise_magnitude;
                sum += perturbations[j];
            }

            // Center perturbations around zero
            const double mean_pert = sum / static_cast<double>(n_tied);
            for (size_t j = 0; j < n_tied; ++j) {
                y[indices[j]] = val + (perturbations[j] - mean_pert);
            }
        }
    }

    if (vl_at_least(verbose_level, verbose_level_t::TRACE)) {
        // Count unique values after tie breaking
        std::unordered_set<double> unique_vals;
        for (Eigen::Index i = 0; i < n; ++i) {
            unique_vals.insert(y[i]);
        }
        Rprintf("  Final unique values: %zu / %lld\n",
                unique_vals.size(), static_cast<long long>(n));
    }
}


/**
 * @brief Prepare binary conditional expectation for extrema detection
 *
 * Standard pipeline for E[Y|X] where Y is binary with values in {0,1}.
 * Enforces [0,1] bounds, applies right winsorization to handle outliers
 * from slow diffusion convergence at isolated vertices, and breaks ties
 * to ensure function values are unique.
 *
 * @param y_hat Estimated conditional expectation (will be modified in place)
 * @param p_right Right tail proportion for winsorization (default 0.01)
 * @param apply_right_winsorization Whether to winsorize right tail
 * @param noise_scale Noise scale for tie breaking (default 1e-10)
 * @param seed Random seed for tie breaking (default 123)
 * @param verbose Print diagnostic information
 *
 * @details
 * Processing pipeline:
 * 1. Floor at 0 (conditional expectation cannot be negative)
 * 2. Right winsorization (if enabled) to handle slow convergence outliers
 * 3. Ceiling at 1 (conditional expectation cannot exceed 1)
 * 4. Tie breaking with adaptive noise
 *
 * Right winsorization is particularly important for k-NN regression where
 * isolated vertices (those with few or distant neighbors) exhibit slow
 * diffusion convergence, leading to overestimated conditional expectations
 * in the right tail. These outliers can dominate extremality score
 * computations and should be winsorized before analysis.
 */
void riem_dcx_t::prepare_binary_cond_exp(
    vec_t& y_hat,
    double p_right,
    bool apply_right_winsorization,
    double noise_scale,
    int seed,
    verbose_level_t verbose_level
) const {

    const Eigen::Index n = y_hat.size();

    if (n == 0) {
        return;
    }

    if (vl_at_least(verbose_level, verbose_level_t::TRACE)) {
        Rprintf("\n=== PREPARING BINARY CONDITIONAL EXPECTATION ===\n");
        Rprintf("Original range: [%.6f, %.6f]\n",
                y_hat.minCoeff(), y_hat.maxCoeff());
    }

    // Step 1: Floor at 0
    if (vl_at_least(verbose_level, verbose_level_t::TRACE)) {
        Rprintf("Step 1: Flooring at 0\n");
    }

    int n_floored = 0;
    for (Eigen::Index i = 0; i < n; ++i) {
        if (y_hat[i] < 0.0) {
            y_hat[i] = 0.0;
            n_floored++;
        }
    }

    if (vl_at_least(verbose_level, verbose_level_t::TRACE) && n_floored > 0) {
        Rprintf("  Floored %d values (%.2f%%)\n",
                n_floored, 100.0 * n_floored / n);
        Rprintf("  Range after floor: [%.6f, %.6f]\n",
                y_hat.minCoeff(), y_hat.maxCoeff());
    }

    // Step 2: Right winsorization (if requested)
    if (apply_right_winsorization) {
        if (vl_at_least(verbose_level, verbose_level_t::TRACE)) {
            Rprintf("Step 2: Right winsorization (p = %.3f)\n", p_right);
        }

        if (p_right <= 0.0 || p_right >= 1.0) {
            Rf_error("prepare_binary_cond_exp: p_right must be in (0, 1), got %.3f",
                     p_right);
        }

        // Compute threshold as (1 - p_right) quantile
        std::vector<double> sorted_vals(n);
        for (Eigen::Index i = 0; i < n; ++i) {
            sorted_vals[i] = y_hat[i];
        }
        std::sort(sorted_vals.begin(), sorted_vals.end());

        const size_t threshold_idx = static_cast<size_t>(
            std::floor(static_cast<double>(n) * (1.0 - p_right))
        );
        const double threshold = sorted_vals[threshold_idx];

        // Winsorize values above threshold
        int n_winsorized = 0;
        for (Eigen::Index i = 0; i < n; ++i) {
            if (y_hat[i] > threshold) {
                y_hat[i] = threshold;
                n_winsorized++;
            }
        }

        if (vl_at_least(verbose_level, verbose_level_t::TRACE)) {
            Rprintf("  Threshold: %.6f\n", threshold);
            Rprintf("  Winsorized %d values (%.2f%%)\n",
                    n_winsorized, 100.0 * n_winsorized / n);
            Rprintf("  Range after winsorization: [%.6f, %.6f]\n",
                    y_hat.minCoeff(), y_hat.maxCoeff());
        }
    }

    // Step 3: Ceiling at 1
    if (vl_at_least(verbose_level, verbose_level_t::TRACE)) {
        Rprintf("Step 3: Ceiling at 1\n");
    }

    int n_ceilinged = 0;
    for (Eigen::Index i = 0; i < n; ++i) {
        if (y_hat[i] > 1.0) {
            y_hat[i] = 1.0;
            n_ceilinged++;
        }
    }

    if (vl_at_least(verbose_level, verbose_level_t::TRACE) && n_ceilinged > 0) {
        Rprintf("  Ceilinged %d values (%.2f%%)\n",
                n_ceilinged, 100.0 * n_ceilinged / n);
        Rprintf("  Range after ceiling: [%.6f, %.6f]\n",
                y_hat.minCoeff(), y_hat.maxCoeff());
    }

    // Step 4: Break ties
    if (vl_at_least(verbose_level, verbose_level_t::TRACE)) {
        Rprintf("Step 4: Breaking ties\n");
    }

    break_ties(y_hat, noise_scale, 1e-12, true, seed, verbose_level);

    if (vl_at_least(verbose_level, verbose_level_t::TRACE)) {
        Rprintf("Final range: [%.6f, %.6f]\n",
                y_hat.minCoeff(), y_hat.maxCoeff());
        Rprintf("=== PREPARATION COMPLETE ===\n\n");
    }
}
