#ifndef MSR2_DIFF_SMOOTHER_H_
#define MSR2_DIFF_SMOOTHER_H_

#include "omp_compat.h"

#include <cstddef>
#include <vector>
#include <memory>
#include <algorithm>
#include <cmath>
#include <random>
#include <chrono>
#include <set>
#include <map>
#include <queue>
#include <limits>
#include <cstdio>  // for sprintf
#include <sstream>  // Add this at the top of your file if not already present

#include <Eigen/Core>

#include <R.h>
#include <Rinternals.h>

using std::size_t;

enum class imputation_method_t {
    LOCAL_MEAN_THRESHOLD,            // the mean of y computed over the training vertices
    NEIGHBORHOOD_MATCHING,           // the matching method
    ITERATIVE_NEIGHBORHOOD_MATCHING, // the iterative matching method
    SUPPLIED_THRESHOLD,              // the supplied threshold,
    GLOBAL_MEAN_THRESHOLD            // the global mean of y (over all vertices)
};

struct iterative_imputation_params_t {
    int max_iterations = 10;
    double convergence_threshold = 1e-6;
};

struct graph_diffusion_smoother_result_t {
    std::vector<std::vector<double>> y_traj;
    std::vector<double> cv_errors;      // An (n_time_steps)-by-(n_CVs) matrix of CV errors for each step
    std::vector<double> mean_cv_errors; // Mean CV errors across iterations for each time step
    std::vector<double> y_optimal;      // Smoothed y values at the optimal time step
    int n_time_steps;
    int n_CVs;
    int optimal_time_step;              // Time step where mean CV error is minimal
    double min_cv_error;                // Minimum mean CV error value
};

struct graph_diffusion_matrix_smoother_result_t {
    std::vector<std::vector<std::vector<double>>> X_traj;
    std::vector<double> mean_cv_error;
};

struct graph_spectral_smoother_result_t {
    Eigen::VectorXd evalues;
    Eigen::MatrixXd evectors;
    int optimal_num_eigenvectors;
    std::vector<double> y_smoothed;
    Eigen::MatrixXd cv_errors;
    std::vector<double> mean_cv_errors;
    Eigen::MatrixXd low_pass_ys;
    int n_filters;
    int min_num_eigenvectors;
    int max_num_eigenvectors;
};

std::unique_ptr<graph_spectral_smoother_result_t>
graph_spectral_smoother(const std::vector<std::vector<int>>& graph,
                        const std::vector<std::vector<double>>& d,
                        const std::vector<double>& weights,
                        const std::vector<double>& y,
                        imputation_method_t imputation_method,
                        iterative_imputation_params_t iterative_params,
                        bool apply_binary_threshold,
                        double binary_threshold,
                        int ikernel,
                        double dist_normalization_factor,
                        int n_CVs,
                        int n_CV_folds,
                        double epsilon,
                        double min_plambda,
                        double max_plambda,
                        unsigned int seed);

// Define convergence criteria types
enum class diffn_sm_convergence_criterion_t {
    AVERAGE_CHANGE,
    MAX_CHANGE
};

/*!
 * @brief Performance metrics container for graph diffusion smoothing process.
 *
 * @details This structure collects and computes various performance metrics during
 * the graph diffusion smoothing process, including trajectory data, energy metrics,
 * and truth-based validation metrics when ground truth is available.
 *
 * The metrics are organized into several categories:
 * - Trajectory and update metrics (vertex values, deltas, step sizes)
 * - Event counting (oscillations, step size adjustments)
 * - Energy metrics (smoothness, fidelity, Laplacian)
 * - Truth-based validation metrics (SNR, MAD)
 *
 * Most vectors are sized according to the number of time steps, while vertex-specific
 * counters are sized according to the number of vertices.
 *
 * @note Truth-based metrics (SNR, MAD) should only be used during algorithm development
 *       and tuning with synthetic data where ground truth is known.
 */
struct graph_diffusion_smoother_performance_t {
    //! Trajectory of vertex values over time: y_trajectory[t][v] is vertex v's value at time t
    std::vector<std::vector<double>> y_trajectory;

    //! Delta values before updates: pre_update_deltas[t][v] is vertex v's delta at time t
    std::vector<std::vector<double>> pre_update_deltas;

    //! Delta values after updates: post_update_deltas[t][v] is vertex v's actual change at time t
    std::vector<std::vector<double>> post_update_deltas;

    //! History of step sizes: step_size_history[t][v] is vertex v's step size at time t
    std::vector<std::vector<double>> step_size_history;

    //! L1 norm of deltas at each iteration
    std::vector<double> global_residual_norm;

    //! Maximum absolute delta at each iteration
    std::vector<double> max_absolute_delta;

    //! Oscillation detection flags: oscillation_events[t][v] true if vertex v oscillated at time t
    std::vector<std::vector<bool>> oscillation_events;

    //! Total oscillation count for each vertex
    std::vector<int> oscillation_count_per_vertex;

    //! Count of step size increases for each vertex
    std::vector<int> increase_events_per_vertex;

    //! Count of step size decreases for each vertex
    std::vector<int> decrease_events_per_vertex;

    //! Count of oscillation-based step size reductions for each vertex
    std::vector<int> oscillation_reductions_per_vertex;

    // Energy metrics
    //! Smoothness energy trajectory: \f$ E_{smoothness} = \sum_{i} \sum_{j \in N(i)} w_{ij} (y_i - y_j)^2 \f$
    std::vector<double> smoothness_energy;

    //! Fidelity energy trajectory: \f$ E_{fidelity} = \sum_{i} (y_i - y^0_i)^2 \f$
    std::vector<double> fidelity_energy;

    //! Laplacian energy trajectory: \f$ E_{\Delta} = \sum_{i} |\sum_{j \in N(i)} w_{ij} (y_i - y_j)| \f$
    std::vector<double> laplacian_energy;

    //! Ratio of smoothness to fidelity energy
    std::vector<double> energy_ratio;

    // Truth-based metrics
    //! Initial signal-to-noise ratio
    double initial_snr;

    //! SNR trajectory over time
    std::vector<double> snr_trajectory;

    //! Mean absolute deviation from truth over time
    std::vector<double> mean_absolute_deviation;

    //! Pointwise curvature error over time
    std::vector<double> pointwise_curvature_error;

    //! Integrated curvature error over time
    std::vector<double> integrated_curvature_error;


    /*!
     * @brief Constructs a performance metrics container.
     * @param n_vertices Number of vertices in the graph
     * @param n_time_steps Number of time steps to be recorded
     */
    graph_diffusion_smoother_performance_t(int n_vertices, int n_time_steps) {
        // Reserve space for time series data
        y_trajectory.reserve(n_time_steps + 1);  // +1 for initial state
        pre_update_deltas.reserve(n_time_steps);
        post_update_deltas.reserve(n_time_steps);
        step_size_history.reserve(n_time_steps);
        global_residual_norm.reserve(n_time_steps);
        max_absolute_delta.reserve(n_time_steps);
        oscillation_events.reserve(n_time_steps);
        smoothness_energy.reserve(n_time_steps);
        fidelity_energy.reserve(n_time_steps);
        laplacian_energy.reserve(n_time_steps);
        energy_ratio.reserve(n_time_steps);
        snr_trajectory.reserve(n_time_steps);
        mean_absolute_deviation.reserve(n_time_steps);
        pointwise_curvature_error.reserve(n_time_steps);
        integrated_curvature_error.reserve(n_time_steps);

        // Initialize per-vertex counters
        oscillation_count_per_vertex.resize(n_vertices, 0);
        increase_events_per_vertex.resize(n_vertices, 0);
        decrease_events_per_vertex.resize(n_vertices, 0);
        oscillation_reductions_per_vertex.resize(n_vertices, 0);
    }

    /*!
     * @brief Computes signal-to-noise ratio given current state and ground truth.
     * @param y_current Current vertex values
     * @param y_true Ground truth vertex values
     * @return SNR in decibels
     * @throws std::invalid_argument if vectors have different sizes or are empty
     */
    double compute_snr(const std::vector<double>& y_current,
                      const std::vector<double>& y_true) const {
        if (y_current.empty() || y_true.empty() || y_current.size() != y_true.size()) {
            throw std::invalid_argument("Invalid input vectors for SNR computation");
        }

        double signal_power = 0.0;
        double noise_power = 0.0;
        for (size_t i = 0; i < y_true.size(); ++i) {
            double true_val = y_true[i];
            double curr_val = y_current[i];

            // Check for NaN/Inf in input values
            #if 0
            if (!std::isfinite(true_val) || !std::isfinite(curr_val)) {
                Rf_error(sprintf("Warning: Non-finite value detected in SNR computation at index %d: true=%f, current=%f", i, true_val, curr_val));
                return std::numeric_limits<double>::lowest();  // Return sentinel value
            }
            #endif

            signal_power += true_val * true_val;
            double noise = true_val - curr_val;
            noise_power += noise * noise;
        }

        if (noise_power <= 0.0) {
            Rf_error("Zero or negative noise power in SNR computation");
        }

        if (noise_power <= std::numeric_limits<double>::epsilon()) {
            Rprintf("Warning: Near-zero noise power in SNR computation");
            return std::numeric_limits<double>::max();  // Perfect SNR
        }

        double snr = 10.0 * std::log10(signal_power / noise_power);
        if (!std::isfinite(snr)) {
            std::ostringstream msg;
            msg << "Warning: Non-finite SNR computed: " << snr
                << " (signal_power=" << signal_power
                << ", noise_power=" << noise_power << ")";
            Rf_warning("%s", msg.str().c_str());
            return std::numeric_limits<double>::lowest();
        }

        return snr;
    }

    /*!
     * @brief Computes mean absolute deviation from ground truth.
     * @param y_current Current vertex values
     * @param y_true Ground truth vertex values
     * @return Mean absolute deviation
     * @throws std::invalid_argument if vectors have different sizes or are empty
     */
    double compute_mad(const std::vector<double>& y_current,
                      const std::vector<double>& y_true) const {
        if (y_current.empty() || y_true.empty() || y_current.size() != y_true.size()) {
            throw std::invalid_argument("Invalid input vectors for MAD computation");
        }

        double sum_abs_dev = 0.0;
        for (size_t i = 0; i < y_true.size(); ++i) {
            sum_abs_dev += std::abs(y_true[i] - y_current[i]);
        }
        return sum_abs_dev / y_true.size();
    }
};

#endif // MSR2_DIFF_SMOOTHER_H_
