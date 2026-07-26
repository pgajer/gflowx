#ifndef ADAPTIVE_LOCAL_LOGISTIC_H
#define ADAPTIVE_LOCAL_LOGISTIC_H

#include <vector>
#include <Eigen/Dense>

/**
 * @brief Parameters controlling adaptive bandwidth selection for local logistic regression
 *
 * @details This structure encapsulates the parameters that control how bandwidths
 * are adapted to local data characteristics in adaptive local logistic regression.
 * The adaptation process combines information about local data density and the
 * complexity of the probability surface to determine appropriate bandwidths for
 * each region.
 *
 * The bandwidth adaptation follows the formula:
 * h(x) = h0 * min{max{(f(x)/f_med)^(-0.2) * (beta2/beta2_med)^(-0.2), c_min}, c_max}
 * where:
 * - h0 is the global_bandwidth
 * - f(x) is the local density
 * - f_med is the median_density
 * - beta2 represents local complexity
 * - beta2_med is the median_complexity
 * - c_min is the min_bandwidth_factor
 * - c_max is the max_bandwidth_factor
 *
 * The power -0.2 in the adaptation formula comes from theoretical arguments about
 * optimal bandwidth rates, providing increased smoothing in regions of low density
 * or low complexity while allowing for more local fitting in regions of high
 * density or rapid probability changes.
 *
 * @param global_bandwidth Initial bandwidth used for pilot estimation.
 *        This serves as a baseline scale for exploring local data structure.
 *        Must be positive. Typical values range from 0.05 to 0.2 for data
 *        scaled to [0,1].
 *
 * @param median_density Reference density for scaling bandwidth adjustments.
 *        This is typically computed as the median of local density estimates
 *        across all points using the global bandwidth. Initially set to 1.0
 *        and updated during the fitting process.
 *
 * @param median_complexity Reference complexity for scaling bandwidth adjustments.
 *        This is typically computed as the median absolute value of quadratic
 *        coefficients from local quadratic fits. Initially set to 1.0 and
 *        updated during the fitting process.
 *
 * @param min_bandwidth_factor Lower bound for bandwidth adjustment factor.
 *        Must be in (0,1]. Default 0.2 allows bandwidth to shrink to 20% of
 *        global bandwidth in regions of high density or complexity. Values
 *        too close to 0 may lead to unstable estimates.
 *
 * @param max_bandwidth_factor Upper bound for bandwidth adjustment factor.
 *        Must be >= 1.0. Default 5.0 allows bandwidth to grow up to 500% of
 *        global bandwidth in regions of low density or complexity. Larger
 *        values increase smoothing in sparse regions but may oversmooth.
 *
 * Example usage:
 * @code
 * // Create parameters with default values
 * adaptive_params_t params;
 *
 * // Create parameters with custom settings
 * adaptive_params_t custom_params(
 *     0.1,    // global_bandwidth
 *     1.0,    // median_density (will be updated)
 *     1.0,    // median_complexity (will be updated)
 *     0.2,    // min_bandwidth_factor
 *     5.0     // max_bandwidth_factor
 * );
 * @endcode
 *
 * @throws std::invalid_argument if parameters violate constraints:
 *         - global_bandwidth <= 0
 *         - min_bandwidth_factor <= 0 or > 1
 *         - max_bandwidth_factor < 1
 *
 * @note The median_density and median_complexity parameters are typically
 * updated automatically during the fitting process. Their initial values
 * are used only for the first iteration of bandwidth selection.
 *
 * @see adaptive_bandwidth_t for the structure storing computed bandwidths
 * and weights for specific points
 */
struct adaptive_params_t {
    double global_bandwidth;
    double median_density;
    double median_complexity;
    double min_bandwidth_factor;
    double max_bandwidth_factor;

    adaptive_params_t(
        double gb = 0.1,
        double md = 1.0,
        double mc = 1.0,
        double min_bf = 0.2,
        double max_bf = 5.0
    );
};

/**
 * @brief Computes adaptive bandwidths for local logistic regression
 *
 * @details This implementation follows a three-step approach:
 * 1. Estimates local data density using a pilot bandwidth
 * 2. Estimates local complexity of the probability surface
 * 3. Adjusts bandwidths based on both density and complexity
 *
 * The final bandwidth at each point adapts to:
 * - Use smaller bandwidths in dense regions with high complexity
 * - Use larger bandwidths in sparse regions or areas of low complexity
 * - Maintain smooth transitions between regions
 */
struct adaptive_bandwidth_t {
    double local_alpha;              // Local bandwidth parameter
    std::vector<double> weights;     // Computed weights for local regression
    double effective_sample_size;    // Sum of weights (for diagnostics)
};

double estimate_local_density(
    const std::vector<double>& x,
    int center_idx,
    double pilot_bandwidth);

double estimate_local_complexity(
    const std::vector<double>& x,
    const std::vector<double>& y,
    int center_idx,
    double pilot_bandwidth);

adaptive_bandwidth_t compute_adaptive_bandwidth(
    const std::vector<double>& x,
    const std::vector<double>& y,
    int center_idx,
    const adaptive_params_t& params);

#endif // ADAPTIVE_LOCAL_LOGISTIC_H
