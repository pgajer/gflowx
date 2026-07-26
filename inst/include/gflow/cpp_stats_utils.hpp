#ifndef CPP_STATS_UTILS_HPP
#define CPP_STATS_UTILS_HPP

#include <vector>
#include <algorithm>
#include <cstddef>
using std::size_t;

std::vector<double> running_window_average(
	const std::vector<double>& values,
	int window_size
	);

double jaccard_index(const std::vector<int>& A,
                     const std::vector<int>& B);

/**
 * @brief Calculates the p-th quantile of a vector of values.
 *
 * This function computes the p-th quantile (percentile) of the input vector.
 * For example, p = 0.25 returns the 25th percentile, p = 0.5 returns the median.
 * The function uses linear interpolation when the quantile position falls between two elements.
 *
 * @tparam T The type of elements in the vector (must support sorting and arithmetic operations)
 * @param x The input vector of values
 * @param p The percentile to calculate (must be between 0 and 1 inclusive)
 * @return T The calculated p-th quantile value
 * @throws std::invalid_argument If the input vector is empty or p is outside [0,1]
 *
 * @note This implementation sorts a copy of the input vector, so the original data is not modified.
 *
 * @example
 * std::vector<double> data = {1.0, 5.0, 3.0, 9.0, 7.0};
 * double q1 = quantile(data, 0.25);  // 25th percentile
 * double median = quantile(data, 0.5);  // median (50th percentile)
 * double q3 = quantile(data, 0.75);  // 75th percentile
 */
template <typename T>
T quantile(
    const std::vector<T>& x,
    double p) {

    if (x.empty()) {
        REPORT_ERROR("Input vector cannot be empty");
    }

    if (p < 0.0 || p > 1.0) {
        REPORT_ERROR("Percentile p must be between 0 and 1");
    }

    // Create a copy of the input vector to avoid modifying the original
    std::vector<T> sorted = x;

    // Sort the values
    std::sort(sorted.begin(), sorted.end());

    // Calculate the position
    double pos = p * (sorted.size() - 1);

    // Get the integer and fractional parts of the position
    size_t idx = static_cast<size_t>(pos);
    double fraction = pos - idx;

    // If the position is an integer, just return the value at that position
    if (fraction == 0.0) {
        return sorted[idx];
    }

    // Otherwise, interpolate between the adjacent values
    return sorted[idx] * (1.0 - fraction) + sorted[idx + 1] * fraction;
}

double calculate_path_evenness(const std::vector<double>& edge_lengths);
double calculate_path_max_threshold(
    const std::vector<double>& edge_lengths,
    double lower_threshold,
    double upper_threshold
	);

#endif // CPP_STATS_UTILS_HPP
