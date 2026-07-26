#ifndef BANDWIDTH_UTILS_HPP
#define BANDWIDTH_UTILS_HPP

#include <vector>
#include <cstddef>
using std::size_t;

/**
 * @brief Generates a set of candidate bandwidth values over the specified range
 */
std::vector<double> get_candidate_bws(
    double min_bw,
    double max_bw,
    size_t n_bws,
    bool log_grid,
    double min_spacing);

#endif // BANDWIDTH_UTILS_HPP
