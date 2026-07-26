#ifndef ADAPTIVE_NBHD_SIZE_H_
#define ADAPTIVE_NBHD_SIZE_H_

#include <vector>      // For std::vector usage throughout
#include "path_graphs.hpp"

struct adaptive_nbhd_size_t {
	// Graph structures
    std::vector<std::pair<std::vector<std::vector<int>>, std::vector<std::vector<double>>>> graphs;
    std::pair<std::vector<std::vector<int>>, std::vector<std::vector<double>>> opt_h_graph;

    // Parameter values
    std::vector<int> h_values;        ///< A vector of h values used to create hHN graphs
    int opt_h;                        ///< The value of h at which cv_errors achieves its miniumum value

    // Error metrics
    std::vector<double> cv_errors;    ///< Mean CV Error for each value of h
    std::vector<double> true_errors;  ///< The mean of true absolute deviation errors

    // Predictions
    std::vector<double> condEy;       ///< conditional expectation of y estimated over the h-hop graph with h corresponding to the smallest CV error
    std::vector<double> local_condEy; ///< conditional expectation of y for each vertex estimated over the h-hop graph with h corresponding to the smallest CV error

    // Bootstrap results
    std::vector<double> bb_condEy;    ///< the central tendency of conditional expectation of y as estimated using Bayesian bootstraps
    std::vector<double> cri_L;        ///< lower limit of credible intervals
    std::vector<double> cri_U;        ///< upper limit of credible intervals
};

#endif // ADAPTIVE_NBHD_SIZE_H_
