#ifndef WASSERSTEIN_PERM_TEST_H_
#define WASSERSTEIN_PERM_TEST_H_

#include "wasserstein_dist.h" // for C_wasserstein_distance_1D()

#include <vector>
#include <cstddef>
using std::size_t;

/**
 * @brief Results structure for vertex-level Wasserstein permutation tests
 */
struct vertex_wasserstein_perm_test_results_t {
    std::vector<double> p_values;           ///< p-value for each vertex
    std::vector<double> effect_sizes;       ///< magnitude of effect at each vertex
    std::vector<bool> significant_vertices; ///< which vertices show significant effects
    std::vector<double> null_distances;     ///< distribution of null Wasserstein distances
};

vertex_wasserstein_perm_test_results_t vertex_wasserstein_perm_test(
    const std::vector<std::vector<double>>& bb_predictions,
    const std::vector<std::vector<double>>& null_predictions,
    size_t n_bootstraps,
    double alpha
	);

#endif // WASSERSTEIN_PERM_TEST_H_
