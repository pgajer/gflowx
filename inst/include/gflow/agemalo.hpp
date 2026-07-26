#ifndef AGEMALO_H_
#define AGEMALO_H_

#include "uniform_grid_graph.hpp"
#include "wasserstein_perm_test.hpp"

#include <cstddef>
using std::size_t;

struct agemalo_result_t {
    double graph_diameter;
    double max_packing_radius;

    std::vector<double> predictions;      ///< predictions[i] is an averaged-model estimate of E(Y|G) at the i-th vertex of the original (before uniform grid construction) graph G
    std::vector<double> errors;           ///< errors[i] is an averaged-model estimate of the prediction error at the i-th vertex of the original (before uniform grid construction) graph
    std::vector<double> scale;            ///< scale[i] is a local scale at the i-th vertex of the original (before uniform grid construction) graph G, which is an approximate radius of a disk in G where predictions is well approximated by a linear model

    std::vector<double> grid_opt_bw;      ///< grid_opt_bw[grid_vertex] = the optimal bandwidth over all models at this vertex; this gives a local scale at each grid vertex
    std::unordered_map<size_t, double> grid_predictions_map; ///< model-averaged predictions at the grid vertices; Note that for the Morse-Smale cells construction we only need estimate of E(y|G) over grid points; this would require reformulation of the notion of tau-gradient flow

    // Bootstrap fields (if n_bb > 0)
    std::vector<std::vector<double>> bb_predictions;
    std::vector<double> cri_lower;
    std::vector<double> cri_upper;

    // Permutation fields (if n_perms > 0)
    std::vector<std::vector<double>> null_predictions;
    std::vector<double> null_predictions_cri_lower;
    std::vector<double> null_predictions_cri_upper;

    // Permutation Tests
    std::optional<vertex_wasserstein_perm_test_results_t> permutation_tests;
};

agemalo_result_t agemalo(
    const uniform_grid_graph_t& grid_graph,
    const std::vector<double>& y,
    // geodesic parameter
    size_t min_path_size,
    // bw parameters
    size_t n_bws,
    bool log_grid,
    double min_bw_factor,
    double max_bw_factor,
    // kernel parameters
    double dist_normalization_factor,
    size_t kernel_type,
    // model parameters
    double model_tolerance,
    double model_blending_coef,
    // Bayesian bootstrap parameters
    size_t n_bb,
    double cri_probability,
    // permutation parameters
    size_t n_perms,
    // verbose
    bool verbose
    );

#endif // AGEMALO_H_
