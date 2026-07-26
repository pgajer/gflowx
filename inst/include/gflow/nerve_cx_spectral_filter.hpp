#ifndef NERVE_CX_SPECTRAL_FILTER_HPP
#define NERVE_CX_SPECTRAL_FILTER_HPP

#include "graph_spectral_filter.hpp" // For enums and structures

#include <vector>
#include <cstddef>

#include <Eigen/Core>

using std::size_t;

class simplex_weight_manager_t {
private:
    std::vector<double> vertex_weights;
    std::vector<double> edge_weights;
    std::vector<double> triangle_weights;
    std::vector<double> tetrahedron_weights;  // if applicable

    // Kernel parameters
    kernel_type_t kernel_type;
    double kernel_bandwidth;
    double weight_epsilon;

public:
    void compute_geometric_weights();
    void apply_kernel_transformation();
    void propagate_weights_downward();
    void validate_compatibility();
};


/**
 * @brief Results of a nerve complex spectral filtering operation
 */
struct nerve_cx_spectral_filter_t {
    // Result data
    std::vector<double>   evalues;        ///< Laplacian eigenvalues
    Eigen::MatrixXd       evectors;       ///< Laplacian eigenvectors (columns)
    std::vector<double>   candidate_ts;   ///< Diffusion times tested
    std::vector<double>   gcv_scores;     ///< GCV score for each t
    size_t                opt_t_idx;      ///< Index of optimal t in candidate_ts
    std::vector<double>   predictions;    ///< y_{t*} (optimal smooth)
    std::vector<std::vector<double>> t_predictions;  ///< Optional: full family y_{t_j}, size = candidate_ts.size()

    // Configuration parameters used
    laplacian_type_t      laplacian_type; ///< Type of Laplacian used
    filter_type_t         filter_type;    ///< Type of spectral filter applied
    size_t                laplacian_power;///< Power to which Laplacian was raised
    std::vector<double>   dim_weights;    ///< Weights for each dimension's contribution
    kernel_params_t       kernel_params;  ///< Kernel parameters used

    // Performance metrics
    double                compute_time_ms = 0.0;  ///< Computation time in milliseconds
    double                gcv_min_score = 0.0;    ///< Minimum GCV score achieved

    // Constructor with default values
    nerve_cx_spectral_filter_t() :
        laplacian_type(laplacian_type_t::STANDARD),
        filter_type(filter_type_t::HEAT),
        laplacian_power(1) {}
};

// Implementation file follows
#endif // NERVE_CX_SPECTRAL_FILTER_HPP
