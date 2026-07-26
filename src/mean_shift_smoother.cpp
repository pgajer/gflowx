#include "mean_shift_smoother.hpp"
#include "cpp_utils.hpp"
#include "SEXP_cpp_conversion_utils.hpp"
#include "graph_diffusion_smoother.hpp"
#include "stats_utils.h"
#include "kernels.h"
#include "kNN.h"

#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <limits>
#include <cmath>
#include <memory>

#include <R.h>
#include <Rinternals.h>

extern "C" {
    SEXP S_mean_shift_data_smoother(SEXP s_X,
                                    SEXP s_k,
                                    SEXP s_density_k,
                                    SEXP s_n_steps,
                                    SEXP s_step_size,
                                    SEXP s_ikernel,
                                    SEXP s_dist_normalization_factor,
                                    SEXP s_method,
                                    SEXP s_momentum,
                                    SEXP s_increase_factor,
                                    SEXP s_decrease_factor);

    SEXP S_mean_shift_data_smoother_with_grad_field_averaging(SEXP s_X,
                                                              SEXP s_k,
                                                              SEXP s_density_k,
                                                              SEXP s_n_steps,
                                                              SEXP s_step_size,
                                                              SEXP s_ikernel,
                                                              SEXP s_dist_normalization_factor,
                                                              SEXP s_average_direction_only);

    SEXP S_mean_shift_data_smoother_adaptive(SEXP s_X,
                                             SEXP s_k,
                                             SEXP s_density_k,
                                             SEXP s_n_steps,
                                             SEXP s_step_size,
                                             SEXP s_ikernel,
                                             SEXP s_dist_normalization_factor,
                                             SEXP s_average_direction_only);
}

/**
 * @brief Performs mean shift smoothing on a matrix of data points.
 *
 * This function implements a mean shift algorithm to smooth a set of data
 * points. It uses k-nearest neighbors (kNN) to estimate local density and
 * compute gradients. The function also tracks the median of distances to the
 * k-th nearest neighbor as a proxy for data entropy.
 *
 * @param X Input matrix of data points (n_points x n_features)
 * @param k Number of nearest neighbors for gradient estimation
 * @param density_k Number of nearest neighbors for density estimation
 * @param n_steps Number of iterations for the mean shift algorithm
 * @param step_size Step size for updating point positions
 * @param ikernel Type of kernel function to use (default: 1 - Epanechnikov).
 *               Available kernels:
 *               - 0-Constant,
 *               - 1-Epanechnikov,
 *               - 2-Triangular,
 *               - 3-TrExponential,
 *               - 4-Laplace,
 *               - 5-Normal,
 *               - 6-Biweight,
 *               - 7-Tricube,
 *               - 8-Cosine
 *               - 9-Hyperbolic
 * @param dist_normalization_factor Factor for normalizing distances (default: 1.01)
 * @return std::unique_ptr<mean_shift_smoother_results_t> Results containing:
 *         - X_traj: Trajectory of data points over iterations
 *         - median_of_kdistances: Median of k-distances for each iteration
 */
std::unique_ptr<mean_shift_smoother_results_t> mean_shift_data_smoother(const std::vector<std::vector<double>>& X,
                                                                          int k,
                                                                          int density_k,
                                                                          int n_steps,
                                                                          double step_size,
                                                                          int ikernel = 1,
                                                                          double dist_normalization_factor = 1.01) {
    const int n_points = X.size();
    const int n_features = X[0].size();

    double scale = 1.0;
    initialize_kernel(ikernel, scale);

    // Precompute KNN for all points
    k++; // Since kNN returns for each point, the point itself among its
         // neighbors, there are only k - 1 neighbors returned that are
         // different from the point; This is why this k++; This also implies
         // that to access distances of neighbors different from the point one
         // has to cycle through knn_res.distances[point * k + j] with j from 1
         // to k - 1
    int k_minus_one = k - 1;
    auto knn_res = kNN(X, k);

    density_k++;
    int density_k_minus_one = density_k - 1;

    std::vector<double> kernel_weights(k_minus_one);
    std::vector<double> local_distances(k_minus_one);
    std::vector<double> kdistances(n_points);
    std::vector<double> median_kdistances(n_steps);

    std::vector<std::vector<std::vector<double>>> X_traj(n_steps, X);

    // Main loop for mean shift steps
    for (int step = 1; step < n_steps; ++step) {
        std::vector<std::vector<double>> gradient_vector_field(n_points, std::vector<double>(n_features));

        // Estimating data density gradient for each point
        for (int point = 0; point < n_points; ++point) {
            std::vector<double> new_point(n_features, 0.0);

            // Finding the maximum distance among k-1 nearest neighbors
          double max_dist = 0.0;
            for (int j = 0; j < k_minus_one; ++j) {
                local_distances[j] = knn_res.distances[point * k + j + 1];
                if (local_distances[j] > max_dist) max_dist = local_distances[j];
            }

            // Normalizing distances
            if (max_dist == 0) max_dist = 1;  // Avoid division by zero
            max_dist *= dist_normalization_factor;

            for (int j = 0; j < k_minus_one; ++j)
                local_distances[j] /= max_dist;

            // Computing kernel weights
            kernel_fn(local_distances.data(), k_minus_one, kernel_weights.data());

            // Computing weighted mean of neighboring points
            double total_weight = 0.0;
            for (int j = 0; j < k_minus_one; ++j) {
                int neighbor_idx = knn_res.indices[point * k + j + 1];
                for (int i = 0; i < n_features; ++i) {
                    new_point[i] += kernel_weights[j] * X_traj[step-1][neighbor_idx][i];
                }
                total_weight += kernel_weights[j];
            }

            if (total_weight <= 0.0) total_weight = 1e-12;

            // Computing and estimate of data density gradient
            double gradient_norm = 0.0;
            for (int i = 0; i < n_features; ++i) {
                new_point[i] /= total_weight;
                gradient_vector_field[point][i] = new_point[i] - X_traj[step-1][point][i];
                gradient_norm += gradient_vector_field[point][i] * gradient_vector_field[point][i];
            }
            gradient_norm = std::sqrt(gradient_norm);

            // Normalize gradient
            if (gradient_norm > 0) {
                for (int i = 0; i < n_features; ++i) {
                    gradient_vector_field[point][i] /= gradient_norm;
                }
            }
        }

        // Updating X_traj
        for (int point = 0; point < n_points; ++point) {
            double gradient_magnitude = 0.0;
            double total_weight = 0.0;

            // Normalizing distances (same as above)
            double max_dist = 0.0;
            for (int j = 0; j < k_minus_one; ++j) {
                local_distances[j] = knn_res.distances[point * k + j + 1];
                if (local_distances[j] > max_dist) max_dist = local_distances[j];
            }

            if (max_dist == 0) max_dist = 1;  // Avoid division by zero
            max_dist *= dist_normalization_factor;

            for (int j = 0; j < k_minus_one; ++j) {
                local_distances[j] /= max_dist;
            }

            // Computing kernel weights
            kernel_fn(local_distances.data(), k_minus_one, kernel_weights.data());

            // Computing gradient magnitude
            for (int j = 0; j < k_minus_one; ++j) {
                int neighbor_idx = knn_res.indices[point * k + j + 1];
                double proj = 0.0;
                for (int i = 0; i < n_features; ++i) {
                    proj += (X_traj[step-1][neighbor_idx][i] - X_traj[step-1][point][i]) * gradient_vector_field[point][i];
                }
                gradient_magnitude += kernel_weights[j] * std::abs(proj);
                total_weight += kernel_weights[j];
            }

            if (total_weight <= 0.0) total_weight = 1e-12;

            gradient_magnitude /= total_weight;

            // Updating point position
            for (int i = 0; i < n_features; ++i) {
                X_traj[step][point][i] = X_traj[step-1][point][i] + step_size * gradient_magnitude * gradient_vector_field[point][i];
            }
        }

        // Computing distances to density_k_minus_one nearest neighbor for all points
        auto step_knn_res = kNN(X_traj[step], density_k);

        for (int point = 0; point < n_points; ++point)
            kdistances[point] = step_knn_res.distances[point * density_k + density_k_minus_one]; // distance to the density_k-th NN

        median_kdistances[step] = median(kdistances.data(), (int)kdistances.size());
    }

    auto results = std::make_unique<mean_shift_smoother_results_t>();
    results->X_traj = std::move(X_traj);
    results->median_kdistances = std::move(median_kdistances);

    return results;
}

/**
 * @brief Performs adaptive mean shift smoothing on a dataset with dynamic k-nearest neighbors computation.
 *
 * This function implements an adaptive mean shift algorithm for smoothing data points in a
 * high-dimensional space. It dynamically computes k-nearest neighbors at each iteration,
 * allowing for more accurate tracking of the evolving data structure during the smoothing process.
 *
 * The algorithm iteratively shifts each point towards the weighted mean of its neighbors,
 * with the neighborhood redefined at each step based on the current state of the data.
 * This approach is particularly useful for datasets where the local structure changes
 * significantly during the smoothing process.
 *
 * @param X A vector of vectors representing the input data points. Each inner vector
 *          represents a single data point, with its elements being the features or
 *          coordinates in the high-dimensional space.
 * @param k The number of nearest neighbors to consider for each point when computing
 *          the mean shift vector. The actual number used will be k-1, as the point
 *          itself is included in the k-NN search.
 * @param density_k The number of nearest neighbors to use when estimating the local
 *                  density. This is used to compute the median k-distances at each step.
 * @param n_steps The number of iterations to perform in the mean shift process.
 * @param step_size A scaling factor for the magnitude of each mean shift step.
 *                  Larger values may lead to faster convergence but might overshoot;
 *                  smaller values provide more accurate but slower shifts.
 * @param ikernel Type of kernel function to use (default: 1 - Epanechnikov).
 *               Available kernels:
 *               - 0-Constant,
 *               - 1-Epanechnikov,
 *               - 2-Triangular,
 *               - 3-TrExponential,
 *               - 4-Laplace,
 *               - 5-Normal,
 *               - 6-Biweight,
 *               - 7-Tricube,
 *               - 8-Cosine
 *               - 9-Hyperbolic
 * @param dist_normalization_factor A factor used to normalize distances when computing
 *                                  kernel weights. Default is 1.01. This helps prevent
 *                                  division by zero and adjusts the influence of distance
 *                                  on the kernel weights.
 *
 * @return A unique pointer to a mean_shift_smoother_results_t structure containing:
 *         - X_traj: A vector of vector of vectors representing the trajectory of each
 *                   point over the iterations. X_traj[i] gives the state of all points
 *                   at iteration i.
 *         - median_kdistances: A vector of median k-distances at each iteration,
 *                              which can be used to estimate the density at each step.
 *
 * @note This function recomputes k-nearest neighbors at each iteration, which can be
 *       computationally expensive but provides more accurate results, especially for
 *       datasets where the local structure changes significantly during smoothing.
 *
 * @Rf_warning The function assumes that the input data X is not empty and that all
 *          points have the same number of features (dimensions). It does not perform
 *          explicit checks for these conditions.
 *
 * @see kNN() function for nearest neighbor computations
 * @see initialize_kernel() function for setting up the kernel function
 * @see kernel_fn() function for computing kernel weights
 * @see median() function for computing median k-distances
 *
 * @todo Explore parallelization opportunities, especially for the k-NN computations
 *       and per-point calculations within each iteration.
 * @todo Implement an early stopping criterion based on the convergence of point positions
 *       or changes in local density estimates.
 */
std::unique_ptr<mean_shift_smoother_results_t> knn_adaptive_mean_shift_smoother(const std::vector<std::vector<double>>& X,
                                                                                int k,
                                                                                int density_k,
                                                                                int n_steps,
                                                                                double step_size,
                                                                                int ikernel = 1,
                                                                                double dist_normalization_factor = 1.01) {
    const int n_points = X.size();
    const int n_features = X[0].size();

    double scale = 1.0;
    initialize_kernel(ikernel, scale);

    k++; // Adjust k to include the point itself
    int k_minus_one = k - 1;
    density_k++;
    int density_k_minus_one = density_k - 1;

    std::vector<double> kernel_weights(k_minus_one);
    std::vector<double> local_distances(k_minus_one);
    std::vector<double> kdistances(n_points);
    std::vector<double> median_kdistances(n_steps);

    std::vector<std::vector<std::vector<double>>> X_traj(n_steps);
    X_traj[0] = X;  // Initialize the first step with the original data

    for (int step = 1; step < n_steps; ++step) {

        int Kmax = std::max(k, density_k);
        auto knn_res = kNN(X_traj[step-1], Kmax);

        X_traj[step].resize(n_points, std::vector<double>(n_features));  // Prepare space for the current step
        std::vector<std::vector<double>> gradient_vector_field(n_points, std::vector<double>(n_features));

        for (int point = 0; point < n_points; ++point) {
            std::vector<double> new_point(n_features, 0.0);
            double max_dist = 0.0;

            for (int j = 0; j < k_minus_one; ++j) {
                local_distances[j] = knn_res.distances[point * Kmax + j + 1];
                if (local_distances[j] > max_dist) max_dist = local_distances[j];
            }

            max_dist = std::max(max_dist, 1.0) * dist_normalization_factor;

            for (int j = 0; j < k_minus_one; ++j)
                local_distances[j] /= max_dist;

            kernel_fn(local_distances.data(), k_minus_one, kernel_weights.data());

            double total_weight = 0.0;
            double gradient_magnitude = 0.0;

            for (int j = 0; j < k_minus_one; ++j) {
                int neighbor_idx = knn_res.indices[point * Kmax + j + 1];
                double proj = 0.0;

                for (int i = 0; i < n_features; ++i) {
                    double diff = X_traj[step-1][neighbor_idx][i] - X_traj[step-1][point][i];
                    new_point[i] += kernel_weights[j] * X_traj[step-1][neighbor_idx][i];
                    proj += diff * gradient_vector_field[point][i];
                }

                total_weight += kernel_weights[j];
                gradient_magnitude += kernel_weights[j] * std::abs(proj);
            }

            if (total_weight <= 0.0) total_weight = 1e-12;

            double gradient_norm = 0.0;
            for (int i = 0; i < n_features; ++i) {
                new_point[i] /= total_weight;
                gradient_vector_field[point][i] = new_point[i] - X_traj[step-1][point][i];
                gradient_norm += gradient_vector_field[point][i] * gradient_vector_field[point][i];
            }
            gradient_norm = std::sqrt(gradient_norm);

            if (gradient_norm > 0) {
                for (int i = 0; i < n_features; ++i) {
                    gradient_vector_field[point][i] /= gradient_norm;
                }
            }

            gradient_magnitude /= total_weight;

            for (int i = 0; i < n_features; ++i) {
                X_traj[step][point][i] = X_traj[step-1][point][i] + step_size * gradient_magnitude * gradient_vector_field[point][i];
            }

            kdistances[point] = knn_res.distances[point * Kmax + density_k_minus_one];
        }

        median_kdistances[step] = median(kdistances.data(), (int)kdistances.size());
    }

    auto results = std::make_unique<mean_shift_smoother_results_t>();
    results->X_traj = std::move(X_traj);
    results->median_kdistances = std::move(median_kdistances);

    return results;
}


/**
 * @brief Performs mean shift smoothing on a dataset with precomputed kernel weights.
 *
 * This function implements a mean shift algorithm for smoothing data points in a
 * high-dimensional space. It precomputes kernel weights for efficiency and performs
 * iterative updates to shift data points towards areas of higher density.
 *
 * @param X A vector of vectors representing the input data points. Each inner vector
 *          represents a single data point, with its elements being the features or
 *          coordinates in the high-dimensional space.
 * @param k The number of nearest neighbors to consider for each point when computing
 *          the mean shift vector. The actual number used will be k-1, as the point
 *          itself is included in the k-NN search.
 * @param density_k The number of nearest neighbors to use when estimating the local
 *                  density. This is used to compute the median k-distances at each step.
 * @param n_steps The number of iterations to perform in the mean shift process.
 * @param step_size A scaling factor for the magnitude of each mean shift step.
 *                  Larger values may lead to faster convergence but might overshoot;
 *                  smaller values provide more accurate but slower shifts.
 * @param ikernel Type of kernel function to use (default: 1 - Epanechnikov).
 *               Available kernels:
 *               - 0-Constant,
 *               - 1-Epanechnikov,
 *               - 2-Triangular,
 *               - 3-TrExponential,
 *               - 4-Laplace,
 *               - 5-Normal,
 *               - 6-Biweight,
 *               - 7-Tricube,
 *               - 8-Cosine
 *               - 9-Hyperbolic
 * @param dist_normalization_factor A factor used to normalize distances when computing
 *                                  kernel weights. Default is 1.01.
 *
 * @return A unique pointer to a mean_shift_smoother_results_t structure containing:
 *         - X_traj: A vector of vector of vectors representing the trajectory of each
 *                   point over the iterations.
 *         - median_kdistances: A vector of median k-distances at each iteration,
 *                              which can be used to estimate the density at each step.
 *
 * @note This function precomputes kernel weights for all points at the beginning,
 *       which can significantly improve performance for large datasets or when
 *       running many iterations.
 *
 * @Rf_warning The function assumes that the input data X is not empty and that all
 *          points have the same number of features (dimensions).
 *
 * @see kNN() function for nearest neighbor computations
 * @see initialize_kernel() function for setting up the kernel function
 * @see kernel_fn() function for computing kernel weights
 * @see median() function for computing median k-distances
 *
 * @todo Explore possibilities for parallelization to further improve performance.
 */
std::unique_ptr<mean_shift_smoother_results_t> mean_shift_data_smoother_precomputed(const std::vector<std::vector<double>>& X,
                                                                                    int k,
                                                                                    int density_k,
                                                                                    int n_steps,
                                                                                    double step_size,
                                                                                    int ikernel = 1,
                                                                                    double dist_normalization_factor = 1.01) {
    const int n_points = X.size();
    const int n_features = X[0].size();

    double scale = 1.0;
    initialize_kernel(ikernel, scale);

    // Precompute KNN for all points
    k++; // Adjust k to include the point itself
    int k_minus_one = k - 1;
    auto knn_res = kNN(X, k);

    density_k++;
    int density_k_minus_one = density_k - 1;

    std::vector<std::vector<double>> kernel_weights(n_points, std::vector<double>(k_minus_one));
    std::vector<double> local_distances(k_minus_one);
    std::vector<double> kdistances(n_points);
    std::vector<double> median_kdistances(n_steps);

    std::vector<std::vector<std::vector<double>>> X_traj(n_steps, X);

    // Precompute kernel weights for all points
    for (int point = 0; point < n_points; ++point) {
        double max_dist = 0.0;
        for (int j = 0; j < k_minus_one; ++j) {
            local_distances[j] = knn_res.distances[point * k + j + 1];
            if (local_distances[j] > max_dist) max_dist = local_distances[j];
        }

        if (max_dist == 0) max_dist = 1;  // Avoid division by zero
        max_dist *= dist_normalization_factor;

        for (int j = 0; j < k_minus_one; ++j)
            local_distances[j] /= max_dist;

        kernel_fn(local_distances.data(), k_minus_one, kernel_weights[point].data());
    }

    // Main loop for mean shift steps
    for (int step = 1; step < n_steps; ++step) {
        std::vector<std::vector<double>> gradient_vector_field(n_points, std::vector<double>(n_features));

        // Estimating data density gradient for each point
        for (int point = 0; point < n_points; ++point) {
            std::vector<double> new_point(n_features, 0.0);
            double total_weight = 0.0;

            for (int j = 0; j < k_minus_one; ++j) {
                int neighbor_idx = knn_res.indices[point * k + j + 1];
                for (int i = 0; i < n_features; ++i) {
                    new_point[i] += kernel_weights[point][j] * X_traj[step-1][neighbor_idx][i];
                }
                total_weight += kernel_weights[point][j];
            }

            if (total_weight <= 0.0) total_weight = 1e-12;

            double gradient_norm = 0.0;
            for (int i = 0; i < n_features; ++i) {
                new_point[i] /= total_weight;
                gradient_vector_field[point][i] = new_point[i] - X_traj[step-1][point][i];
                gradient_norm += gradient_vector_field[point][i] * gradient_vector_field[point][i];
            }
            gradient_norm = std::sqrt(gradient_norm);

            if (gradient_norm > 0) {
                for (int i = 0; i < n_features; ++i) {
                    gradient_vector_field[point][i] /= gradient_norm;
                }
            }
        }

        // Updating X_traj
        for (int point = 0; point < n_points; ++point) {
            double gradient_magnitude = 0.0;
            double total_weight = 0.0;

            for (int j = 0; j < k_minus_one; ++j) {
                int neighbor_idx = knn_res.indices[point * k + j + 1];
                double proj = 0.0;
                for (int i = 0; i < n_features; ++i) {
                    proj += (X_traj[step-1][neighbor_idx][i] - X_traj[step-1][point][i]) * gradient_vector_field[point][i];
                }
                gradient_magnitude += kernel_weights[point][j] * std::abs(proj);
                total_weight += kernel_weights[point][j];
            }

            if (total_weight <= 0.0) total_weight = 1e-12;

            gradient_magnitude /= total_weight;

            for (int i = 0; i < n_features; ++i) {
                X_traj[step][point][i] = X_traj[step-1][point][i] + step_size * gradient_magnitude * gradient_vector_field[point][i];
            }
        }

        // Computing distances to density_k_minus_one nearest neighbor for all points
        auto step_knn_res = kNN(X_traj[step], density_k);

        for (int point = 0; point < n_points; ++point)
            kdistances[point] = step_knn_res.distances[point * density_k + density_k_minus_one];

        median_kdistances[step] = median(kdistances.data(), (int)kdistances.size());
    }

    auto results = std::make_unique<mean_shift_smoother_results_t>();
    results->X_traj = std::move(X_traj);
    results->median_kdistances = std::move(median_kdistances);

    return results;
}

/**
 * @brief Performs mean shift data smoothing with gradient field averaging.
 *
 * This function implements an advanced mean shift algorithm for data smoothing,
 * incorporating gradient field averaging to enhance stability and denoising performance.
 * It can operate in two modes: full gradient vector averaging or direction-only averaging.
 *
 * The algorithm iteratively moves each data point along the estimated gradient of the
 * underlying probability density function. The gradient field is computed and then
 * averaged using kernel weights, which can help stabilize the denoising process.
 *
 * @param X Input data matrix, where each row represents a data point and each column represents a feature.
 * @param k Number of nearest neighbors to consider for gradient estimation (excluding the point itself).
 * @param density_k Number of nearest neighbors to use for density estimation (excluding the point itself).
 * @param n_steps Number of iterations for the mean shift process.
 * @param step_size Step size for updating point positions in each iteration.
 * @param ikernel Type of kernel function to use (default: 1 - Epanechnikov).
 *               Available kernels:
 *               - 0-Constant,
 *               - 1-Epanechnikov,
 *               - 2-Triangular,
 *               - 3-TrExponential,
 *               - 4-Laplace,
 *               - 5-Normal,
 *               - 6-Biweight,
 *               - 7-Tricube,
 *               - 8-Cosine
 *               - 9-Hyperbolic
 * @param dist_normalization_factor Factor for normalizing distances (default is 1.01).
 * @param average_direction_only If true, only the directions of gradients are averaged;
 *                               if false, full gradient vectors are averaged (default is false).
 *
 * @return A unique pointer to a mean_shift_smoother_results_t structure containing:
 *         - X_traj: A vector of data point trajectories across all iterations.
 *         - median_kdistances: A vector of median k-distances for each iteration.
 *
 * @note The function uses k-nearest neighbors (kNN) for both gradient estimation and density estimation.
 * @note The gradient field averaging step helps to reduce noise and increase stability in the mean shift process.
 * @note When average_direction_only is true, the algorithm becomes invariant to local variations in gradient magnitude,
 *       which can be beneficial for data with varying gradient strengths across the feature space.
 * @note The choice between full gradient averaging and direction-only averaging can significantly impact results
 *       and should be chosen based on the characteristics of the data and the desired smoothing properties.
 *
 * @Rf_warning This function assumes that the input data matrix X is non-empty and all rows have the same number of columns.
 * @Rf_warning The function may be computationally intensive for large datasets or high numbers of iterations.
 *
 * @see kNN
 * @see initialize_kernel
 * @see kernel_fn
 * @see mean_shift_smoother_results_t
 *
 * @todo Add parallelization for improved performance on large datasets.
 * @todo Explore adaptive step size mechanisms for potentially faster convergence.
 */
std::unique_ptr<mean_shift_smoother_results_t> mean_shift_data_smoother_with_grad_field_averaging(const std::vector<std::vector<double>>& X,
                                                                                                  int k,
                                                                                                  int density_k,
                                                                                                  int n_steps,
                                                                                                  double step_size,
                                                                                                  int ikernel = 1,
                                                                                                  double dist_normalization_factor = 1.01,
                                                                                                  bool average_direction_only = false) {
    const int n_points = X.size();
    const int n_features = X[0].size();

    double scale = 1.0;
    initialize_kernel(ikernel, scale);

    // Precompute KNN for all points
    k++; // Adjust k to include the point itself
    int k_minus_one = k - 1;
    auto knn_res = kNN(X, k);

    density_k++;
    int density_k_minus_one = density_k - 1;

    std::vector<std::vector<double>> kernel_weights(n_points, std::vector<double>(k_minus_one));
    std::vector<double> local_distances(k_minus_one);
    std::vector<double> kdistances(n_points);
    std::vector<double> median_kdistances(n_steps);

    std::vector<std::vector<std::vector<double>>> X_traj(n_steps, X);

    // Precompute kernel weights for all points
    for (int point = 0; point < n_points; ++point) {
        double max_dist = 0.0;
        for (int j = 0; j < k_minus_one; ++j) {
            local_distances[j] = knn_res.distances[point * k + j + 1];
            if (local_distances[j] > max_dist) max_dist = local_distances[j];
        }

        if (max_dist == 0) max_dist = 1;  // Avoid division by zero
        max_dist *= dist_normalization_factor;

        for (int j = 0; j < k_minus_one; ++j)
            local_distances[j] /= max_dist;

        kernel_fn(local_distances.data(), k_minus_one, kernel_weights[point].data());
    }

    // Main loop for mean shift steps
    for (int step = 1; step < n_steps; ++step) {
        std::vector<std::vector<double>> gradient_vector_field(n_points, std::vector<double>(n_features));

        // Estimating data density gradient for each point
        for (int point = 0; point < n_points; ++point) {
            std::vector<double> new_point(n_features, 0.0);
            double total_weight = 0.0;

            for (int j = 0; j < k_minus_one; ++j) {
                int neighbor_idx = knn_res.indices[point * k + j + 1];
                for (int i = 0; i < n_features; ++i) {
                    new_point[i] += kernel_weights[point][j] * X_traj[step-1][neighbor_idx][i];
                }
                total_weight += kernel_weights[point][j];
            }

            if (total_weight <= 0.0) total_weight = 1e-12;

            double gradient_norm = 0.0;
            for (int i = 0; i < n_features; ++i) {
                new_point[i] /= total_weight;
                gradient_vector_field[point][i] = new_point[i] - X_traj[step-1][point][i];
                gradient_norm += gradient_vector_field[point][i] * gradient_vector_field[point][i];
            }
            gradient_norm = std::sqrt(gradient_norm);

            if (gradient_norm > 0) {
                for (int i = 0; i < n_features; ++i) {
                    gradient_vector_field[point][i] /= gradient_norm;
                }
            }
        }

        // Averaging gradient field using kernel weights
        std::vector<std::vector<double>> averaged_gradient_field(n_points, std::vector<double>(n_features));
        for (int point = 0; point < n_points; ++point) {
            std::vector<double> averaged_gradient(n_features, 0.0);
            double total_weight = 0.0;

            for (int j = 0; j < k_minus_one; ++j) {
                int neighbor_idx = knn_res.indices[point * k + j + 1];
                if (average_direction_only) {
                    double neighbor_grad_norm = 0.0;
                    for (int i = 0; i < n_features; ++i) {
                        neighbor_grad_norm += gradient_vector_field[neighbor_idx][i] * gradient_vector_field[neighbor_idx][i];
                    }
                    neighbor_grad_norm = std::sqrt(neighbor_grad_norm);

                    if (neighbor_grad_norm > 0) {
                        for (int i = 0; i < n_features; ++i) {
                            averaged_gradient[i] += kernel_weights[point][j] * (gradient_vector_field[neighbor_idx][i] / neighbor_grad_norm);
                        }
                    }
                } else {
                    for (int i = 0; i < n_features; ++i) {
                        averaged_gradient[i] += kernel_weights[point][j] * gradient_vector_field[neighbor_idx][i];
                    }
                }
                total_weight += kernel_weights[point][j];
            }

            if (total_weight <= 0.0) total_weight = 1e-12;

            double gradient_norm = 0.0;
            for (int i = 0; i < n_features; ++i) {
                averaged_gradient[i] /= total_weight;
                gradient_norm += averaged_gradient[i] * averaged_gradient[i];
            }
            gradient_norm = std::sqrt(gradient_norm);

            if (gradient_norm > 0) {
                for (int i = 0; i < n_features; ++i) {
                    averaged_gradient_field[point][i] = averaged_gradient[i] / gradient_norm;
                }
            } else {
                averaged_gradient_field[point] = gradient_vector_field[point];
            }
        }

        // Updating X_traj
        for (int point = 0; point < n_points; ++point) {
            double gradient_magnitude = 0.0;
            double total_weight = 0.0;

            for (int j = 0; j < k_minus_one; ++j) {
                int neighbor_idx = knn_res.indices[point * k + j + 1];
                double proj = 0.0;
                for (int i = 0; i < n_features; ++i) {
                    proj += (X_traj[step-1][neighbor_idx][i] - X_traj[step-1][point][i]) * averaged_gradient_field[point][i];
                }
                gradient_magnitude += kernel_weights[point][j] * std::abs(proj);
                total_weight += kernel_weights[point][j];
            }

            if (total_weight <= 0.0) total_weight = 1e-12;

            gradient_magnitude /= total_weight;

            for (int i = 0; i < n_features; ++i) {
                X_traj[step][point][i] = X_traj[step-1][point][i] + step_size * gradient_magnitude * averaged_gradient_field[point][i];
            }
        }

        // Computing distances to density_k_minus_one nearest neighbor for all points
        auto step_knn_res = kNN(X_traj[step], density_k);

        for (int point = 0; point < n_points; ++point)
            kdistances[point] = step_knn_res.distances[point * density_k + density_k_minus_one];

        median_kdistances[step] = median(kdistances.data(), (int)kdistances.size());
    }

    auto results = std::make_unique<mean_shift_smoother_results_t>();
    results->X_traj = std::move(X_traj);
    results->median_kdistances = std::move(median_kdistances);

    return results;
}

/**
 * @brief Adaptive Mean Shift Data Smoother with Gradient Field Averaging
 *
 * This function implements an adaptive version of the mean shift algorithm with gradient field averaging.
 * It iteratively moves data points along the estimated gradient of the underlying probability density function,
 * using nearest neighbor information that is updated at each iteration.
 *
 * @param X Input data matrix, where each row represents a data point and each column represents a feature.
 * @param k Number of nearest neighbors to consider for gradient estimation (excluding the point itself).
 * @param density_k Number of nearest neighbors to use for density estimation (excluding the point itself).
 * @param n_steps Number of iterations for the mean shift process.
 * @param step_size Step size for updating point positions in each iteration.
 * @param ikernel Type of kernel function to use (default: 1 - Epanechnikov).
 *               Available kernels:
 *               - 0-Constant,
 *               - 1-Epanechnikov,
 *               - 2-Triangular,
 *               - 3-TrExponential,
 *               - 4-Laplace,
 *               - 5-Normal,
 *               - 6-Biweight,
 *               - 7-Tricube,
 *               - 8-Cosine
 *               - 9-Hyperbolic
 * @param dist_normalization_factor Factor for normalizing distances (default is 1.01).
 * @param average_direction_only If true, only the directions of gradients are averaged;
 *                               if false, full gradient vectors are averaged (default is false).
 *
 * @return A unique pointer to a mean_shift_smoother_results_t structure containing:
 *         - X_traj: A vector of data point trajectories across all iterations.
 *         - median_kdistances: A vector of median k-distances for each iteration.
 *
 * @note The function is adaptive in the sense that it recomputes nearest neighbors at each iteration,
 *       allowing for more accurate gradient estimation as points move. However, the step size remains constant.
 *
 * @Rf_warning The current implementation does not include a mechanism to adapt the step size dynamically.
 *
 * @pre The input data matrix X must not be empty and all rows must have the same number of columns.
 * @pre k and density_k must be positive integers less than the number of points in X.
 * @pre n_steps must be a positive integer.
 * @pre step_size must be a positive real number.
 *
 * @see initialize_kernel
 * @see kNN
 * @see kernel_fn
 * @see median
 *
 * @example
 * std::vector<std::vector<double>> X = {{1.0, 2.0}, {2.0, 3.0}, {3.0, 4.0}};
 * int k = 2;
 * int density_k = 2;
 * int n_steps = 10;
 * double step_size = 0.1;
 * auto results = adaptive_mean_shift_data_smoother_with_grad_field_averaging(X, k, density_k, n_steps, step_size);
 * // Process results...
 */
std::unique_ptr<mean_shift_smoother_results_t> knn_adaptive_mean_shift_data_smoother_with_grad_field_averaging(const std::vector<std::vector<double>>& X,
                                                                                                               int k,
                                                                                                               int density_k,
                                                                                                               int n_steps,
                                                                                                               double step_size,
                                                                                                               int ikernel = 1,
                                                                                                               double dist_normalization_factor = 1.01,
                                                                                                               bool average_direction_only = false) {
    const int n_points = X.size();
    const int n_features = X[0].size();

    double scale = 1.0;
    initialize_kernel(ikernel, scale);

    // Precompute KNN for all points
    k++; // Since kNN returns for each point, the point itself among its
         // neighbors, there are only k - 1 neighbors returned that are
         // different from the point; This is why this k++; This also implies
         // that to access distances of neighbors different from the point one
         // has to cycle through knn_res.distances[point * k + j] with j from 1
         // to k - 1
    int k_minus_one = k - 1;
    density_k++;
    int density_k_minus_one = density_k - 1;

    std::vector<double> kernel_weights(k_minus_one);
    std::vector<double> local_distances(k_minus_one);
    std::vector<double> kdistances(n_points);
    std::vector<double> median_kdistances(n_steps);

    std::vector<std::vector<std::vector<double>>> X_traj(n_steps);
    X_traj[0] = X;  // Initialize the first step with the original data

    // Main loop for mean shift steps
    for (int step = 1; step < n_steps; ++step) {

        int Kmax = std::max(k, density_k);
        auto knn_res = kNN(X_traj[step-1], Kmax);

        X_traj[step].resize(n_points, std::vector<double>(n_features));  // Prepare space for the current step

        std::vector<std::vector<double>> gradient_vector_field(n_points, std::vector<double>(n_features));

        // Estimating data density gradient for each point
        for (int point = 0; point < n_points; ++point) {

            std::vector<double> new_point(n_features, 0.0);

            // Finding the maximum distance among k-1 nearest neighbors
            double max_dist = 0.0;
            for (int j = 0; j < k_minus_one; ++j) {
                local_distances[j] = knn_res.distances[point * Kmax + j + 1];
                if (local_distances[j] > max_dist) max_dist = local_distances[j];
            }

            // Normalizing distances
            if (max_dist == 0) max_dist = 1;  // Avoid division by zero
            max_dist *= dist_normalization_factor;

            for (int j = 0; j < k_minus_one; ++j)
                local_distances[j] /= max_dist;

            // Computing kernel weights
            kernel_fn(local_distances.data(), k_minus_one, kernel_weights.data());

            // Computing weighted mean of neighboring points
            double total_weight = 0.0;
            for (int j = 0; j < k_minus_one; ++j) {
                int neighbor_idx = knn_res.indices[point * k + j + 1];
                for (int i = 0; i < n_features; ++i) {
                    new_point[i] += kernel_weights[j] * X[neighbor_idx][i];
                }
                total_weight += kernel_weights[j];
            }

            if (total_weight <= 0.0) total_weight = 1e-12;

            // Computing and estimate of data density gradient
            double gradient_norm = 0.0;
            for (int i = 0; i < n_features; ++i) {
                new_point[i] /= total_weight;
                gradient_vector_field[point][i] = new_point[i] - X_traj[step-1][point][i];
                gradient_norm += gradient_vector_field[point][i] * gradient_vector_field[point][i];
            }
            gradient_norm = std::sqrt(gradient_norm);

            // Normalize gradient
            if (gradient_norm > 0) {
                for (int i = 0; i < n_features; ++i) {
                    gradient_vector_field[point][i] /= gradient_norm;
                }
            }
        }

        // Averaging gradient field using kernel weights
        std::vector<std::vector<double>> averaged_gradient_field(n_points, std::vector<double>(n_features));
        for (int point = 0; point < n_points; ++point) {
            std::vector<double> averaged_gradient(n_features, 0.0);
            double total_weight = 0.0;

            // Finding the maximum distance among k-1 nearest neighbors
            double max_dist = 0.0;
            for (int j = 0; j < k_minus_one; ++j) {
                local_distances[j] = knn_res.distances[point * k + j + 1];
                if (local_distances[j] > max_dist) max_dist = local_distances[j];
            }

            // Normalizing distances
            if (max_dist == 0) max_dist = 1;  // Avoid division by zero
            max_dist *= dist_normalization_factor;

            for (int j = 0; j < k_minus_one; ++j)
                local_distances[j] /= max_dist;

            // Computing kernel weights
           kernel_fn(local_distances.data(), k_minus_one, kernel_weights.data());

            // Averaging gradient using kernel weights
            for (int j = 0; j < k_minus_one; ++j) {
                int neighbor_idx = knn_res.indices[point * Kmax + j + 1];
                if (average_direction_only) {
                    // Normalize the neighbor's gradient to a unit vector
                    double neighbor_grad_norm = 0.0;
                    for (int i = 0; i < n_features; ++i) {
                        neighbor_grad_norm += gradient_vector_field[neighbor_idx][i] * gradient_vector_field[neighbor_idx][i];
                    }
                    neighbor_grad_norm = std::sqrt(neighbor_grad_norm);

                    if (neighbor_grad_norm > 0) {
                        for (int i = 0; i < n_features; ++i) {
                            averaged_gradient[i] += kernel_weights[j] * (gradient_vector_field[neighbor_idx][i] / neighbor_grad_norm);
                        }
                    }
                } else {
                    // Original averaging of full gradient vectors
                    for (int i = 0; i < n_features; ++i) {
                        averaged_gradient[i] += kernel_weights[j] * gradient_vector_field[neighbor_idx][i];
                    }
                }
                total_weight += kernel_weights[j];
            }

            if (total_weight <= 0.0) total_weight = 1e-12;

            // Normalizing averaged gradient
            double gradient_norm = 0.0;
            for (int i = 0; i < n_features; ++i) {
                averaged_gradient[i] /= total_weight;
                gradient_norm += averaged_gradient[i] * averaged_gradient[i];
            }
            gradient_norm = std::sqrt(gradient_norm);

            if (gradient_norm > 0) {
                for (int i = 0; i < n_features; ++i) {
                    averaged_gradient_field[point][i] = averaged_gradient[i] / gradient_norm;
                }
            } else {
                averaged_gradient_field[point] = gradient_vector_field[point];
            }
        }

        // Updating X_traj
        for (int point = 0; point < n_points; ++point) {
            double gradient_magnitude = 0.0;
            double total_weight = 0.0;

            // Normalizing distances (same as above)
            double max_dist = 0.0;
            for (int j = 0; j < k_minus_one; ++j) {
                local_distances[j] = knn_res.distances[point * k + j + 1];
                if (local_distances[j] > max_dist) max_dist = local_distances[j];
            }

            if (max_dist == 0) max_dist = 1;  // Avoid division by zero
            max_dist *= dist_normalization_factor;

            for (int j = 0; j < k_minus_one; ++j) {
                local_distances[j] /= max_dist;
            }

            // Computing kernel weights
            kernel_fn(local_distances.data(), k_minus_one, kernel_weights.data());

            // Computing gradient magnitude
            for (int j = 0; j < k_minus_one; ++j) {
                int neighbor_idx = knn_res.indices[point * k + j + 1];
                double proj = 0.0;
                for (int i = 0; i < n_features; ++i) {
                    proj += (X[neighbor_idx][i] - X[point][i]) * gradient_vector_field[point][i];
                }
                gradient_magnitude += kernel_weights[j] * std::abs(proj);
                total_weight += kernel_weights[j];
            }

            if (total_weight <= 0.0) total_weight = 1e-12;
            gradient_magnitude /= total_weight;

            // Updating point position
            for (int i = 0; i < n_features; ++i) {
                X_traj[step][point][i] = X_traj[step-1][point][i] + step_size * gradient_magnitude * averaged_gradient_field[point][i];
            }
        }

        // Computing distances to density_k_minus_one nearest neighbor for all points
        auto step_knn_res = kNN(X_traj[step], density_k);

        for (int point = 0; point < n_points; ++point)
            kdistances[point] = step_knn_res.distances[point * density_k + density_k_minus_one]; // distance to the density_k-th NN

        median_kdistances[step] = median(kdistances.data(), (int)kdistances.size());
    }

    auto results = std::make_unique<mean_shift_smoother_results_t>();
    results->X_traj = std::move(X_traj);
    results->median_kdistances = std::move(median_kdistances);

    return results;
}

/**
 * @brief Fully Adaptive Mean Shift Data Smoother with Gradient Field Averaging
 *
 * This function implements a fully adaptive version of the mean shift algorithm with gradient field averaging.
 * It combines adaptive k-nearest neighbors computation and individual adaptive step sizes for each point.
 * The algorithm iteratively moves data points along the estimated gradient of the underlying probability
 * density function, adapting both the neighborhood and step size at each iteration.
 *
 * @param X Input data matrix, where each row represents a data point and each column represents a feature.
 * @param k Number of nearest neighbors to consider for gradient estimation (excluding the point itself).
 * @param density_k Number of nearest neighbors to use for density estimation (excluding the point itself).
 * @param n_steps Maximum number of iterations for the mean shift process.
 * @param initial_step_size Initial step size for updating point positions. This will be adapted for each point.
 * @param ikernel Type of kernel function to use (default: 1 - Epanechnikov).
 *               Available kernels:
 *               - 0-Constant,
 *               - 1-Epanechnikov,
 *               - 2-Triangular,
 *               - 3-TrExponential,
 *               - 4-Laplace,
 *               - 5-Normal,
 *               - 6-Biweight,
 *               - 7-Tricube,
 *               - 8-Cosine
 *               - 9-Hyperbolic
 * @param dist_normalization_factor Factor for normalizing distances (default is 1.01).
 * @param average_direction_only If true, only the directions of gradients are averaged;
 *                               if false, full gradient vectors are averaged (default is false).
 *
 * @return A unique pointer to a mean_shift_smoother_results_t structure containing:
 *         - X_traj: A vector of data point trajectories across all iterations.
 *         - median_kdistances: A vector of median k-distances for each iteration.
 *
 * @note Adaptive features:
 *       1. K-nearest neighbors are recomputed at each iteration based on updated point positions.
 *       2. Step sizes are individually adapted for each point based on gradient consistency.
 *       3. Momentum is applied to gradients for smoother convergence.
 *
 * @note Adaptive step size parameters (hardcoded in the current implementation):
 *       - momentum: Controls the influence of previous gradients (default: 0.9).
 *       - increase_factor: Factor to increase step size when gradients are consistent (default: 1.2).
 *       - decrease_factor: Factor to decrease step size when gradients are inconsistent (default: 0.5).
 *
 * @Rf_warning The adaptive step size mechanism may lead to very small or very large step sizes in certain scenarios.
 *
 * @pre The input data matrix X must not be empty and all rows must have the same number of columns.
 * @pre k and density_k must be positive integers less than the number of points in X.
 * @pre n_steps must be a positive integer.
 * @pre initial_step_size must be a positive real number.
 *
 * @see initialize_kernel
 * @see kNN
 * @see kernel_fn
 * @see median
 *
 * @example
 * std::vector<std::vector<double>> X = {{1.0, 2.0}, {2.0, 3.0}, {3.0, 4.0}};
 * int k = 2;
 * int density_k = 2;
 * int n_steps = 10;
 * double initial_step_size = 0.1;
 * auto results = adaptive_mean_shift_data_smoother_with_grad_field_averaging(
 *     X, k, density_k, n_steps, initial_step_size);
 * // Process results...
 */
std::unique_ptr<mean_shift_smoother_results_t> adaptive_mean_shift_data_smoother_with_grad_field_averaging(const std::vector<std::vector<double>>& X,
                                                                                                           int k,
                                                                                                           int density_k,
                                                                                                           int n_steps,
                                                                                                           double initial_step_size,
                                                                                                           int ikernel = 1,
                                                                                                           double dist_normalization_factor = 1.01,
                                                                                                           bool average_direction_only = false,
                                                                                                           double momentum = 0.9,
                                                                                                           double increase_factor = 1.2,
                                                                                                           double decrease_factor = 0.5) {
    const int n_points = X.size();
    const int n_features = X[0].size();

    double scale = 1.0;
    initialize_kernel(ikernel, scale);

    // Precompute KNN for all points
    k++; // Since kNN returns for each point, the point itself among its neighbors
    int k_minus_one = k - 1;
    density_k++;
    int density_k_minus_one = density_k - 1;

    std::vector<double> kernel_weights(k_minus_one);
    std::vector<double> local_distances(k_minus_one);
    std::vector<double> kdistances(n_points);
    std::vector<double> median_kdistances(n_steps);

    std::vector<std::vector<std::vector<double>>> X_traj(n_steps);
    X_traj[0] = X;  // Initialize the first step with the original data

    // Adaptive step size parameters
    std::vector<double> step_sizes(n_points, initial_step_size);
    std::vector<std::vector<double>> previous_gradients(n_points, std::vector<double>(n_features, 0.0));

    // Main loop for mean shift steps
    for (int step = 1; step < n_steps; ++step) {
        int Kmax = std::max(k, density_k);
        auto knn_res = kNN(X_traj[step-1], Kmax);

        X_traj[step].resize(n_points, std::vector<double>(n_features));  // Prepare space for the current step

        std::vector<std::vector<double>> gradient_vector_field(n_points, std::vector<double>(n_features));

        // Estimating data density gradient for each point
        for (int point = 0; point < n_points; ++point) {
            std::vector<double> new_point(n_features, 0.0);

            // Finding the maximum distance among k-1 nearest neighbors
            double max_dist = 0.0;
            for (int j = 0; j < k_minus_one; ++j) {
                local_distances[j] = knn_res.distances[point * Kmax + j + 1];
                if (local_distances[j] > max_dist) max_dist = local_distances[j];
            }

            // Normalizing distances
            if (max_dist == 0) max_dist = 1;  // Avoid division by zero
            max_dist *= dist_normalization_factor;

            for (int j = 0; j < k_minus_one; ++j)
                local_distances[j] /= max_dist;

            // Computing kernel weights
            kernel_fn(local_distances.data(), k_minus_one, kernel_weights.data());

            // Computing weighted mean of neighboring points
            double total_weight = 0.0;
            for (int j = 0; j < k_minus_one; ++j) {
                int neighbor_idx = knn_res.indices[point * Kmax + j + 1];
                for (int i = 0; i < n_features; ++i) {
                    new_point[i] += kernel_weights[j] * X_traj[step-1][neighbor_idx][i];
                }
                total_weight += kernel_weights[j];
            }

            if (total_weight <= 0.0) total_weight = 1e-12;

            // Computing and estimate of data density gradient
            double gradient_norm = 0.0;
            for (int i = 0; i < n_features; ++i) {
                new_point[i] /= total_weight;
                gradient_vector_field[point][i] = new_point[i] - X_traj[step-1][point][i];
                gradient_norm += gradient_vector_field[point][i] * gradient_vector_field[point][i];
            }
            gradient_norm = std::sqrt(gradient_norm);

            // Normalize gradient
            if (gradient_norm > 0) {
                for (int i = 0; i < n_features; ++i) {
                    gradient_vector_field[point][i] /= gradient_norm;
                }
            }
        }

        // Averaging gradient field using kernel weights
        std::vector<std::vector<double>> averaged_gradient_field(n_points, std::vector<double>(n_features));
        for (int point = 0; point < n_points; ++point) {
            std::vector<double> averaged_gradient(n_features, 0.0);
            double total_weight = 0.0;

            // Finding the maximum distance among k-1 nearest neighbors
            double max_dist = 0.0;
            for (int j = 0; j < k_minus_one; ++j) {
                local_distances[j] = knn_res.distances[point * k + j + 1];
                if (local_distances[j] > max_dist) max_dist = local_distances[j];
            }

            // Normalizing distances
            if (max_dist == 0) max_dist = 1;  // Avoid division by zero
            max_dist *= dist_normalization_factor;

            for (int j = 0; j < k_minus_one; ++j)
                local_distances[j] /= max_dist;

            // Computing kernel weights
            kernel_fn(local_distances.data(), k_minus_one, kernel_weights.data());

            // Averaging gradient using kernel weights
            for (int j = 0; j < k_minus_one; ++j) {
                int neighbor_idx = knn_res.indices[point * k + j + 1];
                if (average_direction_only) {
                    // Normalize the neighbor's gradient to a unit vector
                    double neighbor_grad_norm = 0.0;
                    for (int i = 0; i < n_features; ++i) {
                        neighbor_grad_norm += gradient_vector_field[neighbor_idx][i] * gradient_vector_field[neighbor_idx][i];
                    }
                    neighbor_grad_norm = std::sqrt(neighbor_grad_norm);

                    if (neighbor_grad_norm > 0) {
                        for (int i = 0; i < n_features; ++i) {
                            averaged_gradient[i] += kernel_weights[j] * (gradient_vector_field[neighbor_idx][i] / neighbor_grad_norm);
                        }
                    }
                } else {
                    // Original averaging of full gradient vectors
                    for (int i = 0; i < n_features; ++i) {
                        averaged_gradient[i] += kernel_weights[j] * gradient_vector_field[neighbor_idx][i];
                    }
                }
                total_weight += kernel_weights[j];
            }

            if (total_weight <= 0.0) total_weight = 1e-12;

            // Normalizing averaged gradient
            double gradient_norm = 0.0;
            for (int i = 0; i < n_features; ++i) {
                averaged_gradient[i] /= total_weight;
                gradient_norm += averaged_gradient[i] * averaged_gradient[i];
            }
            gradient_norm = std::sqrt(gradient_norm);

            if (gradient_norm > 0) {
                for (int i = 0; i < n_features; ++i) {
                    averaged_gradient_field[point][i] = averaged_gradient[i] / gradient_norm;
                }
            } else {
                averaged_gradient_field[point] = gradient_vector_field[point];
            }
        }

        // Updating X_traj with adaptive step size
        for (int point = 0; point < n_points; ++point) {
            double gradient_magnitude = 0.0;
            double total_weight = 0.0;

            // Normalizing distances (same as above)
            double max_dist = 0.0;
            for (int j = 0; j < k_minus_one; ++j) {
                local_distances[j] = knn_res.distances[point * k + j + 1];
                if (local_distances[j] > max_dist) max_dist = local_distances[j];
            }

            if (max_dist == 0) max_dist = 1;  // Avoid division by zero
            max_dist *= dist_normalization_factor;

            for (int j = 0; j < k_minus_one; ++j) {
                local_distances[j] /= max_dist;
            }

            // Computing kernel weights
            kernel_fn(local_distances.data(), k_minus_one, kernel_weights.data());

            // Computing gradient magnitude
            for (int j = 0; j < k_minus_one; ++j) {
                int neighbor_idx = knn_res.indices[point * k + j + 1];
                double proj = 0.0;
                for (int i = 0; i < n_features; ++i) {
                    proj += (X_traj[step-1][neighbor_idx][i] - X_traj[step-1][point][i]) * averaged_gradient_field[point][i];
                }
                gradient_magnitude += kernel_weights[j] * std::abs(proj);
                total_weight += kernel_weights[j];
            }
            if (total_weight <= 0.0) total_weight = 1e-12;
            gradient_magnitude /= total_weight;

            // Compute dot product between current and previous gradient
            double gradient_dot_product = 0.0;
            for (int i = 0; i < n_features; ++i) {
                gradient_dot_product += averaged_gradient_field[point][i] * previous_gradients[point][i];
            }

            // Adjust step size based on gradient direction consistency
            if (gradient_dot_product > 0) {
                // Gradients are pointing in similar directions, increase step size
                step_sizes[point] *= increase_factor;
            } else {
                // Gradients are pointing in different directions, decrease step size
                step_sizes[point] *= decrease_factor;
            }

            // Updating point position with adaptive step size
            for (int i = 0; i < n_features; ++i) {
                X_traj[step][point][i] = X_traj[step-1][point][i] +
                    step_sizes[point] * gradient_magnitude * averaged_gradient_field[point][i];
            }

            // Update previous gradient with momentum
            for (int i = 0; i < n_features; ++i) {
                previous_gradients[point][i] = momentum * previous_gradients[point][i] +
                    (1 - momentum) * averaged_gradient_field[point][i];
            }
        }

        // Computing distances to density_k_minus_one nearest neighbor for all points
        auto step_knn_res = kNN(X_traj[step], density_k);

        for (int point = 0; point < n_points; ++point)
            kdistances[point] = step_knn_res.distances[point * density_k + density_k_minus_one]; // distance to the density_k-th NN

        median_kdistances[step] = median(kdistances.data(), (int)kdistances.size());
    }

    auto results = std::make_unique<mean_shift_smoother_results_t>();
    results->X_traj = std::move(X_traj);
    results->median_kdistances = std::move(median_kdistances);

    return results;
}

/**
 * @brief R wrapper for the mean shift data smoother with gradient field averaging.
 *
 * This function serves as an interface between R and the C++ implementation of the
 * mean shift algorithm with gradient field averaging. It handles the conversion of
 * R objects to C++ types, calls the C++ function, and converts the results back to
 * R objects.
 *
 * @param s_X SEXP object representing the input data matrix. Each row is a data point,
 *            and each column is a feature.
 * @param s_k SEXP object (integer) representing the number of nearest neighbors for
 *            gradient estimation.
 * @param s_density_k SEXP object (integer) representing the number of nearest neighbors
 *                    for density estimation.
 * @param s_n_steps SEXP object (integer) representing the number of iterations for
 *                  the mean shift process.
 * @param s_step_size SEXP object (double) representing the step size for updating
 *                    point positions in each iteration.
 * @param s_ikernel SEXP object (integer) representing the kernel function identifier.
 *                  1: Gaussian, 2: Epanechnikov, 3: Quartic
 * @param s_dist_normalization_factor SEXP object (double) representing the factor
 *                                    for normalizing distances.
 * @param s_average_direction_only SEXP object (logical) indicating whether to average
 *                                 only the directions of gradients (TRUE) or full
 *                                 gradient vectors (FALSE).
 *
 * @return SEXP object (list) containing two elements:
 *         - X_traj: A list of matrices, each representing the positions of all points
 *                   at a given iteration.
 *         - median_kdistances: A vector of median k-distances for each iteration.
 *
 * @throws Rf_error if the mean_shift_smoother returns a null pointer or invalid X_traj.
 *
 * @note The input data is moved rather than copied for efficiency.
 * @note The function performs null checks on the returned pointer to ensure validity.
 *
 * @see mean_shift_data_smoother_with_grad_field_averaging
 * @see create_R_list
 *
 * @example
 * # R code
 * result <- .Call("S_mean_shift_data_smoother_with_grad_field_averaging",
 *                 X, k, density_k, n_steps, step_size, ikernel,
 *                 dist_normalization_factor, average_direction_only)
 * X_traj <- result$X_traj
 * median_kdistances <- result$median_kdistances
 */
SEXP S_mean_shift_data_smoother_with_grad_field_averaging(SEXP s_X,
                                                          SEXP s_k,
                                                          SEXP s_density_k,
                                                          SEXP s_n_steps,
                                                          SEXP s_step_size,
                                                          SEXP s_ikernel,
                                                          SEXP s_dist_normalization_factor,
                                                          SEXP s_average_direction_only) {
    // Convert SEXP inputs to C++ types
    std::vector<std::vector<double>> X = std::move(*Rmatrix_to_cpp(s_X));
    int k = INTEGER(s_k)[0];
    int density_k = INTEGER(s_density_k)[0];
    int n_steps = INTEGER(s_n_steps)[0];
    double step_size = REAL(s_step_size)[0];
    int ikernel = INTEGER(s_ikernel)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    bool average_direction_only = (LOGICAL(s_average_direction_only)[0] == 1);

    // Call the C++ function
    auto results = mean_shift_data_smoother_with_grad_field_averaging(X,
                                                                      k,
                                                                      density_k,
                                                                      n_steps,
                                                                      step_size,
                                                                      ikernel,
                                                                      dist_normalization_factor,
                                                                      average_direction_only);
    // Check if results is valid
    if (!results) {
        Rf_error("mean_shift_smoother returned null pointer");
    }

    // Check if X_traj is valid
    if (results->X_traj.empty() || results->X_traj[0].empty() || results->X_traj[0][0].empty()) {
        Rf_error("mean_shift_smoother returned invalid X_traj");
    }

    // Convert results to R list
    return create_R_list(results->X_traj, results->median_kdistances);
}


/**
 * @brief Performs mean shift data smoothing with an adaptive step size mechanism.
 *
 * This function implements an advanced mean shift algorithm for data smoothing,
 * incorporating gradient field averaging and an adaptive step size mechanism
 * to enhance stability, convergence speed, and denoising performance.
 *
 * The algorithm iteratively moves each data point along the estimated gradient of the
 * underlying probability density function. The gradient field is computed and then
 * averaged using kernel weights. The key feature of this version is the adaptive
 * step size mechanism, which adjusts the step size for each point individually
 * based on the consistency of gradient directions across iterations.
 *
 * @param X Input data matrix, where each row represents a data point and each column represents a feature.
 * @param k Number of nearest neighbors to consider for gradient estimation (excluding the point itself).
 * @param density_k Number of nearest neighbors to use for density estimation (excluding the point itself).
 * @param n_steps Maximum number of iterations for the mean shift process.
 * @param initial_step_size Initial step size for updating point positions.
 * @param ikernel Type of kernel function to use (default: 1 - Epanechnikov).
 *               Available kernels:
 *               - 0-Constant,
 *               - 1-Epanechnikov,
 *               - 2-Triangular,
 *               - 3-TrExponential,
 *               - 4-Laplace,
 *               - 5-Normal,
 *               - 6-Biweight,
 *               - 7-Tricube,
 *               - 8-Cosine
 *               - 9-Hyperbolic
 * @param dist_normalization_factor Factor for normalizing distances (default is 1.01).
 * @param average_direction_only If true, only the directions of gradients are averaged;
 *                               if false, full gradient vectors are averaged (default is false).
 *
 * @return A unique pointer to a mean_shift_smoother_results_t structure containing:
 *         - X_traj: A vector of data point trajectories across all iterations.
 *         - median_kdistances: A vector of median k-distances for each iteration.
 *
 * @note Adaptive Step Size Mechanism:
 *       The function employs an adaptive step size mechanism with the following characteristics:
 *       1. Individual step sizes: Each point has its own step size, allowing for localized adaptation.
 *       2. Gradient consistency check: The dot product between the current and previous gradient is used
 *          to assess the consistency of gradient directions.
 *       3. Step size adjustment:
 *          - If the dot product is positive (consistent direction), the step size is increased.
 *          - If the dot product is negative (inconsistent direction), the step size is decreased.
 *       4. Momentum: A momentum term is applied to smooth out oscillations in gradient directions.
 *
 * @note Adaptive Parameters (hardcoded in the current implementation):
 *       - momentum: Controls the influence of previous gradients (default: 0.9).
 *       - increase_factor: Factor to increase step size when gradients are consistent (default: 1.2).
 *       - decrease_factor: Factor to decrease step size when gradients are inconsistent (default: 0.5).
 *
 * @note The adaptive step size mechanism aims to:
 *       1. Speed up convergence in areas with consistent gradients.
 *       2. Improve stability by reducing step sizes in areas with complex local geometry.
 *       3. Adapt to the local data structure that each point is traversing.
 *
 * @Rf_warning This function assumes that the input data matrix X is non-empty and all rows have the same number of columns.
 * @Rf_warning The adaptive step size mechanism introduces additional parameters that may need tuning for optimal performance on specific datasets.
 * @Rf_warning The function may be computationally intensive for large datasets or high numbers of iterations.
 *
 * @see kNN
 * @see initialize_kernel
 * @see kernel_fn
 * @see mean_shift_smoother_results_t
 *
 * @todo Add user-configurable parameters for momentum, increase_factor, and decrease_factor.
 * @todo Implement minimum and maximum bounds for step sizes to prevent extreme values.
 * @todo Explore alternative convergence criteria based on step size magnitudes or gradient changes.
 * @todo Investigate potential performance optimizations for the adaptive mechanism on large datasets.
 */
std::unique_ptr<mean_shift_smoother_results_t> mean_shift_data_smoother_adaptive(const std::vector<std::vector<double>>& X,
                                                                                 int k,
                                                                                 int density_k,
                                                                                 int n_steps,
                                                                                 double initial_step_size,
                                                                                 int ikernel = 1,
                                                                                 double dist_normalization_factor = 1.01,
                                                                                 bool average_direction_only = false) {
    const int n_points = X.size();
    const int n_features = X[0].size();

    double scale = 1.0;
    initialize_kernel(ikernel, scale);

    // Precompute KNN for all points
    k++; // Since kNN returns for each point, the point itself among its neighbors
    int k_minus_one = k - 1;
    auto knn_res = kNN(X, k);

    density_k++;
    int density_k_minus_one = density_k - 1;

    std::vector<std::vector<double>> kernel_weights(n_points, std::vector<double>(k_minus_one));
    std::vector<double> local_distances(k_minus_one);
    std::vector<double> kdistances(n_points);
    std::vector<double> median_kdistances(n_steps);

    std::vector<std::vector<std::vector<double>>> X_traj(n_steps, X);

    // Adaptive step size parameters
    std::vector<double> step_sizes(n_points, initial_step_size);
    std::vector<std::vector<double>> previous_gradients(n_points, std::vector<double>(n_features, 0.0));
    double momentum = 0.9;
    double increase_factor = 1.2;
    double decrease_factor = 0.5;

    // Precompute kernel weights for all points
    for (int point = 0; point < n_points; ++point) {
        double max_dist = 0.0;
        for (int j = 0; j < k_minus_one; ++j) {
            local_distances[j] = knn_res.distances[point * k + j + 1];
            if (local_distances[j] > max_dist) max_dist = local_distances[j];
        }

        if (max_dist == 0) max_dist = 1;  // Avoid division by zero
        max_dist *= dist_normalization_factor;

        for (int j = 0; j < k_minus_one; ++j)
            local_distances[j] /= max_dist;

        kernel_fn(local_distances.data(), k_minus_one, kernel_weights[point].data());
    }

    // Main loop for mean shift steps
    for (int step = 1; step < n_steps; ++step) {
        std::vector<std::vector<double>> gradient_vector_field(n_points, std::vector<double>(n_features));

        // Estimating data density gradient for each point
        for (int point = 0; point < n_points; ++point) {
            std::vector<double> new_point(n_features, 0.0);
            double total_weight = 0.0;

            for (int j = 0; j < k_minus_one; ++j) {
                int neighbor_idx = knn_res.indices[point * k + j + 1];
                for (int i = 0; i < n_features; ++i) {
                    new_point[i] += kernel_weights[point][j] * X[neighbor_idx][i];
                }
                total_weight += kernel_weights[point][j];
            }

            if (total_weight <= 0.0) total_weight = 1e-12;

            double gradient_norm = 0.0;
            for (int i = 0; i < n_features; ++i) {
                new_point[i] /= total_weight;
                gradient_vector_field[point][i] = new_point[i] - X_traj[step-1][point][i];
                gradient_norm += gradient_vector_field[point][i] * gradient_vector_field[point][i];
            }
            gradient_norm = std::sqrt(gradient_norm);

            if (gradient_norm > 0) {
                for (int i = 0; i < n_features; ++i) {
                    gradient_vector_field[point][i] /= gradient_norm;
                }
            }
        }

        // Averaging gradient field using precomputed kernel weights
        std::vector<std::vector<double>> averaged_gradient_field(n_points, std::vector<double>(n_features));
        for (int point = 0; point < n_points; ++point) {
            std::vector<double> averaged_gradient(n_features, 0.0);
            double total_weight = 0.0;

            for (int j = 0; j < k_minus_one; ++j) {
                int neighbor_idx = knn_res.indices[point * k + j + 1];
                if (average_direction_only) {
                    double neighbor_grad_norm = 0.0;
                    for (int i = 0; i < n_features; ++i) {
                        neighbor_grad_norm += gradient_vector_field[neighbor_idx][i] * gradient_vector_field[neighbor_idx][i];
                    }
                    neighbor_grad_norm = std::sqrt(neighbor_grad_norm);

                    if (neighbor_grad_norm > 0) {
                        for (int i = 0; i < n_features; ++i) {
                            averaged_gradient[i] += kernel_weights[point][j] * (gradient_vector_field[neighbor_idx][i] / neighbor_grad_norm);
                        }
                    }
                } else {
                    for (int i = 0; i < n_features; ++i) {
                        averaged_gradient[i] += kernel_weights[point][j] * gradient_vector_field[neighbor_idx][i];
                    }
                }
                total_weight += kernel_weights[point][j];
            }

            if (total_weight <= 0.0) total_weight = 1e-12;

            double gradient_norm = 0.0;
            for (int i = 0; i < n_features; ++i) {
                averaged_gradient[i] /= total_weight;
                gradient_norm += averaged_gradient[i] * averaged_gradient[i];
            }
            gradient_norm = std::sqrt(gradient_norm);

            if (gradient_norm > 0) {
                for (int i = 0; i < n_features; ++i) {
                    averaged_gradient_field[point][i] = averaged_gradient[i] / gradient_norm;
                }
            } else {
                averaged_gradient_field[point] = gradient_vector_field[point];
            }
        }

        // Updating X_traj with adaptive step size
        for (int point = 0; point < n_points; ++point) {
            double gradient_magnitude = 0.0;
            double total_weight = 0.0;

            for (int j = 0; j < k_minus_one; ++j) {
                int neighbor_idx = knn_res.indices[point * k + j + 1];
                double proj = 0.0;
                for (int i = 0; i < n_features; ++i) {
                    proj += (X[neighbor_idx][i] - X[point][i]) * averaged_gradient_field[point][i];
                }
                gradient_magnitude += kernel_weights[point][j] * std::abs(proj);
                total_weight += kernel_weights[point][j];
            }

            if (total_weight <= 0.0) total_weight = 1e-12;
            gradient_magnitude /= total_weight;

            // Compute dot product between current and previous gradient
            double gradient_dot_product = 0.0;
            for (int i = 0; i < n_features; ++i) {
                gradient_dot_product += averaged_gradient_field[point][i] * previous_gradients[point][i];
            }

            // Adjust step size based on gradient direction consistency
            if (gradient_dot_product > 0) {
                step_sizes[point] *= increase_factor;
            } else {
                step_sizes[point] *= decrease_factor;
            }

            // Updating point position with adaptive step size
            for (int i = 0; i < n_features; ++i) {
                X_traj[step][point][i] = X_traj[step-1][point][i] +
                    step_sizes[point] * gradient_magnitude * averaged_gradient_field[point][i];
            }

            // Update previous gradient with momentum
            for (int i = 0; i < n_features; ++i) {
                previous_gradients[point][i] = momentum * previous_gradients[point][i] +
                    (1 - momentum) * averaged_gradient_field[point][i];
            }
        }

        // Computing distances to density_k_minus_one nearest neighbor for all points
        auto step_knn_res = kNN(X_traj[step], density_k);

        for (int point = 0; point < n_points; ++point)
            kdistances[point] = step_knn_res.distances[point * density_k + density_k_minus_one];

        median_kdistances[step] = median(kdistances.data(), (int)kdistances.size());
    }

    auto results = std::make_unique<mean_shift_smoother_results_t>();
    results->X_traj = std::move(X_traj);
    results->median_kdistances = std::move(median_kdistances);

    return results;
}


/**
 * @brief R wrapper for the adaptive mean shift data smoother.
 *
 * This function serves as an interface between R and the C++ implementation of the
 * adaptive mean shift algorithm. It handles the conversion of R objects to C++ types,
 * calls the C++ function, and converts the results back to R objects. The adaptive
 * version of the algorithm uses variable step sizes for potentially faster convergence.
 *
 * @param s_X SEXP object representing the input data matrix. Each row is a data point,
 *            and each column is a feature.
 * @param s_k SEXP object (integer) representing the number of nearest neighbors for
 *            gradient estimation.
 * @param s_density_k SEXP object (integer) representing the number of nearest neighbors
 *                    for density estimation.
 * @param s_n_steps SEXP object (integer) representing the maximum number of iterations
 *                  for the mean shift process.
 * @param s_step_size SEXP object (double) representing the initial step size for updating
 *                    point positions. This step size will be adapted during the process.
 * @param s_ikernel SEXP object (integer) representing the kernel function identifier.
 *                  1: Gaussian, 2: Epanechnikov, 3: Quartic
 * @param s_dist_normalization_factor SEXP object (double) representing the factor
 *                                    for normalizing distances.
 * @param s_average_direction_only SEXP object (logical) indicating whether to average
 *                                 only the directions of gradients (TRUE) or full
 *                                 gradient vectors (FALSE).
 *
 * @return SEXP object (list) containing two elements:
 *         - X_traj: A list of matrices, each representing the positions of all points
 *                   at a given iteration.
 *         - median_kdistances: A vector of median k-distances for each iteration.
 *
 * @throws Rf_error if the mean_shift_smoother_adaptive returns a null pointer or invalid X_traj.
 *
 * @note The input data is moved rather than copied for efficiency.
 * @note The function performs null checks on the returned pointer to ensure validity.
 * @note The adaptive algorithm adjusts step sizes dynamically based on the local data structure,
 *       which can lead to faster convergence in some cases.
 *
 * @see mean_shift_data_smoother_adaptive
 * @see create_R_list
 *
 * @example
 * # R code
 * result <- .Call("S_mean_shift_data_smoother_adaptive",
 *                 X, k, density_k, n_steps, initial_step_size, ikernel,
 *                 dist_normalization_factor, average_direction_only)
 * X_traj <- result$X_traj
 * median_kdistances <- result$median_kdistances
 */
SEXP S_mean_shift_data_smoother_adaptive(SEXP s_X,
                                         SEXP s_k,
                                         SEXP s_density_k,
                                         SEXP s_n_steps,
                                         SEXP s_step_size,
                                         SEXP s_ikernel,
                                         SEXP s_dist_normalization_factor,
                                         SEXP s_average_direction_only) {
    // Convert SEXP inputs to C++ types
    std::vector<std::vector<double>> X = std::move(*Rmatrix_to_cpp(s_X));
    int k = INTEGER(s_k)[0];
    int density_k = INTEGER(s_density_k)[0];
    int n_steps = INTEGER(s_n_steps)[0];
    double step_size = REAL(s_step_size)[0];
    int ikernel = INTEGER(s_ikernel)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    bool average_direction_only = (LOGICAL(s_average_direction_only)[0] == 1);

    // Call the C++ function
    auto results = mean_shift_data_smoother_adaptive(X,
                                                     k,
                                                     density_k,
                                                     n_steps,
                                                     step_size,
                                                     ikernel,
                                                     dist_normalization_factor,
                                                     average_direction_only);
    // Check if results is valid
    if (!results) {
        Rf_error("mean_shift_smoother returned null pointer");
    }

    // Check if X_traj is valid
    if (results->X_traj.empty() || results->X_traj[0].empty() || results->X_traj[0][0].empty()) {
        Rf_error("mean_shift_smoother returned invalid X_traj");
    }

    // Convert results to R list
    return create_R_list(results->X_traj, results->median_kdistances);
}



/**
 * @brief R-compatible version of mean shift smoothing algorithm.
 *
 * This function performs mean shift smoothing on a dataset and computes median k-distances.
 * It is designed to be called from R using .Call() and uses SEXP types for all parameters and
 * the return value.
 *
 * @param s_X SEXP (numeric matrix)
 *        The input dataset, where each row represents a point and columns represent features.
 *        In R, this should be a numeric matrix.
 *
 * @param s_k SEXP (integer)
 *        The number of nearest neighbors to consider for gradient estimation.
 *        In R, this should be a single integer value.
 *
 * @param s_density_k SEXP (integer)
 *        The number of nearest neighbors to consider for density estimation.
 *        In R, this should be a single integer value.
 *
 * @param s_n_steps SEXP (integer)
 *        The number of smoothing steps to perform.
 *        In R, this should be a single integer value.
 *
 * @param s_step_size SEXP (numeric)
 *        The step size for updating point positions in each iteration.
 *        In R, this should be a single numeric value, typically between 0 and 1.
 *
 * @param s_ikernel SEXP (integer)
 *        An integer specifying the kernel function to use. Valid values are:
 *        1 (Epanechnikov), 2 (Triangular), 3 (Truncated Exponential), 4 (Normal).
 *        In R, this should be a single integer value.
 *
 * @param s_dist_normalization_factor SEXP (numeric)
 *        A scaling factor applied to the maximum distance between a vertex and its neighbors.
 *        This ensures non-zero weights even when all distances are equal.
 *        In R, this should be a single numeric value, typically slightly above 1.0.
 *
 * @param s_method SEXP (integer)
 *        An integer specifying the smoothing method to use. Valid values are:
 *        0 (Basic), 1 (Precomputed), 2-3 (Gradient Field Averaging),
 *        4-5 (Adaptive), 6 (KNN Adaptive), 7-8 (KNN Adaptive with Gradient Field Averaging),
 *        9-10 (Adaptive with Gradient Field Averaging and Momentum).
 *        Odd-numbered methods use average direction only.
 *        In R, this should be a single integer value.
 *
 * @param s_momentum SEXP (numeric)
 *        The momentum factor for methods 9-10.
 *        In R, this should be a single numeric value between 0 and 1.
 *
 * @param s_increase_factor SEXP (numeric)
 *        The factor by which to increase the step size in adaptive methods.
 *        In R, this should be a single numeric value greater than 1.
 *
 * @param s_decrease_factor SEXP (numeric)
 *        The factor by which to decrease the step size in adaptive methods.
 *        In R, this should be a single numeric value between 0 and 1.
 *
 * @return SEXP (list)
 *         Returns an R list containing:
 *         - X_traj: A list of matrices, each representing the smoothed dataset at each step.
 *         - median_kdistances: A numeric vector of median k-distances, one for each step.
 *
 * @note This function assumes the existence of helper functions:
 *       - Rmatrix_to_cpp: Converts R matrix to C++ vector of vectors.
 *       - mean_shift_data_smoother: The main C++ implementation of the algorithm.
 *       - create_R_list: Converts C++ results back to R list format.
 *
 * @Rf_warning This function may be computationally expensive for large datasets or high values of n_steps.
 *          Ensure sufficient memory is available, especially for large datasets.
 */
SEXP S_mean_shift_data_smoother(SEXP s_X,
                                SEXP s_k,
                                SEXP s_density_k,
                                SEXP s_n_steps,
                                SEXP s_step_size,
                                SEXP s_ikernel,
                                SEXP s_dist_normalization_factor,
                                SEXP s_method,
                                SEXP s_momentum,
                                SEXP s_increase_factor,
                                SEXP s_decrease_factor) {
    // Convert SEXP inputs to C++ types
    std::vector<std::vector<double>> X = std::move(*Rmatrix_to_cpp(s_X));
    int k = INTEGER(s_k)[0];
    int density_k = INTEGER(s_density_k)[0];
    int n_steps = INTEGER(s_n_steps)[0];
    double step_size = REAL(s_step_size)[0];
    int ikernel = INTEGER(s_ikernel)[0];
    double dist_normalization_factor = REAL(s_dist_normalization_factor)[0];
    int method = INTEGER(s_method)[0];
    bool average_direction_only = (method % 2 == 1); // Odd methods use average_direction_only = true
    double momentum = REAL(s_momentum)[0];
    double increase_factor = REAL(s_increase_factor)[0];
    double decrease_factor = REAL(s_decrease_factor)[0];

    std::unique_ptr<mean_shift_smoother_results_t> results;

    try {
        switch(method) {
            case 0:
                results = mean_shift_data_smoother(X, k, density_k, n_steps, step_size, ikernel, dist_normalization_factor);
                break;
            case 1:
                results = mean_shift_data_smoother_precomputed(X, k, density_k, n_steps, step_size, ikernel, dist_normalization_factor);
                break;
            case 2:
            case 3:
                results = mean_shift_data_smoother_with_grad_field_averaging(X, k, density_k, n_steps, step_size, ikernel, dist_normalization_factor, average_direction_only);
                break;
            case 4:
            case 5:
                results = mean_shift_data_smoother_adaptive(X, k, density_k, n_steps, step_size, ikernel, dist_normalization_factor, average_direction_only);
                break;
            case 6:
                results = knn_adaptive_mean_shift_smoother(X, k, density_k, n_steps, step_size, ikernel, dist_normalization_factor);
                break;
            case 7:
            case 8:
                results = knn_adaptive_mean_shift_data_smoother_with_grad_field_averaging(X, k, density_k, n_steps, step_size, ikernel, dist_normalization_factor, average_direction_only);
                break;
            case 9:
            case 10:
                results = adaptive_mean_shift_data_smoother_with_grad_field_averaging(X, k, density_k, n_steps, step_size, ikernel, dist_normalization_factor, average_direction_only, momentum, increase_factor, decrease_factor);
                break;
            default:
                Rf_error("Invalid method specified");
        }
    } catch (const std::exception& e) {
        Rf_error("Error in mean_shift_smoother: %s", e.what());
    }

    // Check if results is valid
    if (!results) {
        Rf_error("mean_shift_smoother returned null pointer");
    }
    // Check if X_traj is valid
    if (results->X_traj.empty() || results->X_traj[0].empty() || results->X_traj[0][0].empty()) {
        Rf_error("mean_shift_smoother returned invalid X_traj");
    }
    // Convert results to R list
    return create_R_list(results->X_traj, results->median_kdistances);
}
