#ifndef MEAN_SHIFT_SMOOTHER_H_
#define MEAN_SHIFT_SMOOTHER_H_

#include <vector>
#include <memory>

struct mean_shift_smoother_results_t {
    std::vector<std::vector<std::vector<double>>> X_traj;
    std::vector<double> median_kdistances;
};

std::unique_ptr<mean_shift_smoother_results_t> adaptive_mean_shift_data_smoother_with_grad_field_averaging(const std::vector<std::vector<double>>& X,
                                                                                                           int k,
                                                                                                           int density_k,
                                                                                                           int n_steps,
                                                                                                           double initial_step_size,
                                                                                                           int ikernel,
                                                                                                           double dist_normalization_factor,
                                                                                                           bool average_direction_only,
                                                                                                           double momentum,
                                                                                                           double increase_factor,
                                                                                                           double decrease_factor);

std::unique_ptr<mean_shift_smoother_results_t> knn_adaptive_mean_shift_data_smoother_with_grad_field_averaging(const std::vector<std::vector<double>>& X,
                                                                                                               int k,
                                                                                                               int density_k,
                                                                                                               int n_steps,
                                                                                                               double step_size,
                                                                                                               int ikernel,
                                                                                                               double dist_normalization_factor,
                                                                                                               bool average_direction_only);

#endif // MEAN_SHIFT_SMOOTHER_H_
