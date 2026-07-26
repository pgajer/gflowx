#ifndef MEAN_SHIFT_SMOOTHER_R_H_
#define MEAN_SHIFT_SMOOTHER_R_H_

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

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

#ifdef __cplusplus
}
#endif

#endif // MEAN_SHIFT_SMOOTHER_R_H_
