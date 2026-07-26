#ifndef CENTERED_PATHS_R_H_
#define CENTERED_PATHS_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_get_path_data(
        SEXP adj_list_s,
        SEXP weight_list_s,
        SEXP y_s,
        SEXP ref_vertex_s,
        SEXP bandwidth_s,
        SEXP dist_normalization_factor_s,
        SEXP min_path_size_s,
		SEXP diff_threshold_s,
        SEXP kernel_type_s,
		SEXP verbose_s
        );

	SEXP S_ugg_get_path_data(
        SEXP adj_list_s,
        SEXP weight_list_s,
        SEXP grid_size_s,
        SEXP y_s,
        SEXP ref_vertex_s,
        SEXP bandwidth_s,
        SEXP dist_normalization_factor_s,
        SEXP min_path_size_s,
        SEXP diff_threshold_s,
        SEXP kernel_type_s,
        SEXP verbose_s);

#ifdef __cplusplus
}
#endif
#endif // CENTERED_PATHS_R_H_
