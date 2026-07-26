#ifndef NERVE_CX_R_H_
#define NERVE_CX_R_H_

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_create_nerve_complex(SEXP s_coords, SEXP s_k, SEXP s_max_dim);
	SEXP S_extract_skeleton_graph(SEXP s_complex_ptr);
	SEXP S_set_function_values(SEXP s_complex_ptr, SEXP s_values);
	SEXP S_set_weight_scheme(SEXP s_complex_ptr, SEXP s_weight_type, SEXP s_params);
	SEXP S_solve_full_laplacian(SEXP s_complex_ptr, SEXP s_lambda, SEXP s_dim_weights);
	SEXP S_get_simplex_counts(SEXP s_complex_ptr);
	SEXP S_extract_skeleton_graph(SEXP s_complex_ptr);

	SEXP S_nerve_cx_spectral_filter(
		SEXP s_complex_ptr,
		SEXP s_y,
		SEXP s_laplacian_type,
		SEXP s_filter_type,
		SEXP s_laplacian_power,
		SEXP s_dim_weights,
		SEXP s_kernel_params,
		SEXP s_n_evectors,
		SEXP s_n_candidates,
		SEXP s_log_grid,
		SEXP s_with_t_predictions,
		SEXP s_verbose
		);

#ifdef __cplusplus
}
#endif
#endif // NERVE_CX_R_H_
