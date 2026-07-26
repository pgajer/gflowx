#ifndef GFLOW_CX_R_H_
#define GFLOW_CX_R_H_

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_compute_extrema_hop_nbhds(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y
		);

	SEXP S_create_gflow_cx(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y,
		SEXP s_hop_idx_thld,
		SEXP s_smoother_type,
		SEXP s_max_smoothing_iterations,
		SEXP s_max_inner_iterations,
		SEXP s_smoothing_tolerance,
		SEXP s_sigma,
		SEXP s_process_in_order,
		SEXP s_verbose,
		SEXP s_detailed_recording
		);

	SEXP S_apply_harmonic_extension(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y,
		SEXP s_region_vertices,
		SEXP s_boundary_values,
		SEXP s_smoother_type,
		SEXP s_max_iterations,
		SEXP s_tolerance,
		SEXP s_sigma,
		SEXP s_record_iterations,
		SEXP s_verbose
		);

#ifdef __cplusplus
}
#endif
#endif // GFLOW_CX_R_H_
