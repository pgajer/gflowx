#ifndef GFLOW_BASINS_R_H_
#define GFLOW_BASINS_R_H_

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_find_gflow_basins(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y,
		SEXP s_min_basin_size,
		SEXP s_min_path_size,
		SEXP s_q_edge_thld
		);

	SEXP S_find_local_extrema(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y,
		SEXP s_min_basin_size
		);

	SEXP S_create_basin_cx(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y
		);

	SEXP S_perform_harmonic_smoothing(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_harmonic_predictions,
		SEXP s_region_vertices,
		SEXP s_max_iterations,
		SEXP s_tolerance
		);

#ifdef __cplusplus
}
#endif
#endif // GFLOW_BASINS_R_H_
