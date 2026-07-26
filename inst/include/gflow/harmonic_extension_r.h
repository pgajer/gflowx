#ifndef HARMONIC_EXTENSION_R_H
#define HARMONIC_EXTENSION_R_H

#include <R.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_compute_harmonic_extension(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_trajectory,
		SEXP s_tube_radius,
		SEXP s_tube_type,
		SEXP s_use_edge_weights,
		SEXP s_max_iterations,
		SEXP s_tolerance,
		SEXP s_basin_restriction,
		SEXP s_verbose
		);

	SEXP S_select_max_density_trajectory(
		SEXP s_trajectories,
		SEXP s_density
		);

#ifdef __cplusplus
}
#endif

#endif // HARMONIC_EXTENSION_R_H
