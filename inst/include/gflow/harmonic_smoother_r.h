#ifndef HARMONIC_SMOOTHER_R_H_
#define HARMONIC_SMOOTHER_R_H_

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

SEXP S_perform_harmonic_smoothing(
	SEXP s_adj_list,
	SEXP s_weight_list,
	SEXP s_harmonic_predictions,
	SEXP s_region_vertices,
	SEXP s_max_iterations,
	SEXP s_tolerance
	);

SEXP S_harmonic_smoother(
	SEXP s_adj_list,
	SEXP s_weight_list,
	SEXP s_harmonic_predictions,
	SEXP s_region_vertices,
	SEXP s_max_iterations,
	SEXP s_tolerance,
	SEXP s_record_frequency,
	SEXP s_stability_window,
	SEXP s_stability_threshold
	);

#ifdef __cplusplus
}
#endif
#endif // HARMONIC_SMOOTHER_R_H_
