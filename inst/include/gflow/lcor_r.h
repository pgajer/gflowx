#ifndef LCOR_R_H
#define LCOR_R_H

#include <R.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_lcor_instrumented(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y,
		SEXP s_z,
		SEXP s_type,
		SEXP s_y_diff_type,
		SEXP s_z_diff_type,
		SEXP s_epsilon,
		SEXP s_winsorize_quantile
		);

	SEXP S_lcor(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y,
		SEXP s_z,
		SEXP s_type,
		SEXP s_y_diff_type,
		SEXP s_z_diff_type,
		SEXP s_epsilon,
		SEXP s_winsorize_quantile
		);

	SEXP S_lcor_vector_matrix(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_y,
		SEXP s_Z,
		SEXP s_type,
		SEXP s_y_diff_type,
		SEXP s_z_diff_type,
		SEXP s_epsilon,
		SEXP s_winsorize_quantile,
		SEXP s_instrumented
		);

#ifdef __cplusplus
}
#endif

#endif // LCOR_R_H
