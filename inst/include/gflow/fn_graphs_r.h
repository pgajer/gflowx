#ifndef FN_GRAPHS_R_H_
#define FN_GRAPHS_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

SEXP S_construct_function_aware_graph(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_function_values,
		SEXP s_weight_type,
		SEXP s_epsilon,
		SEXP s_lambda,
		SEXP s_alpha,
		SEXP s_beta,
		SEXP s_tau,
		SEXP s_p,
		SEXP s_q,
		SEXP s_r,
		SEXP s_normalize,
		SEXP s_weight_thld
		);

	SEXP S_analyze_function_aware_weights(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_function_values,
		SEXP s_weight_types,
		SEXP s_epsilon,
		SEXP s_lambda,
		SEXP s_alpha,
		SEXP s_beta,
		SEXP s_tau,
		SEXP s_p,
		SEXP s_q,
		SEXP s_r
		);

#ifdef __cplusplus
}
#endif
#endif // FN_GRAPHS_R_H_
