#ifndef GRAPH_GRADIENT_FLOW_R_H_
#define GRAPH_GRADIENT_FLOW_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

SEXP S_construct_graph_gradient_flow(
	SEXP s_adj_list,
	SEXP s_weight_list,
	SEXP s_y,
	SEXP s_scale,
	SEXP s_quantile_scale_thld,
	SEXP s_with_trajectories
);

#ifdef __cplusplus
}
#endif
#endif // GRAPH_GRADIENT_FLOW_R_H_
