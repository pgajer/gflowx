#ifndef PARAMETERIZE_CIRCULAR_GRAPH_R_H_
#define PARAMETERIZE_CIRCULAR_GRAPH_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_parameterize_circular_graph(
		SEXP s_adj_list,
		SEXP s_weight_list,
		SEXP s_use_edge_lengths
		);

#ifdef __cplusplus
}
#endif
#endif // PARAMETERIZE_CIRCULAR_GRAPH_R_H_
