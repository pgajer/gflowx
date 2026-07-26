#ifndef GRAPH_EDIT_DISTANCE_R_H_
#define GRAPH_EDIT_DISTANCE_R_H_

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

SEXP S_graph_edit_distance(SEXP s_graph1_adj_list,
                           SEXP s_graph1_weights_list,
                           SEXP s_graph2_adj_list,
                           SEXP s_graph2_weights_list,
                           SEXP s_edge_cost,
                           SEXP s_weight_cost_factor);

#ifdef __cplusplus
}
#endif

#endif // GRAPH_EDIT_DISTANCE_R_H_
