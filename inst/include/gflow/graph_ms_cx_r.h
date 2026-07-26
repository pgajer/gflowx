#ifndef GRAPH_MS_CX_R_H_
#define GRAPH_MS_CX_R_H_

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

SEXP S_graph_constrained_gradient_flow_trajectories(SEXP s_graph, SEXP s_core_graph, SEXP s_Ey);
SEXP S_graph_MS_cx_with_path_search(SEXP s_graph, SEXP s_core_graph, SEXP s_Ey);
SEXP S_graph_MS_cx_using_short_h_hops(SEXP s_graph,
                                      SEXP s_hop_list,
                                      SEXP s_core_graph,
                                      SEXP s_Ey);

#ifdef __cplusplus
}
#endif

#endif // GRAPH_MS_CX_R_H_
