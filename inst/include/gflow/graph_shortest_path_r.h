#ifndef GRAPH_SHORTEST_PATH_R_H_
#define GRAPH_SHORTEST_PATH_R_H_

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

SEXP S_shortest_path(SEXP s_graph, SEXP s_edge_lengths, SEXP s_vertices);

#ifdef __cplusplus
}
#endif

#endif // GRAPH_SHORTEST_PATH_R_H_
