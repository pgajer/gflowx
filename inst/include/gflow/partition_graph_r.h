#ifndef PARTITION_GRAPH_R_H
#define PARTITION_GRAPH_R_H

#include <R.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

SEXP S_partition_graph(SEXP s_adj_list, SEXP s_weight_list,
                       SEXP s_partition, SEXP s_weight_type);

#ifdef __cplusplus
}
#endif

#endif // PARTITION_GRAPH_R_H
