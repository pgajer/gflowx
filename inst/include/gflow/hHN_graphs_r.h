#ifndef HHN_GRAPHS_R_H_
#define HHN_GRAPHS_R_H_

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Temporary internal dependency for private .create.hHN.graph(); not public API. */
SEXP S_create_hHN_graph(SEXP s_adj_list, SEXP s_weight_list, SEXP s_h);

#ifdef __cplusplus
}
#endif

#endif // HHN_GRAPHS_R_H_
