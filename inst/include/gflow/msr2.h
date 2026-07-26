#ifndef MSC2_H_
#define MSC2_H_

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h>

/*
 * Legacy compatibility header.
 *
 * This header is kept for backward compatibility with older includes.
 * New code should include scoped module headers directly.
 */

#ifdef __cplusplus
extern "C" {
#endif

void print_2d_double_array_first_n(const double *x, int nr, int nc, int N, int precision);
void print_int_array(const int *x, int n);
void print_double_array(const double *x, int n);
void print_double_array_with_precision(const double *x, int n, int precision);

#ifdef __cplusplus
}
#endif

#include "gflow_macros.h"
#include "grids.h"
#include "mstree.h"
#include "lm.h"
#include "stats_utils.h"

#include "angular_wasserstein_index_r.h"
#include "cpp_mstrees_r.h"
#include "mstree_total_length_r.h"
#include "hHN_graphs_r.h"
#include "graph_shortest_path_r.h"
#include "graph_core_endpoints_r.h"
#include "graph_utils_r.h"
#include "graph_ms_cx_r.h"
#include "fns_over_graphs_r.h"
#include "graph_edit_distance_r.h"
#include "graph_spectrum_r.h"
#include "pruning_long_edges_r.h"
#include "graph_cycles_r.h"
#include "cpp_stats_utils_r.h"
#include "random_sampling_r.h"
#include "graph_conn_components_r.h"

#endif // MSC2_H_
