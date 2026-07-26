#ifndef TEST_MONOTONIC_REACHABILITY_MAP_R_H_
#define TEST_MONOTONIC_REACHABILITY_MAP_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_test_monotonic_reachability_map(
        SEXP s_adj_list,
        SEXP s_weight_list,
        SEXP s_y,
        SEXP s_ref_vertex,
        SEXP s_radius,
        SEXP s_ascending
		);

#ifdef __cplusplus
}
#endif
#endif // TEST_MONOTONIC_REACHABILITY_MAP_R_H_
