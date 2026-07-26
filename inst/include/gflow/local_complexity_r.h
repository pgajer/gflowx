#ifndef LOCAL_COMPLEXITY_H_
#define LOCAL_COMPLEXITY_H_

#include <Rinternals.h> // Needed for SEXP type definition

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

    SEXP S_estimate_local_complexity(
        SEXP x_r,
        SEXP y_r,
        SEXP center_idx_r,
        SEXP pilot_bandwidth_r,
		SEXP kernel_type_r);

	SEXP S_estimate_binary_local_complexity(
        SEXP x_r,
        SEXP y_r,
        SEXP center_idx_r,
        SEXP pilot_bandwidth_r,
        SEXP kernel_type_r,
        SEXP method_r);

	SEXP S_estimate_ma_binary_local_complexity_quadratic(SEXP x_r,
                                                         SEXP y_r,
                                                         SEXP pilot_bandwidth_r,
                                                         SEXP kernel_type_r);
#ifdef __cplusplus
}
#endif

#endif // LOCAL_COMPLEXITY_H_
