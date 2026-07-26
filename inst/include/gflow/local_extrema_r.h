#ifndef LOCAL_EXTREMA_R_H_
#define LOCAL_EXTREMA_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_detect_local_extrema(
        SEXP s_adj_list,
        SEXP s_weight_list,
        SEXP s_y,
        SEXP s_max_radius,
        SEXP s_min_neighborhood_size,
        SEXP s_detect_maxima
		);

#ifdef __cplusplus
}
#endif
#endif // LOCAL_EXTREMA_R_H_
