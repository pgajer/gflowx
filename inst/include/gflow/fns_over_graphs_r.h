#ifndef FNS_OVER_GRAPHS_R_H_
#define FNS_OVER_GRAPHS_R_H_

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

SEXP S_make_response_locally_non_const(SEXP Rgraph,
                                       SEXP Ry,
                                       SEXP Rweights,
                                       SEXP Rstep_factor,
                                       SEXP Rprec,
                                       SEXP Rn_itrs,
                                       SEXP Rmean_adjust);

SEXP S_prop_nbhrs_with_smaller_y(SEXP Rgraph, SEXP Ry);

#ifdef __cplusplus
}
#endif

#endif // FNS_OVER_GRAPHS_R_H_
