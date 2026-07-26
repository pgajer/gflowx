#ifndef RANDOM_SAMPLING_R_H_
#define RANDOM_SAMPLING_R_H_

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

SEXP S_rlaplace(SEXP R_n, SEXP R_location, SEXP R_scale, SEXP R_seed);

#ifdef __cplusplus
}
#endif

#endif // RANDOM_SAMPLING_R_H_
