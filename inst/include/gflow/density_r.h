#ifndef DENSITY_R_H_
#define DENSITY_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

SEXP S_estimate_local_density_over_grid(SEXP s_x,
                                        SEXP s_grid_size,
                                        SEXP s_poffset,
                                        SEXP s_pilot_bandwidth,
                                        SEXP s_kernel_type,
                                        SEXP s_verbose);

#ifdef __cplusplus
}
#endif

#endif // DENSITY_R_H_
