/**
 * @file gfc_r.h
 * @brief SEXP interface declarations for gradient-flow basin utilities.
 */

#ifndef GFC_R_H
#define GFC_R_H

#include <R.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

    SEXP S_compute_gfc_basins(
        SEXP s_adj_list,
        SEXP s_weight_list,
        SEXP s_y,
        SEXP s_modulation_type,
        SEXP s_density,
        SEXP s_edgelen_bandwidth,
        SEXP s_verbose
    );

    SEXP S_compute_vertex_density(
        SEXP s_adj_list,
        SEXP s_weight_list,
        SEXP s_normalize
    );

#ifdef __cplusplus
}
#endif

#endif // GFC_R_H
