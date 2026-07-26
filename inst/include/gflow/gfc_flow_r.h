/**
 * @file gfc_flow_r.h
 * @brief R interface declarations for trajectory-based GFC computation
 *
 * Declares SEXP wrapper functions for compute_gfc_flow() that can be
 * called from R via .Call(). These functions are registered in init.c.
 */

#ifndef GFC_FLOW_R_H
#define GFC_FLOW_R_H

#include <R.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief R interface for compute_gfc_flow()
 *
 * @param s_adj_list Adjacency list (list of integer vectors, 0-based from R after conversion)
 * @param s_weight_list Weight list (list of numeric vectors)
 * @param s_y Numeric vector of function values
 * @param s_density Numeric vector of density values (can be NULL/empty)
 * @param s_params Named list of parameters:
 *   - edge_length_quantile_thld: numeric
 *   - apply_relvalue_filter: logical
 *   - min_rel_value_max: numeric
 *   - max_rel_value_min: numeric
 *   - apply_maxima_clustering: logical
 *   - apply_minima_clustering: logical
 *   - max_overlap_threshold: numeric
 *   - min_overlap_threshold: numeric
 *   - apply_geometric_filter: logical
 *   - p_mean_nbrs_dist_threshold: numeric
 *   - p_mean_hopk_dist_threshold: numeric
 *   - p_deg_threshold: numeric
 *   - min_basin_size: integer
 *   - hop_k: integer
 *   - modulation: character ("NONE", "DENSITY", "EDGELEN", "DENSITY_EDGELEN")
 *   - store_trajectories: logical
 *   - max_trajectory_length: integer
 * @param s_verbose Logical scalar
 * @return Named list containing GFC flow results
 */
SEXP S_compute_gfc_flow(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_density,
    SEXP s_params,
    SEXP s_verbose
);

/**
 * @brief R interface for compute_gfc_flow_matrix()
 *
 * @param s_adj_list Adjacency list
 * @param s_weight_list Weight list
 * @param s_Y Numeric matrix (n_vertices x n_functions)
 * @param s_density Numeric vector of density values
 * @param s_params Named list of parameters
 * @param s_n_cores Integer number of threads
 * @param s_verbose Logical scalar
 * @return List of GFC flow results, one per column
 */
SEXP S_compute_gfc_flow_matrix(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_Y,
    SEXP s_density,
    SEXP s_params,
    SEXP s_n_cores,
    SEXP s_verbose
);

#ifdef __cplusplus
}
#endif

#endif // GFC_FLOW_R_H
