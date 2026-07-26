/**
 * @file madag_r.h
 * @brief R interface declarations for MADAG functions
 */

#ifndef MADAG_R_H
#define MADAG_R_H

#include <R.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Construct MADAG from a source vertex
 *
 * @param s_adj_list Adjacency list (0-based indices)
 * @param s_weight_list Edge weights
 * @param s_y Function values
 * @param s_source_vertex Source vertex (0-based)
 * @param s_params Parameter list
 * @param s_verbose Verbose flag
 * @return R list containing MADAG structure
 */
SEXP S_construct_madag(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_y,
    SEXP s_source_vertex,
    SEXP s_params,
    SEXP s_verbose
);

/**
 * @brief Enumerate trajectories in a specific cell
 *
 * @param s_madag MADAG structure (from S_construct_madag)
 * @param s_y Function values
 * @param s_max_vertex Maximum vertex defining the cell (0-based)
 * @param s_max_trajectories Maximum trajectories to enumerate (0 = unlimited)
 * @return R list containing trajectories
 */
SEXP S_enumerate_cell_trajectories(
    SEXP s_madag,
    SEXP s_y,
    SEXP s_max_vertex,
    SEXP s_max_trajectories
);

/**
 * @brief Sample trajectories from a cell
 *
 * @param s_madag MADAG structure
 * @param s_y Function values
 * @param s_max_vertex Maximum vertex defining the cell (0-based)
 * @param s_n_samples Number of trajectories to sample
 * @param s_seed Random seed
 * @return R list containing sampled trajectories
 */
SEXP S_sample_cell_trajectories(
    SEXP s_madag,
    SEXP s_y,
    SEXP s_max_vertex,
    SEXP s_n_samples,
    SEXP s_seed
);

/**
 * @brief Compute trajectory similarity matrix
 *
 * @param s_trajectories List of trajectories
 * @param s_similarity_type Similarity type ("jaccard", "overlap")
 * @return Numeric matrix of pairwise similarities
 */
SEXP S_trajectory_similarity_matrix(
    SEXP s_trajectories,
    SEXP s_similarity_type
);

/**
 * @brief Identify bottleneck vertices in a cell
 *
 * @param s_madag MADAG structure
 * @param s_max_vertex Maximum vertex defining the cell (0-based)
 * @param s_min_fraction Minimum fraction threshold
 * @return Integer vector of bottleneck vertices (0-based)
 */
SEXP S_identify_bottlenecks(
    SEXP s_madag,
    SEXP s_max_vertex,
    SEXP s_min_fraction
);

#ifdef __cplusplus
}
#endif

#endif // MADAG_R_H
