#ifndef GFASSOC_R_H
#define GFASSOC_R_H

/**
 * @file gfassoc_r.h
 * @brief R C API interface declarations for gradient flow association functions
 *
 * This header declares the SEXP wrapper functions that bridge R and C++
 * implementations. These functions handle conversion between R data types
 * and C++ structures.
 *
 * All functions follow the naming convention:
 *   S_function_name - SEXP wrapper callable from R via .Call()
 */

#include <R.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Compute basin and cell membership from basin structures
 *
 * @param s_lmax_basins List of max basin structures from compute.basins.of.attraction
 * @param s_lmin_basins List of min basin structures from compute.basins.of.attraction
 * @param s_n_vertices Number of vertices in the graph
 *
 * @return R list containing:
 *   - max_basin_indices: list of integer vectors (0-based indices)
 *   - min_basin_indices: list of integer vectors (0-based indices)
 *   - max_membership: list of numeric vectors (weights)
 *   - min_membership: list of numeric vectors (weights)
 *   - max_vertices: integer vector of extremum vertices
 *   - min_vertices: integer vector of extremum vertices
 *   - max_values: numeric vector of extremum values
 *   - min_values: numeric vector of extremum values
 *   - cell_indices: list of 2-column integer matrices (max_idx, min_idx)
 *   - cell_membership: list of numeric vectors (cell weights)
 */
SEXP S_gfassoc_membership(
    SEXP s_lmax_basins,
    SEXP s_lmin_basins,
    SEXP s_n_vertices
);


/**
 * @brief Compute polarity coordinates for a fitted surface
 *
 * @param s_y Numeric vector of function values
 * @param s_membership Membership structure from S_gfassoc_membership
 * @param s_polarity_scale Character: "value" or "rank"
 * @param s_epsilon Numeric: threshold for flat region detection
 *
 * @return R list containing:
 *   - theta: numeric vector [0, 1]
 *   - polarity: numeric vector [-1, 1]
 *   - range: numeric vector (dynamic range at each vertex)
 *   - is_valid: logical vector
 */
SEXP S_gfassoc_polarity(
    SEXP s_y,
    SEXP s_membership,
    SEXP s_polarity_scale,
    SEXP s_epsilon
);


/**
 * @brief Compute vertex-level and global association from two polarity structures
 *
 * @param s_pol_y Polarity structure for y from S_gfassoc_polarity
 * @param s_pol_z Polarity structure for z from S_gfassoc_polarity
 * @param s_vertex_mass Optional numeric vector of vertex weights
 *
 * @return R list containing:
 *   - vertex:
 *     - a_pol: numeric vector [-1, 1]
 *     - sign_pol: numeric vector {-1, 0, 1}
 *     - confidence: numeric vector [0, 1]
 *     - is_valid: logical vector
 *   - global:
 *     - A_pol: scalar
 *     - kappa_pol: scalar
 *     - n_positive: integer
 *     - n_negative: integer
 *     - n_zero: integer
 *     - n_invalid: integer
 */
SEXP S_gfassoc_association(
    SEXP s_pol_y,
    SEXP s_pol_z,
    SEXP s_vertex_mass
);


/**
 * @brief Compute basin association character
 *
 * @param s_y_membership Membership structure for y
 * @param s_z_membership Membership structure for z
 * @param s_pol_y Polarity structure for y
 * @param s_pol_z Polarity structure for z
 * @param s_vertex_mass Optional numeric vector of vertex weights
 *
 * @return R list containing:
 *   - chi_y_max: numeric vector
 *   - chi_y_min: numeric vector
 *   - chi_z_max: numeric vector
 *   - chi_z_min: numeric vector
 *   - mass_y_max: numeric vector
 *   - mass_y_min: numeric vector
 *   - mass_z_max: numeric vector
 *   - mass_z_min: numeric vector
 */
SEXP S_gfassoc_basin_character(
    SEXP s_y_membership,
    SEXP s_z_membership,
    SEXP s_pol_y,
    SEXP s_pol_z,
    SEXP s_vertex_mass
);


/**
 * @brief Compute soft overlap matrices between two basin structures
 *
 * @param s_y_membership Membership structure for y
 * @param s_z_membership Membership structure for z
 * @param s_vertex_mass Optional numeric vector of vertex weights
 *
 * @return R list containing:
 *   - O_pp: matrix (y-max x z-max)
 *   - O_mm: matrix (y-min x z-min)
 *   - O_pm: matrix (y-max x z-min)
 *   - O_mp: matrix (y-min x z-max)
 *   - total_mass: scalar
 */
SEXP S_gfassoc_overlap(
    SEXP s_y_membership,
    SEXP s_z_membership,
    SEXP s_vertex_mass
);


/**
 * @brief Compute deviation from independence for an overlap matrix
 *
 * @param s_overlap_matrix Numeric matrix (an overlap matrix)
 *
 * @return R list containing:
 *   - delta: matrix (raw deviation)
 *   - zeta: matrix (standardized Pearson residuals)
 *   - expected: matrix (expected under independence)
 */
SEXP S_gfassoc_deviation(
    SEXP s_overlap_matrix
);


/**
 * @brief Comprehensive gradient flow correlation analysis
 *
 * This is the main entry point for symmetric association analysis.
 * Combines membership, polarity, association, and overlap computations.
 *
 * @param s_y_hat Numeric vector of fitted values for y
 * @param s_z_hat Numeric vector of fitted values for z
 * @param s_y_lmax_basins List of y max basin structures
 * @param s_y_lmin_basins List of y min basin structures
 * @param s_z_lmax_basins List of z max basin structures
 * @param s_z_lmin_basins List of z min basin structures
 * @param s_vertex_mass Optional numeric vector of vertex weights
 * @param s_options List of options (polarity_scale, epsilon, etc.)
 *
 * @return R list of class "gfcor" containing all results
 */
SEXP S_gfcor(
    SEXP s_y_hat,
    SEXP s_z_hat,
    SEXP s_y_lmax_basins,
    SEXP s_y_lmin_basins,
    SEXP s_z_lmax_basins,
    SEXP s_z_lmin_basins,
    SEXP s_vertex_mass,
    SEXP s_options
);


#ifdef __cplusplus
}
#endif

#endif // GFASSOC_R_H
