#ifndef GFASSOC_MEMBERSHIP_HPP
#define GFASSOC_MEMBERSHIP_HPP

/**
 * @file gfassoc_membership.hpp
 * @brief Functions for computing basin and cell membership from gradient basins
 *
 * This file provides functions to convert the raw gradient basin structures
 * (gradient_basin_t from basin computation) into the soft membership vectors
 * needed for association analysis.
 *
 * The key challenge is handling basin multiplicity: in discrete Morse theory,
 * a vertex may belong to multiple basins simultaneously. The membership
 * functions here normalize this multiplicity into probability-like weights.
 */

#include "gfassoc_types.hpp"
#include "gradient_basin.hpp"

#include <vector>
#include <cstddef>

using std::size_t;

/**
 * @brief Compute soft basin membership from gradient basin structures
 *
 * Converts the discrete gradient basin structures into normalized membership
 * vectors. For each vertex v belonging to basins B_1, ..., B_k, the membership
 * weight for basin B_i is 1/k (uniform weighting).
 *
 * @param max_basins Vector of gradient basins for local maxima (descending flow)
 * @param min_basins Vector of gradient basins for local minima (ascending flow)
 * @param n_vertices Total number of vertices in the graph
 *
 * @return basin_membership_t containing normalized membership vectors
 *
 * @note The gradient_basin_t structures are assumed to use 0-based vertex indices.
 */
basin_membership_t compute_basin_membership(
    const std::vector<gradient_basin_t>& max_basins,
    const std::vector<gradient_basin_t>& min_basins,
    size_t n_vertices
);


/**
 * @brief Compute cell membership from basin membership
 *
 * A Morse-Smale cell C_{ij} is the intersection of max basin B^+_i and min
 * basin B^-_j. The cell membership gamma_{ij}(v) is proportional to the product
 * of basin memberships:
 *
 *   gamma_{ij}(v) = c_{ij}(v) / sum_{k,l} c_{kl}(v)
 *
 * where c_{ij}(v) = 1[v in B^+_i] * 1[v in B^-_j].
 *
 * @param membership Basin membership structure (from compute_basin_membership)
 *
 * @return cell_membership_t containing normalized cell membership vectors
 */
cell_membership_t compute_cell_membership(
    const basin_membership_t& membership
);


/**
 * @brief Compute soft overlap matrices between two basin structures
 *
 * For basin structures from functions y and z, computes four overlap matrices:
 *   O^{++}_{ij} = sum_v m_0(v) * mu^{y,+}_i(v) * mu^{z,+}_j(v)  (max-max)
 *   O^{--}_{ij} = sum_v m_0(v) * mu^{y,-}_i(v) * mu^{z,-}_j(v)  (min-min)
 *   O^{+-}_{ik} = sum_v m_0(v) * mu^{y,+}_i(v) * mu^{z,-}_k(v)  (max-min)
 *   O^{-+}_{jl} = sum_v m_0(v) * mu^{y,-}_j(v) * mu^{z,+}_l(v)  (min-max)
 *
 * @param y_membership Basin membership for function y
 * @param z_membership Basin membership for function z
 * @param vertex_mass Vertex weights m_0(v); if empty, uniform weights are used
 *
 * @return overlap_matrices_t containing all four overlap matrices
 */
overlap_matrices_t compute_soft_overlap(
    const basin_membership_t& y_membership,
    const basin_membership_t& z_membership,
    const std::vector<double>& vertex_mass = {}
);


/**
 * @brief Compute deviation from independence for an overlap matrix
 *
 * For an overlap matrix O, computes:
 *   - Expected overlap E_ij = (row_sum_i * col_sum_j) / total
 *   - Raw deviation delta_ij = O_ij - E_ij
 *   - Standardized deviation zeta_ij = delta_ij / sqrt(E_ij * (1-r_i) * (1-c_j))
 *
 * where r_i and c_j are row and column marginal proportions.
 *
 * @param O Overlap matrix (any of the four types)
 *
 * @return basin_deviation_t containing deviation statistics
 */
basin_deviation_t compute_basin_deviation(
    const Eigen::MatrixXd& O
);


/**
 * @brief Get core vertices (interior) of a basin by trimming boundary
 *
 * Removes vertices within the specified number of hops from the basin boundary.
 * Useful for computing statistics that are robust to boundary effects.
 *
 * @param basin The gradient basin structure
 * @param boundary_hops Number of boundary layers to remove (default 1)
 *
 * @return Vector of core vertex indices
 */
std::vector<size_t> get_basin_core_vertices(
    const gradient_basin_t& basin,
    size_t boundary_hops = 1
);


/**
 * @brief Identify boundary vertices of a cell
 *
 * A vertex v in cell C is on the boundary if it has a neighbor not in C.
 * This function identifies such vertices for boundary trimming.
 *
 * @param cell_indices The (max_idx, min_idx) pair identifying the cell
 * @param y_membership Basin membership structure
 * @param adjacency_list Graph adjacency structure
 *
 * @return Vector of boundary vertex indices
 */
std::vector<size_t> get_cell_boundary_vertices(
    const std::pair<size_t, size_t>& cell_indices,
    const cell_membership_t& cells,
    const std::vector<std::vector<size_t>>& adjacency_list
);


#endif // GFASSOC_MEMBERSHIP_HPP
