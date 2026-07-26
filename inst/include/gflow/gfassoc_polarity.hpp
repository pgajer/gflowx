#ifndef GFASSOC_POLARITY_HPP
#define GFASSOC_POLARITY_HPP

/**
 * @file gfassoc_polarity.hpp
 * @brief Polarity coordinate computation for gradient flow association
 *
 * The polarity coordinate p(v) in [-1, 1] measures where vertex v sits within
 * its accessible dynamic range. This is the fundamental building block for
 * value-aware association measures that distinguish positive from negative
 * association.
 *
 * MATHEMATICAL DEFINITION:
 *
 * For a vertex v in a unique cell C_{ij} (between max M_i and min m_j):
 *
 *   theta(v) = (y(v) - y(m_j)) / (y(M_i) - y(m_j) + epsilon)
 *
 * This normalized height lies in [0, 1] because gradient flow guarantees
 * y(m_j) <= y(v) <= y(M_i) for any vertex reachable from both endpoints.
 *
 * The polarity coordinate rescales to [-1, 1]:
 *
 *   p(v) = 2 * theta(v) - 1
 *
 * So p(v) = +1 at maxima, p(v) = -1 at minima, and p(v) = 0 at midpoints.
 *
 * MULTIPLICITY HANDLING:
 *
 * When v belongs to multiple cells {C_{i1,j1}, ..., C_{ik,jk}}, the polarity
 * is a weighted average:
 *
 *   theta(v) = sum_{a} gamma_{ia,ja}(v) * theta^{ia,ja}(v)
 *
 * where gamma are the normalized cell membership weights.
 *
 * GEOMETRIC INTERPRETATION:
 *
 * The polarity coordinate measures position along the gradient flow from
 * minimum to maximum. Two vertices with the same polarity occupy similar
 * relative positions in their respective gradient flow paths, even if their
 * absolute function values differ.
 */

#include "gfassoc_types.hpp"
#include "gfassoc_membership.hpp"

#include <vector>
#include <cstddef>

using std::size_t;


/**
 * @brief Compute polarity coordinates for a function given basin structure
 *
 * This is the main polarity computation function. For each vertex, it computes:
 *   - theta(v): normalized height in [0, 1]
 *   - p(v): polarity in [-1, 1]
 *   - range(v): dynamic range M(v) - m(v)
 *
 * Vertices in flat regions (range < epsilon) are flagged as invalid.
 *
 * @param y Function values at vertices
 * @param membership Basin membership structure for y
 * @param cells Cell membership structure for y
 * @param epsilon Threshold for flat region detection (default 1e-10)
 *
 * @return polarity_result_t containing polarity coordinates and validity flags
 */
polarity_result_t compute_polarity(
    const std::vector<double>& y,
    const basin_membership_t& membership,
    const cell_membership_t& cells,
    double epsilon = 1e-10
);


/**
 * @brief Compute polarity using rank transformation within cells
 *
 * Instead of using raw function values, this version ranks vertices within
 * each cell and uses the ranks for polarity computation. This provides
 * invariance under monotone transformations of y.
 *
 * For a vertex v in cell C with n_C vertices:
 *   theta_rank(v) = rank(y(v) among cell vertices) / (n_C - 1)
 *
 * The rank transformation is applied before averaging over cells with
 * multiplicity.
 *
 * @param y Function values at vertices
 * @param membership Basin membership structure for y
 * @param cells Cell membership structure for y
 * @param epsilon Threshold for minimum cell size (cells with < 2 vertices invalid)
 *
 * @return polarity_result_t using rank-based theta values
 */
polarity_result_t compute_polarity_rank(
    const std::vector<double>& y,
    const basin_membership_t& membership,
    const cell_membership_t& cells,
    double epsilon = 1e-10
);


/**
 * @brief Compute vertex-level association from two polarity structures
 *
 * The polarity-based association score is:
 *   a_pol(v) = p_y(v) * p_z(v)
 *
 * This is positive when v is high (or low) in both landscapes, and negative
 * when v is high in one and low in the other.
 *
 * @param pol_y Polarity structure for function y
 * @param pol_z Polarity structure for function z
 *
 * @return vertex_association_t containing association scores
 */
vertex_association_t compute_vertex_association(
    const polarity_result_t& pol_y,
    const polarity_result_t& pol_z
);


/**
 * @brief Compute global association statistics from vertex associations
 *
 * Aggregates vertex-level association into global summaries:
 *   - A_pol: mass-weighted mean polarity product
 *   - kappa_pol: mass-weighted sign concordance
 *
 * @param va Vertex association structure
 * @param vertex_mass Vertex weights; if empty, uniform weights are used
 *
 * @return global_association_t containing global summaries
 */
global_association_t compute_global_association(
    const vertex_association_t& va,
    const std::vector<double>& vertex_mass = {}
);


/**
 * @brief Compute basin association character
 *
 * The association character of a y-basin measures the average z-polarity
 * of vertices in that basin:
 *
 *   chi^{y,+}_i = sum_v m_0(v) * mu^{y,+}_i(v) * p_z(v)
 *                 / sum_v m_0(v) * mu^{y,+}_i(v)
 *
 * This is computed for all four basin types: y-max, y-min, z-max, z-min.
 *
 * @param y_membership Basin membership for function y
 * @param z_membership Basin membership for function z
 * @param pol_y Polarity for function y
 * @param pol_z Polarity for function z
 * @param vertex_mass Vertex weights; if empty, uniform weights are used
 *
 * @return basin_character_t containing characters for all basins
 */
basin_character_t compute_basin_character(
    const basin_membership_t& y_membership,
    const basin_membership_t& z_membership,
    const polarity_result_t& pol_y,
    const polarity_result_t& pol_z,
    const std::vector<double>& vertex_mass = {}
);


/**
 * @brief Compute polarity for a single cell (helper function)
 *
 * Computes theta^{ij}(v) for a specific cell (i, j):
 *   theta^{ij}(v) = (y(v) - y(m_j)) / (y(M_i) - y(m_j) + epsilon)
 *
 * @param y_v Function value at vertex v
 * @param y_max Function value at maximum M_i
 * @param y_min Function value at minimum m_j
 * @param epsilon Division safety threshold
 *
 * @return Normalized height in [0, 1], or -1 if invalid (flat cell)
 */
double compute_cell_theta(
    double y_v,
    double y_max,
    double y_min,
    double epsilon = 1e-10
);


/**
 * @brief Get vertices belonging to a specific cell
 *
 * Returns all vertices that have nonzero membership in the specified cell.
 *
 * @param cells Cell membership structure
 * @param max_idx Index of the maximum basin
 * @param min_idx Index of the minimum basin
 *
 * @return Vector of vertex indices in the cell
 */
std::vector<size_t> get_cell_vertices(
    const cell_membership_t& cells,
    size_t max_idx,
    size_t min_idx
);


/**
 * @brief Compute weighted average range for a vertex
 *
 * For vertices with cell multiplicity, computes the weighted average
 * of the dynamic ranges across all containing cells.
 *
 * @param v Vertex index
 * @param membership Basin membership structure
 * @param cells Cell membership structure
 *
 * @return Weighted average dynamic range
 */
double compute_weighted_range(
    size_t v,
    const basin_membership_t& membership,
    const cell_membership_t& cells
);


#endif // GFASSOC_POLARITY_HPP
