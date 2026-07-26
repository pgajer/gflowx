#ifndef GFASSOC_TYPES_HPP
#define GFASSOC_TYPES_HPP

/**
 * @file gfassoc_types.hpp
 * @brief Core data structures for gradient flow association analysis
 *
 * This file defines the fundamental types used throughout the gradient flow
 * association framework. The framework measures association between two fitted
 * surfaces y and z on a graph by analyzing their gradient flow complexes.
 *
 * The key insight is that gradient flow partitions the graph into basins of
 * attraction, and the relationship between these partitions encodes association
 * structure at multiple scales: global, basin-level, cell-level, and trajectory-level.
 *
 * DISCRETE MORSE THEORY CONTEXT:
 * In the discrete setting, gradient trajectories are not unique. A vertex may
 * belong to multiple basins simultaneously when it has multiple neighbors with
 * lower (or higher) function values. This multiplicity is handled throughout
 * via soft membership vectors that distribute unit mass across containing basins.
 */

#include <vector>
#include <utility>
#include <limits>
#include <cstddef>

#include <Eigen/Core>
#include <Eigen/Dense>

using std::size_t;

// ============================================================================
// Membership Structures
// ============================================================================

/**
 * @struct basin_membership_t
 * @brief Soft membership vectors for basin assignment
 *
 * For each vertex v, stores which basins contain v and the normalized membership
 * weights. When a vertex belongs to k basins, each receives weight 1/k under
 * uniform weighting (alternative weighting schemes are supported).
 *
 * The membership vectors satisfy:
 *   sum_i mu_max[v][i] = 1  for all v with at least one max basin
 *   sum_j mu_min[v][j] = 1  for all v with at least one min basin
 */
struct basin_membership_t {
    size_t n_vertices;
    size_t n_max_basins;   ///< Number of local maxima (descending basins)
    size_t n_min_basins;   ///< Number of local minima (ascending basins)

    /// For each vertex: indices of max basins containing it
    /// max_basin_indices[v] = {i1, i2, ...} where v in B^+_{i1}, B^+_{i2}, ...
    std::vector<std::vector<size_t>> max_basin_indices;

    /// For each vertex: indices of min basins containing it
    std::vector<std::vector<size_t>> min_basin_indices;

    /// Normalized membership weights for max basins
    /// max_membership[v][k] = mu^+_{max_basin_indices[v][k]}(v)
    std::vector<std::vector<double>> max_membership;

    /// Normalized membership weights for min basins
    std::vector<std::vector<double>> min_membership;

    /// Vertex indices of the extrema (for value lookup)
    std::vector<size_t> max_vertices;  ///< max_vertices[i] = vertex index of i-th maximum
    std::vector<size_t> min_vertices;  ///< min_vertices[i] = vertex index of i-th minimum

    /// Function values at extrema
    std::vector<double> max_values;    ///< y(M_i) for each maximum M_i
    std::vector<double> min_values;    ///< y(m_j) for each minimum m_j

    /**
     * @brief Get multiplicity of vertex v in max basins
     */
    size_t max_multiplicity(size_t v) const {
        return max_basin_indices[v].size();
    }

    /**
     * @brief Get multiplicity of vertex v in min basins
     */
    size_t min_multiplicity(size_t v) const {
        return min_basin_indices[v].size();
    }

    /**
     * @brief Check if vertex belongs to exactly one max and one min basin
     */
    bool has_unique_cell(size_t v) const {
        return max_multiplicity(v) == 1 && min_multiplicity(v) == 1;
    }
};


/**
 * @struct cell_membership_t
 * @brief Soft membership for Morse-Smale cells
 *
 * A cell C_{ij} = B^+_i cap B^-_j is the intersection of a descending manifold
 * from maximum M_i and an ascending manifold from minimum m_j. With overlapping
 * basins, a vertex may belong to multiple cells.
 *
 * Cell membership is normalized so that sum_{i,j} gamma_{ij}(v) = 1.
 */
struct cell_membership_t {
    size_t n_vertices;
    size_t n_max_basins;
    size_t n_min_basins;

    /// For each vertex: list of (max_idx, min_idx) pairs for cells containing it
    std::vector<std::vector<std::pair<size_t, size_t>>> cell_indices;

    /// Normalized cell membership weights gamma_{ij}(v)
    /// cell_membership[v][k] corresponds to cell cell_indices[v][k]
    std::vector<std::vector<double>> cell_membership;

    /**
     * @brief Get cell multiplicity of vertex v
     */
    size_t cell_multiplicity(size_t v) const {
        return cell_indices[v].size();
    }

    /**
     * @brief Check if vertex belongs to exactly one cell
     */
    bool has_unique_cell(size_t v) const {
        return cell_multiplicity(v) == 1;
    }

    /**
     * @brief Get the unique cell for a vertex (only valid if has_unique_cell(v))
     */
    std::pair<size_t, size_t> get_unique_cell(size_t v) const {
        return cell_indices[v][0];
    }
};


// ============================================================================
// Overlap Structures
// ============================================================================

/**
 * @struct overlap_matrices_t
 * @brief Soft overlap matrices between basin partitions
 *
 * For two functions y and z with their respective basin structures, the overlap
 * matrices quantify how the partitions intersect. Entry O^{++}_{ij} measures the
 * effective overlap between y-max basin i and z-max basin j:
 *
 *   O^{++}_{ij} = sum_v m_0(v) * mu^{y,+}_i(v) * mu^{z,+}_j(v)
 *
 * where m_0(v) is vertex mass and mu are normalized membership weights.
 *
 * The four matrices capture all combinations:
 *   O_pp: y-max with z-max (positive-positive association regions)
 *   O_mm: y-min with z-min (negative-negative association regions)
 *   O_pm: y-max with z-min (opposite association regions)
 *   O_mp: y-min with z-max (opposite association regions)
 */
struct overlap_matrices_t {
    Eigen::MatrixXd O_pp;  ///< [n_y_max x n_z_max] max-max overlap
    Eigen::MatrixXd O_mm;  ///< [n_y_min x n_z_min] min-min overlap
    Eigen::MatrixXd O_pm;  ///< [n_y_max x n_z_min] max-min overlap
    Eigen::MatrixXd O_mp;  ///< [n_y_min x n_z_max] min-max overlap

    double total_mass;     ///< Sum of vertex masses (for normalization)
};


/**
 * @struct basin_deviation_t
 * @brief Deviation from independence for basin pair overlaps
 *
 * Measures whether specific basin pairs have more or less overlap than expected
 * under independence. The standardized deviation zeta follows the form of Pearson
 * residuals in contingency table analysis.
 */
struct basin_deviation_t {
    Eigen::MatrixXd delta;     ///< O_ij - E_ij (raw deviation)
    Eigen::MatrixXd zeta;      ///< Standardized deviation (Pearson residuals)
    Eigen::MatrixXd expected;  ///< E_ij under independence
};


// ============================================================================
// Polarity Structures
// ============================================================================

/**
 * @struct polarity_result_t
 * @brief Polarity coordinates for a single function
 *
 * The polarity coordinate p(v) in [-1, 1] measures where vertex v sits within
 * its accessible dynamic range, defined by the extrema reachable via gradient flow.
 *
 * For a vertex in cell C_{ij} (between max M_i and min m_j):
 *   theta(v) = (y(v) - y(m_j)) / (y(M_i) - y(m_j))  in [0, 1]
 *   p(v) = 2*theta(v) - 1                            in [-1, 1]
 *
 * With cell multiplicity, theta and p are convex combinations over cells.
 */
struct polarity_result_t {
    std::vector<double> theta;     ///< Normalized height in [0, 1]
    std::vector<double> polarity;  ///< p = 2*theta - 1 in [-1, 1]
    std::vector<double> range;     ///< M(v) - m(v), dynamic range at each vertex

    /// Flags for special cases
    std::vector<bool> is_valid;    ///< False if range is too small (flat region)

    double epsilon;                ///< Threshold for flat region detection
};


// ============================================================================
// Association Structures
// ============================================================================

/**
 * @struct vertex_association_t
 * @brief Vertex-level association scores
 *
 * The polarity-based association score a_pol(v) = p_y(v) * p_z(v) measures
 * whether v is simultaneously high or low in both landscapes (positive) versus
 * high in one and low in the other (negative).
 */
struct vertex_association_t {
    std::vector<double> a_pol;       ///< p_y(v) * p_z(v) in [-1, 1]
    std::vector<double> sign_pol;    ///< sign(a_pol(v)) in {-1, 0, +1}
    std::vector<double> confidence;  ///< |a_pol(v)| or custom confidence measure

    /// Validity inherited from polarity computations
    std::vector<bool> is_valid;      ///< Both y and z polarity valid
};


/**
 * @struct global_association_t
 * @brief Global association summary statistics
 *
 * Aggregates vertex-level association into global summaries using vertex masses.
 */
struct global_association_t {
    /// Mass-weighted polarity concordance
    /// A_pol = sum_v m_0(v) * p_y(v) * p_z(v) / sum_v m_0(v)
    double A_pol;

    /// Sign concordance (discrete summary)
    /// kappa_pol = sum_v m_0(v) * sign(p_y(v) * p_z(v)) / sum_v m_0(v)
    double kappa_pol;

    /// Counts for diagnostic purposes
    size_t n_positive;   ///< Vertices with a_pol > 0
    size_t n_negative;   ///< Vertices with a_pol < 0
    size_t n_zero;       ///< Vertices with a_pol = 0
    size_t n_invalid;    ///< Vertices where polarity undefined

    double total_mass;   ///< Sum of vertex masses used
};


/**
 * @struct basin_character_t
 * @brief Association character for each basin
 *
 * The association character chi^{y,+}_i measures the average z-polarity within
 * y-max basin i:
 *
 *   chi^{y,+}_i = sum_v m_0(v) * mu^{y,+}_i(v) * p_z(v) / sum_v m_0(v) * mu^{y,+}_i(v)
 *
 * Values near +1 indicate the y-high region coincides with z-high (direct association).
 * Values near -1 indicate y-high coincides with z-low (opposite association).
 */
struct basin_character_t {
    std::vector<double> chi_y_max;  ///< Character of each y-max basin
    std::vector<double> chi_y_min;  ///< Character of each y-min basin
    std::vector<double> chi_z_max;  ///< Character of each z-max basin
    std::vector<double> chi_z_min;  ///< Character of each z-min basin

    /// Basin masses (denominators in character computation)
    std::vector<double> mass_y_max;
    std::vector<double> mass_y_min;
    std::vector<double> mass_z_max;
    std::vector<double> mass_z_min;
};


// ============================================================================
// Co-monotonicity Structures (Directed Measures)
// ============================================================================

/**
 * @struct path_comonotonicity_t
 * @brief Co-monotonicity statistics for a single trajectory
 *
 * For an ascending y-trajectory gamma = (v_0, ..., v_L), measures whether z
 * changes monotonically along gamma.
 */
struct path_comonotonicity_t {
    double A_path;         ///< Co-monotonicity score in [-1, 1]
    size_t violations;     ///< Number of sign changes in Delta z along path
    size_t path_length;    ///< Number of vertices L+1
    double total_delta_z;  ///< Sum of |Delta z| along path
};


/**
 * @struct cell_comonotonicity_t
 * @brief Co-monotonicity statistics for a Morse-Smale cell
 *
 * Aggregates co-monotonicity over all edges (or trajectories) within a cell.
 */
struct cell_comonotonicity_t {
    size_t max_basin_idx;     ///< Index of the maximum basin
    size_t min_basin_idx;     ///< Index of the minimum basin

    double A_cell;            ///< Edge-based co-monotonicity within cell
    double A_traj;            ///< Trajectory-bundle median (if computed)

    size_t n_edges;           ///< Number of edges in cell
    size_t n_trajectories;    ///< Number of trajectories (if computed)
    size_t n_vertices;        ///< Number of vertices in cell

    bool trajectories_computed;  ///< Whether trajectory statistics are available
};


/**
 * @struct comonotonicity_result_t
 * @brief Complete co-monotonicity analysis results
 *
 * Contains cell-level and optionally trajectory-level co-monotonicity statistics
 * for all cells in the gradient flow complex.
 */
struct comonotonicity_result_t {
    std::vector<cell_comonotonicity_t> cells;

    /// Global summary: median A_cell across all cells
    double median_A_cell;

    /// Global summary: mass-weighted mean A_cell
    double mean_A_cell;

    /// Trajectory-level summaries (if computed)
    double median_A_traj;
    double mean_violations;

    bool trajectories_computed;
};


// ============================================================================
// Configuration Structures
// ============================================================================

/**
 * @enum polarity_scale_t
 * @brief Scale for polarity computation
 */
enum class polarity_scale_t {
    VALUE,  ///< Use raw function values (Pearson-like)
    RANK    ///< Use ranks within each cell (Spearman-like)
};


/**
 * @enum region_mode_t
 * @brief Region selection for statistics computation
 */
enum class region_mode_t {
    ALL,   ///< Use all vertices
    CORE   ///< Use core vertices (boundary trimmed)
};


/**
 * @enum weight_mode_t
 * @brief Vertex weighting scheme
 */
enum class weight_mode_t {
    UNIFORM,     ///< All vertices weighted equally
    MASS,        ///< Weight by sampling density
    CONFIDENCE   ///< Down-weight ambiguous vertices
};


/**
 * @struct gfassoc_options_t
 * @brief Configuration options for association computation
 */
struct gfassoc_options_t {
    polarity_scale_t polarity_scale = polarity_scale_t::VALUE;
    region_mode_t region_mode = region_mode_t::ALL;
    weight_mode_t weight_mode = weight_mode_t::UNIFORM;

    double epsilon = 1e-10;         ///< Threshold for flat region detection
    size_t boundary_hops = 1;       ///< Hops to trim for CORE region mode

    bool compute_trajectories = true;   ///< Whether to compute trajectory statistics
    bool compute_overlap = true;        ///< Whether to compute overlap matrices
    bool compute_deviation = true;      ///< Whether to compute basin deviations

    bool verbose = false;
};


#endif // GFASSOC_TYPES_HPP
