/**
 * @file madag.hpp
 * @brief Monotonic Ascent DAG (MADAG) for trajectory-based gradient flow analysis
 *
 * This header defines data structures and functions for constructing and analyzing
 * Monotonic Ascent DAGs rooted at local minima. The MADAG captures all monotonically
 * ascending paths from a source minimum, enabling identification of reachable maxima
 * (cells) and enumeration of trajectories within each cell.
 *
 * The trajectory-centric approach treats gradient flow lines as fundamental objects
 * for understanding transitions between low-risk and high-risk states in biological
 * systems, particularly microbiome compositional analysis.
 */

#ifndef MADAG_HPP
#define MADAG_HPP

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <string>
#include <cstddef>
#include <utility>
#include <functional>

using std::size_t;

// Forward declaration
struct set_wgraph_t;

/// Sentinel value for invalid vertex
constexpr size_t MADAG_INVALID_VERTEX = static_cast<size_t>(-1);

// ============================================================================
// Trajectory Data Structure
// ============================================================================

/**
 * @brief A single monotonically ascending trajectory from minimum to maximum
 *
 * Represents a directed path through the graph where the function value
 * strictly increases at each step. The trajectory connects a local minimum
 * (source) to a local maximum (sink).
 */
struct madag_trajectory_t {
    std::vector<size_t> vertices;  ///< Ordered vertices from source to sink (0-based)
    size_t source_min;             ///< Local minimum where trajectory starts
    size_t sink_max;               ///< Local maximum where trajectory ends
    double total_ascent;           ///< y[sink] - y[source]
    int cluster_id;                ///< Assigned cluster (-1 if outlier or unassigned)
    
    madag_trajectory_t()
        : source_min(MADAG_INVALID_VERTEX), sink_max(MADAG_INVALID_VERTEX),
          total_ascent(0.0), cluster_id(-1) {}
    
    /// Number of vertices in the trajectory
    size_t length() const { return vertices.size(); }
    
    /// Number of edges in the trajectory
    size_t n_edges() const { return vertices.empty() ? 0 : vertices.size() - 1; }
};

// ============================================================================
// Cell Data Structure
// ============================================================================

/**
 * @brief A Morse-Smale cell defined by a (minimum, maximum) pair
 *
 * The cell contains all vertices that lie on some monotonically ascending
 * trajectory from the source minimum to the sink maximum. In discrete graphs,
 * unlike smooth manifolds, a vertex may belong to multiple cells.
 */
struct ms_cell_t {
    size_t min_vertex;                   ///< Local minimum (source, 0-based)
    size_t max_vertex;                   ///< Local maximum (sink, 0-based)
    double min_value;                    ///< y value at minimum
    double max_value;                    ///< y value at maximum
    
    std::vector<size_t> support;         ///< Vertices in cell support (excluding endpoints)
    std::vector<size_t> bottlenecks;     ///< High-centrality intermediate vertices
    
    // Trajectory information
    std::vector<madag_trajectory_t> trajectories;  ///< Enumerated trajectories
    size_t n_trajectories;               ///< Total count (may exceed trajectories.size() if sampled)
    bool explicitly_enumerated;          ///< Whether all trajectories are stored
    
    // Clustering results (populated by separate clustering step)
    int n_clusters;
    std::vector<std::vector<int>> cluster_members;  ///< cluster_id -> trajectory indices
    
    ms_cell_t()
        : min_vertex(MADAG_INVALID_VERTEX), max_vertex(MADAG_INVALID_VERTEX),
          min_value(0.0), max_value(0.0), n_trajectories(0),
          explicitly_enumerated(false), n_clusters(0) {}
    
    /// Cell identifier string for display
    std::string id_string() const;
    
    /// Total function change across the cell
    double total_change() const { return max_value - min_value; }
};

// ============================================================================
// MADAG Data Structure
// ============================================================================

/**
 * @brief Monotonic Ascent DAG rooted at a single local minimum
 *
 * The MADAG captures all vertices reachable from a source minimum via
 * monotonically ascending paths, along with the multi-parent structure
 * that enables trajectory enumeration.
 *
 * For each vertex v in the MADAG:
 * - predecessors[v] contains all vertices u such that (u,v) is an edge
 *   in the original graph and y[v] > y[u]
 * - successors[v] contains all vertices w such that (v,w) is an edge
 *   in the original graph and y[w] > y[v]
 *
 * The DAG is acyclic because edges only go in the direction of increasing y.
 */
struct madag_t {
    size_t source_vertex;                ///< Root local minimum (0-based)
    double source_value;                 ///< y value at source
    
    /// Vertices reachable from source via ascending paths (includes source)
    std::vector<size_t> reachable_vertices;
    
    /// Set version for O(1) membership queries
    std::unordered_set<size_t> reachable_set;
    
    /// For each vertex, its monotonic predecessors within the MADAG
    /// predecessors[v] = {u : edge (u,v) exists, y[v] > y[u], u reachable from source}
    std::unordered_map<size_t, std::vector<size_t>> predecessors;
    
    /// For each vertex, its monotonic successors within the MADAG
    /// successors[v] = {w : edge (v,w) exists, y[w] > y[v], w reachable from source}
    std::unordered_map<size_t, std::vector<size_t>> successors;
    
    /// Local maxima reachable from source (sinks of the DAG)
    std::vector<size_t> reachable_maxima;
    
    /// Cells: one for each reachable maximum
    std::vector<ms_cell_t> cells;
    
    /// Map from maximum vertex to cell index
    std::unordered_map<size_t, size_t> max_to_cell_idx;
    
    /// Topological ordering of vertices (source first, sinks last)
    std::vector<size_t> topological_order;
    
    /// Reverse topological ordering (sinks first, source last)
    std::vector<size_t> reverse_topological_order;
    
    /// For each vertex, number of paths from source to that vertex
    /// Used for trajectory counting without full enumeration
    std::unordered_map<size_t, size_t> path_count_from_source;
    
    /// For each vertex, number of paths from that vertex to any sink
    std::unordered_map<size_t, size_t> path_count_to_sinks;
    
    madag_t() : source_vertex(MADAG_INVALID_VERTEX), source_value(0.0) {}
    
    /// Number of vertices in the MADAG
    size_t n_vertices() const { return reachable_vertices.size(); }
    
    /// Number of reachable maxima (cells)
    size_t n_cells() const { return reachable_maxima.size(); }
    
    /// Check if a vertex is in the MADAG
    bool contains(size_t v) const { return reachable_set.count(v) > 0; }
    
    /// Check if a vertex is a sink (local maximum) in the MADAG
    bool is_sink(size_t v) const {
        auto it = successors.find(v);
        return it == successors.end() || it->second.empty();
    }
    
    /// Get cell by maximum vertex (returns nullptr if not found)
    const ms_cell_t* get_cell(size_t max_vertex) const {
        auto it = max_to_cell_idx.find(max_vertex);
        if (it == max_to_cell_idx.end()) return nullptr;
        return &cells[it->second];
    }
    
    ms_cell_t* get_cell_mutable(size_t max_vertex) {
        auto it = max_to_cell_idx.find(max_vertex);
        if (it == max_to_cell_idx.end()) return nullptr;
        return &cells[it->second];
    }
};

// ============================================================================
// MADAG Construction Parameters
// ============================================================================

/**
 * @brief Parameters controlling MADAG construction and trajectory enumeration
 */
struct madag_params_t {
    /// Maximum number of trajectories to enumerate per cell
    /// If exceeded, sampling is used instead
    size_t max_trajectories_per_cell = 10000;
    
    /// Whether to compute path counts (for trajectory estimation)
    bool compute_path_counts = true;
    
    /// Whether to enumerate trajectories during construction
    bool enumerate_trajectories = true;
    
    /// Edge length quantile threshold (edges longer than this are excluded)
    double edge_length_quantile_thld = 1.0;  // 1.0 = no filtering
    
    /// Minimum cell support size to retain
    size_t min_cell_support = 1;
    
    madag_params_t() = default;
};

// ============================================================================
// MADAG Construction Functions
// ============================================================================

/**
 * @brief Construct MADAG rooted at a single source vertex
 *
 * Builds the monotonic ascent DAG by exploring all vertices reachable from
 * the source via monotonically ascending paths. For each reachable vertex,
 * records all monotonic predecessors (enabling trajectory enumeration).
 *
 * The algorithm uses BFS-like traversal but allows revisiting vertices
 * through different predecessors, as the strict monotonicity condition
 * y[v] > y[u] prevents cycles.
 *
 * @param graph The weighted graph structure
 * @param y Function values at each vertex
 * @param source_vertex The local minimum to root the MADAG at (0-based)
 * @param params Construction parameters
 * @param verbose Print progress messages
 * @return The constructed MADAG
 */
madag_t construct_madag(
    const set_wgraph_t& graph,
    const std::vector<double>& y,
    size_t source_vertex,
    const madag_params_t& params = madag_params_t(),
    bool verbose = false
);

/**
 * @brief Compute topological ordering of MADAG vertices
 *
 * Returns vertices in an order such that for every edge (u,v),
 * u appears before v. This is always possible since the MADAG is acyclic.
 *
 * @param madag The MADAG structure
 * @return Vector of vertices in topological order
 */
std::vector<size_t> compute_topological_order(const madag_t& madag);

/**
 * @brief Compute path counts from source to all vertices
 *
 * For each vertex v, computes the number of distinct paths from the
 * source to v. Uses dynamic programming in topological order.
 *
 * @param madag The MADAG structure (modified in place to store counts)
 */
void compute_path_counts_from_source(madag_t& madag);

/**
 * @brief Compute path counts from all vertices to sinks
 *
 * For each vertex v, computes the number of distinct paths from v
 * to any sink (local maximum). Uses dynamic programming in reverse
 * topological order.
 *
 * @param madag The MADAG structure (modified in place to store counts)
 */
void compute_path_counts_to_sinks(madag_t& madag);

// ============================================================================
// Cell Analysis Functions
// ============================================================================

/**
 * @brief Compute the support of a cell
 *
 * The support of cell (m, M) is the set of vertices that lie on some
 * trajectory from m to M. This is computed as the intersection of
 * vertices reachable from m and vertices from which M is reachable.
 *
 * @param madag The MADAG structure
 * @param max_vertex The local maximum defining the cell
 * @return Vector of support vertices (excluding m and M)
 */
std::vector<size_t> compute_cell_support(
    const madag_t& madag,
    size_t max_vertex
);

/**
 * @brief Count trajectories in a cell without full enumeration
 *
 * Uses the path count information to compute the exact number of
 * trajectories from source to the specified maximum.
 *
 * @param madag The MADAG structure (must have path counts computed)
 * @param max_vertex The local maximum defining the cell
 * @return Number of trajectories in the cell
 */
size_t count_cell_trajectories(
    const madag_t& madag,
    size_t max_vertex
);

/**
 * @brief Identify bottleneck vertices in a cell
 *
 * Bottlenecks are vertices through which a large fraction of trajectories
 * must pass. These represent obligate intermediate states in the transition
 * from minimum to maximum.
 *
 * @param madag The MADAG structure
 * @param max_vertex The local maximum defining the cell
 * @param min_fraction Minimum fraction of trajectories passing through
 *                     a vertex for it to be considered a bottleneck
 * @return Vector of bottleneck vertices
 */
std::vector<size_t> identify_bottlenecks(
    const madag_t& madag,
    size_t max_vertex,
    double min_fraction = 0.5
);

// ============================================================================
// Trajectory Enumeration Functions
// ============================================================================

/**
 * @brief Enumerate all trajectories in a cell
 *
 * Performs depth-first enumeration of all paths from the source minimum
 * to the specified maximum. Should only be called when the trajectory
 * count is manageable.
 *
 * @param madag The MADAG structure
 * @param y Function values (for computing trajectory statistics)
 * @param max_vertex The local maximum defining the cell
 * @param max_trajectories Maximum number to enumerate (0 = unlimited)
 * @return Vector of trajectories
 */
std::vector<madag_trajectory_t> enumerate_cell_trajectories(
    const madag_t& madag,
    const std::vector<double>& y,
    size_t max_vertex,
    size_t max_trajectories = 0
);

/**
 * @brief Sample trajectories from a cell
 *
 * When the number of trajectories is too large for full enumeration,
 * this function samples trajectories by random walks from source to sink.
 *
 * @param madag The MADAG structure
 * @param y Function values
 * @param max_vertex The local maximum defining the cell
 * @param n_samples Number of trajectories to sample
 * @param seed Random seed (0 = use current time)
 * @return Vector of sampled trajectories
 */
std::vector<madag_trajectory_t> sample_cell_trajectories(
    const madag_t& madag,
    const std::vector<double>& y,
    size_t max_vertex,
    size_t n_samples,
    unsigned int seed = 0
);

// ============================================================================
// Trajectory Similarity Functions
// ============================================================================

/**
 * @brief Compute Jaccard similarity between two trajectories
 *
 * The Jaccard similarity is |V(γ₁) ∩ V(γ₂)| / |V(γ₁) ∪ V(γ₂)|
 * where V(γ) is the vertex set of trajectory γ.
 *
 * @param traj1 First trajectory
 * @param traj2 Second trajectory
 * @return Jaccard similarity in [0, 1]
 */
double trajectory_jaccard_similarity(
    const madag_trajectory_t& traj1,
    const madag_trajectory_t& traj2
);

/**
 * @brief Compute pairwise trajectory similarity matrix
 *
 * @param trajectories Vector of trajectories
 * @param similarity_type Type of similarity: "jaccard", "overlap"
 * @return Symmetric similarity matrix
 */
std::vector<std::vector<double>> compute_trajectory_similarity_matrix(
    const std::vector<madag_trajectory_t>& trajectories,
    const std::string& similarity_type = "jaccard"
);

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * @brief Print MADAG summary statistics
 */
void print_madag_summary(const madag_t& madag, bool verbose = false);

/**
 * @brief Print cell summary statistics
 */
void print_cell_summary(const ms_cell_t& cell, bool verbose = false);

#endif // MADAG_HPP
