#include "riem_dcx.hpp"

#include <unordered_set>
#include <numeric>    // For std::accumulate
#include <queue>
#include <limits>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


/**
 * @brief Extend the complex by one dimension, preserving all existing structure
 *
 * @details
 * Appends a new chain degree to the discrete chain complex, incrementing pmax
 * by one and initializing all data structures for the new dimension. This
 * operation preserves all existing metric structures, boundary operators, and
 * Laplacians at lower dimensions. The function is designed for incremental
 * complex construction where simplices are built dimension by dimension.
 *
 * MOTIVATION:
 * During nerve complex construction from k-NN neighborhoods, we discover
 * simplices dimension by dimension: first vertices, then edges (by checking
 * pairwise neighborhood intersections), then triangles (by checking triple
 * intersections), and so forth. Each dimension depends on having the previous
 * dimension fully constructed to check face existence. This function enables
 * that incremental workflow without destroying lower-dimensional data.
 *
 * OPERATION:
 * Let p_new = pmax + 1 denote the new dimension being added. The function:
 *
 * 1. Increments pmax to p_new
 * 2. Resizes all dimension-indexed containers to accommodate p_new:
 *    - g.M (metric diagonal and matrices)
 *    - g.M_solver (Cholesky factorizations, if maintained)
 *    - L.B (boundary operators)
 *    - L.L (Hodge Laplacians)
 *
 * 3. Initializes the new dimension p_new with default values:
 *    - g.M[p_new]: Sparse diagonal matrix with ones on diagonal
 *    - g.M_solver[p_new]: Null pointer (no factorization yet)
 *    - L.B[p_new]: Empty boundary operator matrix (shape: n_{p-1} × n_new)
 *    - L.L[p_new]: Empty Laplacian (to be assembled after boundary operator is built)
 *
 * USAGE PATTERN:
 * The typical usage in complex construction is:
 *
 * @code
 * // Start with vertices only (dimension 0)
 * riem_dcx_t dcx;
 * dcx.pmax = 0;
 * dcx.vertex_cofaces.resize(n_vertices);
 * // ... populate vertex_cofaces with self-loops and edges ...
 *
 * // Add edge dimension
 * dcx.extend_by_one_dim(n_edges);  // Now pmax = 1
 * dcx.edge_registry.resize(n_edges);
 * // ... populate edge_registry ...
 * dcx.build_boundary_operator_from_edges();
 *
 * // Add triangle dimension
 * dcx.extend_by_one_dim(n_triangles);  // Now pmax = 2
 * dcx.edge_cofaces.resize(n_edges);
 * // ... populate edge_cofaces with triangles ...
 * dcx.build_boundary_operator_from_triangles();
 *
 * // Compute all Laplacians
 * dcx.assemble_operators();
 * @endcode
 *
 * METRIC INITIALIZATION:
 * The new dimension starts with identity metric (diagonal ones). This is a
 * safe default that will be overwritten by the caller based on geometric
 * measurements:
 *   - For vertices: M₀ = diag(ρ₀) set by initialize_metric_from_density()
 *   - For edges: M₁ computed by compute_edge_mass_matrix()
 *   - For triangles: M₂ would follow similar pattern if implemented
 *
 * BOUNDARY OPERATOR:
 * The boundary map B[p_new]: C_{p_new} → C_{p_new-1} is allocated with the
 * correct dimensions but contains no entries initially. The caller must
 * populate this by calling the appropriate boundary construction function:
 *   - B[1]: build_boundary_operator_from_edges()
 *   - B[2]: build_boundary_operator_from_triangles()
 *
 * The boundary operator is built from coface structures (vertex_cofaces,
 * edge_cofaces) rather than simplex tables, making this function independent
 * of the deprecated S and stars structures.
 *
 * DIMENSION COUNTING:
 * Unlike the old version that used S[p-1].size() to determine the number of
 * (p-1)-simplices, this version determines dimensions from coface structures:
 *   - n_0 (vertices) = vertex_cofaces.size()
 *   - n_1 (edges) = edge_registry.size()
 *   - n_2 (triangles) = maximum triangle index in edge_cofaces + 1
 *
 * LAPLACIAN ASSEMBLY:
 * The Hodge Laplacian L[p_new] is left empty and must be computed via
 * assemble_operators() after:
 *   1. All boundary operators B[p_new+1] and B[p_new] are populated
 *   2. The metric g.M[p_new] reflects the actual geometry
 *
 * The Laplacian at dimension p is given by the Hodge decomposition:
 *     L[p] = B[p+1]^T M[p+1]^{-1} B[p+1] + M[p]^{-1} B[p] M[p-1] B[p]^T
 *
 * INVARIANTS MAINTAINED:
 * After calling extend_by_one_dim(n_new):
 *   - pmax has increased by 1
 *   - g.M.size() == L.B.size() == L.L.size() == pmax + 1
 *   - g.M[pmax] is diagonal with n_new ones
 *   - L.B[pmax] is empty sparse matrix of shape (n_{pmax-1}, n_new)
 *   - All structures at dimensions 0, ..., pmax-1 are unchanged
 *
 * EFFICIENCY:
 * This operation is O(n_new) for initializing the new dimension's data
 * structures, plus O(pmax) for resizing dimension-indexed containers.
 * Critically, it does NOT rebuild or copy data from lower dimensions,
 * making it suitable for use during complex construction.
 *
 * THREAD SAFETY:
 * This function modifies the riem_dcx_t object and is NOT thread-safe. Do not
 * call concurrently from multiple threads on the same object.
 *
 * @param n_new Number of simplices at the new dimension p_new = pmax + 1.
 *              Must be non-negative. If zero, the dimension is added but
 *              contains no simplices (useful for maintaining consistent
 *              dimension indexing).
 *
 * @post pmax is incremented by 1
 * @post All dimension-indexed containers have size pmax + 1
 * @post g.M[pmax] is diagonal identity matrix of size n_new × n_new
 * @post L.B[pmax] is empty sparse matrix of appropriate shape
 * @post L.L[pmax] is empty (to be assembled later)
 * @post All data at dimensions 0, ..., pmax-1 (pre-increment) are unchanged
 *
 * @note The caller must subsequently:
 *       1. Populate appropriate coface structures (edge_registry, edge_cofaces, etc.)
 *       2. Set g.M[pmax] entries based on geometric measurements (via metric functions)
 *       3. Build L.B[pmax] by calling boundary construction functions
 *       4. Call assemble_operators() to compute L.L[pmax]
 *
 * @note If n_new = 0, the dimension exists but is empty. This is valid when
 *       no simplices of a given dimension exist (e.g., no triangles in a
 *       sparse graph).
 *
 * @see build_boundary_operator_from_edges() - Builds B[1] from edge_registry
 * @see build_boundary_operator_from_triangles() - Builds B[2] from edge_cofaces
 * @see assemble_operators() - Computes Laplacians after boundary maps are built
 * @see initialize_metric_from_density() - Sets metric from densities in cofaces
 */
void riem_dcx_t::extend_by_one_dim(index_t n_new) {
	const int p_new = pmax + 1;

	// Increment maximum dimension
	pmax = p_new;

	// Resize metric family
	g.M.resize(pmax + 1);
	g.M_solver.resize(pmax + 1);

	// Resize operator family
	L.B.resize(pmax + 1);
	L.L.resize(pmax + 1);

	// Initialize metric at new dimension (identity)
	g.M[p_new] = spmat_t(static_cast<Eigen::Index>(n_new),
						 static_cast<Eigen::Index>(n_new));
	g.M[p_new].reserve(Eigen::VectorXi::Constant(static_cast<int>(n_new), 1));
	for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(n_new); ++i) {
		g.M[p_new].insert(i, i) = 1.0;
	}
	g.M[p_new].makeCompressed();
	g.M_solver[p_new].reset();

	// Initialize boundary operator (empty matrix with correct shape)
	// Determine number of (p-1)-simplices from coface structures
	Eigen::Index n_prev = 0;
	if (p_new == 1) {
		// Boundary of edges: n_0 = number of vertices
		n_prev = static_cast<Eigen::Index>(vertex_cofaces.size());
	} else if (p_new == 2) {
		// Boundary of triangles: n_1 = number of edges
		n_prev = static_cast<Eigen::Index>(edge_registry.size());
	} else if (p_new >= 3) {
		// For higher dimensions, would need additional coface structures
		// For now, this is not implemented as we only handle up to dimension 2
		n_prev = 0;
	}

	const Eigen::Index cols = static_cast<Eigen::Index>(n_new);
	L.B[p_new] = spmat_t(n_prev, cols);

	// Initialize Laplacian (empty, to be assembled later)
	L.L[p_new] = spmat_t();
}

/**
 * @brief Compute number of connected components in the graph
 *
 * This function performs a depth-first search traversal from each unvisited vertex
 * to identify connected components in the graph represented by vertex_cofaces.
 * The graph is defined by the k-nearest neighbor structure, where each vertex i
 * is connected to vertices j appearing in vertex_cofaces[i] (excluding the self-loop
 * at position 0).
 *
 * Two vertices belong to the same connected component if and only if there exists
 * a path of edges connecting them. For the regression framework to function properly,
 * the graph must be connected (have exactly one component), ensuring that geometric
 * information can propagate across the entire domain.
 *
 * The algorithm maintains component labels for each vertex. Starting from each
 * unlabeled vertex, it conducts a depth-first search to label all reachable vertices
 * with the same component identifier. The total number of distinct components equals
 * the number of such searches required to label all vertices.
 *
 * @return Number of connected components in the graph (1 indicates a connected graph)
 *
 * @note This function should be called after vertex_cofaces is populated but before
 *       any geometric computations, to ensure the graph structure is valid.
 *
 * @note Time complexity: O(n + m) where n is the number of vertices and m is the
 *       total number of edges. Space complexity: O(n) for component labels and stack.
 */
int riem_dcx_t::compute_connected_components() {
    const size_t n_vertices = vertex_cofaces.size();

    // Validate that vertex_cofaces is populated
    if (n_vertices == 0) {
        return 0;
    }

    // Component assignment for each vertex (-1 = unvisited)
    std::vector<int> component_id(n_vertices, -1);

    // Stack for depth-first search traversal
    std::vector<index_t> stack;
    stack.reserve(n_vertices);  // Avoid reallocations during traversal

    int num_components = 0;

    // Process each vertex as a potential seed for a new component
    for (size_t seed = 0; seed < n_vertices; ++seed) {
        // Skip vertices already assigned to a component
        if (component_id[seed] >= 0) {
            continue;
        }

        // Start new component from this seed
        stack.clear();
        stack.push_back(seed);
        component_id[seed] = num_components;

        // Depth-first search to find all vertices in this component
        while (!stack.empty()) {
            index_t v = stack.back();
            stack.pop_back();

            // Examine all neighbors of v
            // Note: vertex_cofaces[v][0] is the self-loop, skip it
            for (size_t k = 1; k < vertex_cofaces[v].size(); ++k) {
                index_t neighbor = vertex_cofaces[v][k].vertex_index;

                // If neighbor hasn't been visited, add to current component
                if (component_id[neighbor] < 0) {
                    component_id[neighbor] = num_components;
                    stack.push_back(neighbor);
                }
            }
        }

        // Finished exploring this component
        ++num_components;
    }

    return num_components;
}

/**
 * @brief Assemble all Laplacian operators from current metric
 *
 * Updates all Hodge Laplacians based on the current state of the metric
 * and boundary maps. Should be called after any change to the metric.
 */
void riem_dcx_t::assemble_operators() {

    // Populate edge conductances c1 from edge masses m1.
    // Convention:
    //   vertex_cofaces[i][k].density is the edge mass m1(e) (a measure).
    //   The normalized graph Laplacian uses conductances c_e = 1/m1(e).
    // This makes L0_sym consistent with the diagonal-metric Hodge 0-Laplacian:
    //   L0 = B1 * diag(m1)^{-1} * B1^T * M0.

    L.c1.resize(edge_registry.size());
    for (size_t i = 0; i < vertex_cofaces.size(); ++i) {
        for (size_t k = 1; k < vertex_cofaces[i].size(); ++k) {
            index_t e = vertex_cofaces[i][k].simplex_index;
            double mass = std::max(vertex_cofaces[i][k].density, 1e-10);
            L.c1[e] = 1.0 / mass;
        }
    }

    L.assemble_all(g);
}

/**
 * @brief Compute effective degrees for all vertices based on edge densities
 *
 * The effective degree of a vertex quantifies its connectivity strength in the
 * weighted graph structure. Unlike combinatorial degree (simple neighbor count),
 * effective degree accounts for edge weights through density masses, providing
 * a geometric measure of how well-connected a vertex is to its local neighborhood.
 *
 * MOTIVATION AND INTERPRETATION
 * ==============================
 *
 * In discrete graphs derived from continuous data, not all edges carry equal
 * geometric significance. Edges connecting vertices with large neighborhood
 * intersections represent strong geometric connections, while edges with small
 * intersections may be spurious artifacts of the k-NN construction. The effective
 * degree distinguishes these cases by weighting each edge by its density mass.
 *
 * For vertex i with incident edges e_1, e_2, ..., e_k, the effective degree is:
 *
 *   eff_deg(i) = Σⱼ ρ₁([i,j])
 *
 * where the sum runs over all neighbors j of vertex i, and ρ₁([i,j]) is the
 * edge density computed as the total vertex density in the pairwise neighborhood
 * intersection N̂_k(xᵢ) ∩ N̂_k(xⱼ).
 *
 * Vertices with high effective degree are well-integrated into the graph structure,
 * having many neighbors or neighbors connected through high-density edges. Vertices
 * with low effective degree are poorly connected, either having few neighbors or
 * neighbors connected through sparse intersections. The latter case identifies
 * outlier vertices that are weakly attached to the main data structure.
 *
 * RELATIONSHIP TO GRAPH ISOLATION
 * ================================
 *
 * Effective degree provides a continuous measure of vertex isolation. Classical
 * graph theory uses binary notions: a vertex is either isolated (degree 0) or
 * connected (degree > 0). In weighted graphs derived from statistical data,
 * isolation becomes a matter of degree rather than kind.
 *
 * A vertex with combinatorial degree k but very small effective degree has
 * k geometric connections, but all are weak. Such a vertex is effectively
 * isolated even though technically connected. This occurs for outlier data
 * points far from the main data cloud: their k nearest neighbors may be quite
 * distant, producing small neighborhood intersections and thus low edge densities.
 *
 * The effective degree thus serves as a natural filter for outlier detection
 * in the context of extrema identification. Spurious local extrema often appear
 * at isolated vertices. By filtering vertices with effective degree below some
 * threshold (e.g., 10th percentile), we remove problematic outliers before
 * applying extremum probability scoring.
 *
 * COMPUTATIONAL DETAILS
 * =====================
 *
 * This function computes effective degrees using the *initial* edge densities
 * established during graph initialization, before any iterative refinement.
 * This choice provides several advantages:
 *
 * 1. Stability: Initial densities capture the intrinsic geometry of the data
 *    distribution without response-dependent modulation. Effective degree
 *    becomes a fixed graph property independent of iteration history.
 *
 * 2. Interpretability: The effective degree represents pure geometric connectivity
 *    based on neighborhood overlap, uncontaminated by response coherence effects
 *    that modify edge densities during iteration.
 *
 * 3. Consistency: Initial edge densities are computed via the pairwise intersection
 *    formula ρ₁([i,j]) = Σᵥ∈N̂_k(xᵢ)∩N̂_k(xⱼ) ρ₀(v), matching the construction used
 *    in compute_initial_densities(). Using these values ensures effective degrees
 *    reflect the same geometric quantities used throughout initialization.
 *
 * The computation proceeds by iterating over all vertices and summing the density
 * values of their incident edges. Edge densities are stored in vertex_cofaces[i][k]
 * for k > 0, where vertex_cofaces[i][0] always represents the vertex itself.
 * Since edges appear twice (once in each endpoint's vertex_cofaces), we must
 * avoid double-counting by only accessing each edge once per vertex.
 *
 * COMPLEXITY ANALYSIS
 * ===================
 *
 * Time complexity: O(Σᵢ deg(i)) = O(n·k) where n is the number of vertices
 * and k is the average degree. For k-NN graphs, this is optimal since we must
 * examine every edge incident to every vertex.
 *
 * Space complexity: O(n) for the output vector. No additional space is required
 * beyond the input data structures.
 *
 * The function performs a single pass over vertex_cofaces, making it suitable
 * for large graphs with millions of vertices and edges.
 *
 * USAGE IN EXTREMA DETECTION PIPELINE
 * ====================================
 *
 * Effective degrees support the statistical extrema detection framework:
 *
 * 1. After graph initialization, call compute_effective_degrees() once to
 *    obtain the eff_deg vector.
 *
 * 2. Store eff_deg for use in post-processing extrema detection.
 *
 * 3. When identifying extrema via maxp/minp scoring, apply a pre-filter:
 *      candidates = vertices with maxp > 0.95 AND eff_deg > quantile(eff_deg, 0.10)
 *
 * 4. Compute hop-extremp radius only for candidates passing both filters,
 *    reducing computational cost by eliminating poorly-connected vertices early.
 *
 * This two-stage filtering (probabilistic scoring + connectivity) efficiently
 * identifies genuine extrema while discarding spurious peaks at outlier vertices.
 *
 * @return Vector of length n_vertices where entry i gives the effective degree
 *         of vertex i. All values are non-negative. Sum of all effective degrees
 *         equals twice the sum of all edge densities (since each edge contributes
 *         to two vertices).
 *
 * @pre vertex_cofaces is populated with graph topology
 * @pre vertex_cofaces[i][k].density contains edge densities for k > 0
 * @pre Edge densities have been initialized via compute_initial_densities()
 *
 * @post Returned vector has length equal to number of vertices
 * @post All returned values are non-negative
 * @post Vertices with no neighbors have effective degree 0
 *
 * @note This function uses initial edge densities, not evolved densities from
 *       iterative refinement. This ensures effective degree remains a stable
 *       graph property independent of the regression iteration history.
 *
 * @note The effective degree differs from the weighted degree in classical graph
 *       theory. Weighted degree typically uses edge weights representing costs or
 *       capacities. Effective degree uses edge densities representing geometric
 *       significance of connections in the underlying data distribution.
 *
 * @see compute_initial_densities() for edge density initialization
 * @see update_edge_densities_from_vertices() for evolved edge densities (not used here)
 */
vec_t riem_dcx_t::compute_effective_degrees() const {
    const size_t n_vertices = vertex_cofaces.size();
    vec_t eff_deg = vec_t::Zero(n_vertices);

    for (size_t i = 0; i < n_vertices; ++i) {
        double total_edge_mass = 0.0;

        // Sum density of all incident edges
        // vertex_cofaces[i][0] is the vertex itself, skip it
        // vertex_cofaces[i][k] for k > 0 are the incident edges
        for (size_t k = 1; k < vertex_cofaces[i].size(); ++k) {
            total_edge_mass += vertex_cofaces[i][k].density;
        }

        eff_deg[i] = total_edge_mass;
    }

    return eff_deg;
}

/**
 * @brief Compute hop-extremp radius for a vertex using incremental BFS
 *
 * This function determines the maximum hop distance at which a vertex maintains
 * a specified extremum probability (maxp or minp) threshold. Unlike the naive
 * approach of repeatedly calling BFS for each hop distance, this implementation
 * performs a single BFS traversal, testing the extremp condition incrementally
 * at each hop level.
 *
 * MOTIVATION AND PROBLEM STATEMENT
 * =================================
 *
 * The hop-extremp radius provides a spatial persistence measure for local extrema.
 * A vertex with large hop-extremp radius maintains its extremum character over
 * a large neighborhood, indicating a genuine, stable peak or valley. A vertex
 * with small hop-extremp radius may be a spurious extremum caused by local noise.
 *
 * The hop-extremp radius r_p(i) is defined as the largest hop distance h such that:
 *   extremp(i | N_h(i), y, ρ) ≥ p
 *
 * where N_h(i) is the h-hop neighborhood of vertex i, extremp is the probabilistic
 * extremum score (maxp or minp), and p is a threshold (typically 0.90 or 0.95).
 *
 * ALGORITHMIC APPROACH
 * ====================
 *
 * We perform a single BFS traversal, maintaining:
 * - Current hop neighborhood N_h(i) as we expand level by level
 * - Density-weighted sums for numerator and denominator of extremp score
 * - Cumulative extremp score at each hop level
 *
 * At each hop level h, we:
 * 1. Add newly discovered vertices at distance h to the neighborhood
 * 2. Update the running sums for extremp calculation
 * 3. Test if extremp(i | N_h(i), y, ρ) ≥ p
 * 4. If extremp drops below p, return h-1 as the hop-extremp radius
 *
 * This achieves O(V + E) complexity for the entire computation, where V is the
 * number of vertices and E is the number of edges in the graph.
 *
 * INCREMENTAL EXTREMP CALCULATION
 * ================================
 *
 * For maxima, we track:
 *   numerator   = Σ{ρ(j) : j ∈ N_h(i), y(i) > y(j)}
 *   denominator = Σ{ρ(j) : j ∈ N_h(i)}
 *   maxp(i | N_h(i)) = numerator / denominator
 *
 * As we expand from N_h to N_{h+1}, we add new vertices discovered at hop h+1:
 *   For each new vertex j:
 *     denominator += ρ(j)
 *     if y(i) > y(j): numerator += ρ(j)
 *
 * This incremental update requires only O(degree) work per hop level, avoiding
 * recomputation of the entire neighborhood score.
 *
 * SPECIAL CASES
 * =============
 *
 * Global extrema: If the vertex is a global extremum, it satisfies the extremp
 * condition over the entire graph. We return std::numeric_limits<size_t>::max()
 * to indicate infinite radius.
 *
 * Non-extrema: If the vertex fails the extremp condition even at hop 1, we
 * return 0 to indicate it is not a local extremum at any scale.
 *
 * Maximum hop limit: To prevent excessive computation on large graphs, we impose
 * a maximum hop limit (default 20). If extremp remains above threshold at this
 * limit, we return the limit value (not infinity, to distinguish from global extrema).
 *
 * @param vertex Index of the vertex to evaluate
 * @param y Vector of function values at each vertex (typically fitted values)
 * @param p_threshold Extremp threshold in (0,1], typically 0.90 or 0.95
 * @param detect_maxima If true, compute hop-maxp radius; if false, hop-minp radius
 * @param max_hop Maximum hop distance to explore (default 20)
 *
 * @return Hop-extremp radius: the largest h such that extremp(i | N_h(i)) ≥ p.
 *         Returns 0 if vertex is not a local extremum.
 *         Returns std::numeric_limits<size_t>::max() if vertex is global extremum.
 *         Returns max_hop if threshold still satisfied at maximum hop limit.
 *
 * @pre y.size() == vertex_cofaces.size()
 * @pre vertex < vertex_cofaces.size()
 * @pre 0 < p_threshold ≤ 1.0
 * @pre vertex_cofaces[i][0].density contains vertex densities ρ₀
 *
 * @complexity O(min(V + E, max_hop · average_degree)) where V is vertices, E is edges.
 *             In practice, terminates early when extremp drops below threshold.
 *
 * @note Uses vertex densities ρ₀ from vertex_cofaces[i][0].density for weighting.
 *       This matches the definition of maxp/minp in the R implementation.
 */
size_t riem_dcx_t::compute_hop_extremp_radius(
    size_t vertex,
    const vec_t& y,
    double p_threshold,
    bool detect_maxima,
    size_t max_hop
    ) const {

    const size_t n = vertex_cofaces.size();

    // ================================================================
    // ARGUMENT VALIDATION
    // ================================================================

    if (vertex >= n) {
        Rf_error("compute_hop_extremp_radius: vertex index %zu out of bounds [0, %zu)",
                 vertex, n);
    }

    if (static_cast<size_t>(y.size()) != n) {
        Rf_error("compute_hop_extremp_radius: size mismatch (y: %d, graph: %zu)",
                 static_cast<int>(y.size()), n);
    }

    if (p_threshold <= 0.0 || p_threshold > 1.0) {
        Rf_error("compute_hop_extremp_radius: p_threshold must be in (0,1], got %.3f",
                 p_threshold);
    }

    // ================================================================
    // CHECK FOR GLOBAL EXTREMUM
    // ================================================================

    double global_min = y.minCoeff();
    double global_max = y.maxCoeff();
    const double epsilon = 1e-14 * std::max(std::abs(global_min), std::abs(global_max));

    bool is_global_min = std::abs(y[vertex] - global_min) <= epsilon;
    bool is_global_max = std::abs(y[vertex] - global_max) <= epsilon;

    if ((detect_maxima && is_global_max) || (!detect_maxima && is_global_min)) {
        return std::numeric_limits<size_t>::max();
    }

    // ================================================================
    // CHECK 1-HOP LOCAL EXTREMUM CONDITION
    // ================================================================

    // First verify this is actually a 1-hop local extremum
    // vertex_cofaces[vertex][k] for k > 0 are the neighbors
    for (size_t k = 1; k < vertex_cofaces[vertex].size(); ++k) {
        index_t u = vertex_cofaces[vertex][k].vertex_index;

        if (detect_maxima) {
            if (y[u] >= y[vertex]) {
                return 0;  // Not a local maximum
            }
        } else {
            if (y[u] <= y[vertex]) {
                return 0;  // Not a local minimum
            }
        }
    }

    // ================================================================
    // INCREMENTAL BFS WITH EXTREMP TRACKING
    // ================================================================

    std::vector<bool> visited(n, false);
    std::queue<size_t> current_level;
    std::queue<size_t> next_level;

    // Initialize with start vertex - but don't include it in extremp calculation
    // (extremp is computed over *neighbors*, not including the vertex itself)
    visited[vertex] = true;
    current_level.push(vertex);

    // Running sums for extremp calculation
    double total_rho = 0.0;        // Denominator: Σ ρ(j) over all j in neighborhood
    double satisfying_rho = 0.0;   // Numerator: Σ ρ(j) where j satisfies condition

    size_t current_hop = 0;
    size_t last_valid_hop = 0;

    // Process BFS level by level
    while (!current_level.empty() && current_hop < max_hop) {
        // Expand current level to discover next level
        while (!current_level.empty()) {
            size_t v = current_level.front();
            current_level.pop();

            // Explore all neighbors of v
            // vertex_cofaces[v][k] for k > 0 contains neighbors
            for (size_t k = 1; k < vertex_cofaces[v].size(); ++k) {
                index_t u = vertex_cofaces[v][k].vertex_index;

                if (visited[u]) continue;

                visited[u] = true;
                next_level.push(u);

                // Get vertex density from vertex_cofaces[u][0].density
                double rho_u = vertex_cofaces[u][0].density;

                // Update running sums for extremp calculation
                total_rho += rho_u;

                bool satisfies_condition;
                if (detect_maxima) {
                    satisfies_condition = (y[vertex] > y[u]);
                } else {
                    satisfies_condition = (y[vertex] < y[u]);
                }

                if (satisfies_condition) {
                    satisfying_rho += rho_u;
                }
            }
        }

        // After adding vertices at hop (current_hop + 1), compute extremp
        current_hop++;

        // Check if we have any neighbors
        if (total_rho < std::numeric_limits<double>::epsilon()) {
            // No neighbors with non-zero density - shouldn't happen but be safe
            return last_valid_hop;
        }

        // Compute extremp score over current neighborhood
        double current_extremp = satisfying_rho / total_rho;

        // Check threshold
        if (current_extremp >= p_threshold) {
            last_valid_hop = current_hop;
        } else {
            // Extremp dropped below threshold - return previous hop
            return last_valid_hop;
        }

        // Move to next level
        std::swap(current_level, next_level);
    }

    // If we reached max_hop and still above threshold, return max_hop
    // (not infinity, to distinguish from true global extrema)
    return last_valid_hop;
}


/**
 * @brief Compute hop-extremp radii for multiple vertices
 *
 * Batch computation of hop-extremp radius for a collection of candidate vertices.
 * This is more efficient than calling compute_hop_extremp_radius() individually
 * when processing many vertices, as it allows for potential optimizations and
 * better cache locality.
 *
 * @param vertices Vector of vertex indices to process
 * @param y Vector of function values at each vertex
 * @param p_threshold Extremp threshold in (0,1]
 * @param detect_maxima If true, compute hop-maxp; if false, hop-minp
 * @param max_hop Maximum hop distance to explore
 *
 * @return Vector of hop-extremp radii, parallel to input vertices vector
 *
 * @complexity O(|vertices| · (V + E)) worst case, but typically much faster
 *             due to early termination when extremp drops below threshold
 */
std::vector<size_t> riem_dcx_t::compute_hop_extremp_radii_batch(
    const std::vector<size_t>& vertices,
    const vec_t& y,
    double p_threshold,
    bool detect_maxima,
    size_t max_hop
    ) const {

    std::vector<size_t> radii;
    radii.reserve(vertices.size());

    for (size_t vertex : vertices) {
        size_t radius = compute_hop_extremp_radius(
            vertex, y, p_threshold, detect_maxima, max_hop
        );
        radii.push_back(radius);
    }

    return radii;
}

/**
 * @brief R interface for computing hop-extremp radii for candidate vertices
 *
 * This function provides R access to the C++ hop-extremp radius computation,
 * allowing efficient batch processing of candidate extrema vertices. It accepts
 * the graph structure as returned by fit.knn.riem.graph.regression(), which
 * stores adjacency lists, vertex densities, and edge densities separately rather
 * than in a unified vertex_cofaces structure.
 *
 * INPUT ARGUMENTS (from R)
 * ========================
 *
 * @param s_adj_list List of integer vectors representing graph adjacency.
 *        adj_list[[i]] contains the neighbor indices of vertex i (1-indexed in R).
 *
 * @param s_edge_densities Numeric vector of edge densities in adjacency list order.
 *        For vertex i with neighbors adj_list[[i]], the corresponding edge densities
 *        are the next length(adj_list[[i]]) values in edge_densities.
 *
 * @param s_vertex_densities Numeric vector of vertex densities ρ₀.
 *
 * @param s_candidates Integer vector of candidate vertex indices (1-indexed in R).
 *
 * @param s_y Numeric vector of function values at all vertices.
 *
 * @param s_p_threshold Numeric scalar extremp threshold in (0,1].
 *
 * @param s_detect_maxima Logical scalar: TRUE for maxima, FALSE for minima.
 *
 * @param s_max_hop Integer scalar maximum hop distance to explore.
 */
extern "C" SEXP S_compute_hop_extremp_radii_batch(
    SEXP s_adj_list,
    SEXP s_edge_densities,
    SEXP s_vertex_densities,
    SEXP s_candidates,
    SEXP s_y,
    SEXP s_p_threshold,
    SEXP s_detect_maxima,
    SEXP s_max_hop
    ) {
    // ================================================================
    // INPUT VALIDATION
    // ================================================================

    if (!Rf_isNewList(s_adj_list)) {
        Rf_error("S_compute_hop_extremp_radii_batch: adj_list must be a list");
    }

    if (!Rf_isReal(s_edge_densities)) {
        Rf_error("S_compute_hop_extremp_radii_batch: edge_densities must be numeric");
    }

    if (!Rf_isReal(s_vertex_densities)) {
        Rf_error("S_compute_hop_extremp_radii_batch: vertex_densities must be numeric");
    }

    if (!Rf_isInteger(s_candidates)) {
        Rf_error("S_compute_hop_extremp_radii_batch: candidates must be an integer vector");
    }

    if (!Rf_isReal(s_y)) {
        Rf_error("S_compute_hop_extremp_radii_batch: y must be a numeric vector");
    }

    if (!Rf_isReal(s_p_threshold) || Rf_length(s_p_threshold) != 1) {
        Rf_error("S_compute_hop_extremp_radii_batch: p_threshold must be a numeric scalar");
    }

    if (!Rf_isLogical(s_detect_maxima) || Rf_length(s_detect_maxima) != 1) {
        Rf_error("S_compute_hop_extremp_radii_batch: detect_maxima must be a logical scalar");
    }

    if (!Rf_isInteger(s_max_hop) || Rf_length(s_max_hop) != 1) {
        Rf_error("S_compute_hop_extremp_radii_batch: max_hop must be an integer scalar");
    }

    // ================================================================
    // EXTRACT PARAMETERS
    // ================================================================

    const size_t n_vertices = Rf_length(s_adj_list);
    const size_t n_candidates = Rf_length(s_candidates);

    int* candidates_ptr = INTEGER(s_candidates);
    double* y_ptr = REAL(s_y);
    double* vertex_dens_ptr = REAL(s_vertex_densities);
    double* edge_dens_ptr = REAL(s_edge_densities);
    double p_threshold = REAL(s_p_threshold)[0];
    bool detect_maxima = LOGICAL(s_detect_maxima)[0];
    int max_hop_int = INTEGER(s_max_hop)[0];

    // ================================================================
    // VALIDATE PARAMETER VALUES
    // ================================================================

    if (static_cast<size_t>(Rf_length(s_y)) != n_vertices) {
        Rf_error("S_compute_hop_extremp_radii_batch: length of y (%d) must match number of vertices (%zu)",
                 Rf_length(s_y), n_vertices);
    }

    if (static_cast<size_t>(Rf_length(s_vertex_densities)) != n_vertices) {
        Rf_error("S_compute_hop_extremp_radii_batch: length of vertex_densities (%d) must match number of vertices (%zu)",
                 Rf_length(s_vertex_densities), n_vertices);
    }

    if (p_threshold <= 0.0 || p_threshold > 1.0) {
        Rf_error("S_compute_hop_extremp_radii_batch: p_threshold must be in (0,1], got %.3f",
                 p_threshold);
    }

    if (max_hop_int <= 0) {
        Rf_error("S_compute_hop_extremp_radii_batch: max_hop must be positive, got %d",
                 max_hop_int);
    }
    size_t max_hop = static_cast<size_t>(max_hop_int);

    // Validate candidate indices are in bounds [1, n_vertices]
    for (size_t i = 0; i < n_candidates; ++i) {
        int candidate = candidates_ptr[i];
        if (candidate < 1 || candidate > static_cast<int>(n_vertices)) {
            Rf_error("S_compute_hop_extremp_radii_batch: candidate index %d out of bounds [1, %zu]",
                     candidate, n_vertices);
        }
    }

    // ================================================================
    // RECONSTRUCT VERTEX_COFACES FROM ADJACENCY LIST + DENSITIES
    // ================================================================

    std::vector<std::vector<neighbor_info_t>> vertex_cofaces(n_vertices);
    size_t edge_idx = 0;

    for (size_t i = 0; i < n_vertices; ++i) {
        SEXP s_neighbors_i = VECTOR_ELT(s_adj_list, i);

        if (!Rf_isInteger(s_neighbors_i)) {
            Rf_error("S_compute_hop_extremp_radii_batch: adj_list[[%zu]] must be an integer vector",
                     i + 1);
        }

        int* neighbors_ptr = INTEGER(s_neighbors_i);
        size_t n_neighbors = Rf_length(s_neighbors_i);

        vertex_cofaces[i].reserve(n_neighbors + 1);

        // Add self-entry at position [0]
        neighbor_info_t self_entry;
        self_entry.vertex_index = static_cast<index_t>(i);
        self_entry.simplex_index = static_cast<index_t>(i);
        self_entry.isize = 0;
        self_entry.dist = 0.0;
        self_entry.density = vertex_dens_ptr[i];
        vertex_cofaces[i].push_back(self_entry);

        // Add neighbor entries
        for (size_t k = 0; k < n_neighbors; ++k) {
            int j = neighbors_ptr[k] - 1;  // Convert to 0-based

            if (j < 0 || j >= static_cast<int>(n_vertices)) {
                Rf_error("S_compute_hop_extremp_radii_batch: neighbor index %d out of bounds",
                         neighbors_ptr[k]);
            }

            neighbor_info_t nbr_entry;
            nbr_entry.vertex_index = static_cast<index_t>(j);
            nbr_entry.simplex_index = static_cast<index_t>(edge_idx);  // Edge index
            nbr_entry.isize = 0;
            nbr_entry.dist = 0.0;
            nbr_entry.density = edge_dens_ptr[edge_idx];
            vertex_cofaces[i].push_back(nbr_entry);

            edge_idx++;
        }
    }

    // ================================================================
    // CONVERT Y TO EIGEN VECTOR
    // ================================================================

    vec_t y(n_vertices);
    for (size_t i = 0; i < n_vertices; ++i) {
        y[i] = y_ptr[i];
    }

    // ================================================================
    // CONVERT CANDIDATES TO 0-BASED INDICES
    // ================================================================

    std::vector<size_t> candidates(n_candidates);
    for (size_t i = 0; i < n_candidates; ++i) {
        candidates[i] = static_cast<size_t>(candidates_ptr[i] - 1); // Convert to 0-based
    }

    // ================================================================
    // CREATE MINIMAL RIEM_DCX_T OBJECT
    // ================================================================

    riem_dcx_t rdcx;
    rdcx.vertex_cofaces = std::move(vertex_cofaces);

    // ================================================================
    // COMPUTE HOP-EXTREMP RADII
    // ================================================================

    std::vector<size_t> radii = rdcx.compute_hop_extremp_radii_batch(
        candidates,
        y,
        p_threshold,
        detect_maxima,
        max_hop
    );

    // ================================================================
    // CONVERT RESULTS TO R INTEGER VECTOR
    // ================================================================

    SEXP s_result = PROTECT(Rf_allocVector(INTSXP, n_candidates));
    int* result_ptr = INTEGER(s_result);

    for (size_t i = 0; i < n_candidates; ++i) {
        // Special handling for global extrema (std::numeric_limits<size_t>::max())
        // Encode as -1 since R integer vectors can't represent such large values
        if (radii[i] == std::numeric_limits<size_t>::max()) {
            result_ptr[i] = -1;
        } else {
            result_ptr[i] = static_cast<int>(radii[i]);
        }
    }

    UNPROTECT(1);
    return s_result;
}
