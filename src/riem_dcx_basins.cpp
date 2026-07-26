#include "riem_dcx.hpp"
#include <queue>                 // For std::queue (BFS)
#include <vector>                // For std::vector
#include <unordered_map>         // For basin_t::reachability_map_t::distances, predecessors
#include <map>                   // For basin_t::boundary_vertices_map, boundary_monotonicity_spans_map
#include <algorithm>             // For std::sort
#include <limits>                // For INFINITY
#include <cmath>                 // For std::max (if using it for density clamping)

/**
 * @brief Find the basin of attraction for a rho-local extremum using BFS
 *
 * @details Constructs the basin of attraction around a rho-local extremum by
 *          exploring the graph using breadth-first search. Unlike classical
 *          extrema, rho-extrema use weighted function values that incorporate
 *          vertex densities from the Riemannian structure.
 *
 *          For a rho-local maximum M, the basin includes all vertices from which
 *          there exists an ascending path where weighted function values
 *          w_v * y_hat(v) strictly increase along edges toward M.
 *
 * @param vertex The vertex index of the rho-local extremum
 * @param y_hat Vector of fitted function values at each vertex
 * @param detect_maxima If true, construct maximum basin; if false, minimum basin
 *
 * @return basin_t structure containing the complete basin information
 *
 * @note The weighted function value at vertex v is computed as:
 *       weighted_y[v] = vertex_cofaces[v][0].density * y_hat[v]
 *       where vertex_cofaces[v][0].density is the vertex density rho_0(v)
 */
basin_t riem_dcx_t::find_rho_extremum_basin(
    size_t vertex,
    const vec_t& y_hat,
    bool detect_maxima
) const {

    // ================================================================
    // PRECOMPUTE WEIGHTED FUNCTION VALUES
    // ================================================================

    const size_t n = vertex_cofaces.size();
    std::vector<double> weighted_y(n);

    for (size_t v = 0; v < n; ++v) {
        // vertex_cofaces[v][0] contains the vertex itself with its density
        double rho_v = std::max(vertex_cofaces[v][0].density, 1e-15);
        weighted_y[v] = rho_v * y_hat[static_cast<Eigen::Index>(v)];
    }

    // ================================================================
    // INITIALIZE BASIN STRUCTURE
    // ================================================================

    basin_t basin;
    basin.value = weighted_y[vertex];
    basin.is_maximum = detect_maxima;
    basin.reachability_map.ref_vertex = vertex;
    basin.min_monotonicity_span = INFINITY;
    basin.min_span_vertex = INVALID_VERTEX;
    basin.extremum_vertex = vertex;

    // Seed the basin with the extremum
    basin.reachability_map.distances[vertex] = 0.0;
    basin.reachability_map.predecessors[vertex] = INVALID_VERTEX;
    basin.reachability_map.sorted_vertices.push_back({vertex, 0.0});

    // ================================================================
    // INITIALIZE BFS DATA STRUCTURES
    // ================================================================

    std::vector<int> hop_distance(n, -1);
    std::vector<size_t> prev(n, INVALID_VERTEX);
    std::queue<size_t> q;
    std::vector<bool> in_queue(n, false);
    std::vector<bool> in_basin(n, false);

    hop_distance[vertex] = 0;
    q.push(vertex);
    in_queue[vertex] = true;
    in_basin[vertex] = true;

    // ================================================================
    // BREADTH-FIRST SEARCH WITH MONOTONICITY CONSTRAINTS
    // ================================================================

    while (!q.empty()) {
        size_t u = q.front();
        q.pop();
        in_queue[u] = false;

        // For vertices other than the extremum
        if (u != vertex) {
            // Check monotonicity using weighted function values
            double delta_weighted_y = weighted_y[u] - weighted_y[prev[u]];

            // For maxima: weighted values must decrease (delta < 0)
            // For minima: weighted values must increase (delta > 0)
            bool condition_met = detect_maxima ?
                (delta_weighted_y < 0.0) : (delta_weighted_y > 0.0);

            if (!condition_met) {
                // MONOTONICITY VIOLATION BOUNDARY

                // Calculate monotonicity span using weighted values
                double curr_span = detect_maxima ?
                    (weighted_y[vertex] - weighted_y[prev[u]]) :
                    (weighted_y[prev[u]] - weighted_y[vertex]);

                basin.boundary_monotonicity_spans_map[prev[u]] = curr_span;

                constexpr double eps = 1e-12;
                if (curr_span > eps && curr_span < basin.min_monotonicity_span) {
                    basin.min_monotonicity_span = curr_span;
                    basin.min_span_vertex = prev[u];
                }

                basin.boundary_vertices_map[hop_distance[prev[u]]] = prev[u];

                continue;  // Skip this violating vertex
            }

            // Add vertex to basin
            double dist = static_cast<double>(hop_distance[u]);
            basin.reachability_map.sorted_vertices.push_back({u, dist});
            basin.reachability_map.distances[u] = dist;
            basin.reachability_map.predecessors[u] = prev[u];
            in_basin[u] = true;
        }

        // ================================================================
        // EXPLORE NEIGHBORS USING vertex_cofaces
        // ================================================================

        // vertex_cofaces[u][0] is u itself, skip it
        // vertex_cofaces[u][j] for j > 0 are the incident edges
        for (size_t j = 1; j < vertex_cofaces[u].size(); ++j) {
            size_t v = vertex_cofaces[u][j].vertex_index;

            // If not visited
            if (hop_distance[v] == -1) {
                hop_distance[v] = hop_distance[u] + 1;
                prev[v] = u;

                if (!in_queue[v]) {
                    q.push(v);
                    in_queue[v] = true;
                }
            }
        }
    } // END BFS

    // ================================================================
    // IDENTIFY GRAPH BOUNDARY VERTICES
    // ================================================================

    for (const auto& v_info : basin.reachability_map.sorted_vertices) {
        size_t u = v_info.vertex;

        bool is_boundary = false;

        // Leaf vertices (degree 1) are boundaries
        if (vertex_cofaces[u].size() == 2) {  // [0] is self, [1] is only neighbor
            is_boundary = true;
        } else {
            // Check if any neighbors are outside the basin
            for (size_t j = 1; j < vertex_cofaces[u].size(); ++j) {
                size_t v = vertex_cofaces[u][j].vertex_index;

                if (basin.reachability_map.distances.find(v) ==
                    basin.reachability_map.distances.end()) {
                    is_boundary = true;
                    break;
                }
            }
        }

        if (is_boundary) {
            double dist = basin.reachability_map.distances[u];
            basin.boundary_vertices_map[dist] = u;

            double span = detect_maxima ?
                (weighted_y[vertex] - weighted_y[u]) :
                (weighted_y[u] - weighted_y[vertex]);
            basin.boundary_monotonicity_spans_map[u] = span;
        }
    }

    // ================================================================
    // SORT VERTICES BY DESCENDING DISTANCE
    // ================================================================

    std::sort(
        basin.reachability_map.sorted_vertices.begin(),
        basin.reachability_map.sorted_vertices.end(),
        [](const vertex_info_t& a, const vertex_info_t& b) {
            return a.distance > b.distance;
        }
    );

    return basin;
}
