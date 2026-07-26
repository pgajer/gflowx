#ifndef CENTERED_PATHS_H_
#define CENTERED_PATHS_H_

#include "uniform_grid_graph.hpp"
#include "set_wgraph.hpp" // for vertex_info_t
#include "edge_weights.hpp"

#include <vector>
#include <limits>  // for std::numeric_limits
#include <cstddef>
using std::size_t;

/**
 * @brief A structure containing path data for composite paths constructed from
 * the shortest paths starting at a __grid vertex__ of a graph
 *
 * If there is only one shortest path starting at the given grid vertex (of
 * length bound by the bandwidth parameter) the path itself will be a degenerate
 * composite path. If the set of shortest paths starting at the ref vertex has
 * more than one element, each pair of these paths is combined to form a
 * composite path. In most of the cases the reference vertex that is a graph
 * grid vertex will not be a part of the path as the path consist of only
 * original graph vertices over which the response variable is defined. Thus,
 * ref_vertex is a grid vertex index that may in rare cases be an original
 * vertex. The structure in a non-grid case had a parameter ref_vertex_index
 * that was removed from this structure as in most cases it will not make sense
 * (as the ref vertex is not part of the vertices vector).
 *
 * @details The structure maintains several key components:
 *   - The sequence of vertices in the path
 *   - Information about a reference vertex, including its position in both the
 *     original graph and within this path
 *   - Metrics about the path's structure and the reference vertex's position
 *   - Arrays of computed values along the path for analysis and visualization
 */
struct path_data_t {
    std::vector<size_t> vertices;   ///< Sequence of vertices forming the path
    size_t ref_vertex;              ///< Index of the reference vertex in the original graph
    //size_t ref_vertex_index;        ///< Position of ref_vertex within the vertices vector
    double rel_center_offset;       ///< Relative position of ref_vertex in path (0.0 = center, 0.5 = boundary)
    double total_weight;            ///< Total length of the path using edge weights

    std::vector<double> x_path;     ///< Cumulative distance along path from initial vertex
    std::vector<double> w_path;     ///< Kernel-weighted distances from reference vertex
    std::vector<double> y_path;     ///< Y-values of vertices restricted to the path

    path_data_t()
        : ref_vertex(INVALID_VERTEX),
//          ref_vertex_index(INVALID_VERTEX),
          rel_center_offset(0),
          total_weight(0) {}
};

reachability_map_t precompute_max_bandwidth_paths(
    const uniform_grid_graph_t& uniform_grid_graph,
    size_t ref_vertex,
    double max_bandwidth);

std::vector<path_data_t> ugg_get_path_data_efficient(
    const uniform_grid_graph_t& uniform_grid_graph,
    const std::vector<double>& y,
    size_t ref_vertex,
    double bandwidth,
    const reachability_map_t& reachability_map,
    double dist_normalization_factor,
    size_t min_path_size,
    size_t diff_threshold,
    size_t kernel_type,
    const edge_weights_t& edge_weights);

std::vector<path_data_t> ugg_get_path_data(
    const uniform_grid_graph_t& uniform_grid_graph,
    const std::vector<double>& y,
    size_t ref_vertex,
    double bandwidth,
    double dist_normalization_factor,
    size_t min_path_size,
    size_t diff_threshold,
    size_t kernel_type,
    const edge_weights_t& edge_weights,
    bool verbose
    );

#endif // CENTERED_PATHS_H_
