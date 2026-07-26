/**
 * @file debug_serialization.hpp
 * @brief Debugging utilities for ikNN graph construction verification
 *
 * This file provides functions to serialize intermediate data structures
 * during ikNN graph construction to binary files for detailed comparison
 * and debugging.
 *
 * These functions are only called when DEBUG flags are enabled in source files.
 * When disabled (default), they add zero overhead as the calls are removed
 * by the preprocessor.
 *
 * Usage: Set DEBUG_INITIALIZE_FROM_KNN or DEBUG_CREATE_IKNN_GRAPH to 1
 * in the respective source files to enable debugging output.
 */

#ifndef DEBUG_SERIALIZATION_HPP
#define DEBUG_SERIALIZATION_HPP

#include <string>
#include <vector>
#include <fstream>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <iomanip>
#include "riem_dcx.hpp"

/**
 * @brief Serialization utilities for debugging ikNN graph construction
 *
 * Always available - controlled by local DEBUG flags in source files
 */

namespace debug_serialization {

inline void save_knn_result(
    const std::string& filepath,
    const std::vector<std::vector<index_t>>& knn_indices,
    const std::vector<std::vector<double>>& knn_distances,
    size_t n_points,
    size_t k
) {
    std::ofstream out(filepath, std::ios::binary);
    if (!out) {
        Rf_warning("Cannot open file for writing: %s", filepath.c_str());
        return;
    }

    // Write header
    out.write("KNN_RESULT", 10);
    out.write(reinterpret_cast<const char*>(&n_points), sizeof(size_t));
    out.write(reinterpret_cast<const char*>(&k), sizeof(size_t));

    // Write indices
    for (size_t i = 0; i < n_points; ++i) {
        out.write(reinterpret_cast<const char*>(knn_indices[i].data()), k * sizeof(index_t));
    }

    // Write distances
    for (size_t i = 0; i < n_points; ++i) {
        out.write(reinterpret_cast<const char*>(knn_distances[i].data()), k * sizeof(double));
    }

    out.close();
}

inline void save_neighbor_sets(
    const std::string& filepath,
    const std::vector<std::unordered_set<index_t>>& neighbor_sets
) {
    std::ofstream out(filepath, std::ios::binary);
    if (!out) {
        Rf_warning("Cannot open file for writing: %s", filepath.c_str());
        return;
    }

    out.write("NEIGHBOR_SETS", 13);
    size_t n_points = neighbor_sets.size();
    out.write(reinterpret_cast<const char*>(&n_points), sizeof(size_t));

    for (size_t i = 0; i < n_points; ++i) {
        size_t set_size = neighbor_sets[i].size();
        out.write(reinterpret_cast<const char*>(&set_size), sizeof(size_t));

        // Convert set to sorted vector for deterministic output
        std::vector<index_t> sorted_neighbors(neighbor_sets[i].begin(), neighbor_sets[i].end());
        std::sort(sorted_neighbors.begin(), sorted_neighbors.end());
        out.write(reinterpret_cast<const char*>(sorted_neighbors.data()),
                  set_size * sizeof(index_t));
    }

    out.close();
}

inline void save_vertex_cofaces(
    const std::string& filepath,
    const std::vector<std::vector<neighbor_info_t>>& vertex_cofaces
) {
    std::ofstream out(filepath, std::ios::binary);
    if (!out) {
        Rf_warning("Cannot open file for writing: %s", filepath.c_str());
        return;
    }

    out.write("VERTEX_COFACES", 14);
    size_t n_vertices = vertex_cofaces.size();
    out.write(reinterpret_cast<const char*>(&n_vertices), sizeof(size_t));

    for (size_t i = 0; i < n_vertices; ++i) {
        size_t n_neighbors = vertex_cofaces[i].size();
        out.write(reinterpret_cast<const char*>(&n_neighbors), sizeof(size_t));

        for (const auto& info : vertex_cofaces[i]) {
            out.write(reinterpret_cast<const char*>(&info.vertex_index), sizeof(index_t));
            out.write(reinterpret_cast<const char*>(&info.simplex_index), sizeof(index_t));
            out.write(reinterpret_cast<const char*>(&info.isize), sizeof(size_t));
            out.write(reinterpret_cast<const char*>(&info.dist), sizeof(double));
            out.write(reinterpret_cast<const char*>(&info.density), sizeof(double));
        }
    }

    out.close();
}

inline void save_edge_list(
    const std::string& filepath,
    const std::vector<std::pair<index_t, index_t>>& edges,
    const std::vector<double>& weights,
    const char* phase_name
) {
    std::ofstream out(filepath, std::ios::binary);
    if (!out) {
        Rf_warning("Cannot open file for writing: %s", filepath.c_str());
        return;
    }

    out.write("EDGE_LIST", 9);
    size_t phase_len = std::strlen(phase_name);
    out.write(reinterpret_cast<const char*>(&phase_len), sizeof(size_t));
    out.write(phase_name, phase_len);

    size_t n_edges = edges.size();
    out.write(reinterpret_cast<const char*>(&n_edges), sizeof(size_t));

    for (size_t e = 0; e < n_edges; ++e) {
        out.write(reinterpret_cast<const char*>(&edges[e].first), sizeof(index_t));
        out.write(reinterpret_cast<const char*>(&edges[e].second), sizeof(index_t));
        out.write(reinterpret_cast<const char*>(&weights[e]), sizeof(double));
    }

    out.close();
}

inline void save_pruning_params(
    const std::string& filepath,
    double max_ratio_threshold,
    double threshold_percentile,
    size_t n_edges_before,
    size_t n_edges_after
) {
    std::ofstream out(filepath, std::ios::binary);
    if (!out) {
        Rf_warning("Cannot open file for writing: %s", filepath.c_str());
        return;
    }

    out.write("PRUNING_PARAMS", 14);
    out.write(reinterpret_cast<const char*>(&max_ratio_threshold), sizeof(double));
    out.write(reinterpret_cast<const char*>(&threshold_percentile), sizeof(double));
    out.write(reinterpret_cast<const char*>(&n_edges_before), sizeof(size_t));
    out.write(reinterpret_cast<const char*>(&n_edges_after), sizeof(size_t));

    out.close();
}

inline void save_connectivity(
    const std::string& filepath,
    const std::vector<int>& component_ids,
    size_t n_components
) {
    std::ofstream out(filepath, std::ios::binary);
    if (!out) {
        Rf_warning("Cannot open file for writing: %s", filepath.c_str());
        return;
    }

    out.write("CONNECTIVITY", 12);
    size_t n_vertices = component_ids.size();
    out.write(reinterpret_cast<const char*>(&n_vertices), sizeof(size_t));
    out.write(reinterpret_cast<const char*>(&n_components), sizeof(size_t));
    out.write(reinterpret_cast<const char*>(component_ids.data()),
              n_vertices * sizeof(int));

    out.close();
}

/**
 * @brief Compute connected components from vertex_cofaces structure
 */
inline size_t compute_connected_components_from_vertex_cofaces(
    const std::vector<std::vector<neighbor_info_t>>& vertex_cofaces,
    std::vector<int>& component_ids
) {
    const size_t n_vertices = vertex_cofaces.size();
    component_ids.assign(n_vertices, -1);

    int current_component = 0;
    std::vector<index_t> stack;

    for (size_t seed = 0; seed < n_vertices; ++seed) {
        if (component_ids[seed] >= 0) continue;

        // BFS/DFS from this seed
        stack.clear();
        stack.push_back(seed);
        component_ids[seed] = current_component;

        while (!stack.empty()) {
            index_t v = stack.back();
            stack.pop_back();

            // Visit all neighbors (skip self-loop at position 0)
            for (size_t k = 1; k < vertex_cofaces[v].size(); ++k) {
                index_t neighbor = vertex_cofaces[v][k].vertex_index;
                if (component_ids[neighbor] < 0) {
                    component_ids[neighbor] = current_component;
                    stack.push_back(neighbor);
                }
            }
        }

        ++current_component;
    }

    return current_component;
}

/**
 * @brief Compute connected components from iknn_graph_t structure
 */
inline size_t compute_connected_components_from_iknn_graph(
    const iknn_graph_t& graph,
    std::vector<int>& component_ids
) {
    const size_t n_vertices = graph.graph.size();
    component_ids.assign(n_vertices, -1);

    int current_component = 0;
    std::vector<index_t> stack;

    for (size_t seed = 0; seed < n_vertices; ++seed) {
        if (component_ids[seed] >= 0) continue;

        stack.clear();
        stack.push_back(seed);
        component_ids[seed] = current_component;

        while (!stack.empty()) {
            index_t v = stack.back();
            stack.pop_back();

            for (const auto& neighbor_info : graph.graph[v]) {
                index_t neighbor = neighbor_info.index;
                if (component_ids[neighbor] < 0) {
                    component_ids[neighbor] = current_component;
                    stack.push_back(neighbor);
                }
            }
        }

        ++current_component;
    }

    return current_component;
}

/**
 * @brief Save edge list as CSV for easy import to R
 */
inline void save_edge_list_csv(
    const std::string& filepath,
    const std::vector<std::pair<index_t, index_t>>& edges,
    const std::vector<double>& weights
) {
    std::ofstream out(filepath);
    if (!out) {
        Rf_warning("Cannot open file for writing: %s", filepath.c_str());
        return;
    }

    out << std::fixed << std::setprecision(10);
    out << "i,j,weight\n";
    for (size_t e = 0; e < edges.size(); ++e) {
        out << edges[e].first << "," << edges[e].second << "," << weights[e] << "\n";
    }
    out.close();
}

/**
 * @brief Compute connected components from set_wgraph_t structure
 */
inline size_t compute_connected_components_from_set_wgraph(
    const set_wgraph_t& graph,
    std::vector<int>& component_ids
) {
    const size_t n_vertices = graph.adjacency_list.size();
    component_ids.assign(n_vertices, -1);

    int current_component = 0;
    std::vector<size_t> stack;
    stack.reserve(n_vertices);

    for (size_t seed = 0; seed < n_vertices; ++seed) {
        if (component_ids[seed] >= 0) continue;

        // DFS from this seed
        stack.clear();
        stack.push_back(seed);
        component_ids[seed] = current_component;

        while (!stack.empty()) {
            size_t v = stack.back();
            stack.pop_back();

            // Visit all neighbors
            for (const auto& edge_info : graph.adjacency_list[v]) {
                size_t neighbor = edge_info.vertex;
                if (component_ids[neighbor] < 0) {
                    component_ids[neighbor] = current_component;
                    stack.push_back(neighbor);
                }
            }
        }

        ++current_component;
    }

    return current_component;
}

} // namespace debug_serialization

#endif // DEBUG_SERIALIZATION_HPP
