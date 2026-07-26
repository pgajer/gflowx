// traj_clustering.hpp
//
// Trajectory-graph construction utilities for clustering trajectories
// within a single (m, M) cell.
//
// Leiden clustering is declared as a placeholder for a future C++ backend.

#pragma once

#include <string>
#include <vector>

/**
 * @brief Similarity specification for trajectory overlap measures.
 */
struct traj_similarity_spec_t {
    enum class similarity_t { jaccard, idf_jaccard };
    enum class overlap_mode_t { any, min_shared, min_frac, min_subpath };

    similarity_t similarity_type = similarity_t::idf_jaccard;
    overlap_mode_t overlap_mode = overlap_mode_t::min_shared;

    int min_shared = 2;
    double min_frac = 0.05;
    int min_subpath = 4;

    bool exclude_endpoints = true;
    double idf_smooth = 1.0;
};

/**
 * @brief Sparse weighted graph in edge-list form (undirected).
 */
struct traj_graph_t {
    int n_traj = 0;
    std::vector<int> from;   // 0-based trajectory index
    std::vector<int> to;     // 0-based trajectory index
    std::vector<double> weight;
};

/**
 * @brief Build a kNN-sparsified trajectory similarity graph for a fixed cell (m, M).
 *
 * @details
 * The input is a set of trajectories, each given as a sequence of vertex indices
 * in the underlying data graph (0-based indexing in C++). Similarity between two
 * trajectories is computed from the overlap of intermediate vertices (optionally
 * excluding endpoints) using either unweighted Jaccard or IDF-weighted Jaccard.
 *
 * The full similarity graph is not materialized; instead, per-trajectory top-k
 * neighbors are retained and then symmetrized ("mutual" or "union").
 *
 * @param traj_list Trajectories; each trajectory is a vector of vertex indices (0-based).
 * @param spec Similarity/overlap specification.
 * @param k Number of nearest neighbors retained per trajectory.
 * @param symmetrize One of: "mutual", "union", "none" (treated as "union").
 * @param knn_select One of: "weighted" (use spec similarity for top-k) or "raw"
 *   (use unweighted Jaccard for top-k, still using spec similarity for edge weights).
 * @param n_threads Number of threads used for construction (OpenMP if available).
 *
 * @return A sparse weighted graph over trajectories (edge list).
 */
traj_graph_t build_traj_knn_graph(const std::vector<std::vector<int>>& traj_list,
                                  const traj_similarity_spec_t& spec,
                                  int k,
                                  const std::string& symmetrize,
                                  const std::string& knn_select,
                                  int n_threads);

/**
 * @brief Placeholder for future C++ Leiden implementation.
 *
 * @details
 * This function is declared but intentionally not implemented at this time.
 * The current workflow runs Leiden in R via igraph.
 */
std::vector<int> leiden_cluster(const traj_graph_t& graph,
                                double resolution,
                                int n_iter,
                                int seed);
