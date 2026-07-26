#ifndef PARTITION_GRAPH_HPP
#define PARTITION_GRAPH_HPP

#include <vector>
#include <map>
#include <set>
#include <string>
#include <algorithm>
#include <cmath>

/**
 * @brief Compute aggregated graph from vertex partition
 *
 * Given a graph and a partition of its vertices, constructs a new graph where
 * vertices represent partition cells and edges capture connectivity between cells.
 *
 * @param adj_list Adjacency list representation (0-based indexing)
 * @param weight_list Edge weights corresponding to adjacency list
 * @param partition Partition assignment for each vertex (arbitrary integer labels)
 * @param weight_type Type of edge weight: "count", "jaccard", or "normalized"
 * @param adj_list_out Output adjacency list for partition graph
 * @param weight_list_out Output edge weights for partition graph
 */
void compute_partition_graph(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    const std::vector<int>& partition,
    const std::string& weight_type,
    std::vector<std::vector<int>>& adj_list_out,
    std::vector<std::vector<double>>& weight_list_out
);

#endif // PARTITION_GRAPH_HPP
