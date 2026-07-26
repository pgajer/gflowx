#ifndef PGMALOG_H_
#define PGMALOG_H_

#include <vector>
#include "path_graphs.hpp"

std::pair<std::vector<double>, std::vector<double>> pgmalog(
    const path_graph_plm_t& path_graph,
    const std::vector<double>& y,
    const std::vector<double>& weights,
    int kernel_type,
    int max_distance_deviation,
    double dist_normalization_factor,
    double epsilon,
    bool verbose = false);

#endif // PGMALOG_H_
