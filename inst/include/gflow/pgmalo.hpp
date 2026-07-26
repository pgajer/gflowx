#ifndef PGMALO_HPP_
#define PGMALO_HPP_

#include <vector>
#include "path_graphs.hpp"

struct pgmalo_t {
	// Graph structures
	std::vector<path_graph_plm_t> graphs;
	path_graph_plm_t graph;

	// Parameter values
	std::vector<int> h_values;
	int opt_h;
	int opt_h_idx;

	// Error metrics
	std::vector<double> h_cv_errors;
	std::vector<double> true_errors;

	// Predictions
	std::vector<double> predictions;
	std::vector<double> local_predictions;

	std::vector<std::vector<double>> h_predictions;  ///< predictions for each h value

	// Bootstrap results
	std::vector<double> bb_predictions;
	std::vector<double> ci_lower;
	std::vector<double> ci_upper;

	// Commponent checks
	bool has_bootstrap_results() const { return !bb_predictions.empty(); }
	bool has_true_errors() const { return !true_errors.empty(); }
};

std::pair<std::vector<double>, std::vector<double>> spgmalo(
	const path_graph_plm_t& path_graph,
	const std::vector<double>& y,
	const std::vector<double>& weights,
	int kernel_type,
	int max_distance_deviation,
	double dist_normalization_factor,
	double epsilon,
	bool verbose = false);

#endif // PGMALO_HPP_
