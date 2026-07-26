#ifndef SIMPLEX_WEIGHTS_HPP
#define SIMPLEX_WEIGHTS_HPP

#include "nerve_cx.hpp"

#include <vector>
#include <cmath>
#include <algorithm>
#include <cstddef>

using std::size_t;

namespace simplex_weights {

/**
 * @brief Uniform weights (all simplices have weight 1.0)
 */
	inline double uniform_weight(
		const std::vector<size_t>& vertices,
		const std::vector<std::vector<double>>& coords,
		const std::vector<double>& values
		) {
		(void)vertices;
		(void)coords;
		(void)values;

		return 1.0;
	}

/**
 * @brief Inverse distance weight
 */
	inline double inverse_distance_weight(
		const std::vector<size_t>& vertices,
		const std::vector<std::vector<double>>& coords,
		const std::vector<double>& values
		) {

		(void)values;

		double sum_squared_dist = 0.0;
		size_t count = 0;

		for (size_t i = 0; i < vertices.size(); ++i) {
			for (size_t j = i + 1; j < vertices.size(); ++j) {
				const auto& p1 = coords[vertices[i]];
				const auto& p2 = coords[vertices[j]];

				double dist_sq = 0.0;
				for (size_t d = 0; d < p1.size(); ++d) {
					double diff = p1[d] - p2[d];
					dist_sq += diff * diff;
				}

				sum_squared_dist += dist_sq;
				count++;
			}
		}

		return count > 0 ? 1.0 / (sum_squared_dist / count + 1e-10) : 1.0;
	}

/**
 * @brief Gaussian kernel weight
 */
	inline double gaussian_weight(
		const std::vector<size_t>& vertices,
		const std::vector<std::vector<double>>& coords,
		const std::vector<double>& values,
		double sigma = 1.0
		) {

		(void)values;

		double sum_squared_dist = 0.0;
		size_t count = 0;

		for (size_t i = 0; i < vertices.size(); ++i) {
			for (size_t j = i + 1; j < vertices.size(); ++j) {
				const auto& p1 = coords[vertices[i]];
				const auto& p2 = coords[vertices[j]];

				double dist_sq = 0.0;
				for (size_t d = 0; d < p1.size(); ++d) {
					double diff = p1[d] - p2[d];
					dist_sq += diff * diff;
				}

				sum_squared_dist += dist_sq;
				count++;
			}
		}

		double avg_squared_dist = count > 0 ? sum_squared_dist / count : 0.0;
		return std::exp(-avg_squared_dist / (2.0 * sigma * sigma));
	}

/**
 * @brief Volume-based weight
 */
	inline double volume_weight(
		const std::vector<size_t>& vertices,
		const std::vector<std::vector<double>>& coords,
		const std::vector<double>& values,
		double alpha = 1.0
		) {

		(void)values;

		// For efficiency, approximating volume by average edge length to power of dimension
		double avg_edge_length = 0.0;
		size_t edge_count = 0;

		for (size_t i = 0; i < vertices.size(); ++i) {
			for (size_t j = i + 1; j < vertices.size(); ++j) {
				const auto& p1 = coords[vertices[i]];
				const auto& p2 = coords[vertices[j]];

				double dist = 0.0;
				for (size_t d = 0; d < p1.size(); ++d) {
					double diff = p1[d] - p2[d];
					dist += diff * diff;
				}
				dist = std::sqrt(dist);

				avg_edge_length += dist;
				edge_count++;
			}
		}

		if (edge_count > 0) {
			avg_edge_length /= edge_count;
			double volume = std::pow(avg_edge_length, vertices.size() - 1);
			return std::pow(volume + 1e-10, -alpha);
		}

		return 1.0;
	}

/**
 * @brief Function gradient-based weight
 */
	inline double gradient_weight(
		const std::vector<size_t>& vertices,
		const std::vector<std::vector<double>>& coords,
		const std::vector<double>& values,
		double gamma = 1.0
		) {

		(void)coords;
		(void)values;

		double sum_squared_diff = 0.0;
		size_t count = 0;

		for (size_t i = 0; i < vertices.size(); ++i) {
			for (size_t j = i + 1; j < vertices.size(); ++j) {
				double val_diff = values[vertices[i]] - values[vertices[j]];
				sum_squared_diff += val_diff * val_diff;
				count++;
			}
		}

		double avg_squared_diff = count > 0 ? sum_squared_diff / count : 0.0;
		return std::pow(avg_squared_diff + 1e-10, -gamma);
	}

/**
 * @brief Factory function to create weight function with parameters
 */
	template<typename WeightFunc>
	weight_function_t create_weight_function(WeightFunc func) {
		return [func](
			const std::vector<size_t>& vertices,
			const std::vector<std::vector<double>>& coords,
			const std::vector<double>& values) {
			return func(vertices, coords, values);
		};
	}

/**
 * @brief Factory for gaussian weight with specified sigma
 */
	inline weight_function_t create_gaussian_weight(double sigma) {
		return [sigma](
			const std::vector<size_t>& vertices,
			const std::vector<std::vector<double>>& coords,
			const std::vector<double>& values) {
			return gaussian_weight(vertices, coords, values, sigma);
		};
	}

/**
 * @brief Factory for volume weight with specified alpha
 */
	inline weight_function_t create_volume_weight(double alpha) {
		return [alpha](
			const std::vector<size_t>& vertices,
			const std::vector<std::vector<double>>& coords,
			const std::vector<double>& values) {
			return volume_weight(vertices, coords, values, alpha);
		};
	}

/**
 * @brief Factory for gradient weight with specified gamma
 */
	inline weight_function_t create_gradient_weight(double gamma) {
		return [gamma](
			const std::vector<size_t>& vertices,
			const std::vector<std::vector<double>>& coords,
			const std::vector<double>& values) {
			return gradient_weight(vertices, coords, values, gamma);
		};
	}

} // namespace simplex_weights

#endif // SIMPLEX_WEIGHTS_HPP
