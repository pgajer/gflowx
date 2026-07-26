#ifndef NERVE_COMPLEX_HPP
#define NERVE_COMPLEX_HPP

#include "set_wgraph.hpp"
#include "edge_info.hpp"
#include "nerve_cx_spectral_filter.hpp"

#include <cstddef>
#include <vector>
#include <unordered_map>
#include <set>
#include <memory>
#include <functional>

using std::size_t;

/**
 * @brief Hash function for vectors of size_t (simplex representation)
 */
struct simplex_hash_t {
	std::size_t operator()(const std::vector<size_t>& simplex) const {
		std::size_t hash = 0;
		for (auto v : simplex) {
			hash ^= std::hash<size_t>{}(v) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
		}
		return hash;
	}
};

/**
 * @brief Comparator for simplices (vectors of size_t)
 */
struct simplex_compare_t {
	bool operator()(const std::vector<size_t>& a, const std::vector<size_t>& b) const {
		if (a.size() != b.size()) return a.size() < b.size();
		return a < b;
	}
};

/**
 * @struct simplex_info_t
 * @brief Information associated with a simplex
 */
struct simplex_info_t {
	double weight;   ///< Weight of the simplex
	std::vector<size_t> vertices; ///< Vertices forming the simplex (sorted)
	size_t dimension;///< Dimension of the simplex (|vertices| - 1)

	simplex_info_t() : weight(0.0), dimension(0) {}

	simplex_info_t(const std::vector<size_t>& verts, double w = 1.0)
		: weight(w), vertices(verts), dimension(verts.size() - 1) {
		// Ensure vertices are sorted
		std::sort(vertices.begin(), vertices.end());
	}
};

/**
 * @brief Function type for weight calculation
 * @param vertices The vertices of the simplex
 * @param coords The coordinates of all vertices in the complex
 * @param function_values Function values at vertices (optional)
 */
using weight_function_t = std::function<double(
											const std::vector<size_t>&,
											const std::vector<std::vector<double>>&,
											const std::vector<double>&
											)>;

/**
 * @struct nerve_complex_t
 * @brief Represents a simplicial complex derived from a kNN covering
 */
struct nerve_complex_t {
	// Core data structure: maps dimensions to simplices of that dimension
	std::vector<std::unordered_map<std::vector<size_t>, simplex_info_t, simplex_hash_t>> simplices;

	// Cached k-neighborhood information
	std::vector<std::unordered_set<size_t>> knn_idx;           // Vertex indices in each neighborhood
	std::vector<std::unordered_map<size_t, double>> knn_dist;  // Distances to each neighbor vertex

	// Original point coordinates
	std::vector<std::vector<double>> coordinates;

	// Function values at vertices
	std::vector<double> function_values;

	// The 1-skeleton as a set_wgraph_t
	std::unique_ptr<set_wgraph_t> skeleton;

	// Maximum dimension of simplices in the complex
	size_t max_dimension;

	// Weight calculators for each dimension
	std::vector<weight_function_t> weight_calculators;

	// Default constructor
	nerve_complex_t() : max_dimension(0) {}

	/**
	 * @brief Initialize the complex from point coordinates
	 * @param coords Coordinates of the vertices
	 * @param k Number of nearest neighbors for the covering
	 * @param max_dim Maximum dimension of simplices to compute
	 */
	nerve_complex_t(
		const std::vector<std::vector<double>>& coords,
		size_t k,
		size_t max_dim
		) : coordinates(coords), max_dimension(max_dim) {
		simplices.resize(max_dim + 1);
		initialize_from_knn(k);
	}

	/**
	 * @brief Set function values at vertices
	 */
	void set_function_values(const std::vector<double>& values) {
		function_values = values;
	}

	/**
	 * @brief Initialize the nerve complex from a kNN covering
	 */
	void initialize_from_knn(size_t k);

	/**
	 * @brief Set a weight function for a specific dimension
	 */
	void set_weight_function(size_t dim, weight_function_t weight_func);

	/**
	 * @brief Update weights of all simplices
	 */
	void update_weights();

	/**
	 * @brief Get number of vertices
	 */
	size_t num_vertices() const {
		return coordinates.size();
	}

	/**
	 * @brief Get number of simplices of a given dimension
	 */
	size_t num_simplices(size_t dim) const {
		if (dim >= simplices.size()) return 0;
		return simplices[dim].size();
	}

	/**
	 * @brief Get total number of simplices
	 */
	size_t total_simplices() const {
		size_t total = 0;
		for (const auto& dim_simplices : simplices) {
			total += dim_simplices.size();
		}
		return total;
	}

	/**
	 * @brief Construct the boundary operator matrix for dimension p
	 */
	Eigen::SparseMatrix<double> boundary_operator(size_t p) const;

	/**
	 * @brief Construct the p-th Hodge Laplacian
	 */
	Eigen::SparseMatrix<double> hodge_laplacian(size_t p) const;

	/**
	 * @brief Construct the full Laplacian operator
	 * @param dim_weights Weights for each dimension's contribution
	 */
	Eigen::SparseMatrix<double> full_laplacian(
		const std::vector<double>& dim_weights
		) const;

	/**
	 * @brief Solve the full Laplacian equation
	 * @param lambda Regularization parameter
	 * @param dim_weights Weights for each dimension's contribution
	 * @return Extended function values
	 */
	std::vector<double> solve_full_laplacian(
		double lambda,
		const std::vector<double>& dim_weights
		) const;

	// Helper methods for constructing the complex
	void add_simplex(const std::vector<size_t>& simplex, double weight = 1.0);
	// bool check_nerve_condition(const std::vector<size_t>& simplex, size_t k) const;
	void compute_k_neighborhoods(size_t k);


	/**
	 * @brief Check if two simplices can form a higher-dimensional simplex
	 * @param simplex1 First simplex
	 * @param simplex2 Second simplex
	 * @param target_dim Target dimension plus 1 (number of vertices)
	 * @return True if they can form a higher-dimensional simplex
	 */
	bool can_form_higher_simplex(
		const std::vector<size_t>& simplex1,
		const std::vector<size_t>& simplex2,
		size_t target_dim
		) const;

	/**
	 * @brief Merge two simplices into a higher-dimensional simplex
	 * @param simplex1 First simplex
	 * @param simplex2 Second simplex
	 * @return The merged simplex
	 */
	std::vector<size_t> merge_simplices(
		const std::vector<size_t>& simplex1,
		const std::vector<size_t>& simplex2
		) const;

	Eigen::SparseMatrix<double> compute_B1_principled(double weight) const;

	Eigen::SparseMatrix<double> principled_full_laplacian(
		const std::vector<double>& dim_weights) const;

	nerve_cx_spectral_filter_t
	nerve_cx_spectral_filter(
		const std::vector<double>& y,
		laplacian_type_t laplacian_type,
		filter_type_t filter_type,
		size_t laplacian_power,
		const std::vector<double>& dim_weights,
		kernel_params_t kernel_params,
		size_t n_evectors,
		size_t n_candidates,
		bool log_grid,
		bool with_t_predictions,
		bool verbose
		) const;

};

#endif // NERVE_COMPLEX_HPP
