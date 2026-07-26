#ifndef GRAPH_LAPLACIAN_FOURIER_H_
#define GRAPH_LAPLACIAN_FOURIER_H_

#include <Eigen/Dense>

/**
 * @brief Data structure containing graph representation and its spectral properties
 *
 * Stores the weighted adjacency matrix, graph Laplacian, and its eigendecomposition
 * for a chain graph constructed from ordered data points.
 */
struct graph_data_t {
	Eigen::MatrixXd adjacency_matrix;  ///< Weighted adjacency matrix
	Eigen::MatrixXd laplacian_matrix;  ///< Graph Laplacian matrix (L = D - A)
	Eigen::MatrixXd eigenvectors;      ///< Eigenvectors of the Laplacian matrix
	Eigen::VectorXd eigenvalues;       ///< Eigenvalues of the Laplacian matrix
};

/**
 * @brief Data structure containing both graph data and signal's Fourier coefficients
 *
 * Combines the graph representation with the spectral coefficients of a signal
 * analyzed on this graph structure.
 */
struct spectral_result_t {
	graph_data_t graph_data;           ///< Graph representation and its properties
	Eigen::VectorXd fourier_coeffs;    ///< Fourier coefficients of the analyzed signal
};

graph_data_t create_chain_graph(const Eigen::VectorXd& x);
void compute_laplacian_eigendecomposition(graph_data_t& graph);
Eigen::VectorXd compute_fourier_coefficients(const Eigen::VectorXd& y,const graph_data_t& graph);
spectral_result_t analyze_signal(const Eigen::VectorXd& x,
								 const Eigen::VectorXd& y,
								 const std::optional<graph_data_t>& existing_graph = std::nullopt);

#endif // GRAPH_LAPLACIAN_FOURIER_H_
