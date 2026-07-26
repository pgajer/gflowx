#include "riem_dcx.hpp"

#include <vector>
#include <array>
#include <utility>
#include <algorithm>
#include <memory>
#include <cmath>

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

#include <R_ext/Error.h>

/**
 * @brief Build symmetrized (normalized) vertex Laplacian for spectral analysis
 *
 * @param g Metric family containing vertex masses M[0] and edge information
 *
 * @details
 * This function constructs the mass-symmetrized vertex Laplacian L0_mass_sym:
 *     L0_mass_sym = M[0]^{-1/2} L[0] M[0]^{-1/2}.
 * This is NOT the classical degree-normalized Laplacian.
 *
 * MATHEMATICAL FORMULATION:
 *
 * Given a weighted graph with:
 *   - Vertices V = {1, ..., n}
 *   - Edges E with conductances c_ij > 0
 *   - Vertex masses m_i > 0 (from metric M[0])
 *
 * We construct:
 * 1. Standard graph Laplacian: L_div
 * 2. Mass matrix diagonal: M = diag(m_1, ..., m_n)
 * 3. Normalized Laplacian: L_sym = M^{-1/2} L_div M^{-1/2}
 *
 * The normalized Laplacian has the property that for functions f on vertices:
 *   <f, L_mass_sym f>_{ℓ²} = (1/2) Σ_{edges (i,j)} c_ij (f_i/√m_i - f_j/√m_j)²
 *
 * This measures the "smoothness" of f with respect to the geometry.
 */
void laplacian_bundle_t::build_L0_mass_sym_if_needed(const metric_family_t& g) {

	// Need at least 2 dimensions (vertices and edges)
	if (B.size() < 2) {
		L0_mass_sym.resize(0, 0);
		return;
	}

	// Need a non-empty boundary operator B[1] (vertex-edge incidence)
	if (B[1].nonZeros() == 0) {
		L0_mass_sym.resize(0, 0);
		return;
	}

	// Need vertex metric M[0]
	if (g.M.empty() || g.M[0].rows() == 0) {
		L0_mass_sym.resize(0, 0);
		return;
	}

	// Need conductances c1 for each edge
	// c1[e] represents the "strength" or "weight" of edge e

	// In the DEC interpretation used in this package, c1[e] is derived from the
	// edge 1-form mass as $c_e = 1/m_1(e)$ (up to a small floor), so larger
	// edge mass implies weaker penalization of variation across that edge.

	if (c1.size() != B[1].cols()) {
		L0_mass_sym.resize(0, 0);
		return;
	}

	// ========== EXTRACT DIMENSIONS ==========

	const spmat_t& B1 = B[1];           // Boundary operator ∂₁: C₁ → C₀
	const Eigen::Index n = B1.rows();   // Number of vertices
	const Eigen::Index m = B1.cols();   // Number of edges

	// ========== BUILD STANDARD GRAPH LAPLACIAN L_div ==========

	// The graph Laplacian has the form:
	//   (L_div)_ii = Σ_{j∼i} c_ij  (sum of conductances of incident edges)
	//   (L_div)_ij = -c_ij          (negative conductance if i,j are adjacent)
	//   (L_div)_ij = 0              (zero if i,j not adjacent)
	//
	// This can be written as: L_div = B₁ C B₁ᵀ
	// where C = diag(c₁, ..., c_m) is the diagonal matrix of edge conductances.

	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(4 * m);  // Each edge contributes 4 entries to L_div

	// Iterate over each edge
	for (int e = 0; e < m; ++e) {
		// ===== DECODE EDGE ENDPOINTS FROM BOUNDARY OPERATOR =====

		// The boundary operator B₁ has exactly 2 non-zero entries per column:
		//   B₁[i, e] = -1  (tail vertex of edge e)
		//   B₁[j, e] = +1  (head vertex of edge e)
		//
		// This encodes the oriented edge as: e = [i → j]

		int i = -1;  // Tail vertex (receives -1 in B₁)
		int j = -1;  // Head vertex (receives +1 in B₁)

		for (spmat_t::InnerIterator it(B1, e); it; ++it) {
			if (it.value() > 0) {
				j = it.row();  // Positive coefficient → head vertex
			} else if (it.value() < 0) {
				i = it.row();  // Negative coefficient → tail vertex
			}
		}

		// Skip malformed edges (should not happen in valid complex)
		if (i < 0 || j < 0) continue;

		// ===== GET EDGE CONDUCTANCE =====

		double c = c1[e];

		// Skip edges with non-positive conductance (degenerate or removed)
		if (c <= 0) continue;

		// ===== ADD CONTRIBUTIONS TO LAPLACIAN =====

		// For an edge (i,j) with conductance c, the Laplacian contributions are:
		//
		//   L[i,i] ← L[i,i] + c    (degree term at vertex i)
		//   L[j,j] ← L[j,j] + c    (degree term at vertex j)
		//   L[i,j] ← L[i,j] - c    (off-diagonal coupling)
		//   L[j,i] ← L[j,i] - c    (symmetric)
		//
		// This ensures L_div is symmetric and positive semi-definite.

		triplets.emplace_back(i, i,  c);   // Diagonal: vertex i degree
		triplets.emplace_back(j, j,  c);   // Diagonal: vertex j degree
		triplets.emplace_back(i, j, -c);   // Off-diagonal: coupling
		triplets.emplace_back(j, i, -c);   // Off-diagonal: coupling (symmetric)
	}

	// Assemble sparse matrix from triplets
	// Note: setFromTriplets automatically sums contributions to same (row, col)
	spmat_t L0_div(n, n);
	L0_div.setFromTriplets(triplets.begin(), triplets.end());

	// ========== BUILD NORMALIZATION MATRIX D^{-1/2} ==========

	// The vertex masses m_i are stored on the diagonal of M[0].
	// We construct the diagonal matrix:
	//   D^{-1/2} = diag(1/√m₁, ..., 1/√m_n)
	//
	// This will be used to normalize the Laplacian.

	spmat_t Dm05(n, n);
	Dm05.reserve(Eigen::VectorXi::Constant(n, 1));  // Reserve space for diagonal

	for (Eigen::Index i = 0; i < n; ++i) {
		// Extract vertex mass (with safety threshold to avoid division by zero)
		double m = std::max(g.M[0].coeff(i, i), 1e-15);

		// Insert 1/√m_i on the diagonal
		Dm05.insert(i, i) = 1.0 / std::sqrt(m);
	}
	Dm05.makeCompressed();

	// ========== COMPUTE NORMALIZED LAPLACIAN ==========

	// The normalized (symmetrized) Laplacian is:
	//   L_sym = D^{-1/2} L_div D^{-1/2}
	//
	// This has several nice properties:
	// 1. Symmetric: L_sym = L_symᵀ
	// 2. Eigenvalues in [0, 2]
	// 3. Smallest eigenvalue is 0 (with eigenvector ∝ √masses)
	// 4. Number of 0 eigenvalues = number of connected components
	// 5. Well-conditioned for iterative solvers
	//
	// The normalization makes vertices with large mass contribute less
	// to the energy, balancing the influence of high-degree and low-degree vertices.

	try {
		L0_mass_sym = Dm05 * L0_div * Dm05;
	} catch (const std::bad_alloc&) {
		L0_mass_sym.resize(0, 0);
		Rf_warning("build_L0_mass_sym_if_needed: memory allocation failed while forming "
		           "M^{-1/2} L M^{-1/2}. Falling back to L[0] for spectral computations.");
	}
}

/**
 * @brief Assemble all Hodge Laplacians from boundary maps and metric
 * @param g Metric family defining the inner products
 *
 * Computes L[p] for all dimensions using the Hodge decomposition formula.
 * Each Laplacian combines contributions from boundary maps above and below
 * dimension p, weighted by the metric.
 */
void laplacian_bundle_t::assemble_all(const metric_family_t& g) {
	const int P = static_cast<int>(g.M.size()) - 1;
	L.assign(P + 1, spmat_t());

	// L[1] can be prohibitively large (n_edges x n_edges) on big graphs.
	// Regression uses L[0]; skip higher-order assembly beyond this threshold.
	const Eigen::Index HIGHER_ORDER_LAPLACIAN_MAX_DIM = 50000;

	for (int p = 0; p <= P; ++p) {
		const Eigen::Index n = g.M[p].rows();

		if (p >= 1 && n > HIGHER_ORDER_LAPLACIAN_MAX_DIM) {
			L[p] = spmat_t(n, n);
			Rf_warning("assemble_all: skipping L[%d] assembly for size=%ld to avoid excessive "
			           "memory use; L[0] remains available.", p, static_cast<long>(n));
			continue;
		}

		spmat_t up, down;
		bool has_up = false, has_down = false;

		// Up term: ∂_{p+1} ∂_{p+1}^* = B[p+1] M[p+1]^{-1} B[p+1]^T M[p]
		if (p + 1 <= P && B.size() > static_cast<size_t>(p + 1) && B[p+1].nonZeros() > 0) {
			spmat_t temp = B[p+1] * left_apply_Minv(g, p + 1, B[p+1].transpose());
			up = temp * g.M[p];
			has_up = true;
		}

		// Down term: ∂_p^* ∂_p = M[p]^{-1} B[p]^T M[p-1] B[p]
		if (p >= 1 && B.size() > static_cast<size_t>(p) && B[p].nonZeros() > 0) {
			spmat_t temp = B[p].transpose() * g.M[p - 1] * B[p];
			down = left_apply_Minv(g, p, temp);
			has_down = true;
		}

		// Assemble Laplacian
		if (has_up && has_down) {
			L[p] = up + down;
		} else if (has_up) {
			L[p] = up;
		} else if (has_down) {
			L[p] = down;
		} else {
			L[p] = spmat_t(n, n);
		}
	}

	build_L0_mass_sym_if_needed(g);
}

/**
 * @brief Apply inverse metric from the right to boundary map
 * @param g Metric family
 * @param p Target dimension
 * @param Bp Boundary map B[p]
 * @return B[p] M[p]^{-1} (as a sparse matrix)
 *
 * Helper for assembling Laplacians. Assumes diagonal metric for efficiency.
 */
spmat_t laplacian_bundle_t::right_apply_Minv(const metric_family_t& g, int p, const spmat_t& Bp) const {
	spmat_t A(Bp.rows(), Bp.cols());
	A.reserve(Bp.nonZeros());
	for (int k=0; k<Bp.outerSize(); ++k) {
		for (spmat_t::InnerIterator it(Bp,k); it; ++it) {
			double d = std::max(g.M[p].coeff(it.row(), it.row()), 1e-15);
			A.insert(it.row(), it.col()) = it.value() / d;
		}
	}
	A.makeCompressed();
	return A;
}

/**
 * @brief Apply inverse metric from the left to a matrix
 * @param g Metric family
 * @param p Dimension index
 * @param T Input matrix
 * @return M[p]^{-1} T (as a sparse matrix)
 *
 * Helper for assembling Laplacians. Assumes diagonal metric for efficiency.
 */
spmat_t laplacian_bundle_t::left_apply_Minv(
	const metric_family_t& g,
	int p, const spmat_t& T
	) const {
	spmat_t D(T.rows(), T.cols());
	D.reserve(T.nonZeros());
	for (int k=0; k<T.outerSize(); ++k) {
		for (spmat_t::InnerIterator it(T,k); it; ++it) {
			double d = std::max(g.M[p].coeff(it.row(), it.row()), 1e-15);
			D.insert(it.row(), it.col()) = it.value() / d;
		}
	}
	D.makeCompressed();
	return D;
}

/**
 * @brief Apply heat kernel exp(-t L[p]) to a vector
 * @param p Dimension index
 * @param x Input vector
 * @param t Diffusion time
 * @param m_steps Number of backward Euler steps for approximation
 * @return Approximation to exp(-t L[p]) x
 *
 * Uses m-step backward Euler time stepping to approximate the heat
 * semigroup. Larger m_steps give better approximations at higher cost.
 */
vec_t laplacian_bundle_t::heat_apply(
	int p,
	const vec_t& x,
	double t,
	int m_steps=8
	) const {
	if (p<0 || p>=static_cast<int>(L.size())) Rf_error("heat_apply: bad p");
	if (t<=0) return x;
	spmat_t A = L[p];
	spmat_t I(A.rows(), A.cols()); I.setIdentity();
	double h = t / std::max(1, m_steps);
	A = I + h * A;
	Eigen::ConjugateGradient<spmat_t, Eigen::Lower|Eigen::Upper> cg;
	cg.setTolerance(1e-8); cg.setMaxIterations(1000); cg.compute(A);
	vec_t y = x;
	for (int s=0; s<m_steps; ++s) y = cg.solve(y);
	return y;
}

/**
 * @brief Solve Tikhonov regularization problem (I + eta L[p]) x = b
 * @param p Dimension index
 * @param b Right-hand side vector
 * @param eta Regularization parameter
 * @param tol Conjugate gradient tolerance
 * @param maxit Maximum CG iterations
 * @return Solution vector x
 *
 * Standard ridge regression or Tikhonov regularization against the
 * Laplacian. Used for signal smoothing and fitting with geometric
 * regularization on the complex.
 */
vec_t laplacian_bundle_t::tikhonov_solve(
	int p,
	const vec_t& b,
	double eta,
	double tol=1e-8,
	int maxit=1000
	) const {
	if (p<0 || p>=static_cast<int>(L.size())) Rf_error("tikhonov: bad p");
	spmat_t A = L[p];
	spmat_t I(A.rows(), A.cols()); I.setIdentity();
	A = I + eta * A;
	Eigen::ConjugateGradient<spmat_t, Eigen::Lower|Eigen::Upper> cg;
	cg.setTolerance(tol); cg.setMaxIterations(maxit); cg.compute(A);
	return cg.solve(b);
}
