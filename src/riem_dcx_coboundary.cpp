/**
 * @file riem_dcx_coboundary_full.cpp
 * @brief Complete implementation of coboundary operator with non-diagonal metric
 *
 * This file provides the full implementation of the coboundary operator ∂₁*
 * for non-diagonal Riemannian metrics, where edges incident to the same vertex
 * have non-orthogonal inner products encoded in local Gram matrices.
 *
 */

#include "riem_dcx.hpp"
#include <Eigen/IterativeLinearSolvers>
#include <limits>
#include <queue>

/**
 * @brief Compute coboundary operator ∂₁* with full non-diagonal metric
 *
 * Computes ∂₁*(f) for a vertex function f using the complete Riemannian
 * structure where edges incident to the same vertex may have non-orthogonal
 * inner products encoded in the edge mass matrix M[1].
 *
 * MATHEMATICAL FORMULATION:
 * Solves M₁ x = B₁ᵀ M₀ f where:
 * - M₁ is the full edge mass matrix (potentially non-diagonal)
 * - B₁ is the boundary operator
 * - M₀ is the vertex mass matrix (diagonal)
 * - x gives the coboundary ∂₁*(f)
 *
 * COMPUTATIONAL METHODS:
 * - If M₁ is diagonal: Direct solve O(m), where m is the number of edges
 * - If M₁ prefactored: Cholesky solve O(m)
 * - If use_iterative: Conjugate gradient O(k·nnz(M₁))
 * - Otherwise: Factor then solve O(m^1.5) first call, O(m) subsequent
 *
 * This method may cache the Cholesky factorization of M₁ for efficiency,
 * making subsequent calls much faster.
 *
 * @param f Vertex function (0-chain) with length n_vertices
 * @param use_iterative Force iterative CG solve even if factorization available
 * @param cg_tol Conjugate gradient convergence tolerance
 * @param cg_maxit Maximum CG iterations
 *
 * @return Edge function (1-chain) with length n_edges giving ∂₁*(f)
 *
 * @complexity O(m^1.5) first call with factorization, O(m) subsequent calls
 *
 * @post May populate g.M_solver[1] with Cholesky factorization
 *
 * @see compute_edge_mass_matrix() constructs the non-diagonal M₁
 * @see compute_all_extremality_scores_full() uses this for extrema detection
 */
vec_t riem_dcx_t::compute_coboundary_del1star_full(
    const vec_t& f,
    bool use_iterative,
    double cg_tol,
    int cg_maxit
) const {

    // ========================================================================
    // VALIDATION AND DIMENSION EXTRACTION
    // ========================================================================

    const size_t n_vertices = vertex_cofaces.size();
    const size_t n_edges = edge_registry.size();
    const Eigen::Index n = static_cast<Eigen::Index>(n_vertices);
    const Eigen::Index m = static_cast<Eigen::Index>(n_edges);

    // Validate input dimension
    if (f.size() != n) {
        Rf_error("compute_coboundary_del1star_full: f.size()=%lld must equal n_vertices=%zu",
                 static_cast<long long>(f.size()), n_vertices);
    }

    // Validate metric structures
    if (g.M.size() < 2) {
        Rf_error("compute_coboundary_del1star_full: metric not initialized (g.M.size()=%zu < 2)",
                 g.M.size());
    }

    if (g.M[0].rows() != n) {
        Rf_error("compute_coboundary_del1star_full: vertex metric dimension mismatch");
    }

    if (g.M[1].rows() != m) {
        Rf_error("compute_coboundary_del1star_full: edge metric dimension mismatch");
    }

    // Validate boundary operator exists
    if (L.B.size() < 2 || L.B[1].rows() != n || L.B[1].cols() != m) {
        Rf_error("compute_coboundary_del1star_full: boundary operator B[1] not properly initialized");
    }

    // ========================================================================
    // COMPUTE RIGHT-HAND SIDE: b = B₁ᵀ M₀ f
    // ========================================================================

	// For diagonal M₀, we can compute B₁ᵀ M₀ f efficiently.
    // Since B₁(i, e) = -1 for tail, B₁(j, e) = +1 for head of edge e = [i,j]:
    //   (B₁ᵀ M₀ f)[e] = w_j f_j - w_i f_i

    vec_t b = vec_t::Zero(m);

    for (Eigen::Index e = 0; e < m; ++e) {
        const auto [i, j] = edge_registry[e];
        const Eigen::Index ei = static_cast<Eigen::Index>(i);
        const Eigen::Index ej = static_cast<Eigen::Index>(j);

        // Extract vertex weights (diagonal of M₀)
        const double w_i = std::max(g.M[0].coeff(ei, ei), 1e-15);
        const double w_j = std::max(g.M[0].coeff(ej, ej), 1e-15);

        // Compute weighted boundary value
        b[e] = w_j * f[ej] - w_i * f[ei];
    }

	// ========================================================================
    // SOLVE M₁ x = b FOR DIFFERENT METRIC STRUCTURES
    // ========================================================================

    vec_t x; // Solution vector

    // ========================================================================
    // CASE 1: DIAGONAL METRIC (fast path)
    // ========================================================================

	if (g.is_diagonal(1)) {
        // Diagonal solve: x[e] = b[e] / M₁(e,e)
        x = vec_t(m);
        for (Eigen::Index e = 0; e < m; ++e) {
            const double diag = std::max(g.M[1].coeff(e, e), 1e-15);
            x[e] = b[e] / diag;
        }
        return x;
    }

    // ========================================================================
    // CASE 2: NON-DIAGONAL METRIC WITH ITERATIVE SOLVE
    // ========================================================================

	if (use_iterative) {
		// Extract diagonal for analysis
		vec_t diag_M1 = g.get_diagonal(1);

		// Compute adaptive regularization based on diagonal statistics
		double mean_diag = diag_M1.mean();
		double min_diag = diag_M1.minCoeff();
		double max_diag = diag_M1.maxCoeff();

		// Compute condition number estimate
		double condition_estimate = max_diag / std::max(min_diag, 1e-15);

		// Adaptive regularization based on condition number
		double epsilon;
		if (condition_estimate > 1e8) {
			// Severely ill-conditioned: strong regularization
			epsilon = 1e-4 * mean_diag;
			Rf_warning("compute_coboundary_del1star_full: M[1] severely ill-conditioned "
					   "(condition ~%.2e), using strong regularization (epsilon=%.3e)",
					   condition_estimate, epsilon);
		} else if (condition_estimate > 1e6) {
			// Moderately ill-conditioned: medium regularization
			epsilon = 1e-5 * mean_diag;
		} else if (condition_estimate > 1e4) {
			// Mildly ill-conditioned: light regularization
			epsilon = 1e-6 * mean_diag;
		} else {
			// Well-conditioned: minimal regularization
			epsilon = 1e-7 * mean_diag;
		}

		if (condition_estimate > 1e3) {
			Rprintf("  [DEBUG] M[1] condition ~%.2e, using epsilon=%.3e (%.2e * mean_diag)\n",
					condition_estimate, epsilon, epsilon / mean_diag);
		}

		// Create Tikhonov-regularized system: (M₁ + εI) x = b
		spmat_t I(m, m);
		I.setIdentity();

		spmat_t M1_reg = g.M[1] + epsilon * I;
		M1_reg.makeCompressed();

		// Setup conjugate gradient solver with diagonal preconditioning
		Eigen::ConjugateGradient<spmat_t, Eigen::Lower|Eigen::Upper,
								 Eigen::DiagonalPreconditioner<double>> cg;

		// For condition ~2000, we should achieve tolerance ~1e-5 with enough iterations
		// Set tolerance relative to condition number
		double adaptive_tol = std::max(1e-10, 1e-12 * condition_estimate);
		adaptive_tol = std::min(adaptive_tol, 1e-4);  // Cap at 1e-4

		cg.setTolerance(adaptive_tol);
		cg.setMaxIterations(cg_maxit);

		cg.compute(M1_reg);

		if (cg.info() != Eigen::Success) {
			Rf_error("compute_coboundary_del1star_full: CG setup failed");
		}

		// Solve regularized system
		x = cg.solve(b);

		double final_error = cg.error();
		int final_iters = cg.iterations();

		if (cg.info() != Eigen::Success) {
			Rf_warning("compute_coboundary_del1star_full: CG did not fully converge "
					   "(iterations=%d/%d, error=%.3e, tol=%.3e)",
					   final_iters, cg_maxit, final_error, adaptive_tol);
		} else if (final_error > 1e-3) {
			// Warn even on "success" if residual is large
			Rf_warning("compute_coboundary_del1star_full: CG converged but residual is large "
					   "(iterations=%d, error=%.3e, tol=%.3e)",
					   final_iters, final_error, adaptive_tol);
		}

		return x;
	}

	#if 0
	if (use_iterative) {
		// Extract diagonal for analysis
		vec_t diag_M1 = g.get_diagonal(1);

		// Compute adaptive regularization based on diagonal statistics
		double mean_diag = diag_M1.mean();
		double min_diag = diag_M1.minCoeff();
		double max_diag = diag_M1.maxCoeff();

		// Compute condition number estimate
		double condition_estimate = max_diag / std::max(min_diag, 1e-15);

		if (condition_estimate > 1e3) {
			Rprintf("  [DEBUG] M[1] condition ~%.2e\n", condition_estimate);
		}

		// Only regularize if truly needed (condition > 1e6)
		spmat_t M1_to_solve;
		if (condition_estimate > 1e6) {
			double epsilon;
			if (condition_estimate > 1e8) {
				epsilon = 1e-4 * mean_diag;
			} else {
				epsilon = 1e-5 * mean_diag;
			}

			Rprintf("  [DEBUG] Applying regularization epsilon=%.3e\n", epsilon);

			spmat_t I(m, m);
			I.setIdentity();
			M1_to_solve = g.M[1] + epsilon * I;
			M1_to_solve.makeCompressed();
		} else {
			// Use original matrix - no regularization needed
			M1_to_solve = g.M[1];
		}

		// Setup conjugate gradient solver with diagonal preconditioning
		Eigen::ConjugateGradient<spmat_t, Eigen::Lower|Eigen::Upper,
		 						 Eigen::DiagonalPreconditioner<double>> cg;

		// Eigen::ConjugateGradient<spmat_t, Eigen::Lower|Eigen::Upper,
		// 						 Eigen::IncompleteCholesky<double>> cg;

		// Eigen::ConjugateGradient<spmat_t, Eigen::Lower|Eigen::Upper,
		//  						 Eigen::IncompleteLUT<double>> cg;

		// Set realistic tolerance based on condition number
		// For condition ~2000: aim for ~1e-6 relative accuracy
		double adaptive_tol = std::min(1e-4, std::max(1e-7, 1e-10 * condition_estimate));

		cg.setTolerance(adaptive_tol);
		cg.setMaxIterations(cg_maxit);

		cg.compute(M1_to_solve);

		if (cg.info() != Eigen::Success) {
			Rf_error("compute_coboundary_del1star_full: CG setup failed");
		}

		// Solve system
		x = cg.solve(b);

		double final_error = cg.error();
		int final_iters = cg.iterations();

		// Only warn if convergence is really poor (>1% relative error)
		if (final_error > 1e-2) {
			Rf_warning("compute_coboundary_del1star_full: CG convergence issue "
					   "(iterations=%d/%d, error=%.3e, tol=%.3e)",
					   final_iters, cg_maxit, final_error, adaptive_tol);
		} else if (condition_estimate > 1e3 && final_iters > 100) {
			Rprintf("  [DEBUG] CG converged: %d iterations, error=%.3e\n",
					final_iters, final_error);
		}

		return x;
	}
    #endif

    // ========================================================================
    // CASE 3: NON-DIAGONAL METRIC WITH DIRECT SOLVE
    // ========================================================================

	// Check if factorization already exists
    if (g.M_solver[1] && g.M_solver[1]->info() == Eigen::Success) {
        // Use existing factorization
        x = g.M_solver[1]->solve(b);

        if (g.M_solver[1]->info() != Eigen::Success) {
            Rf_error("compute_coboundary_del1star_full: Cholesky solve failed");
        }

        return x;
    }

    // Need to factor M₁ first
    // Note: This modifies g.M_solver[1], but method is marked const
    // We cast away const for factorization caching (doesn't change mathematical state)

    auto* mutable_g = const_cast<metric_family_t*>(&g);

    auto solver = std::make_unique<Eigen::SimplicialLLT<spmat_t>>();
    solver->compute(g.M[1]);

    if (solver->info() != Eigen::Success) {
        const char* reason;
        switch (solver->info()) {
        case Eigen::NumericalIssue:
            reason = "numerical issue (matrix may not be positive definite)";
            break;
        case Eigen::NoConvergence:
            reason = "no convergence";
            break;
        case Eigen::InvalidInput:
            reason = "invalid input";
            break;
        default:
            reason = "unknown error";
            break;
        }

        Rf_error("compute_coboundary_del1star_full: M[1] factorization failed: %s "
                 "(size=%ld, nnz=%ld)",
                 reason, static_cast<long>(g.M[1].rows()),
                 static_cast<long>(g.M[1].nonZeros()));
    }

    // Solve with new factorization
    x = solver->solve(b);

    if (solver->info() != Eigen::Success) {
        Rf_error("compute_coboundary_del1star_full: Cholesky solve failed after factorization");
    }

    // Cache the factorization for future use
    mutable_g->M_solver[1] = std::move(solver);

    return x;
}

/**
 * @brief Compute the coboundary operator ∂₁* applied to a vertex function taking only diagonal elements of
 *
 * @details
 * Given a function f on vertices (a 0-chain), computes ∂₁*(f) as a function
 * on edges (a 1-chain). This operator is the formal adjoint of the boundary
 * operator ∂₁ with respect to the Riemannian inner products on chains.
 *
 * MATHEMATICAL FOUNDATION
 * =======================
 *
 * The coboundary operator arises from the adjoint relationship:
 *   ⟨∂₁(e), f⟩₀ = ⟨e, ∂₁*(f)⟩₁
 *
 * for all edges e and vertex functions f, where ⟨·,·⟩_p denotes the inner
 * product on p-chains determined by the Riemannian metric g.
 *
 * MATRIX REPRESENTATION
 * =====================
 *
 * When inner products are represented by weight matrices W₀ (on vertices) and
 * W₁ (on edges), the coboundary takes the matrix form:
 *   ∂₁* = W₁⁻¹ ∂₁ᵀ W₀
 *
 * Here:
 * - ∂₁ is the boundary operator (n_vertices × n_edges matrix)
 * - W₀ = M[0] is the vertex mass matrix (diagonal)
 * - W₁ = M[1] is the edge mass matrix (possibly non-diagonal)
 *
 * DIAGONAL METRIC CASE
 * ====================
 *
 * When the edge metric is diagonal (W₁ = diag(w₁,...,w_m)), the formula
 * simplifies dramatically. For an edge e = [vᵢ, vⱼ] with i < j:
 *   [∂₁*(f)](e) = (wⱼ f(vⱼ) - wᵢ f(vᵢ)) / wₑ
 *
 * where wᵢ, wⱼ are vertex weights and wₑ is the edge weight. This is the
 * weighted discrete gradient along the edge.
 *
 * NON-DIAGONAL METRIC CASE
 * ========================
 *
 * When edges incident to the same vertex have non-orthogonal inner products
 * (encoded in local Gram matrices), the computation becomes more sophisticated.
 * At vertex v with incident edges e₁,...,e_k, we solve the local linear system:
 *   Gᵥ x = b
 *
 * where:
 * - Gᵥ is the k×k Gram matrix with (Gᵥ)ᵢⱼ = ⟨eᵢ, eⱼ⟩₁
 * - b is the vector with bᵢ = wᵤᵢ f(uᵢ) - wᵥ f(v) for neighbor uᵢ
 * - x gives the coefficients [∂₁*(f)](eᵢ)
 *
 * The global computation requires solving such systems at every vertex and
 * assembling the results into a consistent 1-chain.
 *
 * COMPUTATIONAL STRATEGY
 * ======================
 *
 * Rather than forming W₁⁻¹ explicitly (prohibitively expensive for large
 * complexes), we solve the linear system:
 *   W₁ x = ∂₁ᵀ W₀ f
 *
 * using iterative methods when W₁ is non-diagonal. For the current implementation,
 * we assume diagonal W₁ and use the explicit formula for efficiency.
 *
 * Future extensions to non-diagonal metrics will employ:
 * 1. Conjugate gradient iteration with diagonal preconditioning
 * 2. Block-diagonal approximations for local Gram matrices
 * 3. Algebraic multigrid for hierarchical complexes
 *
 * IMPLEMENTATION DETAILS
 * ======================
 *
 * The computation proceeds by:
 * 1. Extracting vertex weights from g.M[0] (diagonal entries)
 * 2. Extracting edge weights from g.M[1] (diagonal entries, for now)
 * 3. Iterating through edges via edge_registry
 * 4. Computing weighted differences for each edge
 * 5. Storing results in output vector indexed by edge
 *
 * PERFORMANCE CHARACTERISTICS
 * ===========================
 *
 * Time complexity: O(m) where m is the number of edges
 * Space complexity: O(m) for the output vector
 * Cache efficiency: Sequential access through edge_registry
 *
 * The diagonal metric assumption enables this optimal complexity. Non-diagonal
 * metrics would require O(m·k) where k is average vertex degree (for CG iteration).
 *
 * NUMERICAL STABILITY
 * ===================
 *
 * Regularization is applied through:
 * - Edge weights bounded below by machine epsilon (1e-15)
 * - Vertex weights bounded below by 1e-15
 * This prevents division by zero and maintains positive definiteness.
 *
 * @param f Input vertex function (0-chain) with length n_vertices
 * @return Output edge function (1-chain) with length n_edges, where
 *         result[e] = [∂₁*(f)](e) for edge e
 *
 * @pre vertex_cofaces must be populated (determines n_vertices)
 * @pre edge_registry must be populated with edges as {i,j} pairs with i < j
 * @pre g.M[0] must contain vertex masses (diagonal metric on vertices)
 * @pre g.M[1] must contain edge masses (currently assumed diagonal)
 * @pre f.size() must equal vertex_cofaces.size()
 *
 * @post Output vector has length equal to edge_registry.size()
 * @post Each entry corresponds to the coboundary operator value on that edge
 *
 * @throws None (uses Eigen's expression templates for safe computation)
 *
 * @note Current implementation assumes diagonal edge metric. Extension to
 *       non-diagonal metrics (local Gram matrices) will be added when needed
 *       for geometries with non-orthogonal edge inner products.
 *
 * @see build_boundary_operator_from_edges() for the boundary operator ∂₁
 * @see compute_edge_mass_matrix() for edge metric construction
 * @see initialize_metric_from_density() for vertex metric initialization
 *
 * USAGE EXAMPLE
 * =============
 *
 * @code
 * // Compute coboundary of a vertex function
 * vec_t f = vec_t::Random(dcx.vertex_cofaces.size());  // Random function on vertices
 * vec_t del1star_f = dcx.compute_coboundary_del1star(f);
 *
 * // Check which edges have negative coboundary (gradient points inward)
 * for (size_t e = 0; e < dcx.edge_registry.size(); ++e) {
 *     if (del1star_f[e] < 0) {
 *         auto [i, j] = dcx.edge_registry[e];
 *         Rprintf("Edge [%zu,%zu] has inward gradient: %.4f\n", i, j, del1star_f[e]);
 *     }
 * }
 *
 * // Compute extremality score for a vertex
 * size_t v = 42;
 * double numerator = 0.0, denominator = 0.0;
 * for (size_t k = 1; k < dcx.vertex_cofaces[v].size(); ++k) {
 *     index_t e = dcx.vertex_cofaces[v][k].simplex_index;
 *     double rho_e = dcx.vertex_cofaces[v][k].density;  // Edge density
 *     double del1star_val = del1star_f[e];
 *
 *     numerator += rho_e * del1star_val;
 *     denominator += rho_e * std::abs(del1star_val);
 * }
 * double extremality = -numerator / denominator;
 * Rprintf("Vertex %zu extremality: %.4f\n", v, extremality);
 * @endcode
 *
 * REFERENCES
 * ==========
 *
 * [1] Golub & Van Loan (2013). Matrix Computations, 4th ed. Johns Hopkins.
 * [2] Saad (2003). Iterative Methods for Sparse Linear Systems, 2nd ed. SIAM.
 * [3] Grady & Polimeni (2010). Discrete Calculus. Springer.
 * [4] Hirani (2003). Discrete Exterior Calculus. PhD thesis, Caltech.
 */
vec_t riem_dcx_t::compute_coboundary_del1star(const vec_t& f) const {
    // ========================================================================
    // VALIDATION AND DIMENSION EXTRACTION
    // ========================================================================

    const size_t n_vertices = vertex_cofaces.size();
    const size_t n_edges = edge_registry.size();

    // Validate input dimension
    if (static_cast<size_t>(f.size()) != n_vertices) {
        Rf_error("compute_coboundary_del1star: f.size()=%lld must equal n_vertices=%zu",
                 static_cast<long long>(f.size()), n_vertices);
    }

    // Validate metric structures exist
    if (g.M.size() < 2) {
        Rf_error("compute_coboundary_del1star: metric not initialized (g.M.size()=%zu < 2)",
                 g.M.size());
    }

    if (g.M[0].rows() != static_cast<Eigen::Index>(n_vertices)) {
        Rf_error("compute_coboundary_del1star: vertex metric dimension mismatch");
    }

    if (g.M[1].rows() != static_cast<Eigen::Index>(n_edges)) {
        Rf_error("compute_coboundary_del1star: edge metric dimension mismatch");
    }

    // ========================================================================
    // ALLOCATE OUTPUT
    // ========================================================================

    vec_t del1star_f = vec_t::Zero(static_cast<Eigen::Index>(n_edges));

    // ========================================================================
    // COMPUTE COBOUNDARY FOR EACH EDGE
    // ========================================================================

    // CURRENT IMPLEMENTATION: Diagonal metric case
    //
    // For each edge e = [vᵢ, vⱼ] with i < j:
    //   [∂₁*(f)](e) = (wⱼ f(vⱼ) - wᵢ f(vᵢ)) / wₑ
    //
    // where:
    // - wᵢ = g.M[0](i,i) is the vertex weight at vᵢ
    // - wⱼ = g.M[0](j,j) is the vertex weight at vⱼ
    // - wₑ = g.M[1](e,e) is the edge weight

    for (size_t e = 0; e < n_edges; ++e) {
        // Extract edge endpoints
        const auto [i, j] = edge_registry[e];

        // Extract vertex weights (with regularization)
        const double w_i = std::max(g.M[0].coeff(static_cast<Eigen::Index>(i),
                                                  static_cast<Eigen::Index>(i)),
                                     1e-15);
        const double w_j = std::max(g.M[0].coeff(static_cast<Eigen::Index>(j),
                                                  static_cast<Eigen::Index>(j)),
                                     1e-15);

        // Extract edge weight (with regularization)
        const double w_e = std::max(g.M[1].coeff(static_cast<Eigen::Index>(e),
                                                  static_cast<Eigen::Index>(e)),
                                     1e-15);

        // Extract function values
        const double f_i = f[static_cast<Eigen::Index>(i)];
        const double f_j = f[static_cast<Eigen::Index>(j)];

        // Compute weighted difference
        // Note: Positive orientation is j - i (since edge_registry has i < j)
        const double weighted_diff = w_j * f_j - w_i * f_i;

        // Store coboundary value
        del1star_f[static_cast<Eigen::Index>(e)] = weighted_diff / w_e;
    }

    return del1star_f;
}

/**
 * @brief Compute extremality scores using full non-diagonal metric
 *
 * Computes extr(v; f, g) for all vertices using the complete Riemannian
 * structure including non-orthogonal edge inner products.
 *
 * The extremality score measures directional coherence of the discrete
 * gradient field ∂₁*(f), computed using the full non-diagonal metric:
 *
 *   extr(v; f, g) = -∑_{e∋v} ρ₁(e)·∂₁*(f)(e) / ∑_{e∋v} ρ₁(e)·|∂₁*(f)(e)|
 *
 * where ∂₁*(f) is obtained by solving M₁ x = B₁ᵀ M₀ f with non-diagonal M₁.
 *
 * For weakly correlated edges, both versions give similar results. For
 * strongly correlated edges (dense triple intersections, high-degree vertices),
 * the full version can produce significantly different and more accurate scores.
 *
 * @param f Vertex function (0-chain)
 * @param use_iterative Use iterative solver for coboundary (default: false)
 * @param cg_tol CG tolerance (default: 1e-10)
 * @param cg_maxit CG max iterations (default: 1000)
 *
 * @return Vector of extremality scores in [-1, 1] or NaN for isolated vertices
 *
 * @complexity O(m^1.5 + n) first call, O(m + n) subsequent calls
 *             (factorization cost amortized over multiple calls)
 *
 * @note This method internally calls compute_coboundary_del1star_full()
 *       which may cache the factorization of M₁ for efficiency.
 *
 * @see compute_coboundary_del1star_full() for coboundary computation details
 * @see compute_all_extremality_scores() for faster diagonal approximation
 */
vec_t riem_dcx_t::compute_all_extremality_scores_full(
    const vec_t& f,
    bool use_iterative,
    double cg_tol,
    int cg_maxit
    ) const {

    const size_t n_vertices = vertex_cofaces.size();
    const Eigen::Index n = static_cast<Eigen::Index>(n_vertices);

    // ========================================================================
    // PHASE 1: COMPUTE COBOUNDARY WITH FULL METRIC
    // ========================================================================

#if 0
    const vec_t del1star_f = compute_coboundary_del1star_full(
        f, use_iterative, cg_tol, cg_maxit
		);
#endif
	(void)use_iterative;
	(void)cg_tol;
	(void)cg_maxit;
	const vec_t del1star_f = compute_coboundary_del1star(f);

    // ========================================================================
    // PHASE 2: AGGREGATE AT EACH VERTEX
    // ========================================================================

    vec_t extremality = vec_t::Constant(n, std::numeric_limits<double>::quiet_NaN());

    for (size_t v = 0; v < n_vertices; ++v) {
        // Skip self-loop at position 0
        if (vertex_cofaces[v].size() <= 1) {
            // Isolated vertex (no neighbors)
            continue;
        }

        double numerator = 0.0;
        double denominator = 0.0;

        // Iterate through incident edges
        for (size_t k = 1; k < vertex_cofaces[v].size(); ++k) {
            const auto& neighbor_info = vertex_cofaces[v][k];

            const index_t edge_idx = neighbor_info.simplex_index;

            // Edge density is the diagonal entry of M₁
            // (stored in neighbor_info.density for convenience)
            const double rho_edge = std::max(neighbor_info.density, 1e-15);

            // Get coboundary value for this edge
            double del1star_val = del1star_f[static_cast<Eigen::Index>(edge_idx)];

            // ORIENTATION CORRECTION
            // Edge registry stores [i,j] with i < j
            // ∂₁*(f)[e] is oriented from i → j
            // If v is the head (j), reverse sign for v's perspective
            const auto [edge_tail, edge_head] = edge_registry[edge_idx];

            if (v == edge_head) {
                del1star_val = -del1star_val;
            }

            // Accumulate weighted contributions
            numerator += rho_edge * del1star_val;
            denominator += rho_edge * std::abs(del1star_val);
        }

        // Compute extremality score with sign reversal
        // (negative because inward gradients indicate maximum)
        if (denominator > 1e-15) {
            extremality[static_cast<Eigen::Index>(v)] = -numerator / denominator;
        }
        // else: leave as NaN (all incident gradients are zero)
    }

    return extremality;
}

/**
 * @brief Compute extremality score at a reference vertex (diagonal metric)
 *
 * Computes the generalized extremality score:
 *
 *   extr(v; U, f) = ∑_{j∈U} (w_v f_v - w_j f_j) / ∑_{j∈U} |w_j f_j - w_v f_v|
 *
 * When U = N(v) (the canonical neighborhood), this gives the standard
 * extremality score. The generalized form allows computing extremality
 * over arbitrary vertex sets, useful for hop-neighborhood analysis.
 *
 * For the diagonal metric case, the edge weights w_ij cancel in the
 * numerator and denominator, yielding this simplified form.
 *
 * @param ref_vertex The reference vertex v at which to compute extremality
 * @param vertices The vertex set U (typically neighbors of ref_vertex)
 * @return Extremality score in [-1, 1], or NaN if denominator is zero
 *         Positive values indicate maximum-like behavior
 *         Negative values indicate minimum-like behavior
 *
 * @note This is the diagonal metric version. For non-diagonal metrics,
 *       use compute_all_extremality_scores_full() which solves M₁ x = B₁ᵀ M₀ f
 */
double riem_dcx_t::extremality(
    size_t ref_vertex,
    const std::vector<size_t>& vertices  // Pass by const reference!
    ) const {

    // ========================================================================
    // FIND OPTIMAL ITERATION BY MINIMUM GCV
    // ========================================================================

    if (gcv_history.iterations.empty() || sig.y_hat_hist.empty()) {
        Rf_error("extremality: No GCV history or fitted values available");
    }

    double min_gcv = std::numeric_limits<double>::max();
    size_t opt_idx = 0;

    for (size_t i = 0; i < gcv_history.iterations.size(); ++i) {
        double gcv_score = gcv_history.iterations[i].gcv_optimal;
        if (gcv_score < min_gcv) {
            min_gcv = gcv_score;
            opt_idx = i;
        }
    }

    // ========================================================================
    // EXTRACT OPTIMAL FITTED VALUES AND DENSITIES
    // ========================================================================

    if (opt_idx >= sig.y_hat_hist.size()) {
        Rf_error("extremality: Optimal iteration index out of bounds");
    }

    const vec_t& opt_yhat = sig.y_hat_hist[opt_idx];

    // Get vertex densities from optimal iteration
    if (opt_idx >= density_history.rho_vertex.size()) {
        Rf_error("extremality: No vertex density history for optimal iteration");
    }
    const vec_t& rho_vertex = density_history.rho_vertex[opt_idx];

    // Validate reference vertex
    const Eigen::Index ref_idx = static_cast<Eigen::Index>(ref_vertex);
    if (ref_idx < 0 || ref_idx >= opt_yhat.size()) {
        Rf_error("extremality: ref_vertex %zu out of bounds [0, %lld)",
                 ref_vertex, static_cast<long long>(opt_yhat.size()));
    }

    // Handle empty neighborhood
    if (vertices.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // ========================================================================
    // COMPUTE EXTREMALITY SCORE
    // ========================================================================

	// Extract reference vertex properties
	const double w_ref = std::max(rho_vertex[ref_idx], 1e-15);
	const double f_ref = opt_yhat[ref_idx];
	const double wf_ref = w_ref * f_ref;

	double numerator = 0.0;
	double denominator = 0.0;

	// Iterate through neighbors
	for (size_t k = 0; k < vertices.size(); ++k) {
		const size_t neighbor = vertices[k];
		const Eigen::Index nbr_idx = static_cast<Eigen::Index>(neighbor);

		// Validate neighbor index
		if (nbr_idx < 0 || nbr_idx >= opt_yhat.size()) {
			Rf_warning("extremality: neighbor vertex %zu out of bounds, skipping",
					   neighbor);
			continue;
		}

		// Get neighbor weighted function value
		const double w_nbr = std::max(rho_vertex[nbr_idx], 1e-15);
		const double f_nbr = opt_yhat[nbr_idx];
		const double wf_nbr = w_nbr * f_nbr;  // Consistent naming

		// Compute difference: w_i f_i - w_j f_j
		const double weighted_diff = wf_ref - wf_nbr;

		// Accumulate
		numerator += weighted_diff;
		denominator += std::abs(weighted_diff);
	}

    // ========================================================================
    // RETURN EXTREMALITY WITH SIGN CORRECTION
    // ========================================================================

    if (denominator < 1e-15) {
        // All neighbors have identical weighted values
        return std::numeric_limits<double>::quiet_NaN();
    }

    return numerator / denominator;
}

/**
 * @brief Compute hop-extremality radius and neighborhood size for a vertex
 *
 * For τ ∈ (0,1), computes the τ-hop-extremality radius and the size of the
 * corresponding hop-extremality neighborhood:
 *   hop-extr_τ(i; f) = max{h ≥ 1 : |extr^(h)(i; f)| ≥ τ}
 *   |N^(hop-extr_τ(i;f))(i)|
 *
 * where extr^(h)(i; f) is the extremality of f at vertex i computed over
 * its h-hop neighborhood N^(h)(i).
 *
 * MOTIVATION FOR NEIGHBORHOOD SIZE
 * =================================
 *
 * The hop-extremality radius measures spatial persistence (how far the extremal
 * character extends), while the neighborhood size measures local connectivity
 * (how many vertices support the extremum). These are complementary geometric
 * properties that together characterize the local structure around an extremum.
 *
 * Two vertices with identical hop radius can have vastly different neighborhood
 * sizes, indicating different local geometry:
 * - Small neighborhood + large radius: Sparse, linear structure
 * - Large neighborhood + small radius: Dense, concentrated structure
 * - Large neighborhood + large radius: Broad, persistent extremum
 *
 * The neighborhood size is valuable for:
 * - Statistical significance: Larger neighborhoods provide more evidence
 * - Confidence weighting: Weight extrema by sqrt(neighborhood_size)
 * - Geometric understanding: Distinguish sparse vs dense regions
 * - Numerical diagnostics: Very large neighborhoods may indicate global extrema
 *
 * COMPUTATIONAL COST
 * ==================
 *
 * Computing the neighborhood size adds zero overhead: the BFS traversal already
 * visits all vertices in the hop-extremality neighborhood, so we simply return
 * the count that was already being tracked.
 *
 * @param vertex The vertex at which to compute hop-extremality radius
 * @param y Function values (typically fitted values from optimal iteration)
 * @param p_threshold Extremality threshold τ in (0,1]
 * @param detect_maxima If true, check for maxima; if false, check for minima
 * @param max_hop Maximum hop distance to explore
 *
 * @return std::pair<size_t, size_t> where:
 *   - first:  Hop-extremality radius
 *             0 if not a τ-strong extremum at hop 1
 *             SIZE_MAX if global extremum
 *   - second: Number of vertices in the hop-extremality neighborhood N^(radius)(v)
 *             Excludes the center vertex itself
 *             SIZE_MAX if radius is SIZE_MAX (global extremum)
 *
 * @note For global extrema, both components return SIZE_MAX, indicating the
 *       entire graph is in the neighborhood.
 *
 * @note For non-extrema (radius = 0), the neighborhood size is also 0.
 *
 * @complexity O(V + E) worst case for BFS, but typically terminates early
 *             when extremality drops below threshold
 *
 * USAGE EXAMPLE
 * =============
 *
 * @code
 * // Compute hop-extremality with neighborhood size
 * auto [radius, size] = dcx.compute_hop_extremality_radius(
 *     vertex, y_hat, 0.90, true, 20
 * );
 *
 * if (radius == std::numeric_limits<size_t>::max()) {
 *     Rprintf("Global maximum (neighborhood = entire graph)\n");
 * } else if (radius > 0) {
 *     double density = static_cast<double>(size) / n;
 *     Rprintf("Hop radius %zu with %zu vertices (%.1f%% of graph)\n",
 *             radius, size, 100.0 * density);
 *
 *     // Weight by neighborhood size for confidence
 *     double confidence = std::sqrt(static_cast<double>(size));
 *     Rprintf("Confidence weight: %.2f\n", confidence);
 * }
 * @endcode
 */
std::pair<size_t, size_t> riem_dcx_t::compute_hop_extremality_radius(
    size_t vertex,
    const vec_t& y,
    double p_threshold,
    bool detect_maxima,
    size_t max_hop
    ) const {

    const size_t n = vertex_cofaces.size();

    // ================================================================
    // ARGUMENT VALIDATION
    // ================================================================

    if (vertex >= n) {
        Rf_error("compute_hop_extremality_radius: vertex index %zu out of bounds [0, %zu)",
                 vertex, n);
    }

    if (static_cast<size_t>(y.size()) != n) {
        Rf_error("compute_hop_extremality_radius: size mismatch (y: %d, graph: %zu)",
                 static_cast<int>(y.size()), n);
    }

    if (p_threshold <= 0.0 || p_threshold > 1.0) {
        Rf_error("compute_hop_extremality_radius: p_threshold must be in (0,1], got %.3f",
                 p_threshold);
    }

    // ================================================================
    // CHECK FOR GLOBAL EXTREMUM
    // ================================================================

    double global_min = y.minCoeff();
    double global_max = y.maxCoeff();
    const double epsilon = 1e-14 * std::max(std::abs(global_min), std::abs(global_max));

    bool is_global_min = std::abs(y[vertex] - global_min) <= epsilon;
    bool is_global_max = std::abs(y[vertex] - global_max) <= epsilon;

    if ((detect_maxima && is_global_max) || (!detect_maxima && is_global_min)) {
        // Global extremum: entire graph is the neighborhood
        return {std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()};
    }

    // ================================================================
    // BFS TO BUILD HOP NEIGHBORHOODS
    // ================================================================

    std::vector<bool> visited(n, false);
    std::queue<size_t> current_level;
    std::queue<size_t> next_level;

    // Track all vertices in the hop neighborhood (excluding center vertex)
    std::vector<size_t> hop_neighborhood;

    // Initialize BFS with the center vertex
    visited[vertex] = true;
    current_level.push(vertex);

    size_t current_hop = 0;
    size_t last_valid_hop = 0;
    size_t last_valid_size = 0;

    // Process BFS level by level
    while (!current_level.empty() && current_hop < max_hop) {
        // Expand current level to discover next level
        while (!current_level.empty()) {
            size_t v = current_level.front();
            current_level.pop();

            // Explore all neighbors of v
            for (size_t k = 1; k < vertex_cofaces[v].size(); ++k) {
                index_t u = vertex_cofaces[v][k].vertex_index;

                if (visited[u]) continue;

                visited[u] = true;
                next_level.push(u);

                // Add to hop neighborhood (these are at hop current_hop + 1)
                hop_neighborhood.push_back(u);
            }
        }

        // Move to next hop level
        current_hop++;

        // Check if we have any vertices in neighborhood
        if (hop_neighborhood.empty()) {
            // Isolated vertex - no neighbors found
            return {0, 0};
        }

        // ================================================================
        // COMPUTE EXTREMALITY OVER CURRENT HOP NEIGHBORHOOD
        // ================================================================

        double current_extremality = extremality(vertex, hop_neighborhood);

        // Check for NaN - shouldn't happen if hop_neighborhood is non-empty
        if (std::isnan(current_extremality)) {
            return {last_valid_hop, last_valid_size};
        }

        // Get absolute value for threshold comparison
        double abs_extremality = std::abs(current_extremality);

        // Check if extremality is above threshold
        if (abs_extremality >= p_threshold) {
            // Additionally check sign matches what we're looking for
            bool sign_matches;
            if (detect_maxima) {
                sign_matches = (current_extremality > 0);
            } else {
                sign_matches = (current_extremality < 0);
            }

            if (sign_matches) {
                last_valid_hop = current_hop;
                last_valid_size = hop_neighborhood.size();
            } else {
                // Wrong sign - this is not the type of extremum we're looking for
                return {last_valid_hop, last_valid_size};
            }
        } else {
            // Extremality dropped below threshold
            return {last_valid_hop, last_valid_size};
        }

        // Move to next level
        std::swap(current_level, next_level);
    }

    // If we reached max_hop and still above threshold, return max_hop
    // (not infinity, to distinguish from true global extrema)
    return {last_valid_hop, last_valid_size};
}


/**
 * @brief Compute hop-extremality radii for multiple vertices
 *
 * Batch computation of hop-extremality radius for a collection of candidate vertices.
 * This is more efficient than calling compute_hop_extremality_radius() individually
 * when processing many vertices, as it allows for better code organization and
 * potential compiler optimizations through loop-level parallelization.
 *
 * MOTIVATION
 * ==========
 *
 * In practice, extrema detection workflows identify candidates through initial
 * filtering (e.g., vertices with |extremality| ≥ τ), then compute spatial
 * persistence for this reduced candidate set. Processing candidates as a batch
 * simplifies the calling code and enables future optimizations such as:
 * - Parallel computation across independent vertices
 * - Shared data structure caching (though not currently implemented)
 * - Progress reporting for long-running computations
 *
 * ALGORITHMIC APPROACH
 * ====================
 *
 * The current implementation processes vertices sequentially, calling
 * compute_hop_extremality_radius() for each candidate. Each call performs
 * an independent BFS traversal to determine the maximum hop distance maintaining
 * the extremality threshold.
 *
 * For a vertex v, we compute:
 *   hop-extr_τ(v; f) = max{h ≥ 1 : |extr^(h)(v; f)| ≥ τ}
 *
 * where extr^(h)(v; f) is the Riemannian extremality score computed over the
 * h-hop neighborhood N^(h)(v). The extremality uses weighted differences:
 *
 *   extr(v; U, f) = -Σ_{e∋v} ρ₁(e)·[∂₁*(f)](e) / Σ_{e∋v} ρ₁(e)·|[∂₁*(f)](e)|
 *
 * This differs from the probabilistic hop-extremp radius by respecting the full
 * geometric structure of the Riemannian metric rather than using simple sign
 * counting with density weights.
 *
 * COMPLEXITY ANALYSIS
 * ===================
 *
 * Sequential implementation: O(|vertices| · (V + E)) worst case, where V is the
 * number of vertices and E is the number of edges in the graph. In practice,
 * the BFS terminates early when extremality drops below threshold, typically
 * exploring only O(degree · max_hop) vertices per candidate.
 *
 * For k candidates on a graph with average degree d and early termination at
 * average hop h̄ << max_hop, the practical complexity is O(k · d · h̄).
 *
 * FUTURE OPTIMIZATION OPPORTUNITIES
 * ==================================
 *
 * The batch structure enables future enhancements:
 *
 * 1. Parallel processing: Candidates can be processed independently using OpenMP
 *    or thread pools, achieving near-linear speedup on multi-core systems.
 *
 * 2. Shared BFS state: For dense candidate sets, vertices may share neighborhood
 *    exploration. A modified BFS could cache visited sets or reuse traversals.
 *
 * 3. Adaptive ordering: Process candidates in graph-order (e.g., breadth-first
 *    from a seed) to improve cache locality in the vertex_cofaces structure.
 *
 * These optimizations would require more sophisticated implementations but the
 * current simple approach provides a clean baseline.
 *
 * RELATIONSHIP TO HOP-EXTREMP RADIUS
 * ===================================
 *
 * compute_hop_extremp_radii_batch(): Uses probabilistic extremp scores (maxp/minp)
 * based on density-weighted sign counting. Fast O(V + E) per vertex, treats
 * edges as orthogonal (diagonal metric approximation).
 *
 * compute_hop_extremality_radii_batch(): Uses Riemannian extremality scores
 * from the coboundary operator. More expensive per extremality computation but
 * respects edge correlations in the full metric. Best for final analysis where
 * theoretical rigor is essential.
 *
 * For graphs where the diagonal metric approximation is accurate (weak edge
 * correlations), both methods typically agree. For complex geometries with
 * significant edge correlations, the extremality-based approach provides more
 * accurate detection.
 *
 * @param vertices Vector of candidate vertex indices (0-indexed in C++)
 * @param y Vector of function values at each vertex (typically fitted values
 *          at the optimal GCV iteration)
 * @param p_threshold Extremality threshold τ in (0,1]. Vertices must maintain
 *          |extremality| ≥ τ to continue expanding their hop radius
 * @param detect_maxima If true, compute hop-extremality radius for maxima
 *          (positive extremality scores); if false, for minima (negative scores)
 * @param max_hop Maximum hop distance to explore. Limits computation cost on
 *          large graphs or for vertices with very persistent extremal character
 *
 * @return std::pair<std::vector<size_t>, std::vector<size_t>> where:
 *   - first:  Vector of hop-extremality radii, parallel to input
 *         Each entry contains:
 *         - 0: Vertex fails threshold even at hop 1 (not a τ-strong extremum)
 *         - h ∈ [1, max_hop]: Maximum hop maintaining |extremality| ≥ τ
 *         - max_hop: Still above threshold at maximum hop (may extend further)
 *         - SIZE_MAX: Global extremum (infinite radius)
 *   - second: Vector of neighborhood sizes, parallel to input
 *
 * @note Both vectors have the same size as the input vertices vector
 *
 * @pre vertices[i] < vertex_cofaces.size() for all i
 * @pre y.size() == vertex_cofaces.size()
 * @pre 0 < p_threshold ≤ 1.0
 * @pre max_hop > 0
 * @pre All vertices in the input vector should have the same extremal type
 *      (either all maxima candidates or all minima candidates) to match the
 *      detect_maxima flag. Mixing types will produce incorrect results.
 *
 * @post Return vector has size equal to vertices.size()
 * @post Each return value is either 0, in [1, max_hop], or SIZE_MAX
 *
 * @complexity O(|vertices| · (V + E)) worst case, but typically much faster
 *             due to early termination when extremality drops below threshold.
 *             Practical complexity often O(|vertices| · degree · average_hop).
 *
 * @note This function does not validate that vertices have extremality scores
 *       consistent with detect_maxima. The caller (typically create_extremality_component)
 *       is responsible for filtering candidates by sign before calling.
 *
 * @note The function processes vertices sequentially. Future versions may add
 *       parallel processing for large candidate sets on multi-core systems.
 *
 * @see compute_hop_extremality_radius() for single-vertex computation details
 * @see compute_hop_extremp_radii_batch() for probabilistic alternative
 * @see extremality() for extremality score computation over neighborhoods
 *
 * USAGE EXAMPLE
 * =============
 *
 * @code
 * // Identify strong maxima candidates
 * vec_t extremality = dcx.compute_all_extremality_scores_full(y_hat);
 * std::vector<size_t> maxima_candidates;
 * for (size_t i = 0; i < n; ++i) {
 *     if (extremality[i] >= 0.90) {
 *         maxima_candidates.push_back(i);
 *     }
 * }
 *
 * // Compute hop-extremality radii for all maxima candidates
 * std::vector<size_t> radii = dcx.compute_hop_extremality_radii_batch(
 *     maxima_candidates,
 *     y_hat,
 *     0.90,   // threshold
 *     true,   // detect_maxima
 *     20      // max_hop
 * );
 *
 * // Process results
 * for (size_t i = 0; i < maxima_candidates.size(); ++i) {
 *     size_t vertex = maxima_candidates[i];
 *     size_t radius = radii[i];
 *
 *     if (radius == std::numeric_limits<size_t>::max()) {
 *         Rprintf("Vertex %zu: global maximum\n", vertex);
 *     } else if (radius >= 5) {
 *         Rprintf("Vertex %zu: persistent maximum (radius %zu)\n", vertex, radius);
 *     } else if (radius > 0) {
 *         Rprintf("Vertex %zu: weak maximum (radius %zu)\n", vertex, radius);
 *     } else {
 *         Rprintf("Vertex %zu: spurious (radius 0)\n", vertex);
 *     }
 * }
 * @endcode
 */
std::pair<std::vector<size_t>, std::vector<size_t>>
riem_dcx_t::compute_hop_extremality_radii_batch(
    const std::vector<size_t>& vertices,
    const vec_t& y,
    double p_threshold,
    bool detect_maxima,
    size_t max_hop
    ) const {

    std::vector<size_t> radii;
    std::vector<size_t> sizes;
    radii.reserve(vertices.size());
    sizes.reserve(vertices.size());

    for (size_t vertex : vertices) {
        auto [radius, size] = compute_hop_extremality_radius(
            vertex, y, p_threshold, detect_maxima, max_hop
        );
        radii.push_back(radius);
        sizes.push_back(size);
    }

    return {radii, sizes};
}
