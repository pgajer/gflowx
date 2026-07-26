#include "riem_dcx.hpp"
#include <R.h>
#include <Rinternals.h>
#include <cstring>
#include <set>
#include <vector>

// Forward declare conversion utilities from SEXP_cpp_conversion_utils.cpp
// These should be available if that file is compiled into the package
namespace sexp_utils {
    Eigen::SparseMatrix<double> sexp_to_eigen_sparse(SEXP s_X);
    // Add other utility declarations as needed
}

// Forward declare generic SEXP conversion utilities
std::vector<std::vector<int>> convert_adj_list_from_R(SEXP s_adj_list);
std::vector<std::vector<double>> convert_weight_list_from_R(SEXP s_weight_list);

// ================================================================
// HELPER FUNCTIONS TO BUILD NESTED COMPONENTS
// ================================================================

extern "C" SEXP create_density_history_component(const riem_dcx_t& dcx) {
    const size_t n_iters = dcx.density_history.rho_vertex.size();

    if (n_iters == 0) {
        // Return empty list if no history
        SEXP empty = PROTECT(Rf_allocVector(VECSXP, 0));
        UNPROTECT(1);
        return empty;
    }

    const int n_fields = 2;
    SEXP density = PROTECT(Rf_allocVector(VECSXP, n_fields));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_fields));
    int idx = 0;

    // Field 1: rho.vertex (list of vectors, one per iteration)
    SET_STRING_ELT(names, idx, Rf_mkChar("rho.vertex"));
    SEXP s_rho_list = PROTECT(Rf_allocVector(VECSXP, n_iters));

    for (size_t i = 0; i < n_iters; ++i) {
        const vec_t& rho = dcx.density_history.rho_vertex[i];
        const Eigen::Index n_vertices = rho.size();

        SEXP s_rho_vec = PROTECT(Rf_allocVector(REALSXP, n_vertices));
        for (Eigen::Index j = 0; j < n_vertices; ++j) {
            REAL(s_rho_vec)[j] = rho[j];
        }
        SET_VECTOR_ELT(s_rho_list, i, s_rho_vec);
        UNPROTECT(1);
    }
    SET_VECTOR_ELT(density, idx++, s_rho_list);
    UNPROTECT(1);

    // Field 2: n.iterations (scalar, for convenience)
    SET_STRING_ELT(names, idx, Rf_mkChar("n.iterations"));
    SET_VECTOR_ELT(density, idx++, Rf_ScalarInteger(n_iters));

    Rf_setAttrib(density, R_NamesSymbol, names);
    UNPROTECT(2); // names, density
    return density;
}

extern "C" SEXP create_gcv_component(const riem_dcx_t& dcx) {
    const size_t n_iters = dcx.gcv_history.iterations.size();

    const int n_fields = 4;
    SEXP gcv = PROTECT(Rf_allocVector(VECSXP, n_fields));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_fields));
    int idx = 0;

    // eta.optimal (vector)
    SET_STRING_ELT(names, idx, Rf_mkChar("eta.optimal"));
    SEXP s_eta_opt = PROTECT(Rf_allocVector(REALSXP, n_iters));
    for (size_t i = 0; i < n_iters; ++i) {
        REAL(s_eta_opt)[i] = dcx.gcv_history.iterations[i].eta_optimal;
    }
    SET_VECTOR_ELT(gcv, idx++, s_eta_opt);
    UNPROTECT(1);

    // gcv.optimal (vector)
    SET_STRING_ELT(names, idx, Rf_mkChar("gcv.optimal"));
    SEXP s_gcv_opt = PROTECT(Rf_allocVector(REALSXP, n_iters));
    for (size_t i = 0; i < n_iters; ++i) {
        REAL(s_gcv_opt)[i] = dcx.gcv_history.iterations[i].gcv_optimal;
    }
    SET_VECTOR_ELT(gcv, idx++, s_gcv_opt);
    UNPROTECT(1);

    // eta.grid (list of vectors)
    SET_STRING_ELT(names, idx, Rf_mkChar("eta.grid"));
    SEXP s_eta_grid_list = PROTECT(Rf_allocVector(VECSXP, n_iters));
    for (size_t i = 0; i < n_iters; ++i) {
        const auto& grid = dcx.gcv_history.iterations[i].eta_grid;
        SEXP s_grid = PROTECT(Rf_allocVector(REALSXP, grid.size()));
        for (size_t j = 0; j < grid.size(); ++j) {
            REAL(s_grid)[j] = grid[j];
        }
        SET_VECTOR_ELT(s_eta_grid_list, i, s_grid);
        UNPROTECT(1);
    }
    SET_VECTOR_ELT(gcv, idx++, s_eta_grid_list);
    UNPROTECT(1);

    // gcv.scores (list of vectors)
    SET_STRING_ELT(names, idx, Rf_mkChar("gcv.scores"));
    SEXP s_gcv_scores_list = PROTECT(Rf_allocVector(VECSXP, n_iters));
    for (size_t i = 0; i < n_iters; ++i) {
        const auto& scores = dcx.gcv_history.iterations[i].gcv_scores;
        SEXP s_scores = PROTECT(Rf_allocVector(REALSXP, scores.size()));
        for (size_t j = 0; j < scores.size(); ++j) {
            REAL(s_scores)[j] = scores[j];
        }
        SET_VECTOR_ELT(s_gcv_scores_list, i, s_scores);
        UNPROTECT(1);
    }
    SET_VECTOR_ELT(gcv, idx++, s_gcv_scores_list);
    UNPROTECT(1);

    Rf_setAttrib(gcv, R_NamesSymbol, names);
    UNPROTECT(2); // names, gcv
    return gcv;
}

extern "C" SEXP create_graph_component(const riem_dcx_t& dcx) {
    const int n_fields = 8;  // Expanded from 6 to 8
    SEXP graph = PROTECT(Rf_allocVector(VECSXP, n_fields));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_fields));
    int idx = 0;

    // n.vertices
    SET_STRING_ELT(names, idx, Rf_mkChar("n.vertices"));
    SET_VECTOR_ELT(graph, idx++, Rf_ScalarInteger(dcx.vertex_cofaces.size()));

    // n.edges
    SET_STRING_ELT(names, idx, Rf_mkChar("n.edges"));
    SET_VECTOR_ELT(graph, idx++, Rf_ScalarInteger(dcx.edge_registry.size()));

    // edge.lengths
    SET_STRING_ELT(names, idx, Rf_mkChar("edge.lengths"));
    const size_t n_edge_lengths = dcx.edge_lengths.size();
    SEXP s_edge_lengths = PROTECT(Rf_allocVector(REALSXP, n_edge_lengths));
    for (size_t e = 0; e < n_edge_lengths; ++e) {
        REAL(s_edge_lengths)[e] = dcx.edge_lengths[e];
    }
    SET_VECTOR_ELT(graph, idx++, s_edge_lengths);
    UNPROTECT(1);

    // vertex.densities
    SET_STRING_ELT(names, idx, Rf_mkChar("vertex.densities"));
    const Eigen::Index n_verts = dcx.vertex_cofaces.size();
    SEXP s_vert_dens = PROTECT(Rf_allocVector(REALSXP, n_verts));
    for (Eigen::Index i = 0; i < n_verts; ++i) {
        REAL(s_vert_dens)[i] = dcx.vertex_cofaces[i][0].density;
    }
    SET_VECTOR_ELT(graph, idx++, s_vert_dens);
    UNPROTECT(1);

    // edge.densities
    SET_STRING_ELT(names, idx, Rf_mkChar("edge.densities"));
    const Eigen::Index n_edges_dens = dcx.edge_registry.size();
    SEXP s_edge_dens = PROTECT(Rf_allocVector(REALSXP, n_edges_dens));
    for (Eigen::Index e = 0; e < n_edges_dens; ++e) {
        // Find edge density from first vertex of the edge
        const auto [v0, v1] = dcx.edge_registry[e];
        double edge_density = 0.0;
        for (size_t k = 1; k < dcx.vertex_cofaces[v0].size(); ++k) {
            if (dcx.vertex_cofaces[v0][k].vertex_index == v1) {
                edge_density = dcx.vertex_cofaces[v0][k].density;
                break;
            }
        }
        REAL(s_edge_dens)[e] = edge_density;
    }
    SET_VECTOR_ELT(graph, idx++, s_edge_dens);
    UNPROTECT(1);

    // edge.list (n_edges × 2 matrix)
    SET_STRING_ELT(names, idx, Rf_mkChar("edge.list"));
    const size_t n_edges = dcx.edge_registry.size();
    SEXP s_edge_list = PROTECT(Rf_allocMatrix(INTSXP, n_edges, 2));
    int* edge_data = INTEGER(s_edge_list);
    for (size_t e = 0; e < n_edges; ++e) {
        const auto [v0, v1] = dcx.edge_registry[e];
        edge_data[e] = v0 + 1;  // R uses 1-based indexing
        edge_data[e + n_edges] = v1 + 1;
    }
    SET_VECTOR_ELT(graph, idx++, s_edge_list);
    UNPROTECT(1);

    // ---- adjacency list and edge lengths by vertex ----

    // adj.list (list of integer vectors)
    SET_STRING_ELT(names, idx, Rf_mkChar("adj.list"));
    SEXP adj_list = PROTECT(Rf_allocVector(VECSXP, static_cast<size_t>(n_verts)));

    for (Eigen::Index i = 0; i < n_verts; ++i) {
        const auto& nbhrs = dcx.vertex_cofaces[i];
        // Note: nbhrs[0] corresponds to vertex i itself, so we skip it
        const size_t n_neighbors = nbhrs.size() - 1;

        SEXP RA = PROTECT(Rf_allocVector(INTSXP, (R_len_t)n_neighbors));
        int* A  = INTEGER(RA);
        for (size_t j = 1; j < nbhrs.size(); ++j) {
            *A++ = (int)nbhrs[j].vertex_index + 1;  // +1 for R's 1-based indexing
        }
        SET_VECTOR_ELT(adj_list, i, RA);
        UNPROTECT(1); // RA
    }
    SET_VECTOR_ELT(graph, idx++, adj_list);
    UNPROTECT(1); // adj_list

    // edge.length.list (list of numeric vectors)
    SET_STRING_ELT(names, idx, Rf_mkChar("edge.length.list"));
    SEXP edge_length_list = PROTECT(Rf_allocVector(VECSXP, static_cast<size_t>(n_verts)));

    for (Eigen::Index i = 0; i < n_verts; ++i) {
        const auto& nbhrs = dcx.vertex_cofaces[i];
        const size_t n_neighbors = nbhrs.size() - 1;

        SEXP RD = PROTECT(Rf_allocVector(REALSXP, (R_len_t)n_neighbors));
        double* D = REAL(RD);
        for (size_t j = 1; j < nbhrs.size(); ++j) {
            *D++ = nbhrs[j].dist;
        }
        SET_VECTOR_ELT(edge_length_list, i, RD);
        UNPROTECT(1); // RD
    }
    SET_VECTOR_ELT(graph, idx++, edge_length_list);
    UNPROTECT(1); // edge_length_list

    Rf_setAttrib(graph, R_NamesSymbol, names);
    UNPROTECT(2); // names, graph
    return graph;
}

extern "C" SEXP create_iteration_component(const riem_dcx_t& dcx) {
    const int n_fields = 5;
    SEXP iteration = PROTECT(Rf_allocVector(VECSXP, n_fields));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_fields));
    int idx = 0;

    // converged
    SET_STRING_ELT(names, idx, Rf_mkChar("converged"));
    SET_VECTOR_ELT(iteration, idx++, Rf_ScalarLogical(dcx.converged));

    // n.iterations
    SET_STRING_ELT(names, idx, Rf_mkChar("n.iterations"));
    SET_VECTOR_ELT(iteration, idx++, Rf_ScalarInteger(dcx.n_iterations));

    // response.changes
    SET_STRING_ELT(names, idx, Rf_mkChar("response.changes"));
    const size_t n_changes = dcx.response_changes.size();
    SEXP s_resp_changes = PROTECT(Rf_allocVector(REALSXP, n_changes));
    for (size_t i = 0; i < n_changes; ++i) {
        REAL(s_resp_changes)[i] = dcx.response_changes[i];
    }
    SET_VECTOR_ELT(iteration, idx++, s_resp_changes);
    UNPROTECT(1);

    // density.changes
    SET_STRING_ELT(names, idx, Rf_mkChar("density.changes"));
    const size_t n_dens_changes = dcx.density_changes.size();
    SEXP s_dens_changes = PROTECT(Rf_allocVector(REALSXP, n_dens_changes));
    for (size_t i = 0; i < n_dens_changes; ++i) {
        REAL(s_dens_changes)[i] = dcx.density_changes[i];
    }
    SET_VECTOR_ELT(iteration, idx++, s_dens_changes);
    UNPROTECT(1);

    // fitted.history (list of vectors)
    SET_STRING_ELT(names, idx, Rf_mkChar("fitted.history"));
    const size_t n_hist = dcx.sig.y_hat_hist.size();
    SEXP s_history = PROTECT(Rf_allocVector(VECSXP, n_hist));
    for (size_t iter = 0; iter < n_hist; ++iter) {
        const vec_t& y_hat_iter = dcx.sig.y_hat_hist[iter];
        const Eigen::Index n = y_hat_iter.size();
        SEXP s_y_hat = PROTECT(Rf_allocVector(REALSXP, n));
        for (Eigen::Index i = 0; i < n; ++i) {
            REAL(s_y_hat)[i] = y_hat_iter[i];
        }
        SET_VECTOR_ELT(s_history, iter, s_y_hat);
        UNPROTECT(1);
    }
    SET_VECTOR_ELT(iteration, idx++, s_history);
    UNPROTECT(1);

    Rf_setAttrib(iteration, R_NamesSymbol, names);
    UNPROTECT(2); // names, iteration
    return iteration;
}

extern "C" SEXP create_parameters_component(
    index_t k,
    bool use_counting_measure,
    double density_normalization,
    double t_diffusion,
    double beta_damping,
    double gamma_modulation,
    double t_scale_factor,
    double beta_coefficient_factor,
    int n_eigenpairs,
    rdcx_filter_type_t filter_type,
    double epsilon_y,
    double epsilon_rho,
    int max_iterations,
    double density_alpha,
    double density_epsilon,
    int dense_fallback_mode,
    int triangle_policy_mode
    ) {

    const int n_fields = 17;
    SEXP params = PROTECT(Rf_allocVector(VECSXP, n_fields));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_fields));
    int idx = 0;

    // k
    SET_STRING_ELT(names, idx, Rf_mkChar("k"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarInteger(k));

    // use.counting.measure
    SET_STRING_ELT(names, idx, Rf_mkChar("use.counting.measure"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarLogical(use_counting_measure));

    // density.normalization
    SET_STRING_ELT(names, idx, Rf_mkChar("density.normalization"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(density_normalization));

    // t.diffusion
    SET_STRING_ELT(names, idx, Rf_mkChar("t.diffusion"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(t_diffusion));

    // beta.damping (formerly density.uniform.pull)
    SET_STRING_ELT(names, idx, Rf_mkChar("beta.damping"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(beta_damping));

    // response.penalty.exp (formerly gamma.modulation)
    SET_STRING_ELT(names, idx, Rf_mkChar("response.penalty.exp"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(gamma_modulation));

    // NEW: t.scale.factor
    SET_STRING_ELT(names, idx, Rf_mkChar("t.scale.factor"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(t_scale_factor));

    // NEW: beta.coefficient.factor
    SET_STRING_ELT(names, idx, Rf_mkChar("beta.coefficient.factor"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(beta_coefficient_factor));

    // n.eigenpairs
    SET_STRING_ELT(names, idx, Rf_mkChar("n.eigenpairs"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarInteger(n_eigenpairs));

    // filter.type
    SET_STRING_ELT(names, idx, Rf_mkChar("filter.type"));

    const char* filter_str;
    switch (filter_type) {
    case rdcx_filter_type_t::HEAT_KERNEL:
        filter_str = "heat_kernel";
        break;
    case rdcx_filter_type_t::TIKHONOV:
        filter_str = "tikhonov";
        break;
    case rdcx_filter_type_t::CUBIC_SPLINE:
        filter_str = "cubic_spline";
        break;
    case rdcx_filter_type_t::GAUSSIAN:
        filter_str = "gaussian";
        break;
    case rdcx_filter_type_t::EXPONENTIAL:
        filter_str = "exponential";
        break;
    case rdcx_filter_type_t::BUTTERWORTH:
        filter_str = "butterworth";
        break;
    default:
        filter_str = "unknown";
        break;
    }

    SET_VECTOR_ELT(params, idx++, Rf_mkString(filter_str));

    // epsilon.y
    SET_STRING_ELT(names, idx, Rf_mkChar("epsilon.y"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(epsilon_y));

    // epsilon.rho
    SET_STRING_ELT(names, idx, Rf_mkChar("epsilon.rho"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(epsilon_rho));

    // max.iterations
    SET_STRING_ELT(names, idx, Rf_mkChar("max.iterations"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarInteger(max_iterations));

    // density.alpha
    SET_STRING_ELT(names, idx, Rf_mkChar("density.alpha"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(density_alpha));

    // density.epsilon
    SET_STRING_ELT(names, idx, Rf_mkChar("density.epsilon"));
    SET_VECTOR_ELT(params, idx++, Rf_ScalarReal(density_epsilon));

    // dense.fallback
    SET_STRING_ELT(names, idx, Rf_mkChar("dense.fallback"));
    const char* dense_fallback_str = "unknown";
    switch (dense_fallback_mode) {
    case 0:
        dense_fallback_str = "auto";
        break;
    case 1:
        dense_fallback_str = "never";
        break;
    case 2:
        dense_fallback_str = "always";
        break;
    default:
        break;
    }
    SET_VECTOR_ELT(params, idx++, Rf_mkString(dense_fallback_str));

    // triangle.policy
    SET_STRING_ELT(names, idx, Rf_mkChar("triangle.policy"));
    const char* triangle_policy_str = "unknown";
    switch (triangle_policy_mode) {
    case 0:
        triangle_policy_str = "auto";
        break;
    case 1:
        triangle_policy_str = "never";
        break;
    case 2:
        triangle_policy_str = "always";
        break;
    default:
        break;
    }
    SET_VECTOR_ELT(params, idx++, Rf_mkString(triangle_policy_str));

    Rf_setAttrib(params, R_NamesSymbol, names);
    UNPROTECT(2); // names, params
    return params;
}

extern "C" SEXP create_gamma_selection_component(const riem_dcx_t& dcx) {
    // If gamma selection was not performed, return NULL
    if (!dcx.gamma_was_auto_selected) {
        return R_NilValue;
    }

    const auto& gsr = dcx.gamma_selection_result;

    const int n_fields = 5;
    SEXP gamma_sel = PROTECT(Rf_allocVector(VECSXP, n_fields));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_fields));
    int idx = 0;

    // Field 1: gamma.optimal
    SET_STRING_ELT(names, idx, Rf_mkChar("gamma.optimal"));
    SET_VECTOR_ELT(gamma_sel, idx++, Rf_ScalarReal(gsr.gamma_optimal));

    // Field 2: gcv.optimal (FIRST-ITERATION GCV, not global optimal)
    SET_STRING_ELT(names, idx, Rf_mkChar("gcv.first.iter"));
    SET_VECTOR_ELT(gamma_sel, idx++, Rf_ScalarReal(gsr.gcv_optimal));

    // Field 3: gamma.grid
    SET_STRING_ELT(names, idx, Rf_mkChar("gamma.grid"));
    SEXP s_gamma_grid = PROTECT(Rf_allocVector(REALSXP, gsr.gamma_grid.size()));
    for (size_t i = 0; i < gsr.gamma_grid.size(); ++i) {
        REAL(s_gamma_grid)[i] = gsr.gamma_grid[i];
    }
    SET_VECTOR_ELT(gamma_sel, idx++, s_gamma_grid);
    UNPROTECT(1);

    // Field 4: gcv.scores (GCV at first iteration for each gamma)
    SET_STRING_ELT(names, idx, Rf_mkChar("gcv.scores"));
    SEXP s_gcv_scores = PROTECT(Rf_allocVector(REALSXP, gsr.gcv_scores.size()));
    for (size_t i = 0; i < gsr.gcv_scores.size(); ++i) {
        REAL(s_gcv_scores)[i] = gsr.gcv_scores[i];
    }
    SET_VECTOR_ELT(gamma_sel, idx++, s_gcv_scores);
    UNPROTECT(1);

    // Field 5: y.hat.optimal (fitted values at first iteration with optimal gamma)
    SET_STRING_ELT(names, idx, Rf_mkChar("y.hat.first.iter"));
    const Eigen::Index n = gsr.y_hat_optimal.size();
    SEXP s_y_hat_opt = PROTECT(Rf_allocVector(REALSXP, n));
    for (Eigen::Index i = 0; i < n; ++i) {
        REAL(s_y_hat_opt)[i] = gsr.y_hat_optimal[i];
    }
    SET_VECTOR_ELT(gamma_sel, idx++, s_y_hat_opt);
    UNPROTECT(1);

    Rf_setAttrib(gamma_sel, R_NamesSymbol, names);
    UNPROTECT(2); // names, gamma_sel
    return gamma_sel;
}

/**
 * @brief Create spectral component for R return value
 *
 *  Exports 5 fields:
 *   1. eigenvalues         - raw λ values
 *   2. filtered.eigenvalues - F_η(λ) filter weights
 *   3. eigenvectors        - V matrix
 *   4. eta.optimal         - optimal smoothing parameter
 *   5. filter.type         - filter type string
 *
 * The raw eigenvalues enable per-column GCV in refit.rdgraph.regression()
 */
extern "C" SEXP create_spectral_component(const riem_dcx_t& dcx,
                                          int optimal_iter,
                                          rdcx_filter_type_t filter_type) {
    const int n_fields = 5;  // CHANGED from 2 to 5
    SEXP spectral = PROTECT(Rf_allocVector(VECSXP, n_fields));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_fields));
    int idx = 0;

    if (!dcx.spectral_cache.is_valid) {
        Rf_warning("Spectral cache not valid");
        UNPROTECT(2);
        return R_NilValue;
    }

    const size_t n_eigen = dcx.spectral_cache.eigenvalues.size();
    const Eigen::Index n_verts = dcx.spectral_cache.eigenvectors.rows();

    // Get the optimal eta from GCV history
    double eta_opt = dcx.gcv_history.iterations[optimal_iter].eta_optimal;

    // ================================================================
    // 1. eigenvalues: raw λ values (NEW - enables per-column GCV)
    // ================================================================
    SET_STRING_ELT(names, idx, Rf_mkChar("eigenvalues"));
    SEXP s_lambda = PROTECT(Rf_allocVector(REALSXP, n_eigen));
    for (size_t i = 0; i < n_eigen; ++i) {
        REAL(s_lambda)[i] = dcx.spectral_cache.eigenvalues[i];
    }
    SET_VECTOR_ELT(spectral, idx++, s_lambda);
    UNPROTECT(1);

    // ================================================================
    // 2. filtered.eigenvalues: F_η(λ)
    // ================================================================
    SET_STRING_ELT(names, idx, Rf_mkChar("filtered.eigenvalues"));
    SEXP s_f_lambda = PROTECT(Rf_allocVector(REALSXP, n_eigen));

    // Apply filter based on type
    for (size_t i = 0; i < n_eigen; ++i) {
        double lambda = dcx.spectral_cache.eigenvalues[i];
        double filtered_value;

        switch (filter_type) {
            case rdcx_filter_type_t::HEAT_KERNEL:
                filtered_value = std::exp(-eta_opt * lambda);
                break;
            case rdcx_filter_type_t::TIKHONOV:
                filtered_value = 1.0 / (1.0 + eta_opt * lambda);
                break;
            case rdcx_filter_type_t::CUBIC_SPLINE:
                filtered_value = 1.0 / (1.0 + eta_opt * lambda * lambda);
                break;
            case rdcx_filter_type_t::GAUSSIAN:
                filtered_value = std::exp(-eta_opt * lambda * lambda);
                break;
            case rdcx_filter_type_t::EXPONENTIAL:
                filtered_value = std::exp(-eta_opt * std::sqrt(std::max(lambda, 0.0)));
                break;
            case rdcx_filter_type_t::BUTTERWORTH: {
                // Default to n=2 for Butterworth
                double ratio = lambda / std::max(eta_opt, 1e-15);
                filtered_value = 1.0 / (1.0 + ratio * ratio * ratio * ratio);
                break;
            }
            default:
                filtered_value = 1.0;
                break;
        }

        REAL(s_f_lambda)[i] = filtered_value;
    }
    SET_VECTOR_ELT(spectral, idx++, s_f_lambda);
    UNPROTECT(1);

    // ================================================================
    // 3. eigenvectors: V
    // ================================================================
    SET_STRING_ELT(names, idx, Rf_mkChar("eigenvectors"));
    SEXP s_V = PROTECT(Rf_allocMatrix(REALSXP, n_verts, n_eigen));
    double* V_data = REAL(s_V);
    for (size_t j = 0; j < n_eigen; ++j) {
        for (Eigen::Index i = 0; i < n_verts; ++i) {
            V_data[i + j * n_verts] = dcx.spectral_cache.eigenvectors(i, j);
        }
    }
    SET_VECTOR_ELT(spectral, idx++, s_V);
    UNPROTECT(1);

    // ================================================================
    // 4. eta.optimal: optimal smoothing parameter (NEW)
    // ================================================================
    SET_STRING_ELT(names, idx, Rf_mkChar("eta.optimal"));
    SEXP s_eta = PROTECT(Rf_ScalarReal(eta_opt));
    SET_VECTOR_ELT(spectral, idx++, s_eta);
    UNPROTECT(1);

    // ================================================================
    // 5. filter.type: string name of filter (NEW)
    // ================================================================
    SET_STRING_ELT(names, idx, Rf_mkChar("filter.type"));
    const char* filter_name;
    switch (filter_type) {
        case rdcx_filter_type_t::HEAT_KERNEL:
            filter_name = "heat_kernel";
            break;
        case rdcx_filter_type_t::TIKHONOV:
            filter_name = "tikhonov";
            break;
        case rdcx_filter_type_t::CUBIC_SPLINE:
            filter_name = "cubic_spline";
            break;
        case rdcx_filter_type_t::GAUSSIAN:
            filter_name = "gaussian";
            break;
        case rdcx_filter_type_t::EXPONENTIAL:
            filter_name = "exponential";
            break;
        case rdcx_filter_type_t::BUTTERWORTH:
            filter_name = "butterworth";
            break;
        default:
            filter_name = "unknown";
            break;
    }
    SEXP s_filter_type = PROTECT(Rf_mkString(filter_name));
    SET_VECTOR_ELT(spectral, idx++, s_filter_type);
    UNPROTECT(1);

    Rf_setAttrib(spectral, R_NamesSymbol, names);
    UNPROTECT(2);
    return spectral;
}

/**
 * @brief Create extremality component with hop radii and neighborhood sizes
 *
 * Computes extremality scores and optionally their spatial persistence (hop radii)
 * along with the size of each hop-extremality neighborhood. When hop radius
 * computation is enabled, only vertices exceeding the extremality threshold have
 * their radii and neighborhood sizes computed, while others are assigned NA.
 *
 * MOTIVATION FOR NEIGHBORHOOD SIZES
 * ==================================
 *
 * The neighborhood size provides crucial geometric context for interpreting
 * hop-extremality radii. Two vertices with identical hop radius can have
 * vastly different neighborhood sizes:
 * - Small neighborhood: Sparse region, linear structure, less evidence
 * - Large neighborhood: Dense region, concentrated structure, more evidence
 *
 * This information enables:
 * - Statistical confidence weighting (weight by sqrt(neighborhood_size))
 * - Geometric understanding (distinguish sparse vs dense extrema)
 * - Significance assessment (larger neighborhoods = more robust extrema)
 * - Visualization (color-code by neighborhood size)
 *
 * @param dcx The fitted Riemannian density complex
 * @param y_hat Fitted values at optimal GCV iteration
 * @param p_threshold Extremality threshold in (0,1] for candidate selection
 * @param max_hop Maximum hop distance to explore
 *
 * @return R list with components:
 *   - scores: numeric vector of extremality scores in [-1,1] or NaN
 *   - hop.extremality.radii: numeric vector (hop radius for candidates, NA otherwise)
 *   - hop.neighborhood.sizes: numeric vector (neighborhood size for candidates, NA otherwise)
 *
 * @note The hop_neighborhood_sizes vector is parallel to hop_extremality_radii.
 *       Both are NA for vertices with |extremality| < p_threshold.
 *       Both are -1 (encoded from SIZE_MAX) for global extrema.
 *
 * @complexity O(n·m) for extremality scores + O(k·(V+E)) for hop radii,
 *             where k is number of candidates (typically k << n)
 */
extern "C" SEXP create_extremality_component(
    const riem_dcx_t& dcx,
    const vec_t& y_hat,
    double p_threshold,
    size_t max_hop
    ) {

    const Eigen::Index n = y_hat.size();

    // ================================================================
    // COMPUTE EXTREMALITY SCORES
    // ================================================================

    // Use iterative solver with relaxed tolerance to avoid factorization
    // failures when M[1] has numerical issues
    const bool use_iterative = true;
    const double cg_tol = 1e-6;
    const int cg_maxit = 5000;

    vec_t extremality = dcx.compute_all_extremality_scores_full(
        y_hat,
        use_iterative,
        cg_tol,
        cg_maxit
    );

    // Convert extremality scores to R vector
    SEXP s_extremality = PROTECT(Rf_allocVector(REALSXP, n));
    for (Eigen::Index i = 0; i < n; ++i) {
        REAL(s_extremality)[i] = extremality[i];
    }

    // ================================================================
    // DETERMINE RETURN STRUCTURE
    // ================================================================

    int n_fields = 3;  // scores, hop_extremality_radii, hop_neighborhood_sizes
    SEXP result = PROTECT(Rf_allocVector(VECSXP, n_fields));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_fields));

    // Field 1: scores (always present)
    SET_STRING_ELT(names, 0, Rf_mkChar("scores"));
    SET_VECTOR_ELT(result, 0, s_extremality);

    // ================================================================
    // COMPUTE HOP-EXTREMALITY RADII AND NEIGHBORHOOD SIZES
    // ================================================================

    // Initialize hop radii vector with NA for all vertices
    SEXP s_hop_radii = PROTECT(Rf_allocVector(REALSXP, n));
    double* hop_radii_ptr = REAL(s_hop_radii);
    for (Eigen::Index i = 0; i < n; ++i) {
        hop_radii_ptr[i] = NA_REAL;
    }

    // Initialize neighborhood sizes vector with NA for all vertices
    SEXP s_hop_sizes = PROTECT(Rf_allocVector(REALSXP, n));
    double* hop_sizes_ptr = REAL(s_hop_sizes);
    for (Eigen::Index i = 0; i < n; ++i) {
        hop_sizes_ptr[i] = NA_REAL;
    }

    // ================================================================
    // IDENTIFY CANDIDATES BY SIGN
    // ================================================================

    std::vector<size_t> maxima_candidates;
    std::vector<size_t> minima_candidates;

    for (Eigen::Index i = 0; i < n; ++i) {
        double extr = extremality[i];

        // Skip NaN extremality scores
        if (std::isnan(extr)) {
            continue;
        }

        // Check if above threshold
        if (std::abs(extr) >= p_threshold) {
            if (extr > 0) {
                maxima_candidates.push_back(static_cast<size_t>(i));
            } else {
                minima_candidates.push_back(static_cast<size_t>(i));
            }
        }
    }

    // ================================================================
    // COMPUTE HOP RADII AND SIZES FOR MAXIMA CANDIDATES
    // ================================================================

    if (!maxima_candidates.empty()) {
        auto [max_radii, max_sizes] = dcx.compute_hop_extremality_radii_batch(
            maxima_candidates,
            y_hat,
            p_threshold,
            true,  // detect_maxima
            max_hop
        );

        // Store results
        for (size_t i = 0; i < maxima_candidates.size(); ++i) {
            size_t vertex = maxima_candidates[i];
            size_t radius = max_radii[i];
            size_t size = max_sizes[i];

            if (radius == std::numeric_limits<size_t>::max()) {
                // Global maximum: encode as -1
                hop_radii_ptr[vertex] = -1.0;
                hop_sizes_ptr[vertex] = -1.0;  // SIZE_MAX encoded as -1
            } else {
                hop_radii_ptr[vertex] = static_cast<double>(radius);
                hop_sizes_ptr[vertex] = static_cast<double>(size);
            }
        }
    }

    // ================================================================
    // COMPUTE HOP RADII AND SIZES FOR MINIMA CANDIDATES
    // ================================================================

    if (!minima_candidates.empty()) {
        auto [min_radii, min_sizes] = dcx.compute_hop_extremality_radii_batch(
            minima_candidates,
            y_hat,
            p_threshold,
            false,  // detect_maxima
            max_hop
        );

        // Store results
        for (size_t i = 0; i < minima_candidates.size(); ++i) {
            size_t vertex = minima_candidates[i];
            size_t radius = min_radii[i];
            size_t size = min_sizes[i];

            if (radius == std::numeric_limits<size_t>::max()) {
                // Global minimum: encode as -1
                hop_radii_ptr[vertex] = -1.0;
                hop_sizes_ptr[vertex] = -1.0;  // SIZE_MAX encoded as -1
            } else {
                hop_radii_ptr[vertex] = static_cast<double>(radius);
                hop_sizes_ptr[vertex] = static_cast<double>(size);
            }
        }
    }

    // Field 2: hop_extremality_radii
    SET_STRING_ELT(names, 1, Rf_mkChar("hop.extremality.radii"));
    SET_VECTOR_ELT(result, 1, s_hop_radii);

    // Field 3: hop_neighborhood_sizes
    SET_STRING_ELT(names, 2, Rf_mkChar("hop.neighborhood.sizes"));
    SET_VECTOR_ELT(result, 2, s_hop_sizes);

    // Set names attribute
    Rf_setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(5); // s_hop_sizes, s_hop_radii, names, result, s_extremality
    return result;
}

/**
 * @brief Package posterior summary into R list structure
 *
 * @details
 * Converts C++ posterior_summary_t structure into an R named list with components
 * accessible from R code. The returned list contains credible interval bounds,
 * posterior standard deviations, and metadata about the inference procedure.
 *
 * All numeric vectors are converted to R REALSXP with proper memory protection.
 * The returned SEXP must be PROTECTed by the caller if used in further computations.
 *
 * @param summary C++ posterior_summary_t structure from compute_posterior_summary()
 *
 * @return R named list (SEXP) with components:
 *   - lower: Numeric vector (length n) of lower credible bounds
 *   - upper: Numeric vector (length n) of upper credible bounds
 *   - sd: Numeric vector (length n) of posterior standard deviations
 *   - credible.level: Scalar numeric coverage probability
 *   - sigma: Scalar numeric estimated residual SD
 *
 * @note The returned SEXP is already PROTECTed internally and should be
 *   UNPROTECTed by the caller when no longer needed.
 */
extern "C" SEXP create_posterior_component(const posterior_summary_t& summary) {
    const int n = summary.lower.size();

        // Determine list size based on whether samples are included
    int n_components = summary.has_samples ? 6 : 5;
    SEXP s_posterior = PROTECT(Rf_allocVector(VECSXP, n_components));
    SEXP s_names = PROTECT(Rf_allocVector(STRSXP, n_components));

    // ================================================================
    // Component 1: Lower credible bounds
    // ================================================================

    SEXP s_lower = PROTECT(Rf_allocVector(REALSXP, n));
    double* p_lower = REAL(s_lower);
    for (int i = 0; i < n; ++i) {
        p_lower[i] = summary.lower[i];
    }
    SET_VECTOR_ELT(s_posterior, 0, s_lower);
    SET_STRING_ELT(s_names, 0, Rf_mkChar("lower"));

    // ================================================================
    // Component 2: Upper credible bounds
    // ================================================================

    SEXP s_upper = PROTECT(Rf_allocVector(REALSXP, n));
    double* p_upper = REAL(s_upper);
    for (int i = 0; i < n; ++i) {
        p_upper[i] = summary.upper[i];
    }
    SET_VECTOR_ELT(s_posterior, 1, s_upper);
    SET_STRING_ELT(s_names, 1, Rf_mkChar("upper"));

    // ================================================================
    // Component 3: Posterior standard deviations
    // ================================================================

    SEXP s_sd = PROTECT(Rf_allocVector(REALSXP, n));
    double* p_sd = REAL(s_sd);
    for (int i = 0; i < n; ++i) {
        p_sd[i] = summary.posterior_sd[i];
    }
    SET_VECTOR_ELT(s_posterior, 2, s_sd);
    SET_STRING_ELT(s_names, 2, Rf_mkChar("sd"));

    // ================================================================
    // Component 4: Credible level
    // ================================================================

    SEXP s_level = PROTECT(Rf_allocVector(REALSXP, 1));
    REAL(s_level)[0] = summary.credible_level;
    SET_VECTOR_ELT(s_posterior, 3, s_level);
    SET_STRING_ELT(s_names, 3, Rf_mkChar("credible.level"));

    // ================================================================
    // Component 5: Estimated sigma
    // ================================================================

    SEXP s_sigma = PROTECT(Rf_allocVector(REALSXP, 1));
    REAL(s_sigma)[0] = summary.sigma_hat;
    SET_VECTOR_ELT(s_posterior, 4, s_sigma);
    SET_STRING_ELT(s_names, 4, Rf_mkChar("sigma"));

    // ================================================================
    // Component 6: Posterior samples (optional)
    // ================================================================

    if (summary.has_samples) {
        const int n_samples = summary.samples.cols();

        SEXP s_samples = PROTECT(Rf_allocMatrix(REALSXP, n, n_samples));
        double* p_samples = REAL(s_samples);

        // Copy column-major: R stores matrices in column-major order like Eigen
        for (int j = 0; j < n_samples; ++j) {
            for (int i = 0; i < n; ++i) {
                p_samples[i + j * n] = summary.samples(i, j);
            }
        }

        SET_VECTOR_ELT(s_posterior, 5, s_samples);
        SET_STRING_ELT(s_names, 5, Rf_mkChar("samples"));

        UNPROTECT(1);  // s_samples
    }

    // ================================================================
    // Attach names and return
    // ================================================================

    Rf_setAttrib(s_posterior, R_NamesSymbol, s_names);

    // UNPROTECT all temporary allocations (7 total)
    // s_posterior, s_names, s_lower, s_upper, s_sd, s_level, s_sigma
    UNPROTECT(7);

    // Return the list (caller must PROTECT if storing)
    return s_posterior;
}

// ================================================================
// REFERENCE MEASURE COMPONENT (dk diagnostics + weights)
// ================================================================
extern "C" SEXP create_reference_measure_component(const riem_dcx_t& dcx) {
    const Eigen::Index n = (Eigen::Index)dcx.vertex_cofaces.size();

    const int n_fields = 7;
    SEXP rm = PROTECT(Rf_allocVector(VECSXP, n_fields));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_fields));
    int idx = 0;

    // weights (vertex mass / reference measure)
    SET_STRING_ELT(names, idx, Rf_mkChar("weights"));
    SEXP s_w = PROTECT(Rf_allocVector(REALSXP, n));
    for (Eigen::Index i = 0; i < n; ++i) {
        REAL(s_w)[i] = (i < (Eigen::Index)dcx.reference_measure.size())
            ? dcx.reference_measure[(size_t)i]
            : NA_REAL;
    }
    SET_VECTOR_ELT(rm, idx++, s_w);
    UNPROTECT(1);

    // dk.raw (or NULL)
    SET_STRING_ELT(names, idx, Rf_mkChar("dk.raw"));
    if (dcx.dk_raw.empty()) {
        SET_VECTOR_ELT(rm, idx++, R_NilValue);
    } else {
        SEXP s_dk_raw = PROTECT(Rf_allocVector(REALSXP, n));
        for (Eigen::Index i = 0; i < n; ++i) {
            REAL(s_dk_raw)[i] = dcx.dk_raw[(size_t)i];
        }
        SET_VECTOR_ELT(rm, idx++, s_dk_raw);
        UNPROTECT(1);
    }

    // dk.used (or NULL)
    SET_STRING_ELT(names, idx, Rf_mkChar("dk.used"));
    if (dcx.dk_used.empty()) {
        SET_VECTOR_ELT(rm, idx++, R_NilValue);
    } else {
        SEXP s_dk_used = PROTECT(Rf_allocVector(REALSXP, n));
        for (Eigen::Index i = 0; i < n; ++i) {
            REAL(s_dk_used)[i] = dcx.dk_used[(size_t)i];
        }
        SET_VECTOR_ELT(rm, idx++, s_dk_used);
        UNPROTECT(1);
    }

    // dk.lower
    SET_STRING_ELT(names, idx, Rf_mkChar("dk.lower"));
    SET_VECTOR_ELT(rm, idx++, Rf_ScalarReal(dcx.dk_lower));

    // dk.upper
    SET_STRING_ELT(names, idx, Rf_mkChar("dk.upper"));
    SET_VECTOR_ELT(rm, idx++, Rf_ScalarReal(dcx.dk_upper));

    // dk.outside.low (1-based)
    SET_STRING_ELT(names, idx, Rf_mkChar("dk.outside.low"));
    SEXP s_low = PROTECT(Rf_allocVector(INTSXP, (R_len_t)dcx.dk_used_low.size()));
    for (R_len_t i = 0; i < (R_len_t)dcx.dk_used_low.size(); ++i) {
        INTEGER(s_low)[i] = (int)dcx.dk_used_low[(size_t)i] + 1;
    }
    SET_VECTOR_ELT(rm, idx++, s_low);
    UNPROTECT(1);

    // dk.outside.high (1-based)
    SET_STRING_ELT(names, idx, Rf_mkChar("dk.outside.high"));
    SEXP s_high = PROTECT(Rf_allocVector(INTSXP, (R_len_t)dcx.dk_used_high.size()));
    for (R_len_t i = 0; i < (R_len_t)dcx.dk_used_high.size(); ++i) {
        INTEGER(s_high)[i] = (int)dcx.dk_used_high[(size_t)i] + 1;
    }
    SET_VECTOR_ELT(rm, idx++, s_high);
    UNPROTECT(1);

    Rf_setAttrib(rm, R_NamesSymbol, names);
    UNPROTECT(2); // names, rm
    return rm;
}

/**
 * @brief R interface for kNN Riemannian density graph regression
 *
 * Constructs 1-skeleton kNN complex and iteratively refines geometry to
 * reflect response structure for conditional expectation estimation.
 *
 * PARAMETER ADDITIONS FOR EXTREMALITY ANALYSIS
 * =============================================
 *
 * Three new parameters enable optional extremality analysis during fitting:
 *
 * compute_extremality: When FALSE (default), extremality scores are not computed.
 * This is useful for applications like optimal k selection where only fitted
 * values and GCV scores are needed. When TRUE, computes Riemannian extremality
 * scores for all vertices at the optimal GCV iteration.
 *
 * p_threshold: Extremality threshold in (0,1] used for candidate selection when
 * computing hop-extremality radii. Only vertices with |extremality| ≥ p_threshold
 * have their hop radii computed. Default 0.90 matches compute.gextrema.nbhds().
 *
 * max_hop: Maximum hop distance to explore when computing hop-extremality radii.
 * Default 20 balances computational cost against capturing long-range persistence.
 * Vertices maintaining extremality above threshold beyond max_hop are assigned
 * radius = max_hop (not infinity, to distinguish from global extrema).
 *
 * EXTREMALITY OUTPUT STRUCTURE
 * =============================
 *
 * When compute_extremality=FALSE:
 *   extremality.scores component is not included in the result
 *
 * When compute_extremality=TRUE without hop radii:
 *   extremality.scores = list(scores = numeric vector)
 *
 * When compute_extremality=TRUE with hop radii:
 *   extremality.scores = list(
 *     scores = numeric vector of extremality in [-1,1] or NaN,
 *     hop_extremality_radii = numeric vector with:
 *       - hop radius for vertices with |extremality| ≥ p_threshold
 *       - NA for vertices below threshold
 *       - -1 for global extrema
 *       - 0 for vertices failing threshold even at hop 1
 *   )
 *
 * This structure allows users to identify strong local extrema (high |extremality|)
 * and assess their spatial persistence (large hop radius) in a single fitting call.
 *
 * @param s_X Feature matrix (dense REALSXP or sparse dgCMatrix)
 * @param s_y Response vector (REALSXP)
 * @param s_k Number of nearest neighbors (INTSXP)
 * @param s_use_counting_measure Logical: use counting measure? (LGLSXP)
 * @param s_density_normalization Target density sum (REALSXP)
 * @param s_t_diffusion Heat diffusion time, 0=auto (REALSXP)
 * @param s_beta_damping Damping parameter, 0=auto (REALSXP)
 *
 * @param s_t_update_mode
 * @param s_t_update_max_mult
 *
 * @param s_gamma_modulation Response coherence exponent (REALSXP)
 * @param s_n_eigenpairs Number of eigenpairs for filtering (INTSXP)
 * @param s_filter_type Filter type string (STRSXP)
 * @param s_epsilon_y Response convergence threshold (REALSXP)
 * @param s_epsilon_rho Density convergence threshold (REALSXP)
 * @param s_max_iterations Maximum iteration count (INTSXP)
 *
 * @param s_max_path_edge_ratio_thld SEXP object (double) Maximum acceptable ratio of
 *        alternative path length to edge length for geometric pruning.
 *        Edges with ratio <= this value will be pruned. If <= 0, this pruning stage is skipped.
 *
 * @param s_path_edge_ratio_percentile SEXP object (double) Percentile threshold (0.0-1.0)
 *        for edge lengths to consider in geometric pruning. Only edges with length
 *        greater than this percentile are considered for pruning.
 *
 * @param s_threshold_percentile SEXP object (double) Percentile threshold for quantile-based
 *        edge length pruning. Valid range is [0.0, 0.5]. Value of 0.0 disables quantile pruning.
 *        When > 0, edges in the top (1 - threshold_percentile) quantile by length are
 *        removed if their removal preserves connectivity. For example, 0.9 removes top 10% of edges.
 *
 * @param s_density_alpha Density alpha parameter in [1,2] (REALSXP)
 * @param s_density_epsilon Density regularization epsilon (REALSXP)
 * @param s_compute_extremality Logical: compute extremality scores? (LGLSXP)
 * @param s_p_threshold Extremality threshold for hop radii, 0=skip (REALSXP)
 * @param s_max_hop Maximum hop distance for radii computation (INTSXP)
 * @param s_dense_fallback_mode Integer fallback mode for diffusion linear solve:
 *        0=auto, 1=never, 2=always (INTSXP)
 * @param s_triangle_policy_mode Integer triangle construction policy:
 *        0=auto, 1=never, 2=always (INTSXP)
 * @param s_verbose_level SEXP object (integer) controlling progress reporting during computation; possible values: 0, 1, 2, 3.
 *
 * @return R list of class c("knn.riem.fit", "list") with components:
 *   - fitted.values: Fitted values at optimal GCV iteration
 *   - residuals: Residuals y - fitted.values
 *   - optimal.fit: List with detailed fitting information
 *   - graph: List with graph structure
 *   - iteration: Convergence diagnostics
 *   - parameters: Input parameters
 *   - y: Original response vector
 *   - gcv: GCV history across iterations
 *   - density: Density evolution history
 *   - gamma.selection: Gamma parameter selection info
 *   - spectral: Spectral decomposition at optimal iteration
 *   - extremality: Extremality analysis (only if compute_extremality=TRUE)
 *     R list with components:
 *     - scores: numeric vector of extremality scores in [-1,1] or NaN
 *     - hop.extremality.radii: numeric vector (hop radius for candidates, NA otherwise)
 *     - hop.neighborhood.sizes: numeric vector (neighborhood size for candidates, NA otherwise)
 *
 * @note Input validation is performed on the R side in fit.rdgraph.regression().
 *       Additional defensive checks are included here for robustness.
 *
 * @note This function is called from R via .Call(). Input validation is
 *       performed on the R side in fit.knn.riem.graph.regression().
 *       Additional defensive checks are included here for robustness.
 * @note Setting compute_extremality=FALSE reduces computational cost and memory
 *       usage, making this suitable for large-scale k-selection procedures where
 *       extremality analysis is not needed.
 */
extern "C" SEXP S_fit_rdgraph_regression(
    SEXP s_X,
    SEXP s_y,
    SEXP s_y_vertices,
    SEXP s_k,
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_with_posterior,
    SEXP s_return_posterior_samples,
    SEXP s_credible_level,
    SEXP s_n_posterior_samples,
    SEXP s_posterior_seed,
    SEXP s_use_counting_measure,
    SEXP s_density_normalization,
    SEXP s_t_diffusion,
    SEXP s_beta_damping,
    SEXP s_gamma_modulation,
    SEXP s_t_scale_factor,
    SEXP s_beta_coefficient_factor,
    SEXP s_t_update_mode,
    SEXP s_t_update_max_mult,
    SEXP s_n_eigenpairs,
    SEXP s_filter_type,
    SEXP s_epsilon_y,
    SEXP s_epsilon_rho,
    SEXP s_max_iterations,
    SEXP s_max_ratio_threshold,
    SEXP s_path_edge_ratio_percentile,
    SEXP s_threshold_percentile,
    SEXP s_density_alpha,
    SEXP s_density_epsilon,
    SEXP s_clamp_dk,
    SEXP s_dk_clamp_median_factor,
    SEXP s_target_weight_ratio,
    SEXP s_pathological_ratio_threshold,
    SEXP s_compute_extremality,
    SEXP s_p_threshold,
    SEXP s_max_hop,
    SEXP s_knn_cache_path,
    SEXP s_knn_cache_mode,
    SEXP s_dense_fallback_mode,
    SEXP s_triangle_policy_mode,
    SEXP s_verbose_level
) {
    // ================================================================
    // PART I: INPUT EXTRACTION (same as before)
    // ================================================================

    // -------------------- Feature Matrix X --------------------

    Eigen::SparseMatrix<double> X_sparse;
    Eigen::Index n_points = 0;
    Eigen::Index n_features = 0;

    // Check if dense matrix
    bool is_dense = (Rf_isMatrix(s_X) && TYPEOF(s_X) == REALSXP);
    bool is_sparse = Rf_inherits(s_X, "dgCMatrix");

    if (is_dense) {
        // Dense matrix: convert to sparse
        SEXP s_dim = PROTECT(Rf_getAttrib(s_X, R_DimSymbol));

        if (s_dim == R_NilValue || TYPEOF(s_dim) != INTSXP || Rf_length(s_dim) != 2) {
            UNPROTECT(1);
            Rf_error("X must be a valid matrix with dim attribute");
        }

        n_points = INTEGER(s_dim)[0];
        n_features = INTEGER(s_dim)[1];
        UNPROTECT(1);

        if (n_points < 1 || n_features < 1) {
            Rf_error("X has invalid dimensions: %ld × %ld",
                     (long)n_points, (long)n_features);
        }

        const double* X_data = REAL(s_X);

        // Convert to sparse format
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(n_points * n_features / 2); // Heuristic reserve

        for (Eigen::Index j = 0; j < n_features; ++j) {
            for (Eigen::Index i = 0; i < n_points; ++i) {
                double val = X_data[i + j * n_points];
                if (val != 0.0) {
                    triplets.emplace_back(i, j, val);
                }
            }
        }

        X_sparse.resize(n_points, n_features);
        X_sparse.setFromTriplets(triplets.begin(), triplets.end());

    } else if (is_sparse) {
        // Sparse matrix (dgCMatrix): extract components
        SEXP s_i = PROTECT(Rf_getAttrib(s_X, Rf_install("i")));
        SEXP s_p = PROTECT(Rf_getAttrib(s_X, Rf_install("p")));
        SEXP s_x = PROTECT(Rf_getAttrib(s_X, Rf_install("x")));
        SEXP s_dim = PROTECT(Rf_getAttrib(s_X, Rf_install("Dim")));

        if (s_i == R_NilValue || s_p == R_NilValue ||
            s_x == R_NilValue || s_dim == R_NilValue) {
            UNPROTECT(4);
            Rf_error("Invalid dgCMatrix: missing required slots (i, p, x, or Dim)");
        }

        if (TYPEOF(s_i) != INTSXP || TYPEOF(s_p) != INTSXP ||
            TYPEOF(s_x) != REALSXP || TYPEOF(s_dim) != INTSXP) {
            UNPROTECT(4);
            Rf_error("Invalid dgCMatrix: slots have incorrect types");
        }

        const int* i_data = INTEGER(s_i);
        const int* p_data = INTEGER(s_p);
        const double* x_data = REAL(s_x);
        const int* dim_data = INTEGER(s_dim);

        n_points = dim_data[0];
        n_features = dim_data[1];

        if (n_points < 1 || n_features < 1) {
            UNPROTECT(4);
            Rf_error("dgCMatrix has invalid dimensions: %ld × %ld",
                     (long)n_points, (long)n_features);
        }

        if (Rf_length(s_dim) != 2) {
            UNPROTECT(4);
            Rf_error("dgCMatrix Dim slot must have length 2");
        }

        if (Rf_length(s_p) != n_features + 1) {
            UNPROTECT(4);
            Rf_error("dgCMatrix p slot has incorrect length");
        }

        const int nnz = Rf_length(s_x);

        // Build triplet list
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(nnz);

        for (int j = 0; j < n_features; ++j) {
            int col_start = p_data[j];
            int col_end = p_data[j + 1];

            if (col_start < 0 || col_end > nnz || col_start > col_end) {
                UNPROTECT(4);
                Rf_error("dgCMatrix has invalid p slot (column pointers)");
            }

            for (int idx = col_start; idx < col_end; ++idx) {
                int i = i_data[idx];

                if (i < 0 || i >= n_points) {
                    UNPROTECT(4);
                    Rf_error("dgCMatrix has invalid row index: %d (must be in [0, %ld))",
                             i, (long)n_points);
                }

                double val = x_data[idx];
                triplets.emplace_back(i, j, val);
            }
        }

        X_sparse.resize(n_points, n_features);
        X_sparse.setFromTriplets(triplets.begin(), triplets.end());

        UNPROTECT(4); // s_i, s_p, s_x, s_dim

    } else {
        Rf_error("X must be a numeric matrix or dgCMatrix (from Matrix package)");
    }

    // -------------------- Response Vector y --------------------

    if (TYPEOF(s_y) != REALSXP) {
        Rf_error("y must be a numeric (REAL) vector");
    }

    const R_xlen_t y_len = Rf_xlength(s_y);

    if (y_len != n_points) {
        Rf_error("Length of y (%ld) must equal number of rows in X (%ld)",
                 (long)y_len, (long)n_points);
    }

    vec_t y(y_len);
    const double* y_data = REAL(s_y);

    for (R_xlen_t i = 0; i < y_len; ++i) {
        y[i] = y_data[i];
    }

    // -------------------- Semi-supervised label set y.vertices --------------------
    // If s_y_vertices is non-NULL, it provides 1-based indices of labeled vertices.
    // We convert to 0-based indices and build a boolean mask for downstream use.

    std::vector<index_t> y_vertices;
    std::vector<char> y_mask((size_t)n_points, 1); // default: all labeled
    bool has_y_vertices = (s_y_vertices != R_NilValue);

    if (has_y_vertices) {
        if (TYPEOF(s_y_vertices) != INTSXP) {
            Rf_error("y.vertices must be an integer vector (or NULL)");
        }

        const R_xlen_t m = Rf_xlength(s_y_vertices);
        if (m < 1) {
            Rf_error("y.vertices must have positive length when provided");
        }

        // Build mask: start with all unlabeled
        std::fill(y_mask.begin(), y_mask.end(), (char)0);

        // Copy and validate, convert to 0-based
        std::vector<index_t> idx;
        idx.reserve((size_t)m);

        const int* v_ptr = INTEGER(s_y_vertices);
        for (R_xlen_t j = 0; j < m; ++j) {
            const int v = v_ptr[j];
            if (v == NA_INTEGER) {
                Rf_error("y.vertices cannot contain NA");
            }
            if (v < 1 || v > (int)n_points) {
                Rf_error("y.vertices entries must be in 1..n (got %d, n=%ld)", v, (long)n_points);
            }
            idx.push_back((index_t)(v - 1)); // 0-based
        }

        std::sort(idx.begin(), idx.end());
        for (size_t j = 1; j < idx.size(); ++j) {
            if (idx[j] == idx[j - 1]) {
                Rf_error("y.vertices cannot contain duplicates");
            }
        }

        y_vertices = std::move(idx);
        for (size_t j = 0; j < y_vertices.size(); ++j) {
            y_mask[(size_t)y_vertices[j]] = (char)1;
        }
    }

    // -------------------- Parameter k --------------------

    if (TYPEOF(s_k) != INTSXP || Rf_length(s_k) != 1) {
        Rf_error("k must be a single integer");
    }

    const int k_raw = INTEGER(s_k)[0];

    if (k_raw == NA_INTEGER) {
        Rf_error("k cannot be NA");
    }

    if (k_raw < 2 || k_raw >= n_points) {
        Rf_error("k must satisfy 2 <= k < n (got k=%d, n=%ld)",
                 k_raw, (long)n_points);
    }

    const index_t k = static_cast<index_t>(k_raw);

    // -------------------- Optional precomputed graph (adj.list / weight.list) --------------------

    const bool has_adj_list = !Rf_isNull(s_adj_list);
    const bool has_weight_list = !Rf_isNull(s_weight_list);

    if (has_adj_list != has_weight_list) {
        Rf_error("adj.list and weight.list must be provided together (or both NULL)");
    }

    std::vector<std::vector<index_t>> precomputed_adj_list;
    std::vector<std::vector<double>> precomputed_weight_list;
    const std::vector<std::vector<index_t>>* precomputed_adj_ptr = nullptr;
    const std::vector<std::vector<double>>* precomputed_weight_ptr = nullptr;

    if (has_adj_list) {
        std::vector<std::vector<int>> adj_int = convert_adj_list_from_R(s_adj_list);
        std::vector<std::vector<double>> weight = convert_weight_list_from_R(s_weight_list);

        if (adj_int.size() != static_cast<size_t>(n_points)) {
            Rf_error("adj.list length (%ld) must match nrow(X) (%ld)",
                     (long)adj_int.size(), (long)n_points);
        }
        if (weight.size() != static_cast<size_t>(n_points)) {
            Rf_error("weight.list length (%ld) must match nrow(X) (%ld)",
                     (long)weight.size(), (long)n_points);
        }

        precomputed_adj_list.resize(adj_int.size());
        precomputed_weight_list.resize(weight.size());

        for (size_t i = 0; i < adj_int.size(); ++i) {
            if (adj_int[i].size() != weight[i].size()) {
                Rf_error("Length mismatch at vertex %ld: |adj.list[[i]]|=%ld, |weight.list[[i]]|=%ld",
                         (long)(i + 1),
                         (long)adj_int[i].size(),
                         (long)weight[i].size());
            }

            precomputed_adj_list[i].resize(adj_int[i].size());
            precomputed_weight_list[i] = weight[i];

            for (size_t j = 0; j < adj_int[i].size(); ++j) {
                const int idx = adj_int[i][j];
                if (idx < 0 || idx >= n_points) {
                    Rf_error("adj.list index out of range at vertex %ld (value=%d, expected 0..%ld)",
                             (long)(i + 1), idx, (long)(n_points - 1));
                }
                precomputed_adj_list[i][j] = static_cast<index_t>(idx);
            }
        }

        precomputed_adj_ptr = &precomputed_adj_list;
        precomputed_weight_ptr = &precomputed_weight_list;
    }

    // ==================== Posterior Inference Parameters ====================

    if (!Rf_isLogical(s_with_posterior) || LENGTH(s_with_posterior) != 1) {
        Rf_error("with_posterior must be a single logical value");
    }
    bool with_posterior = static_cast<bool>(LOGICAL(s_with_posterior)[0]);

    if (!Rf_isLogical(s_return_posterior_samples) || LENGTH(s_return_posterior_samples) != 1) {
        Rf_error("return_posterior_samples must be a single logical value");
    }
    bool return_posterior_samples = Rf_asLogical(s_return_posterior_samples);

    if (!Rf_isReal(s_credible_level) || LENGTH(s_credible_level) != 1) {
        Rf_error("credible_level must be a single numeric value");
    }
    double credible_level = REAL(s_credible_level)[0];
    if (with_posterior && (credible_level <= 0.0 || credible_level >= 1.0)) {
        Rf_error("credible_level must be in (0, 1), got %f", credible_level);
    }

    if (!Rf_isInteger(s_n_posterior_samples) || LENGTH(s_n_posterior_samples) != 1) {
        Rf_error("n_posterior_samples must be a single integer value");
    }
    int n_posterior_samples = INTEGER(s_n_posterior_samples)[0];
    if (with_posterior && n_posterior_samples < 100) {
        Rf_error("n_posterior_samples must be at least 100, got %d", n_posterior_samples);
    }

    if (!Rf_isInteger(s_posterior_seed) || LENGTH(s_posterior_seed) != 1) {
        Rf_error("posterior_seed must be a single integer value");
    }
    unsigned int posterior_seed = static_cast<unsigned int>(INTEGER(s_posterior_seed)[0]);

    // -------------------- use.counting.measure --------------------

    if (TYPEOF(s_use_counting_measure) != LGLSXP ||
        Rf_length(s_use_counting_measure) != 1) {
        Rf_error("use.counting.measure must be a single logical value");
    }

    const int ucm_raw = LOGICAL(s_use_counting_measure)[0];

    if (ucm_raw == NA_LOGICAL) {
        Rf_error("use.counting.measure cannot be NA");
    }

    const bool use_counting_measure = (ucm_raw != 0);

    // -------------------- density.normalization --------------------

    if (TYPEOF(s_density_normalization) != REALSXP ||
        Rf_length(s_density_normalization) != 1) {
        Rf_error("density.normalization must be a single numeric value");
    }

    const double density_normalization = REAL(s_density_normalization)[0];

    if (!R_FINITE(density_normalization) || density_normalization < 0.0) {
        Rf_error("density.normalization must be a finite non-negative number (got %.3f)",
                 density_normalization);
    }

    // -------------------- t.diffusion --------------------

    if (TYPEOF(s_t_diffusion) != REALSXP || Rf_length(s_t_diffusion) != 1) {
        Rf_error("t.diffusion must be a single numeric value");
    }

    const double t_diffusion = REAL(s_t_diffusion)[0];

    if (!R_FINITE(t_diffusion) || t_diffusion < 0.0) {
        Rf_error("t.diffusion must be a finite non-negative number (got %.3f)",
                 t_diffusion);
    }

    // -------------------- beta.damping --------------------

    if (TYPEOF(s_beta_damping) != REALSXP || Rf_length(s_beta_damping) != 1) {
        Rf_error("beta.damping must be a single numeric value");
    }

    const double beta_damping = REAL(s_beta_damping)[0];

    if (!R_FINITE(beta_damping) || beta_damping < 0.0) {
        Rf_error("beta.damping must be a finite non-negative number (got %.3f)",
                 beta_damping);
    }

    // -------------------- gamma.modulation --------------------

    if (TYPEOF(s_gamma_modulation) != REALSXP ||
        Rf_length(s_gamma_modulation) != 1) {
        Rf_error("gamma.modulation must be a single numeric value");
    }

    const double gamma_modulation = REAL(s_gamma_modulation)[0];

    if (!R_FINITE(gamma_modulation)) {
        Rf_error("gamma.modulation must be a finite number (got %.3f)",
                 gamma_modulation);
    }

    // -------------------- scale.factor --------------------
    if (TYPEOF(s_t_scale_factor) != REALSXP || Rf_length(s_t_scale_factor) != 1) {
        Rf_error("t.scale.factor must be a single numeric");
    }

    const double t_scale_factor = REAL(s_t_scale_factor)[0];

    if (!R_FINITE(t_scale_factor)) {
        Rf_error("t.scale.factor must be a finite number (got %.3f)",
                 t_scale_factor);
    }

    // -------------------- beta.coefficient.factor --------------------
    if (TYPEOF(s_beta_coefficient_factor) != REALSXP || Rf_length(s_beta_coefficient_factor) != 1) {
        Rf_error("beta.coefficient.factor must be a single numeric");
    }

    const double beta_coefficient_factor = REAL(s_beta_coefficient_factor)[0];

    if (!R_FINITE(beta_coefficient_factor)) {
        Rf_error("beta.coefficient.factor must be a finite number (got %.3f)",
                 beta_coefficient_factor);
    }

    // -------------------- t.update --------------------

    if (TYPEOF(s_t_update_mode) != INTSXP || Rf_length(s_t_update_mode) != 1) {
        Rf_error("t.update.mode must be a single integer (0=fixed, 1=per_iteration)");
    }
    const int t_update_mode = INTEGER(s_t_update_mode)[0];
    if (t_update_mode == NA_INTEGER || (t_update_mode != 0 && t_update_mode != 1)) {
        Rf_error("t.update.mode must be 0 (fixed) or 1 (per_iteration)");
    }

    if (TYPEOF(s_t_update_max_mult) != REALSXP || Rf_length(s_t_update_max_mult) != 1) {
        Rf_error("t.update.max.mult must be a single numeric value >= 1.0");
    }
    const double t_update_max_mult = REAL(s_t_update_max_mult)[0];
    if (!R_FINITE(t_update_max_mult) || t_update_max_mult < 1.0) {
        Rf_error("t.update.max.mult must be finite and >= 1.0 (got %.3f)", t_update_max_mult);
    }

    // -------------------- n.eigenpairs --------------------

    if (TYPEOF(s_n_eigenpairs) != INTSXP || Rf_length(s_n_eigenpairs) != 1) {
        Rf_error("n.eigenpairs must be a single integer");
    }

    const int n_eigenpairs_raw = INTEGER(s_n_eigenpairs)[0];

    if (n_eigenpairs_raw == NA_INTEGER) {
        Rf_error("n.eigenpairs cannot be NA");
    }

    if (n_eigenpairs_raw < 10 || n_eigenpairs_raw > n_points) {
        Rf_error("n.eigenpairs must satisfy 10 <= n.eigenpairs <= n (got %d, n=%ld)",
                 n_eigenpairs_raw, (long)n_points);
    }

    const int n_eigenpairs = n_eigenpairs_raw;

    // -------------------- filter.type --------------------

    if (TYPEOF(s_filter_type) != STRSXP || Rf_length(s_filter_type) != 1) {
        Rf_error("filter.type must be a single string");
    }

    const char* filter_str = CHAR(STRING_ELT(s_filter_type, 0));

    if (filter_str == nullptr || strlen(filter_str) == 0) {
        Rf_error("filter.type cannot be empty string");
    }

    rdcx_filter_type_t filter_type;

    if (strcmp(filter_str, "heat_kernel") == 0) {
        filter_type = rdcx_filter_type_t::HEAT_KERNEL;
    } else if (strcmp(filter_str, "tikhonov") == 0) {
        filter_type = rdcx_filter_type_t::TIKHONOV;
    } else if (strcmp(filter_str, "cubic_spline") == 0) {
        filter_type = rdcx_filter_type_t::CUBIC_SPLINE;
    } else if (strcmp(filter_str, "gaussian") == 0) {
        filter_type = rdcx_filter_type_t::GAUSSIAN;
    } else if (strcmp(filter_str, "exponential") == 0) {
        filter_type = rdcx_filter_type_t::EXPONENTIAL;
    } else if (strcmp(filter_str, "butterworth") == 0) {
        filter_type = rdcx_filter_type_t::BUTTERWORTH;
    } else {
        Rf_error("filter.type must be 'heat_kernel', 'tikhonov', 'cubic_spline', "
                 "'gaussian', 'exponential', or 'butterworth' (got '%s')",
                 filter_str);
    }

    // -------------------- epsilon.y --------------------

    if (TYPEOF(s_epsilon_y) != REALSXP || Rf_length(s_epsilon_y) != 1) {
        Rf_error("epsilon.y must be a single numeric value");
    }

    const double epsilon_y = REAL(s_epsilon_y)[0];

    if (!R_FINITE(epsilon_y) || epsilon_y <= 0.0) {
        Rf_error("epsilon.y must be a finite positive number (got %.3e)",
                 epsilon_y);
    }

    // -------------------- epsilon.rho --------------------

    if (TYPEOF(s_epsilon_rho) != REALSXP || Rf_length(s_epsilon_rho) != 1) {
        Rf_error("epsilon.rho must be a single numeric value");
    }

    const double epsilon_rho = REAL(s_epsilon_rho)[0];

    if (!R_FINITE(epsilon_rho) || epsilon_rho <= 0.0) {
        Rf_error("epsilon.rho must be a finite positive number (got %.3e)",
                 epsilon_rho);
    }

    // -------------------- max.iterations --------------------

    if (TYPEOF(s_max_iterations) != INTSXP || Rf_length(s_max_iterations) != 1) {
        Rf_error("max.iterations must be a single integer");
    }

    const int max_iterations = INTEGER(s_max_iterations)[0];

    if (max_iterations == NA_INTEGER) {
        Rf_error("max.iterations cannot be NA");
    }

    if (max_iterations < 1) {
        Rf_error("max.iterations must be at least 1 (got %d)", max_iterations);
    }

    // -------------------- max.ratio.threshold --------------------

    if (TYPEOF(s_max_ratio_threshold) != REALSXP || Rf_length(s_max_ratio_threshold) != 1) {
        Rf_error("max.ratio.threshold must be a single numeric value");
    }

    const double max_ratio_threshold = REAL(s_max_ratio_threshold)[0];

    if (!R_FINITE(max_ratio_threshold) || max_ratio_threshold < 0.0) {
        Rf_error("max.ratio.threshold must be a finite non-negative number (got %.3e)",
                 max_ratio_threshold);
    }

    // -------------------- path.edge.ratio.percentile --------------------
    if (TYPEOF(s_path_edge_ratio_percentile) != REALSXP || Rf_length(s_path_edge_ratio_percentile) != 1) {
        Rf_error("path.edge.ratio.percentile must be a single numeric value");
    }

    const double path_edge_ratio_percentile = Rf_asReal(s_path_edge_ratio_percentile);

    // -------------------- threshold.percentile --------------------

    if (TYPEOF(s_threshold_percentile) != REALSXP || Rf_length(s_threshold_percentile) != 1) {
        Rf_error("threshold.percentile must be a single numeric value");
    }

    const double threshold_percentile = REAL(s_threshold_percentile)[0];

    if (!R_FINITE(threshold_percentile) || threshold_percentile < 0.0) {
        Rf_error("threshold.percentile must be a non-negative number (got %.3e)",
                 threshold_percentile);
    }

    // -------------------- density.alpha --------------------

    if (TYPEOF(s_density_alpha) != REALSXP || Rf_length(s_density_alpha) != 1) {
        Rf_error("density.alpha must be a single numeric value");
    }

    const double density_alpha = REAL(s_density_alpha)[0];

    if (!R_FINITE(density_alpha) || density_alpha < 0.1 || density_alpha > 2.0) {
        Rf_error("density.alpha must be finite and in [0.1, 2] (got %.3f)", density_alpha);
    }

    // -------------------- s_compute_extremality --------------------

    if (TYPEOF(s_compute_extremality) != LGLSXP ||
        Rf_length(s_compute_extremality) != 1) {
        Rf_error("compute_extremality must be a single logical value");
    }

    const int compute_extremality_int = LOGICAL(s_compute_extremality)[0];

    if (compute_extremality_int == NA_LOGICAL) {
        Rf_error("compute_extremality cannot be NA");
    }

    const bool compute_extremality = (compute_extremality_int != 0);

    // -------------------- s_p_threshold --------------------

    if (TYPEOF(s_p_threshold) != REALSXP || Rf_length(s_p_threshold) != 1) {
        Rf_error("p_threshold must be a single numeric value");
    }

    const double p_threshold = REAL(s_p_threshold)[0];

    // When p_threshold is 0, skip hop radius computation
    // Otherwise validate it's in (0,1]
    if (p_threshold < 0.0 || p_threshold > 1.0) {
        Rf_error("p_threshold must be in [0,1], got %.3f", p_threshold);
    }

    if (compute_extremality && p_threshold > 0.0 && !R_FINITE(p_threshold)) {
        Rf_error("p_threshold must be finite when computing hop radii");
    }

    // -------------------- s_max_hop --------------------

    if (TYPEOF(s_max_hop) != INTSXP || Rf_length(s_max_hop) != 1) {
        Rf_error("max_hop must be a single integer value");
    }

    const int max_hop_int = INTEGER(s_max_hop)[0];

    if (max_hop_int == NA_INTEGER) {
        Rf_error("max_hop cannot be NA");
    }

    if (max_hop_int <= 0) {
        Rf_error("max_hop must be positive, got %d", max_hop_int);
    }

    const size_t max_hop = static_cast<size_t>(max_hop_int);

    // -------------------- s_knn_cache_mode --------------------

    if (TYPEOF(s_knn_cache_mode) != INTSXP || Rf_length(s_knn_cache_mode) != 1) {
        Rf_error("knn.cache.mode must be a single integer in {0,1,2,3}");
    }

    const int knn_cache_mode = INTEGER(s_knn_cache_mode)[0];

    if (knn_cache_mode == NA_INTEGER || knn_cache_mode < 0 || knn_cache_mode > 3) {
        Rf_error("knn.cache.mode must be 0 (none), 1 (read), 2 (write), or 3 (readwrite)");
    }

    // -------------------- s_knn_cache_path --------------------

    std::string knn_cache_path;
    if (!Rf_isNull(s_knn_cache_path)) {
        if (TYPEOF(s_knn_cache_path) != STRSXP ||
            Rf_length(s_knn_cache_path) != 1 ||
            STRING_ELT(s_knn_cache_path, 0) == NA_STRING) {
            Rf_error("knn.cache.path must be NULL or a non-empty character scalar");
        }
        const char* path_cstr = CHAR(STRING_ELT(s_knn_cache_path, 0));
        if (path_cstr == nullptr || path_cstr[0] == '\0') {
            Rf_error("knn.cache.path must be NULL or a non-empty character scalar");
        }
        knn_cache_path.assign(path_cstr);
    }

    if (knn_cache_mode != 0 && knn_cache_path.empty()) {
        Rf_error("knn.cache.path must be provided when knn.cache.mode is not 'none'");
    }

    // -------------------- s_dense_fallback_mode --------------------

    if (TYPEOF(s_dense_fallback_mode) != INTSXP || Rf_length(s_dense_fallback_mode) != 1) {
        Rf_error("dense_fallback_mode must be a single integer (0=auto, 1=never, 2=always)");
    }

    const int dense_fallback_mode = INTEGER(s_dense_fallback_mode)[0];

    if (dense_fallback_mode == NA_INTEGER ||
        (dense_fallback_mode != 0 && dense_fallback_mode != 1 && dense_fallback_mode != 2)) {
        Rf_error("dense_fallback_mode must be 0 (auto), 1 (never), or 2 (always)");
    }

    if (TYPEOF(s_triangle_policy_mode) != INTSXP || Rf_length(s_triangle_policy_mode) != 1) {
        Rf_error("triangle_policy_mode must be a single integer (0=auto, 1=never, 2=always)");
    }

    const int triangle_policy_mode = INTEGER(s_triangle_policy_mode)[0];

    if (triangle_policy_mode == NA_INTEGER ||
        (triangle_policy_mode != 0 && triangle_policy_mode != 1 && triangle_policy_mode != 2)) {
        Rf_error("triangle_policy_mode must be 0 (auto), 1 (never), or 2 (always)");
    }

    // -------------------- density.epsilon --------------------

    if (TYPEOF(s_density_epsilon) != REALSXP || Rf_length(s_density_epsilon) != 1) {
        Rf_error("density.epsilon must be a single numeric value");
    }

    const double density_epsilon = REAL(s_density_epsilon)[0];

    if (!R_FINITE(density_epsilon) || density_epsilon <= 0.0) {
        Rf_error("density.epsilon must be a finite positive number (got %.3e)",
                 density_epsilon);
    }

    // -------------------- s_clamp_dk --------------------

    if (TYPEOF(s_clamp_dk) != LGLSXP ||
        Rf_length(s_clamp_dk) != 1) {
        Rf_error("clamp_dk must be a single logical value");
    }

    const int clamp_dk_int = LOGICAL(s_clamp_dk)[0];

    if (clamp_dk_int == NA_LOGICAL) {
        Rf_error("clamp_dk cannot be NA");
    }

    const bool clamp_dk = (clamp_dk_int != 0);

    
    // -------------------- dk_clamp_median_factor --------------------

    if (!Rf_isNumeric(s_dk_clamp_median_factor) || Rf_length(s_dk_clamp_median_factor) != 1) {
        Rf_error("dk_clamp_median_factor must be a single numeric value");
    }

    const double dk_clamp_median_factor = Rf_asReal(s_dk_clamp_median_factor);

    if (!R_FINITE(dk_clamp_median_factor) || dk_clamp_median_factor <= 1.0) {
        Rf_error("dk_clamp_median_factor must be a finite number greater than 1 (got %.3e)",
                 dk_clamp_median_factor);
    }

    // -------------------- target_weight_ratio --------------------

    if (!Rf_isNumeric(s_target_weight_ratio) || Rf_length(s_target_weight_ratio) != 1) {
        Rf_error("target_weight_ratio must be a single numeric value");
    }

    const double target_weight_ratio = Rf_asReal(s_target_weight_ratio);

    if (!R_FINITE(target_weight_ratio) || target_weight_ratio <= 1.0) {
        Rf_error("target_weight_ratio must be a finite number greater than 1 (got %.3e)",
                 target_weight_ratio);
    }

   // -------------------- pathological_ratio_threshold --------------------

    if (!Rf_isNumeric(s_pathological_ratio_threshold) || Rf_length(s_pathological_ratio_threshold) != 1) {
        Rf_error("pathological_ratio_threshold must be a single numeric value");
    }

    const double pathological_ratio_threshold = Rf_asReal(s_pathological_ratio_threshold);

    if (!R_FINITE(pathological_ratio_threshold) || pathological_ratio_threshold <= 1.0) {
        Rf_error("pathological_ratio_threshold must be a finite number greater than 1 (got %.3e)",
                 pathological_ratio_threshold);
    }

    if (pathological_ratio_threshold <= target_weight_ratio) {
        Rf_error("pathological_ratio_threshold must be greater than target_weight_ratio "
                 "(got %.3e <= %.3e).",
                 pathological_ratio_threshold, target_weight_ratio);
    }
    
    // -------------------- s_verbose_level --------------------

    if (TYPEOF(s_verbose_level) != INTSXP || Rf_length(s_verbose_level) != 1) {
        Rf_error("verbose_level must be a single integer");
    }

    const int verbose_level_int = INTEGER(s_verbose_level)[0];

    if (verbose_level_int == NA_INTEGER) {
        Rf_error("verbose_level cannot be NA");
    }

    if (verbose_level_int < 0) {
        Rf_error("verbose_level must be at least 0 (got %d)", verbose_level_int);
    }

    if (verbose_level_int > 3) {
        Rf_error("verbose_level must be at most 3 (got %d)", verbose_level_int);
    }

    const verbose_level_t verbose_level = vl_from_int(verbose_level_int);

    // ================================================================
    // PART II: CALL MEMBER FUNCTION
    // ================================================================

    riem_dcx_t dcx;  // Stack allocation now! No need for new/delete

    try {
    dcx.fit_rdgraph_regression(
        X_sparse,
        y,
        y_vertices,
        k,
        use_counting_measure,
        density_normalization,
        t_diffusion,
        beta_damping,
        gamma_modulation,
        t_scale_factor,
        beta_coefficient_factor,
        t_update_mode,
        t_update_max_mult,
        n_eigenpairs,
        filter_type,
        epsilon_y,
        epsilon_rho,
        max_iterations,
        max_ratio_threshold,
        path_edge_ratio_percentile,
        threshold_percentile,
        density_alpha,
        density_epsilon,
        clamp_dk,
        dk_clamp_median_factor,
        target_weight_ratio,
        pathological_ratio_threshold,
        knn_cache_path,
        knn_cache_mode,
        dense_fallback_mode,
        triangle_policy_mode,
        verbose_level,
        precomputed_adj_ptr,
        precomputed_weight_ptr
        );

    } catch (const std::exception& e) {
        Rf_error("Regression fitting failed: %s", e.what());
    }

    // ================================================================
    // PART III: BUILD RESULT LIST
    // ================================================================

    const Eigen::Index n = y.size();

    // ---------- Select fitted values based on minimum GCV ----------
    // Instead of simply using the last iteration's fitted values, we
    // select the iteration that minimizes the GCV score across all
    // iterations. This implements proper model selection and can
    // prevent overfitting when the iterative refinement continues
    // beyond the optimal point.

    bool has_fitted_values = !dcx.sig.y_hat_hist.empty();
    bool has_gcv_history = !dcx.gcv_history.iterations.empty();
    vec_t y_hat_raw;
    int optimal_iteration = -1;

    if (has_fitted_values && has_gcv_history) {
        // Verify consistency between fitted values history and GCV history
        if (dcx.sig.y_hat_hist.size() != dcx.gcv_history.iterations.size()) {
            Rf_warning("Inconsistent history sizes: y_hat_hist has %d entries, "
                       "gcv_history has %d entries. Using last iteration.",
                       (int)dcx.sig.y_hat_hist.size(),
                       (int)dcx.gcv_history.iterations.size());
            y_hat_raw = dcx.sig.y_hat_hist.back();
            optimal_iteration = (int)dcx.sig.y_hat_hist.size() - 1;
        } else {
            // Find iteration with minimum GCV score
            double min_gcv = std::numeric_limits<double>::max();
            size_t min_idx = 0;

            for (size_t i = 0; i < dcx.gcv_history.iterations.size(); ++i) {
                double gcv_score = dcx.gcv_history.iterations[i].gcv_optimal;
                if (gcv_score < min_gcv) {
                    min_gcv = gcv_score;
                    min_idx = i;
                }
            }

            // Use fitted values from the iteration with minimum GCV
            y_hat_raw = dcx.sig.y_hat_hist[min_idx];
            optimal_iteration = (int)min_idx;
        }
    } else if (has_fitted_values) {
        // GCV history unavailable (should not happen in normal operation)
        // Fall back to last iteration
        Rf_warning("GCV history unavailable. Using last iteration's fitted values.");
        y_hat_raw = dcx.sig.y_hat_hist.back();
        optimal_iteration = (int)dcx.sig.y_hat_hist.size() - 1;
    } else {
        // Early termination: use observed y as placeholder
        // This allows the result structure to be consistent
        y_hat_raw = y;
        optimal_iteration = -1;
    }

    // ---------- Main result list ----------

    int n_components = 12;
    if (compute_extremality) {
        n_components++;
    }
    if (with_posterior) {
        n_components++;  // Add posterior component
    }

    #if 0
    auto is_binary01 = [](const std::vector<double>& yy, double tol = 1e-12) -> bool {
        for (double v : yy) {
            if (!(std::fabs(v) <= tol || std::fabs(v - 1.0) <= tol)) {
                return false;
            }
        }
        return true;
    };
    const bool y_binary = is_binary01(y);
    #endif

    bool y_binary = (std::set<double>(y.begin(), y.end()) == std::set<double>{0.0, 1.0});

    if (y_binary) {
        n_components++;  // Add posterior component
    }

    SEXP result = PROTECT(Rf_allocVector(VECSXP, n_components));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, n_components));
    int component_idx = 0;

    // Component 1: fitted.values
    SET_STRING_ELT(names, component_idx, Rf_mkChar("fitted.values"));
    SEXP s_raw_fitted = PROTECT(Rf_allocVector(REALSXP, n));
    for (Eigen::Index i = 0; i < n; ++i) {
        REAL(s_raw_fitted)[i] = y_hat_raw[i];
    }
    SET_VECTOR_ELT(result, component_idx++, s_raw_fitted);
    UNPROTECT(1);

    // Component 2: fitted.values (after preprocessing: winsorized + tie-broken)
    if (y_binary) {
        // -------------------- Preprocess Fitted Values --------------------
        // Create a copy for preprocessing (leaves y_hat_raw unchanged)
        vec_t y_hat = y_hat_raw;

        // Apply preprocessing pipeline for binary conditional expectations
        // This handles slow diffusion convergence at isolated vertices
        dcx.prepare_binary_cond_exp(
            y_hat,
            0.02,     // p_right: winsorize top 2%
            true,     // apply_right_winsorization
            1e-10,    // noise_scale: relative to range
            123,      // seed: for reproducibility
            verbose_level   // verbose: match user's verbosity setting
            );

        if (vl_at_least(verbose_level, verbose_level_t::TRACE)) {
            Rprintf("Preprocessing complete:\n");
            Rprintf("  Raw range: [%.6f, %.6f]\n",
                    y_hat_raw.minCoeff(), y_hat_raw.maxCoeff());
            Rprintf("  Preprocessed range: [%.6f, %.6f]\n",
                    y_hat.minCoeff(), y_hat.maxCoeff());
        }

        SET_STRING_ELT(names, component_idx, Rf_mkChar("binary.clipped.fitted.values"));
        SEXP s_fitted = PROTECT(Rf_allocVector(REALSXP, n));
        for (Eigen::Index i = 0; i < n; ++i) {
            REAL(s_fitted)[i] = y_hat[i];
        }
        SET_VECTOR_ELT(result, component_idx++, s_fitted);
        UNPROTECT(1);

        // Component 3: residuals (computed using preprocessed fitted values)
        SET_STRING_ELT(names, component_idx, Rf_mkChar("residuals"));
        SEXP s_resid = PROTECT(Rf_allocVector(REALSXP, n));

        for (Eigen::Index i = 0; i < n; ++i) {
            if (has_y_vertices && y_mask[(size_t)i] == 0) {
                REAL(s_resid)[i] = NA_REAL;
                continue;
            }

            double residual;
            if (dcx.sig.y.size() == n) {
                residual = dcx.sig.y[i] - y_hat[i];
            } else {
                residual = y[i] - y_hat[i];
            }
            REAL(s_resid)[i] = residual;
        }

        SET_VECTOR_ELT(result, component_idx++, s_resid);
        UNPROTECT(1);

    } else {
        // Component 3: residuals (computed using preprocessed fitted values)
        SET_STRING_ELT(names, component_idx, Rf_mkChar("residuals"));
        SEXP s_resid = PROTECT(Rf_allocVector(REALSXP, n));

        for (Eigen::Index i = 0; i < n; ++i) {
            if (has_y_vertices && y_mask[(size_t)i] == 0) {
                REAL(s_resid)[i] = NA_REAL;
                continue;
            }

            double residual;
            if (dcx.sig.y.size() == n) {
                residual = dcx.sig.y[i] - y_hat_raw[i];
            } else {
                residual = y[i] - y_hat_raw[i];
            }
            REAL(s_resid)[i] = residual;
        }

        SET_VECTOR_ELT(result, component_idx++, s_resid);
        UNPROTECT(1);
    }

    // Component 4: optimal.iteration
    SET_STRING_ELT(names, component_idx, Rf_mkChar("optimal.iteration"));
    SET_VECTOR_ELT(result, component_idx++, Rf_ScalarInteger(optimal_iteration + 1));  // R uses 1-based indexing

    // Component 5: graph (nested list)
    SET_STRING_ELT(names, component_idx, Rf_mkChar("graph"));
    SEXP s_graph = PROTECT(create_graph_component(dcx));
    SET_VECTOR_ELT(result, component_idx++, s_graph);
    UNPROTECT(1);

    // Component 6: iteration (nested list)
    SET_STRING_ELT(names, component_idx, Rf_mkChar("iteration"));
    SEXP s_iteration = PROTECT(create_iteration_component(dcx));
    SET_VECTOR_ELT(result, component_idx++, s_iteration);
    UNPROTECT(1);

    // Component 7: parameters (nested list)
    SET_STRING_ELT(names, component_idx, Rf_mkChar("parameters"));
    SEXP s_params = PROTECT(create_parameters_component(
                                k, use_counting_measure, density_normalization,
                                t_diffusion, beta_damping, gamma_modulation,
                                t_scale_factor, beta_coefficient_factor,
                                n_eigenpairs, filter_type, epsilon_y, epsilon_rho, max_iterations,
                                density_alpha, density_epsilon,
                                dense_fallback_mode, triangle_policy_mode
                                ));
    SET_VECTOR_ELT(result, component_idx++, s_params);
    UNPROTECT(1);

    // Component 8: reference.measure (nested list)
    SET_STRING_ELT(names, component_idx, Rf_mkChar("reference.measure"));
    SEXP s_ref = PROTECT(create_reference_measure_component(dcx));
    SET_VECTOR_ELT(result, component_idx++, s_ref);
    UNPROTECT(1);

    // Component 9: y (original response)
    SET_STRING_ELT(names, component_idx, Rf_mkChar("y"));
    SEXP s_y_copy = PROTECT(Rf_allocVector(REALSXP, n));

    // CRITICAL FIX: Check if sig.y was populated during fitting
    if (dcx.sig.y.size() == n) {
        // Use the stored response from dcx
        for (Eigen::Index i = 0; i < n; ++i) {
            REAL(s_y_copy)[i] = dcx.sig.y[i];
        }
    } else {
        // Early termination: sig.y not populated, use input y
        for (Eigen::Index i = 0; i < n; ++i) {
            REAL(s_y_copy)[i] = y[i];
        }
    }
    SET_VECTOR_ELT(result, component_idx++, s_y_copy);
    UNPROTECT(1);

    // Component 10: gcv (nested list)
    SET_STRING_ELT(names, component_idx, Rf_mkChar("gcv"));
    SEXP s_gcv = PROTECT(create_gcv_component(dcx));
    SET_VECTOR_ELT(result, component_idx++, s_gcv);
    UNPROTECT(1);

    // Component 11: density
    SET_STRING_ELT(names, component_idx, Rf_mkChar("density"));
    SEXP s_density = PROTECT(create_density_history_component(dcx));
    SET_VECTOR_ELT(result, component_idx++, s_density);
    UNPROTECT(1);

    // Component 12: gamma.selection
    SET_STRING_ELT(names, component_idx, Rf_mkChar("gamma.selection"));
    SEXP s_gamma_sel = PROTECT(create_gamma_selection_component(dcx));
    SET_VECTOR_ELT(result, component_idx++, s_gamma_sel);
    UNPROTECT(1);

    // Component 13: spectral
    SET_STRING_ELT(names, component_idx, Rf_mkChar("spectral"));
    SEXP s_spectral = PROTECT(create_spectral_component(dcx, optimal_iteration, filter_type));
    SET_VECTOR_ELT(result, component_idx++, s_spectral);
    UNPROTECT(1);

    // Component 14: extremality.scores
    if (compute_extremality) {
        SET_STRING_ELT(names, component_idx, Rf_mkChar("extremality"));

        SEXP s_extremality = PROTECT(create_extremality_component(
                                         dcx,
                                         y_hat_raw,
                                         p_threshold,
                                         max_hop
                                         ));
        SET_VECTOR_ELT(result, component_idx++, s_extremality);
        UNPROTECT(1);
    }

    // ================================================================
    // POSTERIOR INFERENCE (OPTIONAL)
    // ================================================================

    if (with_posterior) {
        // Extract necessary components from dcx
        // Assumes dcx has public access to spectral_cache and optimal iteration results
        const Eigen::MatrixXd& V = dcx.spectral_cache.eigenvectors;
        const vec_t& eigenvalues = dcx.spectral_cache.eigenvalues;
        const vec_t& filtered_eigenvalues = dcx.spectral_cache.filtered_eigenvalues;
        double eta_opt = dcx.gcv_history.iterations[optimal_iteration].eta_optimal;

        try {
            posterior_summary_t posterior = dcx.compute_posterior_summary(
                V,
                eigenvalues,
                filtered_eigenvalues,
                y,
                y_hat_raw,
                eta_opt,
                credible_level,
                n_posterior_samples,
                posterior_seed,
                return_posterior_samples
                );

            // Package into SEXP and add to result list
            SEXP s_posterior_component = PROTECT(create_posterior_component(posterior));

            // Add as new field to result list
            // Assumes result list was created earlier with space for posterior component
            // You'll need to adjust list size and indexing accordingly
            SET_VECTOR_ELT(result, component_idx, s_posterior_component);
            SET_STRING_ELT(names, component_idx, Rf_mkChar("posterior"));

            UNPROTECT(1);  // s_posterior_component

        } catch (const std::exception& e) {
            Rf_warning("Posterior computation failed: %s. Continuing without posterior inference.", e.what());
            // Don't fail entire computation, just skip posterior
        }
    }

    // Set names attribute
    Rf_setAttrib(result, R_NamesSymbol, names);

    // Set class attribute
    SEXP class_attr = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(class_attr, 0, Rf_mkChar("knn.riem.fit"));
    SET_STRING_ELT(class_attr, 1, Rf_mkChar("list"));
    Rf_setAttrib(result, R_ClassSymbol, class_attr);

    UNPROTECT(3); // class_attr, names, result

    return result;
}
