/**
 * @file riem_dcx_initialization.cpp
 *
 * @brief Initialization routines for Riemannian simplicial complex regression
 *
 * This file implements the geometric initialization phase that constructs
 * the simplicial complex from k-NN data, computes initial densities, and
 * assembles the initial Hodge Laplacian operators.
 *
 * DEBUG INSTRUMENTATION
 * =====================
 * This file contains debug serialization code controlled by DEBUG_INITIALIZE_FROM_KNN.
 *
 * To enable debugging:
 * 1. Set DEBUG_INITIALIZE_FROM_KNN to 1 in this file
 * 2. Create directory: mkdir -p /tmp/gflow_debug/initialize_from_knn
 * 3. Recompile and run
 * 4. Debug files will be written to /tmp/gflow_debug/initialize_from_knn/
 *
 * Debug files include:
 * - phase_1a_knn_result.bin: k-NN computation results
 * - phase_1b_edges_pre_pruning.bin: Edge list before geometric pruning
 * - phase_1c_edges_post_pruning.bin: Edge list after geometric pruning
 * - phase_1e_connectivity.bin: Final connectivity information
 *
 * Use R/debug_comparison.R::compare_debug_outputs() to compare with create_iknn_graph
 *
 * IMPORTANT: Debug code must remain disabled (flag = 0) for CRAN submission.
 */

// ============================================================
// DEBUG CONTROL
// ============================================================
// Creates files in /tmp/gflow_debug/initialize_from_knn/
#define DEBUG_INITIALIZE_FROM_KNN 0

#include "riem_dcx.hpp"
#include "omp_compat.h"
#include "progress_utils.hpp" // for elapsed_time
#include <algorithm>
#include <limits>
#include <unordered_set>
#include <numeric>
#include <cmath>
#include <array>

// ================================================================
// WEIGHT COMPRESSION HELPERS (reference measure stabilization)
// ================================================================
static inline double median_of_copy(std::vector<double> x, bool average_even = false) {
    const size_t n = x.size();
    if (n == 0) return std::numeric_limits<double>::quiet_NaN();

    auto mid = x.begin() + (n / 2);
    std::nth_element(x.begin(), mid, x.end());
    double m_hi = *mid;

    if (!average_even || (n % 2 == 1)) {
        return m_hi; // odd n -> exact median; even n -> upper median
    }

    auto mid_lo = x.begin() + (n / 2 - 1);
    std::nth_element(x.begin(), mid_lo, x.end());
    double m_lo = *mid_lo;

    return 0.5 * (m_lo + m_hi);
}

static inline double positive_floor(double scale_hint) {
    const double base = 1e-12;
    if (!std::isfinite(scale_hint) || scale_hint <= 0.0) scale_hint = 1.0;
    return base * std::max(1.0, scale_hint);
}

// log compression (robust, soft)
static inline void validate_and_log_compress_weights_in_place(
    std::vector<double>& w,
    bool average_even_median = false,
    bool clamp_nonpositive_to_zero = false
) {
    const size_t n = w.size();
    if (n == 0) return;

    for (size_t i = 0; i < n; ++i) {
        if (!std::isfinite(w[i])) {
            Rf_error("reference_measure: vertex_weights contains non-finite values.");
        }
        if (w[i] < 0.0) {
            if (clamp_nonpositive_to_zero) w[i] = 0.0;
            else Rf_error("reference_measure: vertex_weights contains negative values.");
        }
    }

    double w_median = median_of_copy(w, average_even_median);
    if (!std::isfinite(w_median) || w_median <= 0.0) {
        w_median = positive_floor(1.0);
    } else {
        w_median = std::max(w_median, positive_floor(w_median));
    }

    for (size_t i = 0; i < n; ++i) {
        w[i] = std::log1p(w[i] / w_median);
    }
}

// power compression to enforce exact ratio cap
static inline void power_compress_weights_to_ratio_in_place(
    std::vector<double>& w,
    double target_ratio
) {
    if (!(target_ratio > 1.0) || w.empty()) return;

    const double w_min = *std::min_element(w.begin(), w.end());
    const double w_max = *std::max_element(w.begin(), w.end());

    if (!std::isfinite(w_min) || !std::isfinite(w_max) || w_min <= 0.0) {
        Rf_error("reference_measure: non-finite or non-positive weights encountered.");
    }

    const double ratio = w_max / w_min;
    if (!(ratio > target_ratio)) return;

    const double q = std::log(target_ratio) / std::log(ratio); // in (0,1)
    for (size_t i = 0; i < w.size(); ++i) {
        w[i] = std::pow(w[i], q);
    }
}

static inline size_t intersection_size_of_sets(
    const std::unordered_set<index_t>& a,
    const std::unordered_set<index_t>& b
) {
    const auto* smaller = &a;
    const auto* larger = &b;
    if (a.size() > b.size()) {
        smaller = &b;
        larger = &a;
    }

    size_t count = 0;
    for (index_t v : *smaller) {
        if (larger->find(v) != larger->end()) {
            ++count;
        }
    }
    return count;
}

void riem_dcx_t::initialize_from_graph(
    const std::vector<std::vector<index_t>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    bool use_counting_measure,
    double density_normalization,
    double density_alpha,
    double density_epsilon,
    bool clamp_dk,
    double dk_clamp_median_factor,
    double target_weight_ratio,
    double pathological_ratio_threshold,
    index_t neighborhood_k,
    triangle_policy_t triangle_policy,
    double gamma_modulation,
    verbose_level_t verbose_level
) {
    const size_t n_points = adj_list.size();
    auto phase_time = std::chrono::steady_clock::now();

    if (n_points == 0) {
        Rf_error("Precomputed graph is empty.");
    }
    if (weight_list.size() != n_points) {
        Rf_error("Precomputed graph adjacency/weight list size mismatch.");
    }
    if (neighborhood_k < 1) {
        Rf_error("neighborhood_k must be at least 1 in precomputed graph initialization.");
    }

    const size_t neighborhood_k_sz = static_cast<size_t>(neighborhood_k);

    // ================================================================
    // PHASE 0: INITIALIZE DIMENSION STRUCTURE
    // ================================================================

    pmax = -1;
    pmax = 0;
    g.M.resize(1);
    g.M_solver.resize(1);
    L.B.resize(1);
    L.L.resize(1);
    g.M[0] = spmat_t();
    g.M_solver[0].reset();
    L.B[0] = spmat_t();
    L.L[0] = spmat_t();

    // ================================================================
    // PHASE 1: PRECOMPUTED GRAPH INGESTION
    // ================================================================
    if (vl_at_least(verbose_level, verbose_level_t::PROGRESS)) {
        progress_log_inline("  Phase 1: ingest precomputed graph ... ");
        phase_time = std::chrono::steady_clock::now();
    }

    neighbor_sets.clear();
    neighbor_sets.resize(n_points);

    for (size_t i = 0; i < n_points; ++i) {
        const auto& nbrs = adj_list[i];
        const auto& wts = weight_list[i];

        if (nbrs.size() != wts.size()) {
            Rf_error("Precomputed graph row %zu has mismatched adjacency/weight lengths.", i + 1);
        }

        std::unordered_set<index_t> seen;
        std::vector<std::pair<double, index_t>> by_dist;
        by_dist.reserve(nbrs.size() + 1);
        by_dist.emplace_back(0.0, static_cast<index_t>(i)); // include self to mirror kNN path

        for (size_t j = 0; j < nbrs.size(); ++j) {
            const index_t v = nbrs[j];
            const double w = wts[j];

            if (v >= n_points) {
                Rf_error("Precomputed graph row %zu has invalid neighbor index %zu (n=%zu).",
                         i + 1, static_cast<size_t>(v), n_points);
            }
            if (v == i) {
                Rf_error("Precomputed graph row %zu contains self-loop.", i + 1);
            }
            if (!std::isfinite(w) || w <= 0.0) {
                Rf_error("Precomputed graph row %zu contains non-finite or non-positive edge weights.",
                         i + 1);
            }
            if (!seen.insert(v).second) {
                Rf_error("Precomputed graph row %zu contains duplicate neighbor %zu.",
                         i + 1, static_cast<size_t>(v));
            }

            by_dist.emplace_back(w, v);
        }

        std::sort(
            by_dist.begin(),
            by_dist.end(),
            [](const std::pair<double, index_t>& lhs, const std::pair<double, index_t>& rhs) {
                if (lhs.first < rhs.first) return true;
                if (lhs.first > rhs.first) return false;
                return lhs.second < rhs.second;
            }
        );

        const size_t k_effective = std::min(neighborhood_k_sz, by_dist.size());
        std::unordered_set<index_t> local_neighbors;
        local_neighbors.reserve(k_effective);
        for (size_t t = 0; t < k_effective; ++t) {
            local_neighbors.insert(by_dist[t].second);
        }
        neighbor_sets[i] = std::move(local_neighbors);
    }

    // Validate undirected symmetry with reciprocal-weight consistency
    constexpr double symmetry_tol = 1e-10;
    for (size_t i = 0; i < n_points; ++i) {
        for (size_t j_idx = 0; j_idx < adj_list[i].size(); ++j_idx) {
            const index_t j = adj_list[i][j_idx];
            const double w_ij = weight_list[i][j_idx];

            const auto& nbrs_j = adj_list[j];
            const auto it = std::find(nbrs_j.begin(), nbrs_j.end(), static_cast<index_t>(i));
            if (it == nbrs_j.end()) {
                Rf_error("Precomputed graph must be undirected: edge %zu -> %zu missing reciprocal edge.",
                         i + 1, static_cast<size_t>(j) + 1);
            }
            const size_t rev_idx = static_cast<size_t>(std::distance(nbrs_j.begin(), it));
            const double w_ji = weight_list[j][rev_idx];
            const double scale = std::max({1.0, std::abs(w_ij), std::abs(w_ji)});
            if (std::abs(w_ij - w_ji) > symmetry_tol * scale) {
                Rf_error("Precomputed graph reciprocal weight mismatch for edge (%zu,%zu): %.12g vs %.12g.",
                         i + 1, static_cast<size_t>(j) + 1, w_ij, w_ji);
            }
        }
    }

    vertex_cofaces.clear();
    vertex_cofaces.resize(n_points);
    for (index_t i = 0; i < n_points; ++i) {
        const auto& nbrs = adj_list[i];
        const auto& wts = weight_list[i];

        vertex_cofaces[i].reserve(nbrs.size() + 1);
        vertex_cofaces[i].push_back(neighbor_info_t{
                i,  // self-loop
                i,
                neighbor_sets[i].size(),
                0.0,
                0.0
            });

        for (size_t j_idx = 0; j_idx < nbrs.size(); ++j_idx) {
            const index_t j = nbrs[j_idx];
            const size_t isize = intersection_size_of_sets(neighbor_sets[i], neighbor_sets[j]);
            vertex_cofaces[i].push_back(neighbor_info_t{
                    j,
                    0,
                    isize,
                    wts[j_idx],
                    0.0
                });
        }
    }

    // Connectivity validation
    const int n_conn_comps = compute_connected_components();
    if (n_conn_comps > 1) {
        Rf_error("The provided precomputed graph has %d connected components. "
                 "rdgraph regression requires a connected graph.",
                 n_conn_comps);
    }

    // Build edge registry from undirected adjacency
    std::vector<std::pair<index_t, index_t>> edge_list;
    edge_list.reserve(n_points * 4);
    for (index_t i = 0; i < n_points; ++i) {
        for (index_t j : adj_list[i]) {
            if (j > i) {
                edge_list.push_back({i, j});
            }
        }
    }

    const index_t n_edges = static_cast<index_t>(edge_list.size());
    edge_registry.resize(n_edges);
    for (index_t e = 0; e < n_edges; ++e) {
        edge_registry[e] = {edge_list[e].first, edge_list[e].second};
    }
    extend_by_one_dim(n_edges);

    for (index_t e = 0; e < n_edges; ++e) {
        const auto [i, j] = edge_registry[e];

        for (size_t k_idx = 1; k_idx < vertex_cofaces[i].size(); ++k_idx) {
            if (vertex_cofaces[i][k_idx].vertex_index == j) {
                vertex_cofaces[i][k_idx].simplex_index = e;
                break;
            }
        }
        for (size_t k_idx = 1; k_idx < vertex_cofaces[j].size(); ++k_idx) {
            if (vertex_cofaces[j][k_idx].vertex_index == i) {
                vertex_cofaces[j][k_idx].simplex_index = e;
                break;
            }
        }
    }

    // Triangle policy decision
    const bool offdiag_edge_mass_enabled =
        (static_cast<size_t>(n_edges) <= EDGE_MASS_OFFDIAG_MAX_EDGES);
    bool build_triangles = false;
    const char* triangle_policy_label = "auto";
    const char* triangle_reason = "";

    switch (triangle_policy) {
    case triangle_policy_t::ALWAYS:
        triangle_policy_label = "always";
        build_triangles = true;
        triangle_reason = "forced by triangle.policy='always'";
        break;
    case triangle_policy_t::NEVER:
        triangle_policy_label = "never";
        build_triangles = false;
        triangle_reason = "disabled by triangle.policy='never'";
        break;
    case triangle_policy_t::AUTO:
    default:
        triangle_policy_label = "auto";
        build_triangles = (gamma_modulation > 0.0) && offdiag_edge_mass_enabled;
        if (build_triangles) {
            triangle_reason = "gamma_modulation > 0 and off-diagonal edge mass is enabled";
        } else if (gamma_modulation <= 0.0) {
            triangle_reason = "gamma_modulation <= 0 (response modulation disabled)";
        } else {
            triangle_reason = "edge count exceeds off-diagonal edge-mass threshold";
        }
        break;
    }

    if (vl_at_least(verbose_level, verbose_level_t::PROGRESS)) {
        progress_log(
            "  [triangles] policy=%s -> %s (%s; n_edges=%ld, offdiag_threshold=%ld, gamma=%.3g)",
            triangle_policy_label,
            build_triangles ? "build" : "skip",
            triangle_reason,
            static_cast<long>(n_edges),
            static_cast<long>(EDGE_MASS_OFFDIAG_MAX_EDGES),
            gamma_modulation
            );
    }

    // Build triangles and edge cofaces
    if (build_triangles) {
        edge_cofaces.resize(n_edges);
        for (index_t e = 0; e < n_edges; ++e) {
            const auto [i, j] = edge_registry[e];
            (void)j;
            edge_cofaces[e].push_back(neighbor_info_t{
                    static_cast<index_t>(i),
                    static_cast<index_t>(e),
                    0,
                    0.0,
                    0.0
                });
        }

        std::vector<std::array<index_t, 3>> triangle_list;
        triangle_list.reserve(n_points * 4);

        for (size_t i = 0; i < n_points; ++i) {
            for (size_t a = 1; a < vertex_cofaces[i].size(); ++a) {
                const index_t j = vertex_cofaces[i][a].vertex_index;
                for (size_t b = a + 1; b < vertex_cofaces[i].size(); ++b) {
                    const index_t s = vertex_cofaces[i][b].vertex_index;

                    index_t edge_js_idx = NO_EDGE;
                    for (size_t k_idx = 1; k_idx < vertex_cofaces[j].size(); ++k_idx) {
                        if (vertex_cofaces[j][k_idx].vertex_index == s) {
                            edge_js_idx = vertex_cofaces[j][k_idx].simplex_index;
                            break;
                        }
                    }
                    if (edge_js_idx == NO_EDGE) continue;

                    bool has_intersection = false;
                    for (index_t v : neighbor_sets[i]) {
                        if (neighbor_sets[j].find(v) != neighbor_sets[j].end() &&
                            neighbor_sets[s].find(v) != neighbor_sets[s].end()) {
                            has_intersection = true;
                            break;
                        }
                    }
                    if (!has_intersection) continue;

                    std::array<index_t, 3> tri = {
                        static_cast<index_t>(i), j, s
                    };
                    std::sort(tri.begin(), tri.end());
                    triangle_list.push_back(tri);
                }
            }
        }

        std::sort(triangle_list.begin(), triangle_list.end());
        triangle_list.erase(std::unique(triangle_list.begin(), triangle_list.end()),
                            triangle_list.end());

        const index_t n_triangles = static_cast<index_t>(triangle_list.size());
        for (index_t t = 0; t < n_triangles; ++t) {
            const auto& tri = triangle_list[t];
            const index_t v0 = tri[0];
            const index_t v1 = tri[1];
            const index_t v2 = tri[2];

            index_t e01 = NO_EDGE;
            for (size_t k = 1; k < vertex_cofaces[v0].size(); ++k) {
                if (vertex_cofaces[v0][k].vertex_index == v1) {
                    e01 = vertex_cofaces[v0][k].simplex_index;
                    break;
                }
            }
            index_t e02 = NO_EDGE;
            for (size_t k = 1; k < vertex_cofaces[v0].size(); ++k) {
                if (vertex_cofaces[v0][k].vertex_index == v2) {
                    e02 = vertex_cofaces[v0][k].simplex_index;
                    break;
                }
            }
            index_t e12 = NO_EDGE;
            for (size_t k = 1; k < vertex_cofaces[v1].size(); ++k) {
                if (vertex_cofaces[v1][k].vertex_index == v2) {
                    e12 = vertex_cofaces[v1][k].simplex_index;
                    break;
                }
            }

            size_t triple_isize = 0;
            for (index_t v : neighbor_sets[v0]) {
                if (neighbor_sets[v1].find(v) != neighbor_sets[v1].end() &&
                    neighbor_sets[v2].find(v) != neighbor_sets[v2].end()) {
                    triple_isize++;
                }
            }

            edge_cofaces[e01].push_back(neighbor_info_t{v2, t, triple_isize, 0.0, 0.0});
            edge_cofaces[e02].push_back(neighbor_info_t{v1, t, triple_isize, 0.0, 0.0});
            edge_cofaces[e12].push_back(neighbor_info_t{v0, t, triple_isize, 0.0, 0.0});
        }
    } else {
        edge_cofaces.clear();
    }

    if (vl_at_least(verbose_level, verbose_level_t::PROGRESS)) {
        elapsed_time(phase_time, "DONE", true, true);
    }

    // ================================================================
    // PHASE 2: DENSITY INITIALIZATION
    // ================================================================
    if (vl_at_least(verbose_level, verbose_level_t::PROGRESS)) {
        progress_log_inline("  Phase 2: reference measure + densities ... ");
        phase_time = std::chrono::steady_clock::now();
    }

    // Build pseudo-kNN vectors from local weighted neighborhoods.
    // Distances are sorted increasingly so d_k is the tail element.
    // Mirror kNN-branch semantics by including self (distance 0)
    // and retaining only the first neighborhood_k entries.
    std::vector<std::vector<index_t>> pseudo_knn_indices(n_points);
    std::vector<std::vector<double>> pseudo_knn_distances(n_points);
    for (size_t i = 0; i < n_points; ++i) {
        std::vector<std::pair<double, index_t>> by_dist;
        by_dist.reserve(adj_list[i].size() + 1);
        by_dist.emplace_back(0.0, static_cast<index_t>(i));
        for (size_t j = 0; j < adj_list[i].size(); ++j) {
            by_dist.emplace_back(weight_list[i][j], adj_list[i][j]);
        }
        std::sort(by_dist.begin(), by_dist.end(),
                  [](const std::pair<double, index_t>& lhs, const std::pair<double, index_t>& rhs) {
                      if (lhs.first < rhs.first) return true;
                      if (lhs.first > rhs.first) return false;
                      return lhs.second < rhs.second;
                  });

        const size_t k_effective = std::min(neighborhood_k_sz, by_dist.size());
        pseudo_knn_indices[i].reserve(k_effective);
        pseudo_knn_distances[i].reserve(k_effective);
        for (size_t t = 0; t < k_effective; ++t) {
            const double dist = by_dist[t].first;
            const index_t nbr = by_dist[t].second;
            pseudo_knn_indices[i].push_back(nbr);
            pseudo_knn_distances[i].push_back(dist);
        }
    }

    initialize_reference_measure(
        pseudo_knn_indices,
        pseudo_knn_distances,
        use_counting_measure,
        density_normalization,
        density_alpha,
        density_epsilon,
        clamp_dk,
        dk_clamp_median_factor,
        target_weight_ratio,
        pathological_ratio_threshold,
        verbose_level
    );

    compute_initial_densities(verbose_level);

    if (vl_at_least(verbose_level, verbose_level_t::PROGRESS)) {
        elapsed_time(phase_time, "DONE", true, true);
    }

    // ================================================================
    // PHASE 3: METRIC + OPERATORS
    // ================================================================
    if (vl_at_least(verbose_level, verbose_level_t::PROGRESS)) {
        progress_log_inline("  Phase 3: metric + operators assembly ... ");
        phase_time = std::chrono::steady_clock::now();
    }

    initialize_metric_from_density();
    build_boundary_operator_from_edges();
    if (!edge_cofaces.empty()) {
        build_boundary_operator_from_triangles();
    }
    assemble_operators();

    if (vl_at_least(verbose_level, verbose_level_t::PROGRESS)) {
        elapsed_time(phase_time, "DONE", true, true);
    }
}

/**
 * @brief Initialize Riemannian simplicial complex from k-NN structure
 *
 * This function constructs the complete geometric structure needed for
 * regression from k-nearest neighbor data. The initialization proceeds
 * through several phases that build the simplicial complex, compute
 * probability densities, construct the Riemannian metric, and assemble
 * the Hodge Laplacian operators.
 *
 * The construction follows van der Waerden's pedagogical approach: we begin
 * with the motivating observation that the k-NN graph captures local
 * geometric structure in the data. From this graph, we build a simplicial
 * complex whose topology and geometry encode both the intrinsic structure
 * of the feature space and the distribution of data points.
 *
 * PHASE 1: GEOMETRIC CONSTRUCTION
 *
 * We start by computing k-nearest neighborhoods for each data point. These
 * neighborhoods define a covering of the data space, with each neighborhood
 * serving as a local coordinate chart. The nerve theorem ensures that the
 * topology of this covering is captured by the simplicial complex whose
 * vertices correspond to data points and whose edges connect points with
 * overlapping neighborhoods.
 *
 * The edge construction uses an efficient O(nk) algorithm that examines
 * neighborhood intersections only between neighbors. For each pair of
 * vertices i,j that share neighborhood overlap, we create an edge weighted
 * by the intersection size and the minimum path length through common
 * neighbors. This weighting encodes both combinatorial (intersection size)
 * and metric (path length) information.
 *
 * Geometric edge pruning removes spurious connections that arise from
 * ambient space geometry but do not reflect intrinsic structure. We apply
 * ratio-based pruning that compares each edge length to the typical scale
 * of edges in the local neighborhood, removing edges whose length ratio
 * exceeds a threshold. This cleaning step produces a 1-skeleton that better
 * represents the intrinsic geometry.
 *
 * PHASE 2: DENSITY INITIALIZATION
 *
 * We compute initial probability densities from a reference measure on the
 * data points. The reference measure can be either uniform (counting measure)
 * or weighted by local density estimates. From this reference measure, we
 * derive vertex densities by summing measure over neighborhoods, and edge
 * densities by summing measure over neighborhood intersections.
 *
 * These densities provide the initial probability distribution that will
 * evolve during iterative refinement. The normalization ensures that vertex
 * densities sum to n and edge densities sum to n_edges, maintaining
 * interpretability as probability masses.
 *
 * PHASE 3: METRIC CONSTRUCTION
 *
 * The Riemannian metric is constructed from the density distributions. At
 * vertices (dimension 0), the metric is always diagonal with M_0 = diag(ρ_0).
 * At edges (dimension 1), the full mass matrix M_1 captures geometric
 * interactions through triple neighborhood intersections.
 *
 * For two edges sharing a vertex, their inner product is determined by the
 * density mass in the triple intersection of their endpoints' neighborhoods.
 * This construction ensures positive semidefiniteness while encoding the
 * geometric relationship between edges through their shared vertex structure.
 *
 * PHASE 4: LAPLACIAN ASSEMBLY
 *
 * The Hodge Laplacian L_0 is assembled from the boundary operator B_1 and
 * the metric tensors M_0, M_1 according to the formula L_0 = B_1 M_1^{-1} B_1^T M_0.
 * This operator encodes diffusion dynamics on the complex, with diffusion
 * rates determined by the Riemannian structure.
 *
 * The normalized Laplacian L_sym = M_0^{-1/2} L_0 M_0^{-1/2} is also computed
 * for spectral analysis and filtering operations.
 *
 * @param X Data matrix (n × d) where rows are observations and columns are features
 * @param k Number of nearest neighbors for graph construction
 * @param use_counting_measure If true, use uniform reference measure; if false,
 *                             use density-weighted measure based on local scales
 * @param density_normalization Power for density weighting (when not using counting measure).
 *                              Typical value: 1.0 gives 1/d_local weighting
 * @param max_ratio_threshold Maximum edge length ratio for geometric pruning.
 *                            Edges longer than this ratio times the local scale are removed.
 *                            Typical value: 2.0 to 3.0
 * @param threshold_percentile Percentile for computing local scale in geometric pruning.
 *                             Typical value: 0.5 (median)
 *
 * @post Simplicial complex S[0], S[1] populated with vertices and edges
 * @post Densities rho.rho[0], rho.rho[1] computed and normalized
 * @post Metric tensors g.M[0], g.M[1] constructed
 * @post Hodge Laplacian L.L[0] and mass-symmetrized Laplacian L.L0_mass_sym assembled
 * @post neighbor_sets member variable populated for subsequent use
 * @post reference_measure member variable populated
 *
 * @note This function performs substantial computation: O(nk^2) for edge
 *       construction, O(nk^3) for metric assembly. For large n or k, this
 *       initialization can take significant time.
 */
/**
 * @file riem_dcx_initialization.cpp
 * @brief Initialization routines for Riemannian simplicial complex regression
 */

#if 0
#include "riem_dcx.hpp"
#include "set_wgraph.hpp"
#include "kNN.h"
#include "debug_serialization.hpp"
#include <algorithm>
#include <limits>
#include <unordered_set>

/**
 * @brief Initialize geometric complex from k-nearest neighbor structure
 *
 * Constructs the 1-skeleton simplicial complex from feature matrix X by building
 * a k-nearest neighbor graph and initializing all geometric structures required
 * for the regression framework. This function establishes the topological foundation,
 * validates graph connectivity, initializes density measures, and prepares the
 * geometric operators for subsequent iterative refinement.
 *
 * The initialization proceeds through five phases:
 *
 * **Phase 1A: k-NN Construction** computes the k nearest neighbors for each point
 * using Euclidean distance in the feature space. The neighborhood information
 * defines local covering sets that will determine edge connectivity.
 *
 * **Phase 1B: Edge Construction** forms edges between vertices i and j when their
 * k-neighborhoods overlap (share at least one common neighbor). For each edge,
 * the function records the intersection size and computes the minimum path length
 * through common neighbors as an initial edge weight estimate.
 *
 * **Connectivity Validation** verifies that the resulting graph forms a single
 * connected component. Disconnected graphs cannot support the geometric diffusion
 * and spectral methods underlying the regression framework. If connectivity fails,
 * the function terminates with an error message suggesting parameter adjustments.
 *
 * **Phase 1C: Geometric Edge Pruning** (optional) removes edges with anomalously
 * large weights relative to the local weight distribution. This pruning step
 * eliminates spurious long-range connections that could distort the geometry,
 * typically arising from the k-NN construction in regions of varying density.
 *
 * **Phase 2: Density Initialization** constructs the reference measure μ and
 * computes initial vertex and edge densities by aggregating μ over appropriate
 * neighborhoods. The reference measure can be uniform (counting measure) or
 * density-weighted based on k-nearest neighbor distances.
 *
 * **Phase 3: Metric Construction** builds the initial Riemannian metric from
 * the initialized densities, constructing mass matrices M_0 (vertices) and M_1
 * (edges) that encode the geometric structure.
 *
 * **Phase 4: Build Boundary Operators** constructs the discrete differential
 * operators B_1 and B_2 encoding the simplicial complex topology.
 *
 * **Phase 5: Laplacian Assembly** assembles the Hodge Laplacian operators L_0
 * and L_1 from the metric and boundary operators, preparing them for use in
 * spectral filtering and density evolution.
 *
 * @param X Feature matrix (n × d) in sparse format, where n is the number of
 *          observations and d is the number of features. Each row represents
 *          one observation in feature space. Features should be appropriately
 *          scaled if they have different units or ranges.
 *
 * @param k Number of nearest neighbors to use for graph construction. Must
 *          satisfy 2 ≤ k < n. Larger k creates more edges and makes connectivity
 *          more likely but may include geometrically irrelevant long-range
 *          connections. Typical values: 10-50. The minimum k for connectivity
 *          depends on data geometry but is often around log(n) + 5.
 *
 * @param use_counting_measure If true, use uniform reference measure μ(x) = 1
 *                             for all vertices (counting measure). If false,
 *                             use distance-based reference measure
 *                             μ(x) = (ε + d_k(x))^(-α) where d_k(x) is the
 *                             distance to the k-th nearest neighbor. Counting
 *                             measure is appropriate for uniformly sampled data;
 *                             distance-based measure adapts to varying sampling
 *                             density.
 *
 * @param density_normalization Target sum for normalized density measures.
 *                              If 0 (default), densities are normalized to sum
 *                              to n, giving average vertex density = 1. If positive,
 *                              densities sum to the specified value. Most users
 *                              should use the default (0 for automatic n-normalization).
 *
 * @param max_ratio_threshold Maximum allowed ratio between edge weight and local
 *                           median edge weight during geometric pruning. Edges
 *                           with weights exceeding this ratio times the local
 *                           median are removed as anomalous long-range connections.
 *                           Set to Inf to disable pruning. Typical value: 3-10.
 *                           Only used when pruning is enabled.
 *
 * @param threshold_percentile Percentile (in [0,1]) of local edge weight
 *                            distribution used as reference for pruning threshold.
 *                            Value of 0.5 uses local median; 0.75 uses 75th
 *                            percentile. Higher percentiles create more aggressive
 *                            pruning. Only used when pruning is enabled.
 *
 * @param density_alpha Exponent α in [1,2] for distance-based reference measure
 *                     formula μ(x) = (ε + d_k(x))^(-α). Larger values (near 2)
 *                     create stronger density-adaptive weighting, concentrating
 *                     measure on densely sampled regions. Smaller values (near 1)
 *                     produce more uniform weighting. Only used when
 *                     use_counting_measure = false. Default: 1.5.
 *
 * @param density_epsilon Regularization parameter ε > 0 in reference measure
 *                       formula μ(x) = (ε + d_k(x))^(-α). Prevents numerical
 *                       issues when nearest neighbor distances are very small.
 *                       Should be small relative to typical k-NN distances but
 *                       large enough for stability. Only used when
 *                       use_counting_measure = false. Default: 1e-10.
 *
 * @throws Rf_error If the constructed k-NN graph is disconnected (has more than
 *                  one connected component). Error message suggests increasing k
 *                  or removing outlier samples.
 *
 * @pre The feature matrix X must have at least 2 rows (n ≥ 2)
 * @pre Parameter k must satisfy 1 ≤ k < n
 * @pre If use_counting_measure = false, density_alpha must be in [1, 2]
 * @pre If use_counting_measure = false, density_epsilon must be positive
 *
 * @post vertex_cofaces is populated with graph topology
 * @post reference_measure contains initialized vertex weights
 * @post Initial densities rho.rho[0] and rho.rho[1] are computed and normalized
 * @post Metric matrices g.M[0] and g.M[1] are constructed
 * @post Boundary operators B[1] and B[2] are assembled
 * @post Hodge Laplacians L.L[0] and L.L[1] are assembled
 * @post The object is ready for fit_rdgraph_regression() to proceed
 *       with iterative refinement
 *
 * @note To explore graph connectivity before calling this function, use
 *       dgraphs::create.iknn.graphs() or dgraphs::create.single.iknn.graph()
 *       from R to examine connectivity structure across different k values.
 *
 * @note Edge pruning (controlled by max_ratio_threshold) is optional but
 *       recommended for data with highly variable local density, where k-NN
 *       construction may create spurious long-range edges in sparse regions.
 *
 * @see fit_rdgraph_regression() for the main regression function
 * @see compute_connected_components() for connectivity validation
 * @see initialize_reference_measure() for density weighting details
 */
void riem_dcx_t::initialize_from_knn(
    const spmat_t& X,
    index_t k,
    bool use_counting_measure,
    double density_normalization,
    double max_ratio_threshold,
	double path_edge_ratio_percentile,
    double threshold_percentile,
	double density_alpha,
	double density_epsilon,
	bool clamp_dk,
	double dk_clamp_median_factor,
	double target_weight_ratio,
	double pathological_ratio_threshold,
	const std::string& knn_cache_path,
	int knn_cache_mode,
	triangle_policy_t triangle_policy,
	double gamma_modulation,
	verbose_level_t verbose_level
) {
    const size_t n_points = static_cast<size_t>(X.rows());
    auto phase_time = std::chrono::steady_clock::now();

    // Setup debug directory if enabled
    std::string debug_dir;
    #if DEBUG_INITIALIZE_FROM_KNN
        debug_dir = "/tmp/gflow_debug/initialize_from_knn/";
        // Note: In production code, you'd want to ensure this directory exists
        // For debugging, create it manually: mkdir -p /tmp/gflow_debug/initialize_from_knn/
    #endif

    // ================================================================
    // PHASE 0: INITIALIZE DIMENSION STRUCTURE
    // ================================================================

    pmax = -1;
    pmax = 0;
    g.M.resize(1);
    g.M_solver.resize(1);
    L.B.resize(1);
    L.L.resize(1);
    g.M[0] = spmat_t();
    g.M_solver[0].reset();
    L.B[0] = spmat_t();
    L.L[0] = spmat_t();

    // ================================================================
    // PHASE 1A: K-NN COMPUTATION AND NEIGHBORHOOD CONSTRUCTION
    // ================================================================
    if (vl_at_least(verbose_level, verbose_level_t::PROGRESS)) {
        progress_log_inline("  Phase 1: kNN neighborhoods ... ");
        phase_time = std::chrono::steady_clock::now();
    }

    bool knn_cache_hit = false;
    bool knn_cache_written = false;
    knn_result_t knn_result = compute_knn_from_eigen(
        X,
        static_cast<int>(k),
        knn_cache_path,
        knn_cache_mode,
        &knn_cache_hit,
        &knn_cache_written
    );

    if (knn_result.k < static_cast<int>(k)) {
        Rf_error("kNN results contain fewer neighbors than required (%d < %d)",
                 knn_result.k,
                 static_cast<int>(k));
    }

    if (vl_at_least(verbose_level, verbose_level_t::DEBUG) && knn_cache_mode != 0) {
        if (knn_cache_hit) {
            progress_log("  [kNN cache] hit: %s (cached k=%d)",
                         knn_cache_path.c_str(),
                         knn_result.k);
        } else if (knn_cache_written) {
            progress_log("  [kNN cache] wrote: %s (k=%d)",
                         knn_cache_path.c_str(),
                         knn_result.k);
        } else {
            progress_log("  [kNN cache] mode=%d, no cache IO performed", knn_cache_mode);
        }
    }

    std::vector<std::vector<index_t>> knn_indices(n_points, std::vector<index_t>(k));
    std::vector<std::vector<double>> knn_distances(n_points, std::vector<double>(k));
    const size_t knn_stride = static_cast<size_t>(knn_result.k);

    for (size_t i = 0; i < n_points; ++i) {
        for (size_t j = 0; j < k; ++j) {
            const size_t offset = i * knn_stride + j;
            knn_indices[i][j] = static_cast<index_t>(knn_result.indices[offset]);
            knn_distances[i][j] = knn_result.distances[offset];
        }
    }

    #if DEBUG_INITIALIZE_FROM_KNN
        debug_serialization::save_knn_result(
            debug_dir + "phase_1a_knn_result.bin",
            knn_indices, knn_distances, n_points, k
        );
    #endif

    neighbor_sets.resize(n_points);
    for (size_t i = 0; i < n_points; ++i) {
        for (index_t neighbor_idx : knn_indices[i]) {
            neighbor_sets[i].insert(neighbor_idx);
        }
    }

    #if DEBUG_INITIALIZE_FROM_KNN
        debug_serialization::save_neighbor_sets(
            debug_dir + "phase_1a_neighbor_sets.bin",
            neighbor_sets
        );
    #endif

    if (vl_at_least(verbose_level, verbose_level_t::PROGRESS)) {
        elapsed_time(phase_time, "DONE", true, true);
    }

    // ================================================================
    // PHASE 1B: EDGE CONSTRUCTION VIA NEIGHBORHOOD INTERSECTIONS
    // ================================================================
    if (vl_at_least(verbose_level, verbose_level_t::PROGRESS)) {
        progress_log_inline("  Phase 2: simplicial graph build/pruning ... ");
        phase_time = std::chrono::steady_clock::now();
    }

    const int iknn_num_threads = std::max(1, gflow_get_num_procs());
    // Temporary split dependency: rdgraph initialization still reuses gflow's
    // internal ikNN builder until the native backend is fully owned by dgraphs.
    iknn_graph_t iknn_graph = create_iknn_graph_from_knn_result(
        knn_result,
        static_cast<int>(k),
        iknn_num_threads > 1,
        iknn_num_threads
    );

	// Initialize vertex_cofaces with self-loops and shared-backend ikNN edges
	vertex_cofaces.resize(n_points);
	for (index_t i = 0; i < n_points; ++i) {
        const auto& neighbors = iknn_graph.graph[i];
		vertex_cofaces[i].reserve(neighbors.size() + 1);

		// Add self-loop at position [0]
		// Convention: vertex_cofaces[i][0] is always the vertex itself
		vertex_cofaces[i].push_back(neighbor_info_t{
				i,  // vertex_index = i (self-reference)
				i,  // simplex_index = i (vertex index)
				neighbor_sets[i].size(),
				0.0,
				0.0
			});

        for (const auto& nbr : neighbors) {
            vertex_cofaces[i].push_back(neighbor_info_t{
                    static_cast<index_t>(nbr.index),
                    0,
                    nbr.isize,
                    nbr.dist,
                    0.0
                });
        }
	}


    #if DEBUG_INITIALIZE_FROM_KNN
	debug_serialization::save_vertex_cofaces(
		debug_dir + "phase_1b_vertex_cofaces_pre_pruning.bin",
		vertex_cofaces
        );

	// Extract edge list for easier comparison
	std::vector<std::pair<index_t, index_t>> edges_pre_pruning;
	std::vector<double> weights_pre_pruning;
	for (size_t i = 0; i < n_points; ++i) {
		for (size_t k_idx = 1; k_idx < vertex_cofaces[i].size(); ++k_idx) {
			index_t j = vertex_cofaces[i][k_idx].vertex_index;
			if (j > i) {
				edges_pre_pruning.push_back({static_cast<index_t>(i), j});
				weights_pre_pruning.push_back(vertex_cofaces[i][k_idx].dist);
			}
		}
	}
	debug_serialization::save_edge_list(
		debug_dir + "phase_1b_edges_pre_pruning.bin",
		edges_pre_pruning, weights_pre_pruning, "PHASE_1B_PRE_PRUNING"
        );
    #endif


	// ================================================================
	// CONNECTIVITY VALIDATION
	// ================================================================

	/**
	 * Verify that the constructed k-nearest neighbor graph is connected.
	 *
	 * The graph must form a single connected component for the geometric
	 * diffusion and spectral methods to function properly. If the graph has
	 * multiple components, geometric information cannot propagate between
	 * components, and the Laplacian eigendecomposition will have multiple
	 * zero eigenvalues corresponding to the disconnected pieces.
	 *
	 * Common causes of disconnected graphs include:
	 * - k value too small relative to data geometry
	 * - Isolated outlier points far from main data cloud
	 * - Data naturally divided into separate clusters
	 *
	 * Solutions:
	 * - Increase k to create more connections
	 * - Remove outlier observations before fitting
	 * - Process each connected component separately if appropriate
	 */
	int n_conn_comps = compute_connected_components();

	if (n_conn_comps > 1) {
		Rf_error("The constructed k-NN graph (k=%d) has %d connected components.\n"
				 "The regression requires a fully connected graph.\n"
				 "Solutions:\n"
				 "  1. Increase k to create more edges\n"
				 "  2. Remove outlier samples that are isolated\n"
				 "  3. Use dgraphs::create.iknn.graphs() to explore connectivity across k values",
				 (int)k, n_conn_comps);
	}

    // ================================================================
    // PHASE 1C: GEOMETRIC EDGE PRUNING
    // ================================================================

    set_wgraph_t temp_graph(iknn_graph);


    #if DEBUG_INITIALIZE_FROM_KNN
	size_t n_edges_before_pruning = 0;
    for (size_t i = 0; i < n_points; ++i) {
		n_edges_before_pruning += vertex_cofaces[i].size() - 1;
	}
	n_edges_before_pruning /= 2;

	Rprintf("=== Before geometric pruning ===\n");
	Rprintf("Edges before: %zu\n", n_edges_before_pruning);
    #endif

    const bool do_geometric_prune = (max_ratio_threshold > 1.0);
    set_wgraph_t pruned_graph = temp_graph;
    if (do_geometric_prune) {
        pruned_graph = temp_graph.prune_edges_geometrically(
            max_ratio_threshold,
            path_edge_ratio_percentile
        );
    } else if (vl_at_least(verbose_level, verbose_level_t::PROGRESS)) {
        progress_log("[geometric_prune] skipped (ratio threshold <= 1.0)");
    }

    #if DEBUG_INITIALIZE_FROM_KNN
	size_t n_edges_after_pruning = 0;
	for (const auto& nbrs : pruned_graph.adjacency_list) {
		n_edges_after_pruning += nbrs.size();
	}
	n_edges_after_pruning /= 2;

	Rprintf("=== After geometric pruning ===\n");
	Rprintf("Edges after: %zu\n", n_edges_after_pruning);
	Rprintf("Edges removed: %zu\n", n_edges_before_pruning - n_edges_after_pruning);

	std::vector<std::pair<index_t, index_t>> edges_post_pruning;
	std::vector<double> weights_post_pruning;
	for (size_t i = 0; i < n_points; ++i) {
		for (const auto& edge_info : pruned_graph.adjacency_list[i]) {
			index_t j = edge_info.vertex;
			if (j > i) {
				edges_post_pruning.push_back({static_cast<index_t>(i), j});
				weights_post_pruning.push_back(edge_info.weight);
			}
		}
	}

	debug_serialization::save_edge_list(
		debug_dir + "phase_1c_edges_post_pruning.bin",
		edges_post_pruning, weights_post_pruning, "PHASE_1C_POST_PRUNING"
        );

	debug_serialization::save_pruning_params(
		debug_dir + "phase_1c_pruning_params.bin",
		max_ratio_threshold, threshold_percentile,
		n_edges_before_pruning, edges_post_pruning.size()
        );
    #endif

	// ---- Quantile-based edge length pruning
	if (threshold_percentile > 0.0) {
        #if DEBUG_INITIALIZE_FROM_KNN
		Rprintf("=== Applying quantile-based edge length pruning ===\n");
		Rprintf("threshold_percentile=%.3f (removing top %.1f%% of edges by length)\n",
				threshold_percentile, 100.0 * (1.0 - threshold_percentile));
        #endif

		pruned_graph = pruned_graph.prune_long_edges(threshold_percentile);

		// Recount edges after quantile pruning
		size_t n_edges_in_pruned_graph_sz = 0;
		for (const auto& nbrs : pruned_graph.adjacency_list) {
			n_edges_in_pruned_graph_sz += nbrs.size();
		}
		n_edges_in_pruned_graph_sz /= 2;

        #if DEBUG_INITIALIZE_FROM_KNN
		Rprintf("=== After quantile pruning ===\n");
		Rprintf("Edges after quantile pruning: %zu\n", n_edges_in_pruned_graph_sz);
		Rprintf("Additional edges removed: %zu\n",
				n_edges_after_geometric_sz - n_edges_in_pruned_graph_sz);
        #endif
	}



    // ================================================================
    // PHASE 1D: FILTER vertex_cofaces IN PLACE
    // ================================================================

    for (size_t i = 0; i < n_points; ++i) {
        const auto& pruned_neighbors = pruned_graph.adjacency_list[i];

        std::unordered_set<index_t> survivors;
        for (const auto& edge_info : pruned_neighbors) {
            survivors.insert(edge_info.vertex);
        }

        auto new_end = std::remove_if(
            vertex_cofaces[i].begin() + 1,
            vertex_cofaces[i].end(),
            [&survivors](const neighbor_info_t& info) {
                return !survivors.count(info.vertex_index);
            }
        );
        vertex_cofaces[i].erase(new_end, vertex_cofaces[i].end());
    }

    #if DEBUG_INITIALIZE_FROM_KNN
        debug_serialization::save_vertex_cofaces(
            debug_dir + "phase_1d_vertex_cofaces_post_pruning.bin",
            vertex_cofaces
        );
    #endif

    // ================================================================
    // PHASE 1E: BUILD edge_registry AND ASSIGN FINAL INDICES
    // ================================================================

    std::vector<std::pair<index_t, index_t>> edge_list;
    edge_list.reserve(pruned_graph.adjacency_list.size() * k / 2);

    for (size_t i = 0; i < n_points; ++i) {
        for (const auto& edge_info : pruned_graph.adjacency_list[i]) {
            index_t j = edge_info.vertex;
            if (j > i) {
                edge_list.push_back({static_cast<index_t>(i), j});
            }
        }
    }

    const index_t n_edges = edge_list.size();
    edge_registry.resize(n_edges);
    for (index_t e = 0; e < n_edges; ++e) {
        edge_registry[e] = {edge_list[e].first, edge_list[e].second};
    }

    extend_by_one_dim(n_edges);

    for (index_t e = 0; e < n_edges; ++e) {
        const auto [i, j] = edge_registry[e];

        for (size_t k_idx = 1; k_idx < vertex_cofaces[i].size(); ++k_idx) {
            if (vertex_cofaces[i][k_idx].vertex_index == j) {
                vertex_cofaces[i][k_idx].simplex_index = e;
                break;
            }
        }

        for (size_t k_idx = 1; k_idx < vertex_cofaces[j].size(); ++k_idx) {
            if (vertex_cofaces[j][k_idx].vertex_index == i) {
                vertex_cofaces[j][k_idx].simplex_index = e;
                break;
            }
        }
    }

    #if DEBUG_INITIALIZE_FROM_KNN
	debug_serialization::save_vertex_cofaces(
		debug_dir + "phase_1e_final_vertex_cofaces.bin",
		vertex_cofaces
        );

	std::vector<double> final_weights;
	for (const auto& [i, j] : edge_list) {
		for (size_t k_idx = 1; k_idx < vertex_cofaces[i].size(); ++k_idx) {
			if (vertex_cofaces[i][k_idx].vertex_index == j) {
				final_weights.push_back(vertex_cofaces[i][k_idx].dist);
				break;
			}
		}
	}
	debug_serialization::save_edge_list(
		debug_dir + "phase_1e_final_edge_registry.bin",
		edge_list, final_weights, "PHASE_1E_FINAL"
        );

	// Compute and save connectivity
	std::vector<int> component_ids;
	size_t n_components = debug_serialization::compute_connected_components_from_vertex_cofaces(
		vertex_cofaces, component_ids
        );
	debug_serialization::save_connectivity(
		debug_dir + "phase_1e_connectivity.bin",
		component_ids, n_components
        );
    #endif

	// ================================================================
	// PHASE 1F: TRIANGLE POLICY DECISION
	// ================================================================
	const bool offdiag_edge_mass_enabled =
		(static_cast<size_t>(n_edges) <= EDGE_MASS_OFFDIAG_MAX_EDGES);
	bool build_triangles = false;
	const char* triangle_policy_label = "auto";
	const char* triangle_reason = "";

	switch (triangle_policy) {
	case triangle_policy_t::ALWAYS:
		triangle_policy_label = "always";
		build_triangles = true;
		triangle_reason = "forced by triangle.policy='always'";
		break;
	case triangle_policy_t::NEVER:
		triangle_policy_label = "never";
		build_triangles = false;
		triangle_reason = "disabled by triangle.policy='never'";
		break;
	case triangle_policy_t::AUTO:
	default:
		triangle_policy_label = "auto";
		build_triangles = (gamma_modulation > 0.0) && offdiag_edge_mass_enabled;
		if (build_triangles) {
			triangle_reason = "gamma_modulation > 0 and off-diagonal edge mass is enabled";
		} else if (gamma_modulation <= 0.0) {
			triangle_reason = "gamma_modulation <= 0 (response modulation disabled)";
		} else {
			triangle_reason = "edge count exceeds off-diagonal edge-mass threshold";
		}
		break;
	}

	if (vl_at_least(verbose_level, verbose_level_t::PROGRESS)) {
		progress_log(
			"  [triangles] policy=%s -> %s (%s; n_edges=%ld, offdiag_threshold=%ld, gamma=%.3g)",
			triangle_policy_label,
			build_triangles ? "build" : "skip",
			triangle_reason,
			static_cast<long>(n_edges),
			static_cast<long>(EDGE_MASS_OFFDIAG_MAX_EDGES),
			gamma_modulation
			);
	}

	// ================================================================
	// PHASE 1G: BUILD TRIANGLES AND POPULATE edge_cofaces (optional)
	// ================================================================
	if (build_triangles) {
		// Keep self-loop convention at [0] for downstream triangle iterators.
		edge_cofaces.resize(n_edges);

		for (size_t e = 0; e < n_edges; ++e) {
			edge_cofaces[e].reserve(1);

			const auto [i, j] = edge_registry[e];
			(void)j;

			edge_cofaces[e].push_back(neighbor_info_t{
					static_cast<index_t>(i),
					static_cast<index_t>(e),
					0,
					0.0,
					0.0
				});
		}

		std::vector<std::array<index_t, 3>> triangle_list;
		triangle_list.reserve(n_points * k);

		for (size_t i = 0; i < n_points; ++i) {
			for (size_t a = 1; a < vertex_cofaces[i].size(); ++a) {
				index_t j = vertex_cofaces[i][a].vertex_index;

				for (size_t b = a + 1; b < vertex_cofaces[i].size(); ++b) {
					index_t s = vertex_cofaces[i][b].vertex_index;

					index_t edge_js_idx = NO_EDGE;
					for (size_t k_idx = 1; k_idx < vertex_cofaces[j].size(); ++k_idx) {
						if (vertex_cofaces[j][k_idx].vertex_index == s) {
							edge_js_idx = vertex_cofaces[j][k_idx].simplex_index;
							break;
						}
					}

					if (edge_js_idx == NO_EDGE) continue;

					bool has_intersection = false;
					for (index_t v : neighbor_sets[i]) {
						if (neighbor_sets[j].find(v) != neighbor_sets[j].end() &&
							neighbor_sets[s].find(v) != neighbor_sets[s].end()) {
							has_intersection = true;
							break;
						}
					}

					if (!has_intersection) continue;

					std::array<index_t, 3> tri = {i, j, s};
					std::sort(tri.begin(), tri.end());
					triangle_list.push_back(tri);
				}
			}
		}

		std::sort(triangle_list.begin(), triangle_list.end());
		triangle_list.erase(std::unique(triangle_list.begin(), triangle_list.end()),
							triangle_list.end());

		const index_t n_triangles = triangle_list.size();

		for (index_t t = 0; t < n_triangles; ++t) {
			const auto& tri = triangle_list[t];
			const index_t v0 = tri[0];
			const index_t v1 = tri[1];
			const index_t v2 = tri[2];

			index_t e01 = NO_EDGE;
			for (size_t k = 1; k < vertex_cofaces[v0].size(); ++k) {
				if (vertex_cofaces[v0][k].vertex_index == v1) {
					e01 = vertex_cofaces[v0][k].simplex_index;
					break;
				}
			}

			index_t e02 = NO_EDGE;
			for (size_t k = 1; k < vertex_cofaces[v0].size(); ++k) {
				if (vertex_cofaces[v0][k].vertex_index == v2) {
					e02 = vertex_cofaces[v0][k].simplex_index;
					break;
				}
			}

			index_t e12 = NO_EDGE;
			for (size_t k = 1; k < vertex_cofaces[v1].size(); ++k) {
				if (vertex_cofaces[v1][k].vertex_index == v2) {
					e12 = vertex_cofaces[v1][k].simplex_index;
					break;
				}
			}

			size_t triple_isize = 0;
			for (index_t v : neighbor_sets[v0]) {
				if (neighbor_sets[v1].find(v) != neighbor_sets[v1].end() &&
					neighbor_sets[v2].find(v) != neighbor_sets[v2].end()) {
					triple_isize++;
				}
			}

			edge_cofaces[e01].push_back(neighbor_info_t{
					v2,
					t,
					triple_isize,
					0.0,
					0.0
				});

			edge_cofaces[e02].push_back(neighbor_info_t{
					v1,
					t,
					triple_isize,
					0.0,
					0.0
				});

			edge_cofaces[e12].push_back(neighbor_info_t{
					v0,
					t,
					triple_isize,
					0.0,
					0.0
				});
		}
	} else {
		edge_cofaces.clear();
	}

    if (vl_at_least(verbose_level, verbose_level_t::PROGRESS)) {
        elapsed_time(phase_time, "DONE", true, true);
    }

    // ================================================================
    // PHASE 2: DENSITY INITIALIZATION
    // ================================================================
    if (vl_at_least(verbose_level, verbose_level_t::PROGRESS)) {
        progress_log_inline("  Phase 3: reference measure + densities ... ");
        phase_time = std::chrono::steady_clock::now();
    }

    // Initialize reference measure
    initialize_reference_measure(
        knn_indices,
        knn_distances,
        use_counting_measure,
        density_normalization,
		density_alpha,
		density_epsilon,
		clamp_dk,
		dk_clamp_median_factor,
		target_weight_ratio,
		pathological_ratio_threshold,
		verbose_level
    );

    // Compute initial densities from reference measure
    compute_initial_densities(verbose_level);

    if (vl_at_least(verbose_level, verbose_level_t::PROGRESS)) {
        elapsed_time(phase_time, "DONE", true, true);
    }

    // ================================================================
    // PHASE 3: METRIC CONSTRUCTION
    // ================================================================
    if (vl_at_least(verbose_level, verbose_level_t::PROGRESS)) {
        progress_log_inline("  Phase 4: metric + operators assembly ... ");
        phase_time = std::chrono::steady_clock::now();
    }

    initialize_metric_from_density();

	// ================================================================
    // PHASE 4: BUILD BOUNDARY OPERATORS
    // ================================================================

	build_boundary_operator_from_edges();      // B[1]: edges → vertices

	if (edge_cofaces.size() > 0) {
		build_boundary_operator_from_triangles();  // B[2]: triangles → edges
	}

    // ================================================================
    // PHASE 5: LAPLACIAN ASSEMBLY
    // ================================================================

    assemble_operators();

    if (vl_at_least(verbose_level, verbose_level_t::PROGRESS)) {
        elapsed_time(phase_time, "DONE", true, true);
    }
}

/**
 * @brief Initialize the reference measure (vertex weights) used by the Riemannian DCX solver.
 *
 * @details
 * This method defines an initial per-vertex reference measure (a positive weight vector)
 * used by downstream diffusion/regularization steps. Two modes are supported:
 *
 * 1) Counting measure (uniform):
 *    - If use_counting_measure == true, sets all weights to 1 and normalizes to the target sum.
 *
 * 2) Distance-based measure (robust density surrogate):
 *    - Let d_k(i) be the k-th neighbor distance for vertex i (the last entry of knn_distances[i]).
 *    - Raw weights are computed as:
 *        w(i) = (density_epsilon + d_k_used(i))^{-density_alpha},
 *      with input validation to prevent non-finite/invalid values.
 *
 *    - Optional d_k clamping (new):
 *      If clamp_dk == true, d_k_used(i) is clamped to:
 *        [median(d_k) / dk_clamp_median_factor, median(d_k) * dk_clamp_median_factor]
 *      before converting to weights.
 *      If clamp_dk == false, NO clamping is performed; the same bounds are computed only
 *      for diagnostics, and dk_used == dk_raw.
 *
 *    - Weight dynamic range control:
 *      After raw weights are computed, their ratio is capped to a fixed policy target
 *      (currently 1000:1). If the ratio is pathological (non-finite or extremely large),
 *      an additional log-compression step is applied first, followed by an exact ratio cap.
 *
 * Diagnostics:
 *  - dk_raw: raw d_k values
 *  - dk_used: d_k values used to compute weights (equals dk_raw if clamp_dk == false)
 *  - dk_lower, dk_upper: the median-factor bounds used (for reporting/diagnostics)
 *  - dk_used_low/high: indices that fall outside bounds; if clamp_dk==true these were clamped,
 *    otherwise these are "would-have-been clamped" indices.
 *
 * Normalization:
 *  - The weights are scaled so that sum(weights) == target, where
 *      target = density_normalization if > 0, else target = n.
 *
 * @param knn_neighbors List of kNN neighbor indices for each vertex (used for input validation only).
 * @param knn_distances List of kNN neighbor distances for each vertex; must be same shape as knn_neighbors.
 * @param use_counting_measure If true, use uniform weights (counting measure).
 * @param density_normalization Target sum for the weights. If <= 0, uses n (number of vertices).
 * @param density_alpha Exponent alpha in (epsilon + d_k)^{-alpha}. Must be finite and >= 0.
 * @param density_epsilon Regularization epsilon added to d_k. Must be finite and >= 0.
 * @param clamp_dk If true, clamp d_k prior to computing weights using median-factor bounds.
 * @param dk_clamp_median_factor Median-factor for clamping bounds; must be finite and > 1.
 * @param verbose_level Verbosity level for diagnostics.
 *
 * @throws Rf_error on invalid inputs or non-finite intermediate values.
 */
#endif

void riem_dcx_t::initialize_reference_measure(
    const std::vector<std::vector<index_t>>& knn_neighbors,
    const std::vector<std::vector<double>>& knn_distances,
    bool use_counting_measure,
    double density_normalization,
    double density_alpha,
    double density_epsilon,
    bool clamp_dk,
    double dk_clamp_median_factor,
	double target_weight_ratio,
	double pathological_ratio_threshold,
    verbose_level_t verbose_level
) {
    const size_t n = knn_neighbors.size();
    if (n == 0) {
        Rf_error("reference_measure: empty kNN inputs (n=0).");
    }
    if (knn_distances.size() != n) {
        Rf_error("reference_measure: knn_distances size mismatch (got %zu, expected %zu).",
                 (size_t)knn_distances.size(), n);
    }
    if (!std::isfinite(density_alpha) || density_alpha < 0.0) {
        Rf_error("reference_measure: density_alpha must be finite and >= 0.");
    }
    if (!std::isfinite(density_epsilon) || density_epsilon < 0.0) {
        Rf_error("reference_measure: density_epsilon must be finite and >= 0.");
    }
    if (!std::isfinite(density_normalization)) {
        Rf_error("reference_measure: density_normalization must be finite.");
    }
    if (!std::isfinite(dk_clamp_median_factor) || dk_clamp_median_factor <= 1.0) {
        Rf_error("reference_measure: dk_clamp_median_factor must be finite and > 1.");
    }

    std::vector<double> vertex_weights(n);

    // Reset dk diagnostics each time
    dk_raw.clear();
    dk_used.clear();
    dk_used_low.clear();
    dk_used_high.clear();
    dk_lower = NA_REAL;
    dk_upper = NA_REAL;

    if (use_counting_measure) {
        // Counting measure: uniform weights
        std::fill(vertex_weights.begin(), vertex_weights.end(), 1.0);

        // For counting measure we leave dk diagnostics empty/NA.
    } else {
        // ================================================================
        // DISTANCE-BASED MEASURE (OPTIONAL d_k CLAMPING; WEIGHT RATIO CAPPED)
        // ================================================================

        // Policy constants (kept internal for now)
        // const double target_weight_ratio = 1000.0;        // enforce w_max / w_min <= this
        // const double pathological_ratio_threshold = 1e12; // trigger log compression first

        // Step 1: Collect all k-th neighbor distances (raw d_k)
        std::vector<double> all_dk(n);
        for (size_t i = 0; i < n; ++i) {
            const auto& nbrs = knn_neighbors[i];
            const auto& dists = knn_distances[i];

            if (nbrs.size() != dists.size()) {
                Rf_error("reference_measure: knn_neighbors/knn_distances row %zu size mismatch.", i + 1);
            }
            if (dists.empty()) {
                Rf_error("reference_measure: knn_distances row %zu is empty (need at least 1 neighbor).", i + 1);
            }

            const double d_k = dists.back();
            if (!std::isfinite(d_k) || d_k < 0.0) {
                Rf_error("reference_measure: non-finite or negative d_k at vertex %zu.", i + 1);
            }
            all_dk[i] = d_k;
        }

        // Store dk diagnostics
        dk_raw = all_dk;
        dk_used = all_dk; // may be modified if clamping enabled

        // Step 2: Compute median-factor bounds
        std::vector<double> sorted_dk = all_dk;
        std::sort(sorted_dk.begin(), sorted_dk.end());
        const double d_median = sorted_dk[n / 2];

        // Note: if d_median == 0, bounds become [0,0]; that's OK as long as (epsilon + d) > 0
        const double d_min_allowed = d_median / dk_clamp_median_factor;
        const double d_max_allowed = d_median * dk_clamp_median_factor;

        dk_lower = d_min_allowed;
        dk_upper = d_max_allowed;

        // Step 3: Determine out-of-bounds indices; clamp if requested
        for (size_t i = 0; i < n; ++i) {
            double d = all_dk[i];

            const bool is_low  = (d < d_min_allowed);
            const bool is_high = (d > d_max_allowed);

            if (is_low) {
                dk_used_low.push_back((index_t)i);
                if (clamp_dk) d = d_min_allowed;
            } else if (is_high) {
                dk_used_high.push_back((index_t)i);
                if (clamp_dk) d = d_max_allowed;
            }

            if (clamp_dk) {
                dk_used[i] = d;
            }
        }

        if (vl_at_least(verbose_level, verbose_level_t::DEBUG)) {
            Rprintf("\n\tReference measure: d_k range [%.3e, %.3e], median=%.3e\n",
                    sorted_dk.front(), sorted_dk.back(), d_median);

            if (!dk_used_low.empty() || !dk_used_high.empty()) {
                if (clamp_dk) {
                    Rprintf("\tClamping enabled (factor=%.3g): bounds [%.3e, %.3e]; clamped: low=%d, high=%d\n",
                            dk_clamp_median_factor, d_min_allowed, d_max_allowed,
                            (int)dk_used_low.size(), (int)dk_used_high.size());
                } else {
                    Rprintf("\tDiagnostic (no-clamp) bounds [%.3e, %.3e]; outside: low=%d, high=%d\n",
                            d_min_allowed, d_max_allowed,
                            (int)dk_used_low.size(), (int)dk_used_high.size());
                }
            } else {
                if (clamp_dk) {
                    Rprintf("\tClamping enabled (factor=%.3g): no vertices outside bounds.\n",
                            dk_clamp_median_factor);
                }
            }
        }

        // Step 4: Raw weights from d_k (clamped or raw depending on clamp_dk)
        // w(i) = (epsilon + d_k_used(i))^{-alpha}
        for (size_t i = 0; i < n; ++i) {
            const double base = density_epsilon + dk_used[i];

            if (!(base > 0.0)) {
                // If alpha == 0, weight is 1 regardless (but still prefer a strict base check)
                if (density_alpha == 0.0) {
                    vertex_weights[i] = 1.0;
                } else {
                    Rf_error("reference_measure: density_epsilon + d_k must be > 0 (vertex %zu).", i + 1);
                }
            } else {
                vertex_weights[i] = std::pow(base, -density_alpha);
            }

            if (!std::isfinite(vertex_weights[i]) || vertex_weights[i] < 0.0) {
                Rf_error("reference_measure: non-finite or negative weight produced at vertex %zu.", i + 1);
            }
        }

        // Step 5: Cap weight ratio to target_weight_ratio (if needed)
        double w_min = *std::min_element(vertex_weights.begin(), vertex_weights.end());
        double w_max = *std::max_element(vertex_weights.begin(), vertex_weights.end());

        if (!(w_min > 0.0) || !std::isfinite(w_min) || !std::isfinite(w_max)) {
            Rf_error("reference_measure: invalid weight range before compression (min=%.3e, max=%.3e).",
                     w_min, w_max);
        }

        double ratio = w_max / w_min;

        if (!std::isfinite(ratio) || ratio > target_weight_ratio) {

            if (vl_at_least(verbose_level, verbose_level_t::DEBUG)) {
                Rprintf("\tInitial weight ratio = %.2e (target <= %.2e)\n",
                        ratio, target_weight_ratio);
            }

            // Robust mode for pathological ratios
            if (!std::isfinite(ratio) || ratio > pathological_ratio_threshold) {

                if (vl_at_least(verbose_level, verbose_level_t::DEBUG)) {
                    Rprintf("\tPathological ratio detected (%.2e). Applying log compression...\n", ratio);
                }

                validate_and_log_compress_weights_in_place(
                    vertex_weights,
                    /*average_even_median=*/false,
                    /*clamp_nonpositive_to_zero=*/false
                );

                // After log compression, still enforce exact cap via power compression if needed
                power_compress_weights_to_ratio_in_place(vertex_weights, target_weight_ratio);

            } else {
                // Normal case: enforce exact cap via power compression
                power_compress_weights_to_ratio_in_place(vertex_weights, target_weight_ratio);
            }

            w_min = *std::min_element(vertex_weights.begin(), vertex_weights.end());
            w_max = *std::max_element(vertex_weights.begin(), vertex_weights.end());

            if (!(w_min > 0.0) || !std::isfinite(w_min) || !std::isfinite(w_max)) {
                Rf_error("reference_measure: invalid weight range after compression (min=%.3e, max=%.3e).",
                         w_min, w_max);
            }

            ratio = w_max / w_min;

            if (vl_at_least(verbose_level, verbose_level_t::DEBUG)) {
                Rprintf("\tAfter compression: weight ratio = %.2e\n", ratio);
            }
        }
    }

    // Normalize to target sum
    double total = std::accumulate(vertex_weights.begin(), vertex_weights.end(), 0.0);
    const double target = (density_normalization > 0.0) ?
        density_normalization : static_cast<double>(n);

    if (!std::isfinite(total) || total <= 0.0) {
        Rf_error("reference_measure: non-positive or non-finite total weight (total=%.3e).", total);
    }
    if (!std::isfinite(target) || target <= 0.0) {
        Rf_error("reference_measure: invalid target normalization (target=%.3e).", target);
    }

    {
        const double scale = target / total;
        for (double& w : vertex_weights) {
            w *= scale;
        }
    }

    // Final safety check
    const double final_min = *std::min_element(vertex_weights.begin(), vertex_weights.end());
    const double final_max = *std::max_element(vertex_weights.begin(), vertex_weights.end());

    if (!(final_min > 0.0) || !std::isfinite(final_min) || !std::isfinite(final_max)) {
        Rf_error("reference_measure: invalid final weight range (min=%.3e, max=%.3e).",
                 final_min, final_max);
    }

    const double final_ratio = final_max / final_min;

    if (vl_at_least(verbose_level, verbose_level_t::DEBUG)) {
        Rprintf("\tFinal reference measure: range [%.3e, %.3e], ratio=%.2e\n",
                final_min, final_max, final_ratio);
    }

    if (final_ratio > 1e6) {
        Rf_warning("Reference measure has extreme ratio %.2e. "
                   "Consider use_counting_measure=TRUE for more uniform initialization.",
                   final_ratio);
    }

    this->reference_measure = std::move(vertex_weights);
}
