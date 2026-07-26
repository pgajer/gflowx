#ifndef RIEM_DCX_R_H_
#define RIEM_DCX_R_H_

#include <Rinternals.h>

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_fit_rdgraph_regression(
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
			);

	SEXP S_compute_hop_extremp_radii_batch(
		SEXP s_adj_list,
		SEXP s_edge_densities,
		SEXP s_vertex_densities,
		SEXP s_candidates,
		SEXP s_y,
		SEXP s_p_threshold,
		SEXP s_detect_maxima,
		SEXP s_max_hop
		);

		SEXP S_compute_basins_of_attraction(
			SEXP s_adj_list,
			SEXP s_weight_list,
			SEXP s_y,
			SEXP s_edge_length_thld,
			SEXP s_with_trajectories
			//SEXP s_k_paths
			);

		SEXP S_compute_basins_of_attraction_rtcb(
			SEXP s_adj_list,
			SEXP s_weight_list,
			SEXP s_y,
			SEXP s_edge_length_thld,
			SEXP s_with_trajectories,
			SEXP s_params
			);

#ifdef __cplusplus
}
#endif

#endif // RIEM_DCX_R_H_
