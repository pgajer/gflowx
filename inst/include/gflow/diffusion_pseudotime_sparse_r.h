#ifndef DIFFUSION_PSEUDOTIME_SPARSE_R_H
#define DIFFUSION_PSEUDOTIME_SPARSE_R_H

#include <R.h>
#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

SEXP S_build_sparse_transition(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_weight_mode,
    SEXP s_weight_param,
    SEXP s_lazy
);

SEXP S_compute_diffusion_pseudotime_sparse(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_root_vertices,
    SEXP s_root_weights,
    SEXP s_t_steps,
    SEXP s_n_probes,
    SEXP s_seed,
    SEXP s_weight_mode,
    SEXP s_weight_param,
    SEXP s_lazy,
    SEXP s_normalize,
    SEXP s_return_transition
);

SEXP S_compute_potential_pseudotime_sparse(
    SEXP s_adj_list,
    SEXP s_weight_list,
    SEXP s_root_vertices,
    SEXP s_root_weights,
    SEXP s_t_steps,
    SEXP s_potential_eps,
    SEXP s_landmark_vertices,
    SEXP s_weight_mode,
    SEXP s_weight_param,
    SEXP s_lazy,
    SEXP s_normalize,
    SEXP s_return_transition
);

#ifdef __cplusplus
}
#endif

#endif // DIFFUSION_PSEUDOTIME_SPARSE_R_H
