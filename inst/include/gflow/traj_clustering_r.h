// traj_clustering_r.h
//
// R interface declarations for trajectory clustering utilities.

#pragma once

#include <R.h>
#include <Rinternals.h>

SEXP S_cluster_cell_trajectories(SEXP s_traj_list,
								 SEXP s_similarity_type,
								 SEXP s_overlap_mode,
								 SEXP s_min_shared,
								 SEXP s_min_frac,
								 SEXP s_min_subpath,
								 SEXP s_exclude_endpoints,
								 SEXP s_idf_smooth,
								 SEXP s_k,
								 SEXP s_symmetrize,
								 SEXP s_knn_select,
								 SEXP s_n_threads);
