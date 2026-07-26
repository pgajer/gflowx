#ifndef IKNN_GRAPHS_R_H_
#define IKNN_GRAPHS_R_H_

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_verify_pruning(SEXP s_X,
						  SEXP s_k,
						  SEXP s_max_alt_path_length);

#ifdef __cplusplus
}
#endif
#endif // IKNN_GRAPHS_R_H_
