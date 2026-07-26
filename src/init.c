#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#include "riem_dcx_r.h"
#include "riem_dcx_posterior_r.h"
#include "ulogit_r.h"
#include "mean_shift_smoother_r.h"

extern SEXP _gflowx_rcpp_adaptive_mean_shift_gfa(
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP
);
extern SEXP _gflowx_rcpp_knn_adaptive_mean_shift_gfa(
    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP
);

static const R_CallMethodDef call_methods[] = {
    {"S_fit_rdgraph_regression", (DL_FUNC) &S_fit_rdgraph_regression, 42},
    {"S_compute_hop_extremp_radii_batch", (DL_FUNC) &S_compute_hop_extremp_radii_batch, 8},
    {"S_compute_posterior_summary", (DL_FUNC) &S_compute_posterior_summary, 10},
    {"S_ulogit", (DL_FUNC) &S_ulogit, 8},
    {"S_eigen_ulogit", (DL_FUNC) &S_eigen_ulogit, 8},
    {"S_mean_shift_data_smoother", (DL_FUNC) &S_mean_shift_data_smoother, 11},
    {"S_mean_shift_data_smoother_with_grad_field_averaging",
     (DL_FUNC) &S_mean_shift_data_smoother_with_grad_field_averaging, 8},
    {"S_mean_shift_data_smoother_adaptive",
     (DL_FUNC) &S_mean_shift_data_smoother_adaptive, 8},
    {"_gflowx_rcpp_adaptive_mean_shift_gfa",
     (DL_FUNC) &_gflowx_rcpp_adaptive_mean_shift_gfa, 11},
    {"_gflowx_rcpp_knn_adaptive_mean_shift_gfa",
     (DL_FUNC) &_gflowx_rcpp_knn_adaptive_mean_shift_gfa, 8},
    {NULL, NULL, 0}
};

void attribute_visible R_init_gflowx(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, FALSE);
}
