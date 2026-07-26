#ifndef ULOGIT_R_H
#define ULOGIT_R_H

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

	SEXP S_ulogit(SEXP x_sexp,
				  SEXP y_sexp,
				  SEXP w_sexp,
				  SEXP max_iterations_sexp,
				  SEXP ridge_lambda_sexp,
				  SEXP max_beta_sexp,
				  SEXP tolerance_sexp,
				  SEXP verbose_sexp);

	SEXP S_eigen_ulogit(SEXP x_sexp,
						SEXP y_sexp,
						SEXP w_sexp,
						SEXP fit_quadratic_sexp,
						SEXP with_errors_sexp,
						SEXP max_iterations_sexp,
						SEXP ridge_lambda_sexp,
						SEXP tolerance_sexp);
#ifdef __cplusplus
}
#endif

#endif // ULOGIT_R_H
