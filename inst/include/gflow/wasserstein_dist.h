#ifndef WASSERSTEIN_DIST_H_
#define WASSERSTEIN_DIST_H_

// C interface declarations
#ifdef __cplusplus
extern "C" {
#endif

    void C_wasserstein_distance_1D(const double *p,
								   const double *q,
								   const int    *rn,
								   double *d);

#ifdef __cplusplus
}
#endif
#endif // WASSERSTEIN_DIST_H_
