#ifndef SAMPLING_H_
#define SAMPLING_H_

#ifdef __cplusplus
extern "C" {
#endif

    void C_runif_simplex( const int *rK, double *lambda);
    void C_rsimplex(const double *w, const int *rn, double *lambda);

#ifdef __cplusplus
}
#endif

#endif // SAMPLING_H_
