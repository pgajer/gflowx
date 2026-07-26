#ifndef GFLOW_LM_H_
#define GFLOW_LM_H_

#ifdef __cplusplus
extern "C" {
#endif

void C_flm(const double *x,
           const double *y,
           const int    *rnr,
           const int    *rnc,
           double *beta);

#ifdef __cplusplus
}
#endif

#endif // GFLOW_LM_H_
