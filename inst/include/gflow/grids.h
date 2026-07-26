#ifndef GFLOW_GRIDS_H_
#define GFLOW_GRIDS_H_

#ifdef __cplusplus
extern "C" {
#endif

void C_create_ED_grid_2D(const double *rdx,
                         const double *rx1L,
                         const int    *rn1,
                         const double *rx2L,
                         const int    *rn2,
                         double *grid);

void C_create_ED_grid_3D(const double *rdx,
                         const double *rx1L,
                         const int    *rn1,
                         const double *rx2L,
                         const int    *rn2,
                         const double *rx3L,
                         const int    *rn3,
                         double *grid);

void C_create_ED_grid_xD(const double *rw,
                         const double *L,
                         const int    *rdim,
                         const int    *size_d,
                         const int    *rTotElts,
                         double *grid);

void C_create_ENPs_grid_2D(const int *rn,
                           const double *rx1L,
                           const double *rx1R,
                           const double *rx2L,
                           const double *rx2R,
                           const double *rf,
                           double *grid);

void C_create_ENPs_grid_3D(const int *rn,
                           const double *rx1L,
                           const double *rx1R,
                           const double *rx2L,
                           const double *rx2R,
                           const double *rx3L,
                           const double *rx3R,
                           const double *rf,
                           double *grid);

#ifdef __cplusplus
}
#endif

#endif // GFLOW_GRIDS_H_
