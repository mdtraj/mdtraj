#ifndef _THEOBALD_RMSD_H_
#define _THEOBALD_RMSD_H_
#ifdef __cplusplus
extern "C" {
#endif

int solve_cubic_equation(double  c3, double  c2,  double c1, double c0,
                         double *x1, double *x2, double *x3);

int quartic_equation_solve_exact(double *r1, double *r2, double *r3, double *r4,
				 int *nr12, int *nr34,double d0,double d1,double d2, double d3, double d4);

float msd_axis_major(const int nrealatoms, const int npaddedatoms, const int rowstride,
                     const float* aT, const float* bT, const float G_a, const float G_b);

float msd_atom_major(const int nrealatoms, const int npaddedatoms,
                     const float* a, const float* b, const float G_a, const float G_b,
                     int computeRot, float rot[9]);

#ifdef __cplusplus
}
#endif
#endif
