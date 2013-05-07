#ifndef _THEOBALD_RMSD_H_
#define _THEOBALD_RMSD_H_

#ifdef USE_DOUBLE
#define real double
#else
#define real float
#endif
typedef real rvec[3];

int solve_cubic_equation(double  c3, double  c2,  double c1, double c0,
                         double *x1, double *x2, double *x3);

int quartic_equation_solve_exact(double *r1, double *r2, double *r3, double *r4,
				 int *nr12, int *nr34,double d0,double d1,double d2, double d3, double d4);

real ls_rmsd2_aligned_T_g(const int nrealatoms, const int npaddedatoms, const int rowstride,
                          const float* aT, const float* bT, const real G_a, const real G_b);

#endif
