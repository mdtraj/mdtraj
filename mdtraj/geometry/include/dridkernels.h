#ifndef _MDTRAJ_DRIDKERNELS_H_
#define _MDTRAJ_DRIDKERNELS_H_

#ifdef __cplusplus
extern "C" {
#endif

int drid_moments(float* coords, int index, int* partners, int n_partners, double* moments);

#ifdef __cplusplus
}
#endif

#endif
