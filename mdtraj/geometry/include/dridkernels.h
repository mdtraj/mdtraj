#ifndef _MDTRAJ_DRIDKERNELS_H_
#define _MDTRAJ_DRIDKERNELS_H_

#ifdef _MSC_VER
typedef __int32 int32_t;
#else
#include <stdint.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

void drid_moments(float* coords, int32_t index, int32_t* partners, int32_t n_partners, double* moments);

#ifdef __cplusplus
}
#endif

#endif
