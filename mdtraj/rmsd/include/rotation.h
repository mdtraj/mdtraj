#ifndef _ROTATION_H_
#define _ROTATION_H_
#ifdef __cplusplus
extern "C" {
#endif

void rot_atom_major(const int n_atoms, float* a, const float rot[9]);

float rot_msd_atom_major(const int n_real_atoms, const int n_padded_atoms,
                         const float* a, const float* b, const float rot[9]);

void sgemm33(const float A[9], const float B[9], float out[9]);

#ifdef __cplusplus
}
#endif
#endif /* _ROTATION_H_ */
