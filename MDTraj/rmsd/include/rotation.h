#ifndef _ROTATION_H_
#define _ROTATION_H_

float rot_msd_atom_major(const int n_real_atoms, const int n_padded_atoms,
                         const float* a, const float* b, const float rot[9]);

#endif /* _ROTATION_H_ */