#ifndef _MDTRAJ_SASA_H_
#define _MDTRAJ_SASA_H_
#ifdef __cplusplus
extern "C" {
#endif

void sasa(const int n_frames, const int n_atoms, const float* xyzlist,
          const float* atom_radii, const int n_sphere_points,
          const int* atom_mapping, const int n_groups, float* out);


#ifdef __cplusplus
}
#endif
#endif
