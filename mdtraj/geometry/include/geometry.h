#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_
#ifdef __cplusplus
extern "C" {
#endif

int dist_mic(const float* xyz, const int* pairs, const float* box_matrix,
             float* distance_out, float* displacement_out,
             const int n_frames, const int n_atoms, const int n_pairs);

int dist(const float* xyz, const int* pairs, float* distance_out,
         float* displacement_out, const int n_frames, const int n_atoms,
         const int n_pairs);

int angle(const float* xyz, const int* triplets, float* out,
          const int n_frames, const int n_atoms, const int n_angles);

int angle_mic(const float* xyz, const int* triplets,
              const float* box_matrix, float* out,
              const int n_frames, const int n_atoms, const int n_angles);

int dihedral(const float* xyz, const int* quartets, float* out,
             const int n_frames, const int n_atoms, const int n_quartets);

int dihedral_mic(const float* xyz, const int* quartets,
                 const float* box_matrix, float* out,
                 const int n_frames, const int n_atoms, const int n_quartets);

int kabsch_sander(const float* xyz, const int* nco_indices, const int* ca_indices,
                  const int* is_proline, const int n_frames, const int n_atoms,
                  const int n_residues, int* hbonds, float* henergies);

int dssp(const float* xyz, const int* nco_indices, const int* ca_indices,
         const int* is_proline, const int* chains_ids, const int n_frames,
         const int n_atoms, const int n_residues, char* secondary);

#ifdef __cplusplus
}
#endif
#endif
