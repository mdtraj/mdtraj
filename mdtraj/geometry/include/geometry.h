#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_
#ifdef __cplusplus
extern "C" {
#endif

void dist_mic(const float* xyz, const int* pairs, const float* box_matrix,
              float* distance_out, float* displacement_out,
              const int n_frames, const int n_atoms, const int n_pairs);

void dist_mic_t(const float* xyz, const int* pairs, const int* times,
                const float* box_matrix, float* distance_out,
                float* displacement_out, const int n_frames, const int n_atoms,
                const int n_pairs);

void dist(const float* xyz, const int* pairs, float* distance_out,
          float* displacement_out, const int n_frames, const int n_atoms,
          const int n_pairs);

void dist_t(const float* xyz, const int* pairs, const int* times,
        float* distance_out, float* displacement_out, const int n_frames, 
        const int n_atoms, const int n_pairs);

void dist_mic_triclinic(const float* xyz, const int* pairs, const float* box_matrix,
                        float* distance_out, float* displacement_out,
                        const int n_frames, const int n_atoms, const int n_pairs);

void dist_mic_triclinic_t(const float* xyz, const int* pairs, const int* times, 
                        const float* box_matrix,
                        float* distance_out, float* displacement_out,
                        const int n_frames, const int n_atoms, const int n_pairs);

void angle(const float* xyz, const int* triplets, float* out,
           const int n_frames, const int n_atoms, const int n_angles);

void angle_mic(const float* xyz, const int* triplets,
               const float* box_matrix, float* out,
               const int n_frames, const int n_atoms, const int n_angles);

void angle_mic_triclinic(const float* xyz, const int* triplets,
                         const float* box_matrix, float* out,
                         const int n_frames, const int n_atoms, const int n_angles);

void dihedral(const float* xyz, const int* quartets, float* out,
              const int n_frames, const int n_atoms, const int n_quartets);

void dihedral_mic(const float* xyz, const int* quartets,
                  const float* box_matrix, float* out,
                  const int n_frames, const int n_atoms, const int n_quartets);

void dihedral_mic_triclinic(const float* xyz, const int* quartets,
                            const float* box_matrix, float* out,
                            const int n_frames, const int n_atoms, const int n_quartets);

void kabsch_sander(const float* xyz, const int* nco_indices, const int* ca_indices,
                   const int* is_proline, const int n_frames, const int n_atoms,
                   const int n_residues, int* hbonds, float* henergies);

void dssp(const float* xyz, const int* nco_indices, const int* ca_indices,
          const int* is_proline, const int* chains_ids, const int n_frames,
          const int n_atoms, const int n_residues, char* secondary);

void find_closest_contact(const float* positions, const int* group1, const int* group2, int n_group1, int n_group2,
                          const float* box_vectors_pointer, int* atom1, int* atom2, float* distance);

#ifdef __cplusplus
}
#endif
#endif
