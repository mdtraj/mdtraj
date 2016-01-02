#include "geometry.h"
#include "unitcell.h"

#ifndef MIN
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#endif
#ifndef MAX
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#endif
#define CLIP(X, X_min, X_max) (MIN(MAX(X, X_min), X_max))


int compute_distances(
    const float* xyz, const int* pairs, float* distance_out,
    float* displacement_out, const int n_frames, const int n_atoms,
    const int n_pairs)
{
    int store_displacement = displacement_out == NULL ? 0 : 1;
    int store_distance = distance_out == NULL ? 0 : 1;
    fvec4 ones(1,1,1,1);

    for (int i = 0; i < n_frames; i++) {
        for (int j = 0; j < n_pairs; j++) {
            fvec4 x0 = fvec4(load_float3(&xyz[i*n_atoms*3 + 3*pairs[2*j + 0]]));
            fvec4 x1 = fvec4(load_float3(&xyz[i*n_atoms*3 + 3*pairs[2*j + 1]]));
            fvec4 r12 = x1 - x0;

            if (store_displacement) {
                float f_r12[4];
                r12.store(f_r12);

                displacement_out[i*n_pairs*3 + 3*j + 0] = f_r12[0];
                displacement_out[i*n_pairs*3 + 3*j + 1] = f_r12[1];
                displacement_out[i*n_pairs*3 + 3*j + 2] = f_r12[2];
            }

            if (store_distance) {
                distance_out[i*n_pairs + j] = sqrt(dot3(r12, r12));
            }
        }
    }

    return 0;
}


int compute_distances_orthorhombic(
    const float* xyz, const int* pairs, const float* box_matrix,
    float* distance_out, float* displacement_out,
    const int n_frames, const int n_atoms, const int n_pairs)
{
    int store_displacement = displacement_out == NULL ? 0 : 1;
    int store_distance = distance_out == NULL ? 0 : 1;

    for (int i = 0; i < n_frames; i++) {
        fvec4 h(box_matrix[9*i + 0], box_matrix[9*i + 4], box_matrix[9*i + 8], 0);

        for (int j = 0; j < n_pairs; j++) {
            fvec4 x0 = fvec4(load_float3(&xyz[i*n_atoms*3 + 3*pairs[2*j + 0]]));
            fvec4 x1 = fvec4(load_float3(&xyz[i*n_atoms*3 + 3*pairs[2*j + 1]]));
            fvec4 r12 = x1 - x0;

            fvec4 s = r12 / h;
            r12 = h * (s - round(s));

            if (store_displacement) {
                float f_r12[4];
                r12.store(f_r12);

                displacement_out[i*n_pairs*3 + 3*j + 0] = f_r12[0];
                displacement_out[i*n_pairs*3 + 3*j + 1] = f_r12[1];
                displacement_out[i*n_pairs*3 + 3*j + 2] = f_r12[2];
            }

            if (store_distance) {
                distance_out[i*n_pairs + j] = sqrt(dot3(r12, r12));
            }
        }
    }

    return 0;
}


int compute_distances_triclinic(
    const float* xyz, const int* pairs, const float* box_matrix,
    float* distance_out, float* displacement_out,
    const int n_frames, const int n_atoms, const int n_pairs)
{
    int store_displacement = displacement_out == NULL ? 0 : 1;
    int store_distance = distance_out == NULL ? 0 : 1;

    for (int i = 0; i < n_frames; i++) {
        UnitCell u(&box_matrix[9*i]);

        for (int j = 0; j < n_pairs; j++) {
            fvec4 x0 = fvec4(load_float3(&xyz[i*n_atoms*3 + 3*pairs[2*j + 0]]));
            fvec4 x1 = fvec4(load_float3(&xyz[i*n_atoms*3 + 3*pairs[2*j + 1]]));
            fvec4 r = x1 - x0;
            r = u.minimum_image(r);

            if (store_displacement) {
                float fr[4];
                r.store(fr);

                displacement_out[i*n_pairs*3 + 3*j + 0] = fr[0];
                displacement_out[i*n_pairs*3 + 3*j + 1] = fr[1];
                displacement_out[i*n_pairs*3 + 3*j + 2] = fr[2];
            }

            if (store_distance) {
                distance_out[i*n_pairs + j] = sqrt(dot3(r, r));
            }

        }
    }

    return 0;
}


int compute_angles(
    const float* xyz, const int* triplets, float* out,
    const int n_frames, const int n_atoms, const int n_angles)
{
    for (int i = 0; i < n_frames; i++) {
        for (int j = 0; j < n_angles; j++) {
            fvec4 x0 = fvec4(load_float3(&xyz[i*n_atoms*3 + 3*triplets[3*j + 0]]));
            fvec4 x1 = fvec4(load_float3(&xyz[i*n_atoms*3 + 3*triplets[3*j + 1]]));
            fvec4 x2 = fvec4(load_float3(&xyz[i*n_atoms*3 + 3*triplets[3*j + 2]]));

            fvec4 u_prime = x0 - x1;
            fvec4 v_prime = x2 - x1;

            fvec4 u = u_prime / sqrt(dot3(u_prime, u_prime));
            fvec4 v = v_prime / sqrt(dot3(v_prime, v_prime));

            out[i*n_angles + j] = static_cast<float>(acos(CLIP(dot3(u, v), -1, 1)));
        }
    }

    return 0;
}


int compute_angles_orthorhombic(
    const float* xyz, const int* triplets,
    const float* box_matrix, float* out,
    const int n_frames, const int n_atoms, const int n_angles)
{
    for (int i = 0; i < n_angles; i++) {
        fvec4 h(box_matrix[9*i + 0], box_matrix[9*i + 4], box_matrix[9*i + 8], 0);

        for (int j = 0; j < n_angles; j++) {
            fvec4 x0 = fvec4(load_float3(&xyz[i*n_atoms*3 + 3*triplets[3*j + 0]]));
            fvec4 x1 = fvec4(load_float3(&xyz[i*n_atoms*3 + 3*triplets[3*j + 1]]));
            fvec4 x2 = fvec4(load_float3(&xyz[i*n_atoms*3 + 3*triplets[3*j + 2]]));

            fvec4 u_prime = x0 - x1;
            fvec4 v_prime = x2 - x1;

            fvec4 s1 = u_prime / h;
            fvec4 s2 = v_prime / h;
            u_prime = h * (s1 - round(s1));
            v_prime = h * (s2 - round(s2));

            fvec4 u = u_prime / sqrt(dot3(u_prime, u_prime));
            fvec4 v = v_prime / sqrt(dot3(v_prime, v_prime));

            out[i*n_angles + j] = static_cast<float>(acos(CLIP(dot3(u, v), -1, 1)));
        }
    }

    return 0;
}


int compute_dihedrals(
    const float* xyz, const int* quartets, float* out,
    const int n_frames, const int n_atoms, const int n_quartets)
{
    for (int i = 0; i < n_frames; i++) {
        for (int j = 0; j < n_quartets; j++) {
            fvec4 x0 = fvec4(load_float3(&xyz[i*n_atoms*3 + 3*quartets[4*j + 0]]));
            fvec4 x1 = fvec4(load_float3(&xyz[i*n_atoms*3 + 3*quartets[4*j + 1]]));
            fvec4 x2 = fvec4(load_float3(&xyz[i*n_atoms*3 + 3*quartets[4*j + 2]]));
            fvec4 x3 = fvec4(load_float3(&xyz[i*n_atoms*3 + 3*quartets[4*j + 3]]));

            fvec4 b1 = x1 - x0;
            fvec4 b2 = x2 - x1;
            fvec4 b3 = x3 - x2;

            fvec4 c1 = cross(b2, b3);
            fvec4 c2 = cross(b1, b2);

            float p1 = dot3(b1, c1) * sqrt(dot3(b2, b2));
            float p2 = dot3(c1, c2);

            out[i*n_quartets + j] = static_cast<float>(atan2(p1, p2));
        }
    }

    return 0;
}


int compute_dihedrals_orthorhombic(
    const float* xyz, const int* quartets,
    const float* box_matrix, float* out,
    const int n_frames, const int n_atoms, const int n_quartets)
{
    for (int i = 0; i < n_frames; i++) {
        fvec4 h(box_matrix[9*i + 0], box_matrix[9*i + 4], box_matrix[9*i + 8], 0);

        for (int j = 0; j < n_quartets; j++) {
            fvec4 x0 = fvec4(load_float3(&xyz[i*n_atoms*3 + 3*quartets[4*j + 0]]));
            fvec4 x1 = fvec4(load_float3(&xyz[i*n_atoms*3 + 3*quartets[4*j + 1]]));
            fvec4 x2 = fvec4(load_float3(&xyz[i*n_atoms*3 + 3*quartets[4*j + 2]]));
            fvec4 x3 = fvec4(load_float3(&xyz[i*n_atoms*3 + 3*quartets[4*j + 3]]));

            fvec4 b1 = x1 - x0;
            fvec4 b2 = x2 - x1;
            fvec4 b3 = x3 - x2;

            fvec4 s1 = b1 / h;
            fvec4 s2 = b2 / h;
            fvec4 s3 = b3 / h;
            b1 = h * (s1 - round(s1));
            b1 = h * (s2 - round(s2));
            b1 = h * (s3 - round(s3));

            fvec4 c1 = cross(b2, b3);
            fvec4 c2 = cross(b1, b2);

            float p1 = dot3(b1, c1) * sqrt(dot3(b2, b2));
            float p2 = dot3(c1, c2);

            out[i*n_quartets + j] = static_cast<float>(atan2(p1, p2));
        }
    }

    return 0;
}
