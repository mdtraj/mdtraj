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
    const float* xyz, const Py_ssize_t* pairs, const float* box_matrix,
    float* distance_out, float* displacement_out,
    const Py_ssize_t n_frames, const Py_ssize_t n_atoms, const Py_ssize_t n_pairs)
{
    const bool store_displacement = displacement_out == NULL ? 0 : 1;
    const bool store_distance = distance_out == NULL ? 0 : 1;
    UnitCell u;

    for (Py_ssize_t i = 0; i < n_frames; i++) {
        if (box_matrix != NULL) {
             u.setBoxMatrix(&box_matrix[9*i]);
        }

        for (Py_ssize_t j = 0; j < n_pairs; j++) {
            fvec4 x0 = load3(&xyz[i*n_atoms*3 + 3*pairs[2*j + 0]]);
            fvec4 x1 = load3(&xyz[i*n_atoms*3 + 3*pairs[2*j + 1]]);
            fvec4 r = x1 - x0;

            if (box_matrix != NULL) {
                r = u.minimumImage(r);
            }

            if (store_displacement) {
                store3(r, &displacement_out[i*n_pairs*3 + 3*j]);
            }

            if (store_distance) {
                distance_out[i*n_pairs + j] = sqrt(dot3(r, r));
            }

        }
    }

    return 0;
}


int compute_angles(
    const float* xyz, const Py_ssize_t* triplets,
    const float* box_matrix, float* out,
    const Py_ssize_t n_frames, const Py_ssize_t n_atoms, const Py_ssize_t n_angles)
{
    UnitCell u;

    for (Py_ssize_t i = 0; i < n_frames; i++) {
        if (box_matrix != NULL) {
             u.setBoxMatrix(&box_matrix[9*i]);
        }

        for (Py_ssize_t j = 0; j < n_angles; j++) {
            fvec4 x0 = load3(&xyz[i*n_atoms*3 + 3*triplets[3*j + 0]]);
            fvec4 x1 = load3(&xyz[i*n_atoms*3 + 3*triplets[3*j + 1]]);
            fvec4 x2 = load3(&xyz[i*n_atoms*3 + 3*triplets[3*j + 2]]);

            fvec4 u_prime = x0 - x1;
            fvec4 v_prime = x2 - x1;

            if (box_matrix != NULL) {
                u_prime = u.minimumImage(u_prime);
                v_prime = u.minimumImage(v_prime);
            }

            fvec4 u = u_prime / sqrt(dot3(u_prime, u_prime));
            fvec4 v = v_prime / sqrt(dot3(v_prime, v_prime));

            out[i*n_angles + j] = static_cast<float>(acos(CLIP(dot3(u, v), -1, 1)));
        }
    }

    return 0;
}

int compute_dihedrals(
    const float* xyz, const Py_ssize_t* quartets,
    const float* box_matrix, float* out,
    const Py_ssize_t n_frames, const Py_ssize_t n_atoms, const Py_ssize_t n_quartets)
{
    UnitCell u;

    for (Py_ssize_t i = 0; i < n_frames; i++) {
        if (box_matrix != NULL) {
             u.setBoxMatrix(&box_matrix[9*i]);
        }

        for (Py_ssize_t j = 0; j < n_quartets; j++) {
            fvec4 x0 = load3(&xyz[i*n_atoms*3 + 3*quartets[4*j + 0]]);
            fvec4 x1 = load3(&xyz[i*n_atoms*3 + 3*quartets[4*j + 1]]);
            fvec4 x2 = load3(&xyz[i*n_atoms*3 + 3*quartets[4*j + 2]]);
            fvec4 x3 = load3(&xyz[i*n_atoms*3 + 3*quartets[4*j + 3]]);

            fvec4 b1 = x1 - x0;
            fvec4 b2 = x2 - x1;
            fvec4 b3 = x3 - x2;

            if (box_matrix != NULL) {
                b1 = u.minimumImage(b1);
                b2 = u.minimumImage(b2);
                b3 = u.minimumImage(b3);
            }

            fvec4 c1 = cross(b2, b3);
            fvec4 c2 = cross(b1, b2);

            float p1 = dot3(b1, c1) * sqrt(dot3(b2, b2));
            float p2 = dot3(c1, c2);

            out[i*n_quartets + j] = static_cast<float>(atan2(p1, p2));
        }
    }

    return 0;
}
