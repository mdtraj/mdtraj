#include "ssetools.h"
#include <mmintrin.h>
#include "stdio.h"
#include "stdlib.h"

int main(int argc, char** argv) {
    int i;
    int n_frames = 10;
    int n_atoms = 5;
    int n_sphere_points = 960;
    int pairs[2] = {0, 1};
    int triplets[3] = {0, 1, 2};
    int quartets[4] = {0, 1, 2, 3};
    float box_matrix[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    float * atom_radii = calloc(n_atoms, sizeof(float));
    float * xyzlist = calloc(n_frames*n_atoms*3, sizeof(float));
    float * geom_out = calloc(n_frames, sizeof(float));
    float * array_of_areas = calloc(n_frames*n_atoms, sizeof(float));
    for (i = 0; i < n_frames*n_atoms*3; i++)
        xyzlist[i] = i;

    sasa(n_frames, n_atoms, xyzlist, atom_radii, n_sphere_points, array_of_areas);
    dihedral(xyzlist, quartets, geom_out, n_frames, n_atoms, 1);
    dist(xyzlist, pairs, geom_out, NULL, n_frames, n_atoms, 1);
    dist(xyzlist, pairs, NULL, geom_out, n_frames, n_atoms, 1);
    dist_mic(xyzlist, pairs, box_matrix, geom_out, NULL, n_frames, n_atoms, 1);
    dist_mic(xyzlist, pairs, box_matrix, NULL, geom_out, n_frames, n_atoms, 1);
    angle(xyzlist, triplets, geom_out, n_frames, n_atoms, 1);




    // example: this is an invalid read that should trigger valgrid
    // or the clang address tool
    // __m128 f = _mm_loadu_ps(atom_radii+(n_atoms-3));
    
    free(atom_radii);
    free(xyzlist);
    free(array_of_areas);
    free(geom_out);
}
         
