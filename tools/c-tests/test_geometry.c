#include "ssetools.h"
#include <mmintrin.h>
#include "stdio.h"
#include "stdlib.h"

int dist_mic(const float* xyz, const int* pairs, const float* box_matrix,
             float* distance_out, float* displacement_out,
             const int n_frames, const int n_atoms, const int n_pairs);

int main(int argc, char** argv) {
    int i, t;
    int n_frames = 10;
    int n_atoms = 5;
    int n_sphere_points = 960;
    int* pairs = calloc(2, sizeof(int));
    pairs[0] = 0;
    pairs[1] = 1;
    int triplets[3] = {0, 1, 2};
    int quartets[4] = {0, 1, 2, 3};
    float * box_matrix = calloc(n_frames*9, sizeof(float));
    for (t = 0; t < n_frames; t++)
        for (i = 0; i < 3; i++)
            box_matrix[t*9 + i*3+i] = 1;
    float * atom_radii = calloc(n_atoms, sizeof(float));
    float * xyzlist = calloc(n_frames*n_atoms*3, sizeof(float));
    float * geom_out = calloc(n_frames, sizeof(float));
    float * displacement_out = calloc(n_frames*1*3, sizeof(float));
    float * array_of_areas = calloc(n_frames*n_atoms, sizeof(float));
    for (i = 0; i < n_frames*n_atoms*3; i++)
        xyzlist[i] = i;

    sasa(n_frames, n_atoms, xyzlist, atom_radii, n_sphere_points, array_of_areas);
    printf("Finished sasa\n");
    dihedral(xyzlist, quartets, geom_out, n_frames, n_atoms, 1);
    printf("finished dihedral\n");
    dist(xyzlist, pairs, geom_out, NULL, n_frames, n_atoms, 1);
    printf("finished dist 1\n");
    dist(xyzlist, pairs, NULL, displacement_out, n_frames, n_atoms, 1);
    printf("finished dist 2\n");
    dist_mic(xyzlist, pairs, box_matrix, geom_out, NULL, n_frames, n_atoms, 1);
    printf("finished mic 1\n");
    dist_mic(xyzlist, pairs, box_matrix, NULL, displacement_out, n_frames, n_atoms, 1);
    printf("finished mic 2\n");
    angle(xyzlist, triplets, geom_out, n_frames, n_atoms, 1);
    printf("finished angle\n");




    // example: this is an invalid read that should trigger valgrid
    // or the clang address tool
    // __m128 f = _mm_loadu_ps(atom_radii+(n_atoms-3));
    
    free(atom_radii);
    free(xyzlist);
    free(array_of_areas);
    free(geom_out);
    free(displacement_out);
    free(pairs);
    free(box_matrix);
}
         
