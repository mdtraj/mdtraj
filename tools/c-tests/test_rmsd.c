#include "theobald_rmsd.h"
#include "stdlib.h"
#include "stdio.h"

int main(int argc, char** argv) {
    int natoms = 10;
    int rowstride = natoms;
    float* a = calloc(natoms*3, sizeof(float));
    float* b = calloc(natoms*3, sizeof(float));
    float* rot = calloc(9, sizeof(float));
    float G_a = 0;
    float G_b = 0;

    msd_axis_major(natoms, natoms, natoms, a, b, G_a, G_b);
    msd_atom_major(natoms, natoms, a, b, G_a, G_b, 0, NULL);
    msd_atom_major(natoms, natoms, a, b, G_a, G_b, 1, rot);

    free(a);
    free(b);
    free(rot);

    return 0;
}

