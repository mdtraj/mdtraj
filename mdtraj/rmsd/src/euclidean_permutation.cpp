#include "Munkres.h"
#include "euclidean_permutation.hpp"

#include "stdio.h"
#include <vector>
#include <limits>       // std::numeric_limits

template <class T> inline T square(T x) { return x * x; };

/*
 *
 *
 *
 */
std::vector<int> euclidean_permutation(
float* target,
float* reference,
int n_atoms,
int n_dims,
std::vector<std::vector<int> >& permute_groups)
{
    // the cost matrix A[i,j] of matching target[i] to reference[j]
    std::vector<double> A(n_atoms * n_atoms, std::numeric_limits<double>::max());
    // is atom [i] in a permute group?
    std::vector<int> in_permute_group(n_atoms, 0);

    // fill in the cost matrix with the distance between all pairs that are in
    // the same permute group
    for (int g = 0; g < (int) permute_groups.size(); g++) {
        for (int i = 0; i < (int) permute_groups[g].size(); i++) {
            const int ii = permute_groups[g][i];
            in_permute_group[ii] = 1;
            for (int j = 0; j < (int) permute_groups[g].size(); j++) {
                const int jj = permute_groups[g][j];
                double sq_euclidean_ii_jj = 0;
                for (int d = 0; d < n_dims; d++)
                    sq_euclidean_ii_jj +=  square(target[ii*n_dims + d] - reference[jj*n_dims + d]);
                A[ii*n_atoms + jj] = sq_euclidean_ii_jj;
            }
        }
    }

    // set the diagonal entries of the cost matrix for elements that are not
    // in a permute group
    for (int i = 0; i < n_atoms; i++) {
        if (!in_permute_group[i]) {
            double sq_euclidean_i_i = 0;
            for (int d = 0; d < n_dims; d++)
                sq_euclidean_i_i += square(target[i*n_dims + d] - reference[i*n_dims + d]);
            A[i*n_atoms + i] = sq_euclidean_i_i;
        }
    }

    // solve the assignment problem with this cost matrix
    Munkres munk;
    std::vector<int> mask(n_atoms * n_atoms);
    // printf("running solve...\n");
    munk.solve(&A[0], &mask[0], n_atoms, n_atoms);
    // printf("done\n");

    std::vector<int> mapping(n_atoms);

    for (int i = 0; i < n_atoms; i++) {
        for (int j = 0; j < n_atoms; j++) {
            if (mask[i*n_atoms + j]) {
                mapping[i] = j;
                break;
            }
        }
    }

    return mapping;
}

