#include <vector>
#ifndef __PERMUTATION_MSD__
#define __PERMUTATION_MSD__

std::vector<int> euclidean_permutation(
float* target,
float* reference,
int n_atoms,
int n_dims,
std::vector<std::vector<int> >& permute_groups);


#endif // __PERMUTATION_MSD__
