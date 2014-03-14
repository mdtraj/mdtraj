#include <vector>
#ifndef __PERMUTATION_MSD__
#define __PERMUTATION_MSD__

std::vector<long> euclidean_permutation(
float* target,
float* reference,
long n_atoms,
long n_dims,
std::vector<std::vector<long> >& permute_groups);


#endif // __PERMUTATION_MSD__
