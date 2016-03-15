#ifndef MDTRAJ_NEIGHBORLIST_H_
#define MDTRAJ_NEIGHBORLIST_H_

#include <vector>

std::vector<std::vector<int> > _compute_neighborlist(const float* atomLocations, int numAtoms, float maxDistance, const float* boxVectors);

#endif
