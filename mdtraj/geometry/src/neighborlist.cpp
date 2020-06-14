/**
 * MDTraj: A Python Library for Loading, Saving, and Manipulating
 *         Molecular Dynamics Trajectories.
 * Copyright 2013-2016 Stanford University and the Authors
 *
 * Authors: Peter Eastman
 * Contributors:
 *
 * MDTraj is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 2.1
 * of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.
 */

#include "vectorize.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include "math_patch.h"

using namespace std;

class VoxelIndex
{
public:
    VoxelIndex() : y(0), z(0) {
    }
    VoxelIndex(int y, int z) : y(y), z(z) {
    }
    int y;
    int z;
};

/**
 * This data structure organizes the particles spatially.  It divides them into bins along the x and y axes,
 * then sorts each bin along the z axis so ranges can be identified quickly with a binary search.
 */
class Voxels {
public:
    Voxels(float vsy, float vsz, float miny, float maxy, float minz, float maxz, const float (&periodicBoxVectors)[3][3], bool usePeriodic) :
            voxelSizeY(vsy), voxelSizeZ(vsz), miny(miny), minz(minz), periodicBoxVectors(periodicBoxVectors), usePeriodic(usePeriodic) {
        periodicBoxSize[0] = periodicBoxVectors[0][0];
        periodicBoxSize[1] = periodicBoxVectors[1][1];
        periodicBoxSize[2] = periodicBoxVectors[2][2];
        recipBoxSize[0] = 1/periodicBoxVectors[0][0];
        recipBoxSize[1] = 1/periodicBoxVectors[1][1];
        recipBoxSize[2] = 1/periodicBoxVectors[2][2];
        triclinic = (periodicBoxVectors[0][1] != 0.0 || periodicBoxVectors[0][2] != 0.0 ||
                     periodicBoxVectors[1][0] != 0.0 || periodicBoxVectors[1][2] != 0.0 ||
                     periodicBoxVectors[2][0] != 0.0 || periodicBoxVectors[2][1] != 0.0);
        if (usePeriodic) {
            ny = max(1, (int) floorf(periodicBoxVectors[1][1]/voxelSizeY+0.5f));
            nz = max(1, (int) floorf(periodicBoxVectors[2][2]/voxelSizeZ+0.5f));
            voxelSizeY = periodicBoxVectors[1][1]/ny;
            voxelSizeZ = periodicBoxVectors[2][2]/nz;
        } else {
            ny = max(1, (int) floorf((maxy-miny)/voxelSizeY+0.5f));
            nz = max(1, (int) floorf((maxz-minz)/voxelSizeZ+0.5f));
            if (maxy > miny)
                voxelSizeY = (maxy-miny)/ny;
            if (maxz > minz)
                voxelSizeZ = (maxz-minz)/nz;
        }
        bins.resize(ny);
        for (int i = 0; i < ny; i++) {
            bins[i].resize(nz);
            for (int j = 0; j < nz; j++)
                bins[i][j].resize(0);
        }
    }

    /**
     * Insert a particle into the voxel data structure.
     */
    void insert(const int& atom, const float* location) {
        VoxelIndex voxelIndex = getVoxelIndex(location);
        bins[voxelIndex.y][voxelIndex.z].push_back(make_pair(location[0], atom));
    }

    /**
     * Sort the particles in each voxel by x coordinate.
     */
    void sortItems() {
        for (int i = 0; i < ny; i++)
            for (int j = 0; j < nz; j++)
                sort(bins[i][j].begin(), bins[i][j].end());
    }

    /**
     * Find the index of the first particle in voxel (y,z) whose x coordinate is >= the specified value.
     */
    int findLowerBound(int y, int z, double x, int lower, int upper) const {
        const vector<pair<float, int> >& bin = bins[y][z];
        while (lower < upper) {
            int middle = (lower+upper)/2;
            if (bin[middle].first < x)
                lower = middle+1;
            else
                upper = middle;
        }
        return lower;
    }

    /**
     * Find the index of the first particle in voxel (y,z) whose x coordinate is greater than the specified value.
     */
    int findUpperBound(int y, int z, double x, int lower, int upper) const {
        const vector<pair<float, int> >& bin = bins[y][z];
        while (lower < upper) {
            int middle = (lower+upper)/2;
            if (bin[middle].first > x)
                upper = middle;
            else
                lower = middle+1;
        }
        return upper;
    }

    /**
     * Get the voxel index containing a particular location.
     */
    VoxelIndex getVoxelIndex(const float* location) const {
        float yperiodic, zperiodic;
        if (!usePeriodic) {
            yperiodic = location[1]-miny;
            zperiodic = location[2]-minz;
        }
        else {
            float scale2 = floorf(location[2]*recipBoxSize[2]);
            yperiodic = location[1]-periodicBoxVectors[2][1]*scale2;
            zperiodic = location[2]-periodicBoxVectors[2][2]*scale2;
            float scale1 = floorf(yperiodic*recipBoxSize[1]);
            yperiodic -= periodicBoxVectors[1][0]*scale1;
        }
        int y = max(0, min(ny-1, int(floorf(yperiodic / voxelSizeY))));
        int z = max(0, min(nz-1, int(floorf(zperiodic / voxelSizeZ))));

        return VoxelIndex(y, z);
    }

    void getNeighbors(vector<int>& neighbors, int atomIndex, float maxDistance, const float* atomLocations, VoxelIndex atomVoxelIndex) const {
        neighbors.resize(0);
        fvec4 boxSize(periodicBoxSize[0], periodicBoxSize[1], periodicBoxSize[2], 0);
        fvec4 invBoxSize(recipBoxSize[0], recipBoxSize[1], recipBoxSize[2], 0);
        fvec4 periodicBoxVec4[3];
        periodicBoxVec4[0] = fvec4(periodicBoxVectors[0][0], periodicBoxVectors[0][1], periodicBoxVectors[0][2], 0);
        periodicBoxVec4[1] = fvec4(periodicBoxVectors[1][0], periodicBoxVectors[1][1], periodicBoxVectors[1][2], 0);
        periodicBoxVec4[2] = fvec4(periodicBoxVectors[2][0], periodicBoxVectors[2][1], periodicBoxVectors[2][2], 0);

        float maxDistanceSquared = maxDistance * maxDistance;

        int dIndexY = int(maxDistance/voxelSizeY)+1; // How may voxels away do we have to look?
        int dIndexZ = int(maxDistance/voxelSizeZ)+1;
        if (usePeriodic) {
            dIndexY = min(ny/2, dIndexY);
            dIndexZ = min(nz/2, dIndexZ);
        }
        const float* centerAtomPos = &atomLocations[3*atomIndex];
        fvec4 centerPosVec(centerAtomPos[0], centerAtomPos[1], centerAtomPos[2], 0);

        // Loop over voxels along the z axis.

        int startz = atomVoxelIndex.z-dIndexZ;
        int endz = atomVoxelIndex.z+dIndexZ;
        if (usePeriodic)
            endz = min(endz, startz+nz-1);
        else {
            startz = max(startz, 0);
            endz = min(endz, nz-1);
        }
        VoxelIndex voxelIndex(0, 0);
        for (int z = startz; z <= endz; ++z) {
            voxelIndex.z = z;
            if (usePeriodic)
                voxelIndex.z = (z < 0 ? z+nz : (z >= nz ? z-nz : z));

            // Loop over voxels along the y axis.

            float boxz = floor((float) z/nz);
            int starty = atomVoxelIndex.y-dIndexY;
            int endy = atomVoxelIndex.y+dIndexY;
            float yoffset = (float) (usePeriodic ? boxz*periodicBoxVectors[2][1] : 0);
            if (usePeriodic) {
                starty -= (int) ceil(yoffset/voxelSizeY);
                endy -= (int) floor(yoffset/voxelSizeY);
                endy = min(endy, starty+ny-1);
            }
            else {
                starty = max(starty, 0);
                endy = min(endy, ny-1);
            }
            for (int y = starty; y <= endy; ++y) {
                voxelIndex.y = y;
                if (usePeriodic)
                    voxelIndex.y = (y < 0 ? y+ny : (y >= ny ? y-ny : y));
                float boxy = floor((float) y/ny);

                // Identify the range of atoms within this bin we need to search.  When using periodic boundary
                // conditions, there may be two separate ranges.

                float minx = centerAtomPos[0];
                float maxx = centerAtomPos[0];
                if (usePeriodic && triclinic) {
                    fvec4 delta1(0, voxelSizeY*voxelIndex.y-centerPosVec[1], voxelSizeZ*voxelIndex.z-centerPosVec[2], 0);
                    fvec4 delta2 = delta1+fvec4(0, voxelSizeY, 0, 0);
                    fvec4 delta3 = delta1+fvec4(0, 0, voxelSizeZ, 0);
                    fvec4 delta4 = delta1+fvec4(0, voxelSizeY, voxelSizeZ, 0);
                    delta1 -= periodicBoxVec4[2]*floorf(delta1[2]*recipBoxSize[2]+0.5f);
                    delta1 -= periodicBoxVec4[1]*floorf(delta1[1]*recipBoxSize[1]+0.5f);
                    delta1 -= periodicBoxVec4[0]*floorf(delta1[0]*recipBoxSize[0]+0.5f);
                    delta2 -= periodicBoxVec4[2]*floorf(delta2[2]*recipBoxSize[2]+0.5f);
                    delta2 -= periodicBoxVec4[1]*floorf(delta2[1]*recipBoxSize[1]+0.5f);
                    delta2 -= periodicBoxVec4[0]*floorf(delta2[0]*recipBoxSize[0]+0.5f);
                    delta3 -= periodicBoxVec4[2]*floorf(delta3[2]*recipBoxSize[2]+0.5f);
                    delta3 -= periodicBoxVec4[1]*floorf(delta3[1]*recipBoxSize[1]+0.5f);
                    delta3 -= periodicBoxVec4[0]*floorf(delta3[0]*recipBoxSize[0]+0.5f);
                    delta4 -= periodicBoxVec4[2]*floorf(delta4[2]*recipBoxSize[2]+0.5f);
                    delta4 -= periodicBoxVec4[1]*floorf(delta4[1]*recipBoxSize[1]+0.5f);
                    delta4 -= periodicBoxVec4[0]*floorf(delta4[0]*recipBoxSize[0]+0.5f);
                    if (delta1[1] < 0 && delta1[1]+voxelSizeY > 0)
                        delta1 = fvec4(delta1[0], 0, delta1[2], 0);
                    if (delta1[2] < 0 && delta1[2]+voxelSizeZ > 0)
                        delta1 = fvec4(delta1[0], delta1[1], 0, 0);
                    if (delta3[1] < 0 && delta3[1]+voxelSizeY > 0)
                        delta3 = fvec4(delta3[0], 0, delta3[2], 0);
                    if (delta2[2] < 0 && delta2[2]+voxelSizeZ > 0)
                        delta2 = fvec4(delta2[0], delta2[1], 0, 0);
                    fvec4 delta = min(min(min(abs(delta1), abs(delta2)), abs(delta3)), abs(delta4));
                    float dy = (voxelIndex.y == atomVoxelIndex.y ? 0.0f : delta[1]);
                    float dz = (voxelIndex.z == atomVoxelIndex.z ? 0.0f : delta[2]);
                    float dist2 = maxDistanceSquared-dy*dy-dz*dz;
                    if (dist2 > 0) {
                        float dist = sqrtf(dist2);
                        minx = min(minx, centerAtomPos[0]-dist-max(max(max(delta1[0], delta2[0]), delta3[0]), delta4[0]));
                        maxx = max(maxx, centerAtomPos[0]+dist-min(min(min(delta1[0], delta2[0]), delta3[0]), delta4[0]));
                    }
                }
                else {
                    float xoffset = (float) (usePeriodic ? boxy*periodicBoxVectors[1][0]+boxz*periodicBoxVectors[2][0] : 0);
                    fvec4 offset(-xoffset, -yoffset+voxelSizeY*y+(usePeriodic ? 0.0f : miny), voxelSizeZ*z+(usePeriodic ? 0.0f : minz), 0);
                    fvec4 delta1 = offset-centerPosVec;
                    fvec4 delta2 = delta1+fvec4(0, voxelSizeY, voxelSizeZ, 0);
                    if (usePeriodic) {
                        delta1 -= round(delta1*invBoxSize)*boxSize;
                        delta2 -= round(delta2*invBoxSize)*boxSize;
                    }
                    fvec4 delta = min(abs(delta1), abs(delta2));
                    float dy = (y == atomVoxelIndex.y ? 0.0f : delta[1]);
                    float dz = (z == atomVoxelIndex.z ? 0.0f : delta[2]);
                    float dist2 = maxDistanceSquared-dy*dy-dz*dz;
                    if (dist2 > 0) {
                        float dist = sqrtf(dist2);
                        minx = min(minx, centerAtomPos[0]-dist-xoffset);
                        maxx = max(maxx, centerAtomPos[0]+dist-xoffset);
                    }
                }
                if (minx == maxx)
                    continue;
                bool needPeriodic = usePeriodic && (centerAtomPos[1] < maxDistance || centerAtomPos[1] > periodicBoxSize[1]-maxDistance ||
                                                    centerAtomPos[2] < maxDistance || centerAtomPos[2] > periodicBoxSize[2]-maxDistance ||
                                                    minx < 0.0f || maxx > periodicBoxVectors[0][0]);
                int numRanges;
                int rangeStart[2];
                int rangeEnd[2];
                int binSize = bins[voxelIndex.y][voxelIndex.z].size();
                rangeStart[0] = findLowerBound(voxelIndex.y, voxelIndex.z, minx, 0, binSize);
                if (needPeriodic) {
                    numRanges = 2;
                    rangeEnd[0] = findUpperBound(voxelIndex.y, voxelIndex.z, maxx, rangeStart[0], binSize);
                    if (rangeStart[0] > 0 && rangeEnd[0] < binSize)
                        numRanges = 1;
                    else if (rangeStart[0] > 0) {
                        rangeStart[1] = 0;
                        rangeEnd[1] = min(findUpperBound(voxelIndex.y, voxelIndex.z, maxx-periodicBoxSize[0], 0, rangeStart[0]), rangeStart[0]);
                    }
                    else {
                        rangeStart[1] = max(findLowerBound(voxelIndex.y, voxelIndex.z, minx+periodicBoxSize[0], rangeEnd[0], binSize), rangeEnd[0]);
                        rangeEnd[1] = bins[voxelIndex.y][voxelIndex.z].size();
                    }
                }
                else {
                    numRanges = 1;
                    rangeEnd[0] = findUpperBound(voxelIndex.y, voxelIndex.z, maxx, rangeStart[0], binSize);
                }
                bool periodicRectangular = (needPeriodic && !triclinic);

                // Loop over atoms and check to see if they are neighbors of this one.

                const vector<pair<float, int> >& voxelBins = bins[voxelIndex.y][voxelIndex.z];
                for (int range = 0; range < numRanges; range++) {
                    for (int item = rangeStart[range]; item < rangeEnd[range]; item++) {
                        const int index = voxelBins[item].second;

                        // Avoid duplicate entries.
                        if (index >= atomIndex)
                            continue;

                        fvec4 atomPos(atomLocations[3*index], atomLocations[3*index+1], atomLocations[3*index+2], 0);
                        fvec4 delta = atomPos-centerPosVec;
                        if (periodicRectangular)
                            delta -= round(delta*invBoxSize)*boxSize;
                        else if (needPeriodic) {
                            delta -= periodicBoxVec4[2]*floorf(delta[2]*recipBoxSize[2]+0.5f);
                            delta -= periodicBoxVec4[1]*floorf(delta[1]*recipBoxSize[1]+0.5f);
                            delta -= periodicBoxVec4[0]*floorf(delta[0]*recipBoxSize[0]+0.5f);
                        }
                        delta = abs(delta);
                        float dSquared = dot3(delta, delta);
                        if (dSquared > maxDistanceSquared)
                            continue;

                        // Add this atom to the list of neighbors.

                        neighbors.push_back(index);
                    }
                }
            }
        }
    }

private:
    float voxelSizeY, voxelSizeZ;
    float miny, minz;
    int ny, nz;
    float periodicBoxSize[3], recipBoxSize[3];
    bool triclinic;
    const float (&periodicBoxVectors)[3][3];
    const bool usePeriodic;
    vector<vector<vector<pair<float, int> > > > bins;
};

vector<vector<int> > _compute_neighborlist(const float* atomLocations, int numAtoms, float maxDistance, const float* boxVectors) {
    float periodicBoxVectors[3][3];
    bool usePeriodic = (boxVectors != NULL);
    if (usePeriodic) {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                periodicBoxVectors[i][j] = boxVectors[3*i+j];

        // Make sure box vectors are in reduced form.

        for (int i = 0; i < 3; i++)
            periodicBoxVectors[2][i] -= periodicBoxVectors[1][i]*roundf(periodicBoxVectors[2][1]/periodicBoxVectors[1][1]);
        for (int i = 0; i < 3; i++)
            periodicBoxVectors[2][i] -= periodicBoxVectors[0][i]*roundf(periodicBoxVectors[2][0]/periodicBoxVectors[0][0]);
        for (int i = 0; i < 3; i++)
            periodicBoxVectors[1][i] -= periodicBoxVectors[0][i]*roundf(periodicBoxVectors[1][0]/periodicBoxVectors[0][0]);
    }

    // Identify the range of atom positions along each axis.

    fvec4 minPos(atomLocations[0], atomLocations[1], atomLocations[2], 0);
    fvec4 maxPos = minPos;
    for (int i = 0; i < numAtoms; i++) {
        fvec4 pos(atomLocations[3*i], atomLocations[3*i+1], atomLocations[3*i+2], 0);
        minPos = min(minPos, pos);
        maxPos = max(maxPos, pos);
    }

    // Build the voxel hash.

    float edgeSizeY, edgeSizeZ;
    if (!usePeriodic)
        edgeSizeY = edgeSizeZ = maxDistance; // TODO - adjust this as needed
    else {
        edgeSizeY = 0.6f*periodicBoxVectors[1][1]/floorf(periodicBoxVectors[1][1]/maxDistance);
        edgeSizeZ = 0.6f*periodicBoxVectors[2][2]/floorf(periodicBoxVectors[2][2]/maxDistance);
    }
    Voxels voxels(edgeSizeY, edgeSizeZ, minPos[1], maxPos[1], minPos[2], maxPos[2], periodicBoxVectors, usePeriodic);
    for (int i = 0; i < numAtoms; i++)
        voxels.insert(i, &atomLocations[3*i]);
    voxels.sortItems();

    // Compute neighbors.

    vector<vector<int> > neighbors(numAtoms);
#pragma omp parallel for default(shared)
    for (int i = 0; i < numAtoms; i++)
        voxels.getNeighbors(neighbors[i], i, maxDistance, atomLocations, voxels.getVoxelIndex(&atomLocations[3*i]));

    // Add in the symmetric entries.

    for (int i = 0; i < numAtoms; i++)
        for (int j = 0; j < neighbors[i].size(); j++)
            neighbors[neighbors[i][j]].push_back(i);
    return neighbors;
}
