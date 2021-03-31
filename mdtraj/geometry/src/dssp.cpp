/**
 * DSSP secondary structure assignment
 * Copyright [2014] Stanford University and the Authors
 * Authors: Robert T. McGibbon, Maarten L. Hekkelman
 *
 * This code is adapted from DSSP-2.2.0, written by Maarten L. Hekkelman,
 * and ported to MDTraj by Robert T. McGibbon. DSSP-2.2.0 is distributed
 * under the Boost Software License, Version 1.0. This code, as part of
 * MDTraj, is distributed under the GNU LGPL.
 */
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <deque>
#include <set>
#include <map>
#include <algorithm>
#include "geometry.h"
#include "vectorize.h"

#define MIN(X, Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X, Y) ((X) > (Y) ? (X) : (Y))
#define CLIP(X, X_min, X_max) (MIN(MAX(X, X_min), X_max))

//#include <pmmintrin.h>
#include "msvccompat.h"
#include <math.h>

enum helix_flag_t {HELIX_NONE, HELIX_START, HELIX_END, HELIX_START_AND_END, HELIX_MIDDLE};
enum bridge_t {BRIDGE_NONE, BRIDGE_PARALLEL, BRIDGE_ANTIPARALLEL};
enum ss_t {SS_LOOP, SS_ALPHAHELIX, SS_BETABRIDGE, SS_STRAND, SS_HELIX_3, SS_HELIX_5,
          SS_TURN, SS_BEND};

/* This struct tracks information about beta bridges and sheets */
struct Bridge {
    bridge_t type;
    std::deque<int> i, j;
    int chain_i, chain_j;

    Bridge(bridge_t type, int chain_i, int chain_j, int first_i, int first_j):
        type(type),
        chain_i(chain_i),
        chain_j(chain_j)
    {
        i.push_back(first_i);
        j.push_back(first_j);
    };

    bool operator<(const Bridge& b) const {
        return chain_i < b.chain_i || (chain_i == b.chain_i && i.front() < b.i.front());
    }
};


/**
 * Is there an h-bond from donor to acceptor
 */
static bool _test_bond(int donor, int acceptor, const int* hbonds)
{
    return (hbonds[donor*2 + 0] == acceptor) || (hbonds[donor*2+1] == acceptor);
}


/**
 * Test whether two residues are engaged in a beta-bridge
 *
 * Equivalent to MBridgeType MResidue::TestBridge(MResidue* test)
 * from dssp-2.2.0/structure.cpp:687
 */
static bridge_t _residue_test_bridge(int i, int j, int n_residues,
    const int* chain_ids, const int* hbonds)
{
    // printf("testing bridge %d, %d\n", i, j);
    const int a=i-1, b=i, c=i+1;
    const int d=j-1, e=j, f=j+1;
    // printf("e=%d, b=%d: %d\n", e, b, _test_bond(e, b, hbonds));

    if (a >= 0 && c < n_residues && chain_ids[a] == chain_ids[c] &&
        d >= 0 && f < n_residues && chain_ids[d] == chain_ids[f]) {

        if ((_test_bond(c, e, hbonds) && _test_bond(e, a, hbonds)) ||
            (_test_bond(f, b, hbonds) && _test_bond(b, d, hbonds)))
                return BRIDGE_PARALLEL;
        if ((_test_bond(c, d, hbonds) && _test_bond(f, a, hbonds)) ||
            (_test_bond(e, b, hbonds) && _test_bond(b, e, hbonds)))
                return BRIDGE_ANTIPARALLEL;
    }

    return BRIDGE_NONE;
}


/**
 * Identify the beta secondary structure elements. This modifies the
 * vector `secondary` in place.
 *
 * Equivalent to MProtein::CalculateBetaSheets(const vector<MResidue*>& inResidues)
 * from dssp-2.2.0/structure.cpp:1793
 */
static void calculate_beta_sheets(const int* chain_ids, const int* hbonds,
    const std::vector<int>& skip, const int n_residues, std::vector<ss_t>& secondary)
{
    std::vector<Bridge> bridges;

    // Calculate bridges
    for (int i = 1; i < n_residues - 4; i++) {
        for (int j = i+3; j < n_residues - 1; j++) {
            bridge_t type = _residue_test_bridge(j, i, n_residues, chain_ids, hbonds);
            if (type == BRIDGE_NONE || skip[i] || skip[j]) {
                continue;
            }

            // printf("Initial bridge between %d and %d\n", i, j);
            bool found = false;
            for (std::vector<Bridge>::iterator bridge = bridges.begin();
                 bridge != bridges.end(); ++bridge) {

                if (type != bridge->type || i != bridge->i.back() + 1)
                    continue;

                if (type == BRIDGE_PARALLEL && bridge->j.back() + 1 == j) {
                    bridge->i.push_back(i);
                    bridge->j.push_back(j);
                    found = true;
                    break;
                }

                if (type == BRIDGE_ANTIPARALLEL && bridge->j.front() - 1 == j) {
                    bridge->i.push_back(i);
                    bridge->j.push_front(j);
                    found = true;
                    break;
                }
            }
            if (!found) {
	        Bridge bridge(type, chain_ids[i], chain_ids[j], i, j);
                bridges.push_back(bridge);
            }
        }
    }

    // Extend ladders
    sort(bridges.begin(), bridges.end());
    for (int i = 0; i < (int) bridges.size(); ++i) {
        for (int j = i + 1; j < (int) bridges.size(); ++j) {
            int ibi = bridges[i].i.front();
            int iei = bridges[i].i.back();
            int jbi = bridges[i].j.front();
            int jei = bridges[i].j.back();
            int ibj = bridges[j].i.front();
            int iej = bridges[j].i.back();
            int jbj = bridges[j].j.front();
            int jej = bridges[j].j.back();

            if ((bridges[i].type != bridges[j].type) ||
                chain_ids[std::min(ibi, ibj)] != chain_ids[std::max(iei, iej)] ||
                chain_ids[std::min(jbi, jbj)] != chain_ids[std::max(jei, jej)] ||
                ibj - iei >= 6 || (iei >= ibj && ibi <= iej)) {
                    continue;
            }

            bool bulge;
            if (bridges[i].type == BRIDGE_PARALLEL)
                bulge = (jbj > jbi) && ((jbj - jei < 6 && ibj - iei < 3) || (jbj - jei < 3));
            else
                bulge = (jbj < jbi) && ((jbi - jej < 6 && ibj - iei < 3) || (jbi - jej < 3));

            if (bulge) {
                bridges[i].i.insert(bridges[i].i.end(), bridges[j].i.begin(), bridges[j].i.end());
                if (bridges[i].type == BRIDGE_PARALLEL)
                    bridges[i].j.insert(bridges[i].j.end(), bridges[j].j.begin(), bridges[j].j.end());
                else
                    bridges[i].j.insert(bridges[i].j.begin(), bridges[j].j.begin(), bridges[j].j.end());
                bridges.erase(bridges.begin() + j);
                --j;
            }
        }
    }

    for (std::vector<Bridge>::iterator bridge = bridges.begin();
         bridge != bridges.end(); ++bridge) {
         // printf("Bridge from i in (%d, %d)    j in (%d, %d)\n", bridge->i.front(), bridge->i.back(), bridge->j.front(), bridge->j.back());

        ss_t ss = SS_BETABRIDGE;
        if (bridge->i.size() > 1)
            ss = SS_STRAND;

        for (int i = bridge->i.front(); i <= bridge->i.back(); ++i) {
            if (secondary[i] != SS_STRAND)
                secondary[i] = ss;
         }

        for (int i = bridge->j.front(); i <= bridge->j.back(); ++i) {
            if (secondary[i] != SS_STRAND)
                secondary[i] = ss;
        }
    }
}


/**
 * Identify bends in the chain, where the kappa angle (virtual bond angle from
 * c-alpha i-2, to i, to i+2) is greater than 70 degrees
 * dssp-2.2.0/structure.cpp:1729
 */
static std::vector<int> calculate_bends(const float* xyz, const int* ca_indices,
    const int* chain_ids, const int n_residues, std::vector<int>& skip)
{
    std::vector<int> is_bend(n_residues, 0);
    for (int i = 2; i < n_residues-2; i++) {
        if (chain_ids[i-2] == chain_ids[i+2] && !skip[i-2] && !skip[i] && !skip[i+2]) {
            fvec4 prev_ca(xyz[3*ca_indices[i-2]], xyz[3*ca_indices[i-2]+1], xyz[3*ca_indices[i-2]+2], 0);
            fvec4 this_ca(xyz[3*ca_indices[i]], xyz[3*ca_indices[i]+1], xyz[3*ca_indices[i]+2], 0);
            fvec4 next_ca(xyz[3*ca_indices[i+2]], xyz[3*ca_indices[i+2]+1], xyz[3*ca_indices[i+2]+2], 0);
            fvec4 u_prime = prev_ca-this_ca;
            fvec4 v_prime = this_ca-next_ca;
            float cosangle = dot3(u_prime, v_prime)/sqrtf(dot3(u_prime, u_prime)*dot3(v_prime, v_prime));
            float kappa = acosf(CLIP(cosangle, -1, 1));
            is_bend[i] = kappa > (70 * (M_PI / 180.0));
        }
    }
    return is_bend;
}


/**
 * Identify the helical secondary structure elements. This modifies the `secondary`
 * vector inplace.
 *
 * Corresponds to MProtein::CalculateAlphaHelices(const vector<MResidue*>& inResidues, bool inPreferPiHelices)
 * dssp-2.2.0/structure.cpp:1693. Note that `inPreferHelices` is set to true, in this code.
 */
static void calculate_alpha_helices(const float* xyz,
    const int* ca_indices, const int* chain_ids,
    const int* hbonds, std::vector<int>& skip, const int n_atoms, const int n_residues,
    std::vector<ss_t>& secondary)
{
    std::map<int, std::vector<int> > chains;
    for (int i = 0; i < n_residues; i++)
        chains[chain_ids[i]].push_back(i);
    std::vector< std::vector< helix_flag_t> > helix_flags(n_residues, std::vector<helix_flag_t>(6, HELIX_NONE));

    // Helix and Turn
    for (std::map<int, std::vector<int> >::iterator it = chains.begin(); it != chains.end(); it++) {
        std::vector<int> residues = it->second;

        for (int stride = 3; stride <= 5; stride++) {
            for (int ii = 0; ii < (int) residues.size(); ii++) {
                int i = residues[ii];

                if ((i+stride) < n_residues && _test_bond(i+stride, i, hbonds) && (chain_ids[i] == chain_ids[i+stride])) {
                    // printf("%d->%d\n", i+stride, i);

                    helix_flags[i+stride][stride] = HELIX_END;
                    for (int j = i+1; j < i+stride; j++) {
                        if (helix_flags[j][stride] == HELIX_NONE)
                            helix_flags[j][stride] = HELIX_MIDDLE;
                    }

                    if (helix_flags[i][stride] == HELIX_END)
                        helix_flags[i][stride] = HELIX_START_AND_END;
                    else
                        helix_flags[i][stride] = HELIX_START;
                }
            }
        }
    }

    for (int i = 1; i < n_residues-4; i++)
        if ((helix_flags[i][4] == HELIX_START || helix_flags[i][4] == HELIX_START_AND_END) &&
            (helix_flags[i-1][4] == HELIX_START || helix_flags[i-1][4] == HELIX_START_AND_END)) {

            for (int j = i; j <= i + 3; j++)
                secondary[j] = SS_ALPHAHELIX;
        }

    for (int i = 1; i < n_residues - 3; i++)
        if ((helix_flags[i][3] == HELIX_START || helix_flags[i][3] == HELIX_START_AND_END) &&
            (helix_flags[i-1][3] == HELIX_START || helix_flags[i-1][3] == HELIX_START_AND_END)) {

            bool empty = true;
            for (int j = i; empty && j <= i + 2; ++j)
                empty = (secondary[j] == SS_LOOP || secondary[j] == SS_HELIX_3);
            if (empty)
                for (int j = i; j <= i + 2; ++j)
                    secondary[j] = SS_HELIX_3;
        }

    for (int i = 1; i < n_residues - 5; i++)
        if ((helix_flags[i][5] == HELIX_START || helix_flags[i][5] == HELIX_START_AND_END) &&
            (helix_flags[i-1][5] == HELIX_START || helix_flags[i-1][5] == HELIX_START_AND_END)) {

            bool empty = true;
            for (int j = i; empty && j <= i + 4; ++j)
                empty = (secondary[j] == SS_LOOP || secondary[j] == SS_HELIX_5 || secondary[j] == SS_ALPHAHELIX);
            if (empty)
                for (int j = i; j <= i + 4; ++j)
                    secondary[j] = SS_HELIX_5;
        }

    const std::vector<int> is_bend = calculate_bends(xyz, ca_indices, chain_ids, n_residues, skip);
    for (int i = 1; i < n_residues-1; i++)
        if (secondary[i] == SS_LOOP && !skip[i]) {
            bool isTurn = false;
            for (int stride = 3; stride <= 5 && !isTurn; ++stride)
                for (int k = 1; k < stride && !isTurn; ++k)
                    isTurn = (i >= k) && (helix_flags[i-k][stride] == HELIX_START || helix_flags[i-k][stride] == HELIX_START_AND_END);

            if (isTurn)
                secondary[i] = SS_TURN;
            else if (is_bend[i])
                secondary[i] = SS_BEND;
        }
}

/**
 * Calculate DSSP secondary structure assignments
 *
 * Parameters
 * ----------
 * xyz : array, shape=(n_frames, n_atoms, 3)
 *     The cartesian coordinates of all of the atoms in each frame.
 * nco_indices : array, shape=(n_residues, 3)
 *     The indices of the backbone N, C, and O atoms for each residue. The value
 *     should be -1 for any missing residues.
 * ca_indices : array, shape=(n_residues,)
 *     The index of the CA atom of each residue. Should be -1 for any missing
 *     residues
 * is_proline : array, shape=(n_residue,)
 *     Is the residue a proline. These need to be handled slightly differently.
 * chain_ids : array, shape=(n_residue,)
 *     The index of the chain each residue is in. Various parts of this code
 *     require continuity of different secondary structure elements along a chain.
 *
 * Returns
 * -------
 * secondary : array, shape=(n_frames, n_residues)
 *     The DSSP assignment codes for each residue, in each frame. These are
 *     chars, with one of the 8 DSSP codes per residue. Note that the char
 *     array is _not_ null-terminated (at least, this function doesn't add any
 *     null bytes) so you should be careful using it a input to libc string
 *     functions.
 */
void dssp(const float* xyz, const int* nco_indices, const int* ca_indices,
          const int* is_proline, const int* chain_ids, const int n_frames,
          const int n_atoms, const int n_residues, char* secondary)
{
    std::vector<int> skip(n_residues, 0);
    for (int i = 0; i < n_residues; i++)
        if ((nco_indices[i*3] == -1) || (nco_indices[i*3+1] == -1) ||
             (nco_indices[i*3+2] == -1) || ca_indices[i] == -1) {
             skip[i] = 1;
        }

    for (int i = 0; i < n_frames; i++) {
        const float* framexyz = xyz + (i * n_atoms * 3);
        std::vector<int> hbonds(n_residues*2, -1);
        std::vector<float> henergies(n_residues*2, 0);
        // loop is the 'default' secondary structure, which applies
        // when nothing else matches.
        std::vector<ss_t> framesecondary(n_residues, SS_LOOP);

        // call kabsch_sander to calculate the hydrogen bonds
        kabsch_sander(framexyz, nco_indices, ca_indices, is_proline,
                      1, n_atoms, n_residues, &hbonds[0], &henergies[0]);
        // identify the secndary structure elements
        calculate_beta_sheets(chain_ids, &hbonds[0], skip, n_residues, framesecondary);
        calculate_alpha_helices(framexyz, ca_indices, chain_ids,
            &hbonds[0], skip, n_atoms, n_residues, framesecondary);

        // replace the enums with the character codes
        for (int j = 0; j < n_residues; j++) {
            char ss = ' ';
            switch (framesecondary[j]) {
                case SS_ALPHAHELIX:  ss='H'; break;
                case SS_BETABRIDGE:  ss='B'; break;
                case SS_STRAND:      ss='E'; break;
                case SS_HELIX_3:     ss='G'; break;
                case SS_HELIX_5:     ss='I'; break;
                case SS_TURN:        ss='T'; break;
                case SS_BEND:        ss='S'; break;
                case SS_LOOP:        ss=' '; break;
            }
            secondary[i*n_residues+j] = ss;
        }
    }
}
