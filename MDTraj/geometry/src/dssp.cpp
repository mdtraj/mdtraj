/**
 * DSSP secondary structure assignment
 * Copyright [2014] Stanford University and the Authors
 * Authors: Robert T. McGibbon, Maarten L. Hekkelman
<<<<<<< HEAD
 *
=======
 * 
>>>>>>> bc0920cd44886717e7a6b398e1ff0015a32c41f6
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

#define MIN(X, Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X, Y) ((X) > (Y) ? (X) : (Y))
#define CLIP(X, X_min, X_max) (MIN(MAX(X, X_min), X_max))

#ifndef __SSE4_1__
#error
#else
#include <pmmintrin.h>
#include <smmintrin.h>
#include "ssetools.h"
#include "msvccompat.h"
#include <math.h>

enum helix_flag_t {HELIX_NONE, HELIX_START, HELIX_END, HELIX_START_AND_END, HELIX_MIDDLE};
enum bridge_t {BRIDGE_NONE, BRIDGE_PARALLEL, BRIDGE_ANTIPARALLEL};
enum ss_t {SS_LOOP, SS_ALPHAHELIX, SS_BETABRIDGE, SS_STRAND, SS_HELIX_3, SS_HELIX_5,
          SS_TURN, SS_BEND};

struct MBridge {
    bridge_t type;
    int sheet, ladder;
    std::set<MBridge*> link;
    std::deque<int> i, j;
    int chain_i, chain_j;

    bool operator<(const MBridge& b) const {
        return chain_i < b.chain_i || (chain_i == b.chain_i && i.front() < b.i.front());
    }
};


static bool _test_bond(int i, int j, const int* hbonds)
{
    return (hbonds[j*2 + 0] == i) || (hbonds[j*2+1] == i);
}


static bridge_t _residue_test_bridge(int i, int j, int n_residues,
    const int* chain_ids, const int* hbonds)
{
    // printf("testing bridge %d, %d\n", i, j);
    const int a=i-1, b=i, c=i+1;
    const int d=j-1, e=j, f=j+1;
    // printf("e=%d, b=%d: %d\n", e, b, _test_bond(e, b, hbonds));

    if (a >= 0 && c < n_residues && chain_ids[a] == chain_ids[c] &&
        d >= 0 && f < n_residues && chain_ids[d] == chain_ids[f]) {

        if ((_test_bond(e, c, hbonds) && _test_bond(a, e, hbonds)) ||
            (_test_bond(b, f, hbonds) && _test_bond(d, b, hbonds)))
                return BRIDGE_PARALLEL;
        if ((_test_bond(d, c, hbonds) && _test_bond(a, f, hbonds)) ||
            (_test_bond(b, e, hbonds) && _test_bond(e, b, hbonds)))
                return BRIDGE_ANTIPARALLEL;
    }

    return BRIDGE_NONE;
}

static void calculate_beta_sheets(const int* chain_ids, const int* hbonds,
    const std::vector<int>& skip, const int n_residues, std::vector<ss_t>& secondary)
{
    std::vector<MBridge> bridges;

    // Calculate bridges
    for (int i = 1; i < n_residues - 4; i++) {
        for (int j = i+3; j < n_residues - 1; j++) {
            bridge_t type = _residue_test_bridge(j, i, n_residues, chain_ids, hbonds);
            if (type == BRIDGE_NONE || skip[i] || skip[j]) {
                continue;
            }

            // printf("Initial bridge between %d and %d\n", i, j);
            bool found = false;
            for (std::vector<MBridge>::iterator bridge = bridges.begin();
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
                MBridge bridge = {};
                bridge.type = type;
                bridge.i.push_back(i);
                bridge.chain_i = chain_ids[i];
                bridge.j.push_back(j);
                bridge.chain_i = chain_ids[j];
                bridges.push_back(bridge);
            }
        }
    }

    // Extend ladders
    sort(bridges.begin(), bridges.end());
    for (unsigned int i = 0; i < bridges.size(); ++i) {
        for (unsigned int j = i + 1; j < bridges.size(); ++j) {
            unsigned int ibi = bridges[i].i.front();
            unsigned int iei = bridges[i].i.back();
            unsigned int jbi = bridges[i].j.front();
            unsigned int jei = bridges[i].j.back();
            unsigned int ibj = bridges[j].i.front();
            unsigned int iej = bridges[j].i.back();
            unsigned int jbj = bridges[j].j.front();
            unsigned int jej = bridges[j].j.back();

            if ((bridges[i].type != bridges[j].type) ||
                chain_ids[std::min(ibi, ibj)] != chain_ids[std::max(iei, iej)] ||
                chain_ids[std::min(jbi, jbj)] != chain_ids[std::max(jei, jej)] ||
                ibj - iei >= 6 || (iei >= ibj && ibi <= iej)) {
                    continue;
            }

            bool bulge;
            if (bridges[i].type == BRIDGE_PARALLEL)
                bulge = ((jbj - jei < 6 && ibj - iei < 3) || (jbj - jei < 3));
            else
                bulge = ((jbi - jej < 6 && ibj - iei < 3) || (jbi - jej < 3));

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

    for (std::vector<MBridge>::iterator bridge = bridges.begin();
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

static std::vector<int> calculate_bends(const float* xyz, const int* ca_indices,
    const int* chain_ids, const int n_residues, std::vector<int>& skip)
{
    __m128 prev_ca, this_ca, next_ca, u_prime, v_prime, u, v;
    float kappa;
    std::vector<int> is_bend(n_residues, 0);
    for (int i = 2; i < n_residues-2; i++) {
        if (chain_ids[i-2] == chain_ids[i+2] && !skip[i-2] && !skip[i+2]) {
            prev_ca = load_float3(xyz + 3*ca_indices[i-2]);
            this_ca = load_float3(xyz + 3*ca_indices[i]);
            next_ca = load_float3(xyz + 3*ca_indices[i+2]);
            u_prime = _mm_sub_ps(prev_ca, this_ca);
            v_prime = _mm_sub_ps(this_ca, next_ca);
            /* normalize the vectors u_prime and v_prime */
            u = _mm_div_ps(u_prime, _mm_sqrt_ps(_mm_dp_ps(u_prime, u_prime, 0x7F)));
            v = _mm_div_ps(v_prime, _mm_sqrt_ps(_mm_dp_ps(v_prime, v_prime, 0x7F)));
            /* compute the arccos of the dot product, && store the result. */
            kappa = (float) acos(CLIP(_mm_cvtss_f32(_mm_dp_ps(u, v, 0x71)), -1, 1));
            is_bend[i] = kappa > (70 * (M_PI / 180.0));
        }
    }
    return is_bend;
}


static void calculate_alpha_helicies(const float* xyz,
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
            for (unsigned int ii = 0; ii < residues.size(); ii++) {
                int i = residues[ii];

                if ((i+stride) < n_residues && _test_bond(i, i+stride, hbonds) && (chain_ids[i] == chain_ids[i+stride])) {
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

    for (int i = 0; i < n_residues-4; i++)
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
 *     The indices of the backbone N, C, and O atoms for each residue.
 * ca_indices : array, shape=(n_residues,)
 *     The index of the CA atom of each residue.
 * is_proline : array, shape=(n_residue,)
 *     If a particular residue does not contain a CA atom, or you want to skip
 *     the residue for another reason, this value should evaluate to True.
 * chain_ids : array, shape=(n_residue,)
 *     The index of the chain each residue is in
 *
 * Returns
 * -------
 * secondary : array, shape=(n_frames, n_residues)
 *     The DSSP assignment codes for each residue, in each frame.
 */
int dssp(const float* xyz, const int* nco_indices, const int* ca_indices,
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
        std::vector<ss_t> framesecondary(n_residues, SS_LOOP);

        // call kabsch_sander to calculate the hydrogen bonds
        kabsch_sander(framexyz, nco_indices, ca_indices, is_proline,
                      1, n_atoms, n_residues, &hbonds[0], &henergies[0]);

        calculate_beta_sheets(chain_ids, &hbonds[0], skip, n_residues, framesecondary);
        calculate_alpha_helicies(framexyz, ca_indices, chain_ids,
            &hbonds[0], skip, n_atoms, n_residues, framesecondary);

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

    return 1;
}
#endif
