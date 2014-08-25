#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <deque>
#include <set>
#include <map>
#include <algorithm>
#include "geometry.h"

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
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
    return (hbonds[i*2 + 0] == j) || (hbonds[i*2 + 1] == j);
}


static bridge_t _residue_test_bridge(int i, int j, int n_residues,
    const int* chain_ids, const int* hbonds)
{
    const int a=i-1, b=i, c=i+1;
    const int d=j-1, e=j, f=j+1;
    // printf("Testing bridge from %d to %d\n", i, j);

    if (a >= 0 && c < n_residues && chain_ids[a] == chain_ids[c] &&
        d >= 0 && f < n_residues && chain_ids[d] == chain_ids[f]) {

        if ((_test_bond(c, e, hbonds) && _test_bond(a, e, hbonds)) ||
            (_test_bond(f, b, hbonds) && _test_bond(b, d, hbonds)))
                return BRIDGE_PARALLEL;

        // printf("cd: %d, fa:%d, eb:%d, be:%d\n", _test_bond(c, d, hbonds), _test_bond(f, a, hbonds), _test_bond(e, b, hbonds), _test_bond(b, e, hbonds));
        if ((_test_bond(d, c, hbonds) && _test_bond(a, f, hbonds)) ||
            (_test_bond(e, b, hbonds) && _test_bond(b, e, hbonds)))
                return BRIDGE_ANTIPARALLEL;
    }

    return BRIDGE_NONE;
}


static void calculate_beta_sheets(const float* xyz, const int* nco_indices,
    const int* ca_indices, const int* chain_ids, const int* hbonds,
    const int n_atoms, const int n_residues, std::vector<ss_t>& secondary)
{
    std::vector<MBridge> bridges;

    // Calculate bridges
    for (int i = 1; i < n_residues - 4; i++) {
        for (int j = i+3; j < n_residues - 1; j++) {
            bridge_t type = _residue_test_bridge(j, i, n_residues, chain_ids, hbonds);
            if (type == BRIDGE_NONE) {
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
                ibj - iei >= 6 || (iei >= ibj && ibi <= iej))
                    continue;

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
    const int* chain_ids, const int n_residues)
{
    __m128 prev_ca, this_ca, next_ca, u_prime, v_prime, u, v;
    float kappa;
    std::vector<int> is_bend(n_residues, 0);
    for (int i = 2; i < n_residues-2; i++) {
        if (chain_ids[i-2] == chain_ids[i+2]) {
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

static void calculate_alpha_helicies(const float* xyz, const int* nco_indices,
    const int* ca_indices, const int* chain_ids,
    const int* hbonds, const int n_atoms, const int n_residues,
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

    // printf("HelixFlags\n");
    // for (int stride = 3; stride <= 5; stride++) {
    //     printf("%d: ", stride);
    //     for (int i = 0; i < n_residues; i++)
    //         printf("%d", helix_flags[i][stride]);
    //     printf("\n");
    // }


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

    const std::vector<int> is_bend = calculate_bends(xyz, ca_indices, chain_ids, n_residues);
    for (int i = 1; i < n_residues-1; i++)
        if (secondary[i] == SS_LOOP) {
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



int dssp(const float* xyz, const int* nco_indices, const int* ca_indices,
         const int* is_proline, const int* chain_ids, const int n_frames,
         const int n_atoms, const int n_residues, char* secondary)
{
    const float *framexyz;
    const int *framehbonds;

    int* hbonds = (int*) malloc(n_frames*n_residues*2*sizeof(int));
    float* henergies = (float*) calloc(n_frames*n_residues*2, sizeof(float));
    if ((hbonds == NULL) || (henergies == NULL)) {
        fprintf(stderr, "Memory Error\n");
        exit(1);
    }
    for (int i = 0; i < n_frames*n_residues*2; i++)
        hbonds[i] = -1;

    // call kabsch_sander to calculate the hydrogen bonds
    kabsch_sander(xyz, nco_indices, ca_indices, is_proline, n_frames, n_atoms,
                  n_residues, hbonds, henergies);

    for (int i = 0; i < n_frames; i++) {
        framexyz = xyz + (i * n_atoms * 3);
        framehbonds = hbonds + (i * n_residues * 2);
        std::vector<ss_t> framesecondary(n_residues, SS_LOOP);

        calculate_beta_sheets(framexyz, nco_indices, ca_indices, chain_ids,
            framehbonds, n_atoms, n_residues, framesecondary);
        calculate_alpha_helicies(framexyz, nco_indices, ca_indices, chain_ids,
            framehbonds, n_atoms, n_residues, framesecondary);

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
