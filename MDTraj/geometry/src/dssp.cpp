#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <deque>
#include <set>
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



struct MBridge
{
	char type;
	int sheet, ladder;
	std::set<MBridge*> link;
	std::deque<int> i, j;
	int chain_i, chain_j;

	bool operator<(const MBridge& b) const {
	    return chain_i < b.chain_i || (chain_i == b.chain_i and i.front() < b.i.front());
    }    
};


static bool _test_bond(int i, int j, const int* hbonds)
{
    return (hbonds[i*2 + 0] == j) || (hbonds[i*2 + 1] == j);
}
    

static char _residue_test_bridge(int i, int j, int n_residues,
    const int* chain_ids, const int* hbonds)
{
    const int a=i-1, b=i, c=i+1;
    const int d=j-1, e=j, f=j+1;

    if (a >= 0 && c < n_residues && chain_ids[a] == chain_ids[c] &&
        d >= 0 && f < n_residues && chain_ids[d] == chain_ids[f]) {

        if ((_test_bond(c, e, hbonds) && _test_bond(a, e, hbonds)) ||
            (_test_bond(f, b, hbonds) && _test_bond(b, d, hbonds)))
                return 'P';
        if ((_test_bond(c, d, hbonds) && _test_bond(f, a, hbonds)) ||
            (_test_bond(e, b, hbonds) && _test_bond(b, e, hbonds)))
                return 'A';
    }
    return 'N';
}


static int calculate_beta_sheets(const float* xyz, const int* nco_indices,
    const int* ca_indices, const int* chain_ids, const int* hbonds,
    const int n_atoms, const int n_residues, char* secondary)
{
    int i, j;
    char type;
    std::vector<MBridge> bridges;
    
    
    // Calculate bridges
    for (i = 0; i < n_residues; i++) {
        for (j = i+3; j < n_residues; j++) {
            type = _residue_test_bridge(i, j, n_residues, chain_ids, hbonds);
            if (type == 'N') {
                continue;
            }
            
            printf("Initial bridge between %d and %d\n", i, j);
            
            bool found = false;
            for (std::vector<MBridge>::iterator bridge = bridges.begin();
                 bridge != bridges.end(); ++bridge) {

    			if (type != bridge->type or i != bridge->i.back() + 1)
    				continue;
                
				if (type == 'P' and bridge->j.back() + 1 == j) {
					bridge->i.push_back(i);
					bridge->j.push_back(j);
					found = true;
					break;
				}
                
				if (type == 'A' and bridge->j.front() - 1 == j) {
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
	for (int i = 0; i < bridges.size(); ++i) {
		for (int j = i + 1; j < bridges.size(); ++j) {
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
                ibj - iei >= 6 || (iei >= ibj and ibi <= iej))
                    continue;
			
			bool bulge;
			if (bridges[i].type == 'P')
				bulge = ((jbj - jei < 6 and ibj - iei < 3) or (jbj - jei < 3));
			else
				bulge = ((jbi - jej < 6 and ibj - iei < 3) or (jbi - jej < 3));

			if (bulge) {
				bridges[i].i.insert(bridges[i].i.end(), bridges[j].i.begin(), bridges[j].i.end());
				if (bridges[i].type == 'P')
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
        
        char ss = 'B';
 		if (bridge->i.size() > 1)
 			ss = 'E';

 		for (int i = bridge->i.front(); i <= bridge->i.back(); ++i) {
            if (secondary[i] != 'E')
                secondary[i] = ss;
     	}
        
 		for (int i = bridge->j.front(); i <= bridge->j.back(); ++i) {
            if (secondary[i] != 'E')
                secondary[i] = ss;
     	}
    }
    
}

static int calculate_alpha_helicies(const float* xyz, const int* nco_indices,
    const int* ca_indices, const int* chain_ids,
    const int* hbonds, const int n_atoms, const int n_residues, char* secondary)
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
            /* compute the arccos of the dot product, and store the result. */
            kappa = (float) acos(CLIP(_mm_cvtss_f32(_mm_dp_ps(u, v, 0x71)), -1, 1));
            printf("%d %f\n", i, kappa * 180 / M_PI);
            is_bend[i] = kappa > (70 * (M_PI / 180.0));
        }
    }
    
    printf("is_bend\n");
    for (int i = 0; i < n_residues; i++)
        printf("%d", is_bend[i]);
    printf("\n");
        
    
}


int dssp(const float* xyz, const int* nco_indices, const int* ca_indices,
         const int* is_proline, const int* chain_ids, const int n_frames,
         const int n_atoms, const int n_residues)         
{
    int i;
    const float *framexyz;
    const int *framehbonds;
    char* framesecondary;

    int* hbonds = (int*) malloc(n_frames*n_residues*2*sizeof(int));
    float* henergies = (float*) calloc(n_frames*n_residues*2, sizeof(float));
    if ((hbonds == NULL) || (henergies == NULL)) {
        fprintf(stderr, "Memory Error\n");
        exit(1);
    }
    for (i = 0; i < n_frames*n_residues*2; i++)
        hbonds[i] = -1;
    
    char* secondary = (char*) malloc(n_frames*n_residues*sizeof(char));
    memset(secondary, 'C', n_frames*n_residues);
    
    // call kabsch_sander to calculate the hydrogen bonds
    kabsch_sander(xyz, nco_indices, ca_indices, is_proline, n_frames, n_atoms,
                  n_residues, hbonds, henergies);
    
    printf("Here I am!\n");

    for (i = 0; i < n_frames; i++) {
        framexyz = xyz + (i * n_atoms * 3);
        framehbonds = hbonds + (i * n_residues * 2);
        framesecondary = secondary + (i * n_residues);
        calculate_beta_sheets(framexyz, nco_indices, ca_indices, chain_ids,
            framehbonds, n_atoms, n_residues, framesecondary);
        calculate_alpha_helicies(framexyz, nco_indices, ca_indices, chain_ids,
            framehbonds, n_atoms, n_residues, framesecondary);
    
    
        for (int j = 0; j < n_residues; j++) {
            printf("%c ", secondary[j]);
        }
        
        // printf("\nTestBond matrix\n");
        // for (int i = 0; i < n_residues; i++) {
        //     for (int j = 0; j < n_residues; j++) {
        //         if (_test_bond(j, i, framehbonds))
        //             printf("%d, %d\n", i, j);
        //         if (_test_bond(i, j, framehbonds))
        //             printf(" %d, %d\n", i, j);
        //     }
        // }
    }

    printf("\ndone\n");
    return 1;
}
#endif