#include "unitcell.h"
#include <cfloat>

UnitCell::UnitCell(const float M[9])
    : v0(M[0], M[1], M[2], 0.0f)
    , v1(M[3], M[4], M[5], 0.0f)
    , v2(M[6], M[7], M[8], 0.0f)
    , det(M[0] * (M[4] * M[8] - M[7] * M[5])
        + M[1] * (M[5] * M[6] - M[8] * M[1])
        + M[2] * (M[3] * M[7] - M[6] * M[4]))
    , vi0((1.0 / det) * fvec4(
        M[4]*M[8] - M[5]*M[7],
      -(M[1]*M[8] - M[7]*M[2]),
        M[1]*M[5] - M[2]*M[4], 0.0f))
    , vi1((1.0 / det) * fvec4(
       -(M[3]*M[8] - M[6]*M[5]),
         M[0]*M[8] - M[6]*M[2],
       -(M[0]*M[5] - M[2]*M[3]), 0.0f))
    , vi2((1.0 / det) * fvec4(
         M[3]*M[7] - M[6]*M[4],
       -(M[0]*M[7] - M[1]*M[6]),
         M[0]*M[4] - M[1]*M[3] , 0.0f))
 { }


fvec4 UnitCell::minimum_image(fvec4 r) {
    fvec4 s = vi0*r[0] + vi1*r[1] + vi2*r[2];
    s = s - round(s);
    r = v0*s[0] + v1*s[1] + v2*s[2];

    fvec4 min_disp;
    float min_d2 = FLT_MAX;

    for (int k0 = -1; k0 < 2; k0++) {
        for (int k1 = -1; k1 < 2; k1++) {
            for (int k2 = -1; k2 < 2; k2++) {
                fvec4 r_prime = r + k0*v0 + k1*v1 + k2*v2;

                float d2 = dot3(r_prime, r_prime);
                if (d2 < min_d2) {
                    min_d2 = d2;
                    min_disp = r_prime;
                }
            }
        }
    }

    return min_disp;
};
