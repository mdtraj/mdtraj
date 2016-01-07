#include "unitcell.h"
#include <cfloat>

void UnitCell::setBoxMatrix(const float M[9])
{
    isOrthorhombic = (M[1] == 0 && M[2] == 0 && M[3] == 0 &&
                      M[5] == 0 && M[6] == 0 && M[7] == 0);

    if (isOrthorhombic) {
        boxLengths = fvec4(M[0], M[4], M[8], 0.0f);
        invBoxLengths = 1.0f / boxLengths;
    } else {
        v0 = fvec4(M[0], M[1], M[2], 0.0f);
        v1 = fvec4(M[3], M[4], M[5], 0.0f);
        v2 = fvec4(M[6], M[7], M[8], 0.0f);
        det = (M[0] * (M[4] * M[8] - M[7] * M[5])
             + M[1] * (M[5] * M[6] - M[8] * M[1])
             + M[2] * (M[3] * M[7] - M[6] * M[4]));

        // inverse of the box matrix
        vi0 = fvec4((1.0 / det) * fvec4(
          M[4]*M[8] - M[5]*M[7],
        -(M[1]*M[8] - M[7]*M[2]),
          M[1]*M[5] - M[2]*M[4], 0.0f));
        vi1 = fvec4((1.0 / det) * fvec4(
         -(M[3]*M[8] - M[6]*M[5]),
           M[0]*M[8] - M[6]*M[2],
         -(M[0]*M[5] - M[2]*M[3]), 0.0f));
        vi2 = fvec4((1.0 / det) * fvec4(
           M[3]*M[7] - M[6]*M[4],
         -(M[0]*M[7] - M[1]*M[6]),
          M[0]*M[4] - M[1]*M[3] , 0.0f));
    }

}


fvec4 UnitCell::minimumImage(const fvec4& r) {
    if (isOrthorhombic) {
        fvec4 s = invBoxLengths*r;
        s = s - round(s);
        return boxLengths*s;
    } else {
        fvec4 s = vi0*r[0] + vi1*r[1] + vi2*r[2];
        s = s - round(s);
        fvec4 r_prime = v0*s[0] + v1*s[1] + v2*s[2];

        fvec4 min_disp;
        float min_d2 = FLT_MAX;

        for (int k0 = -1; k0 < 2; k0++) {
            for (int k1 = -1; k1 < 2; k1++) {
                for (int k2 = -1; k2 < 2; k2++) {
                    fvec4 r_dprime = r_prime + k0*v0 + k1*v1 + k2*v2;

                    float d2 = dot3(r_dprime, r_dprime);
                    if (d2 < min_d2) {
                        min_d2 = d2;
                        min_disp = r_dprime;
                    }
                }
            }
        }

        return min_disp;
    }
};
