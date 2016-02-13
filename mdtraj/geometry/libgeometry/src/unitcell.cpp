#include "unitcell.h"
#include <cfloat>
#include <cmath>
#include <cstdio>

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

        //  Make sure they're in reduced form (from OpenMM)
        v2 = v2 - v1*round(v2[1]/v1[1]);
        v2 = v2 - v0*round(v2[0]/v0[0]);
        v1 = v1 - v0*round(v1[0]/v0[0]);

        float M2[9];
        store3(v0, M2);
        store3(v1, M2+3);
        store3(v2, M2+6);

        det = (M2[0] * (M2[4] * M2[8] - M2[7] * M2[5])
             + M2[1] * (M2[5] * M2[6] - M2[8] * M2[1])
             + M2[2] * (M2[3] * M2[7] - M2[6] * M2[4]));

        // inverse of the box matrix
        vi0 = fvec4((1.0 / det) * fvec4(
          M2[4]*M2[8] - M2[5]*M2[7],
        -(M2[1]*M2[8] - M2[7]*M2[2]),
          M2[1]*M2[5] - M2[2]*M2[4], 0.0f));
        vi1 = fvec4((1.0 / det) * fvec4(
         -(M2[3]*M2[8] - M2[6]*M2[5]),
           M2[0]*M2[8] - M2[6]*M2[2],
         -(M2[0]*M2[5] - M2[2]*M2[3]), 0.0f));
        vi2 = fvec4((1.0 / det) * fvec4(
           M2[3]*M2[7] - M2[6]*M2[4],
         -(M2[0]*M2[7] - M2[1]*M2[6]),
          M2[0]*M2[4] - M2[1]*M2[3] , 0.0f));
    }

}


fvec4 UnitCell::minimumImage(const fvec4& x1, const fvec4& x2) {
// fvec4 UnitCell::minimumImage(const fvec4& r) {
    if (isOrthorhombic) {
        fvec4 s = invBoxLengths * (x1-x2);
        s = s - round(s);
        // the unused last element, r.w might be NaN or something,
        // which can break things downsteam, so we clear it.
        return clearw(boxLengths*s);
    } else {
    //     fvec4 f1 = vi0*x1[0] + vi1*x1[1] + vi2*x1[2];
    //     fvec4 f2 = vi0*x2[0] + vi1*x2[1] + vi2*x2[2];
    //     f1 = f1 - floor(f1);
    //     f2 = f2 - floor(f2);
    //
    //     fvec4 s = f1 - f2;
    //     fvec4 r2 = v0*s[0] + v1*s[1] + v2*s[2];
    //
    //     // x1 = v0*f1[0] + v1*f1[1] + v2*f1[2];
    //     // x2 = v0*f2[0] + v1*f2[1] + v2*f2[2];
    //     //
    //     //
    //     // fvec4 r = x2 - x1;
    //     r2 = r2 - v2*floor(r2[2]/v2[2]+0.5);
    //     r2 = r2 - v1*floor(r2[1]/v1[1]+0.5);
    //     r2 = r2 - v0*floor(r2[0]/v0[0]+0.5);
    //     return r2;
    // } if (0) {
        printf("v0: [%f, %f, %f]\n", v0[0], v0[1], v0[2]);
        printf("v1: [%f, %f, %f]\n", v1[0], v1[1], v1[2]);
        printf("v2: [%f, %f, %f]\n", v2[0], v2[1], v2[2]);



        fvec4 f1 = vi0*x1[0] + vi1*x1[1] + vi2*x1[2];
        fvec4 f2 = vi0*x2[0] + vi1*x2[1] + vi2*x2[2];

        f1 = f1 - floor(f1);
        f2 = f2 - floor(f2);
        fvec4 s = f1 - f2;

        fvec4 r_prime = v0*s[0] + v1*s[1] + v2*s[2];

        fvec4 min_disp;
        float min_d2 = FLT_MAX;

        const int minq = -2;
        const int maxq = 3;

        for (int k0 = minq; k0 < maxq; k0++) {
            for (int k1 = minq; k1 < maxq; k1++) {
                for (int k2 = minq; k2 < maxq; k2++) {
                    fvec4 r_dprime = r_prime + k0*v0 + k1*v1 + k2*v2;

                    float d2 = dot3(r_dprime, r_dprime);
                    if (d2 < min_d2) {
                        min_d2 = d2;
                        min_disp = r_dprime;
                        printf("%d %d %d\n", k0, k1, k2);
                    }
                }
            }
        }

        return min_disp;
    }
};
