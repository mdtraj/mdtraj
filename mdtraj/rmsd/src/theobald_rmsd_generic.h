#include "msvccompat.h"
#include <math.h>
#include <stdio.h>
#include "theobald_rmsd.h"


float msd_axis_major(const int nrealatoms, const int npaddedatoms, const int rowstride,
                     const float* aT, const float* bT, const float G_a, const float G_b)
{

    int k;
    float ax,ay,az,b;
    float M[9];
    const float* aTx = aT;
    const float* aTy = aT+rowstride;
    const float* aTz = aT+2*rowstride;
    const float* bTx = bT;
    const float* bTy = bT+rowstride;
    const float* bTz = bT+2*rowstride;


    if (aT==bT && G_a==G_b)
        return 0.0;

    float xx = 0, xy = 0, xz = 0,  yx = 0,  yy = 0, yz = 0, zx = 0, zy = 0, zz = 0;
    for (k = 0; k < nrealatoms; k++) {
        ax = aTx[k];
        ay = aTy[k];
        az = aTz[k];
        b =  bTx[k];

        xx += (ax * b);
        yx += (ay * b);
        zx += (az * b);

        b = bTy[k];

        xy += (ax * b);
        yy += (ay * b);
        zy += (az * b);

        b = bTz[k];

        xz += (ax * b);
        yz += (ay * b);
        zz += (az * b);
    }

    M[0]=xx;
    M[1]=xy;
    M[2]=xz;
    M[3]=yx;
    M[4]=yy;
    M[5]=yz;
    M[6]=zx;
    M[7]=zy;
    M[8]=zz;

    return msdFromMandG(M, G_a, G_b, nrealatoms, 0, NULL);
}


float msd_atom_major(const int nrealatoms, const int npaddedatoms,
                     const float* a, const float* b, const float G_a, const float G_b,
                     int computeRot, float rot[9])
{

    int k;
    float M[9];

    float ax,ay, az, bx, by, bz;


    if (a==b && G_a==G_b) {
        if (computeRot) {
            rot[0] = rot[4] = rot[8] = 1.0;
            rot[1] = rot[2] = rot[3] = rot[5] = rot[6] = rot[7] = 0.0;
        }
        return 0.0f;
    }

    float xx = 0, xy = 0, xz = 0,  yx = 0,  yy = 0, yz = 0, zx = 0, zy = 0, zz = 0;
    for (k = 0; k < nrealatoms; k++) {

        bx = b[k*3];
        by = b[k*3+1];
        bz = b[k*3+2];

        ax = a[k*3];
        ay = a[k*3+1];
        az = a[k*3+2];

        xx += (bx*ax);
        xy += (by*ax);
        xz += (bz*ax);

        yx += (bx*ay);
        yy += (by*ay);
        yz += (bz*ay);

        zx += (bx*az);
        zy += (by*az);
        zz += (bz*az);

    }


    M[0]=xx;
    M[1]=xy;
    M[2]=xz;
    M[3]=yx;
    M[4]=yy;
    M[5]=yz;
    M[6]=zx;
    M[7]=zy;
    M[8]=zz;

    return msdFromMandG(M, G_a, G_b, nrealatoms, computeRot, rot);
}
