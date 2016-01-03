#ifndef MDTRAJ_UNITCELL_H
#define MDTRAJ_UNITCELL_H

#include "vectorize_sse.h"

class UnitCell {
public:
    void setBoxMatrix(const float M[9]);
    fvec4 minimumImage(fvec4 r);

private:
    bool isOrthorhombic;
    // for orthorhombic
    fvec4 boxLengths, invBoxLengths;
    // for triclinic
    fvec4 v0, v1, v2;
    double det;
    fvec4 vi0, vi1, vi2;
};

#endif
