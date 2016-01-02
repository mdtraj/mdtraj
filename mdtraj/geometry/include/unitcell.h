#ifndef MDTRAJ_UNITCELL_H
#define MDTRAJ_UNITCELL_H

#include "vectorize_sse.h"

class UnitCell {
public:
    UnitCell(const float M[9]);
    fvec4 minimum_image(fvec4 r);
    const fvec4 v0, v1, v2;
    const double det;
    const fvec4 vi0, vi1, vi2;
};

#endif
