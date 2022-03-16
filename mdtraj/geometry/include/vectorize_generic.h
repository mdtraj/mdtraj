#ifndef VECTORIZE_GENERIC_H_
#define VECTORIZE_GENERIC_H_

#include <algorithm>
#include <cstdlib>
#include <cmath>


class ivec4;

class fvec4 {
public:
    float val[4];

    fvec4() = default;
    fvec4(float v) : val {v,v,v,v} {}
    fvec4(float v1, float v2, float v3, float v4) : val {v1, v2, v3, v4} {}
    fvec4(const float* v) : val {v[0], v[1], v[2], v[3]} {}
    operator float() const
    {
        return *val;
    }

    fvec4(const float* table, const int idx[4])
        : fvec4(table[idx[0]], table[idx[1]], table[idx[2]], table[idx[3]]) { }

    float operator[](int i) const
    {
        switch (i) {
        case 0:
            return val[0];
        case 1:
            return val[1];
        case 2:
            return val[2];
        case 3:
            return val[3];
        }
        return 0.0f;
    }

    void store(float* v) const
    {
        int i;
        for(i = 0; i<4; i++) {
            v[i]=val[i];
        }
    }

    void storeVec3(float* v) const
    {
        int i;
        for(i = 0; i<3; i++) {
            v[i] = val[i];
        }
    }

    fvec4 operator+(fvec4 other) const
    {
        float new_fvec4[4];
        int i;
        for(i = 0; i<4; i++) {
            new_fvec4[i] = val[i] + other[i];

        }
        return fvec4(new_fvec4);

    }

    fvec4 operator-(fvec4 other) const
    {
        float new_fvec4[4];
        int i;
        for(i = 0; i<4; i++) {
            new_fvec4[i] = val[i] - other[i];
        }
        return fvec4(new_fvec4);
    }

    fvec4 operator*(fvec4 other) const
    {
        float new_fvec4[4];
        int i;
        for(i = 0; i<4; i++) {
            new_fvec4[i] = val[i] * other[i];
        }
        return fvec4(new_fvec4);
    }

    fvec4 operator/(fvec4 other) const
    {
        float new_fvec4[4];
        int i;
        for(i = 0; i<4; i++) {
            new_fvec4[i] = val[i] / other[i];
        }
        return fvec4(new_fvec4);
    }

    void operator+=(fvec4 other)
    {
        int i;
        for(i = 0; i<4; i++) {
            val[i] += other[i];
        }
    }

    void operator-=(fvec4 other)
    {
        int i;
        for(i = 0; i<4; i++) {
            val[i] = val[i]-other[i];
        }
    }

    void operator*=(fvec4 other)
    {
        int i;
        for(i = 0; i<4; i++) {
            val[i] *= other[i];
        }
    }

    void operator/=(fvec4 other)
    {
        int i;
        for(i = 0; i<4; i++) {
            val[i] = val[i]/other[i];
        }
    }

    fvec4 operator-() const
    {
        float new_fvec4[4];
        int i;
        for(i = 0; i<4; i++) {
            new_fvec4[i] = 0-val[i];
        }
        return fvec4(new_fvec4);
    }

    fvec4 operator&(fvec4 other) const
    {
        float new_fvec4[4];
        int i;
        for(i = 0; i<4; i++) {
            new_fvec4[i] = (float)((unsigned int)val[i] && (unsigned int)other[i]);
        }
        return fvec4(new_fvec4);
    }

    fvec4 operator|(fvec4 other) const
    {
        float new_fvec4[4];
        int i;
        for(i = 0; i<4; i++) {
            new_fvec4[i] = (float)((unsigned int)val[i] || (unsigned int)other[i]);
        }
        return fvec4(new_fvec4);
    }

    ivec4 operator==(fvec4 other) const;
    ivec4 operator!=(fvec4 other) const;
    ivec4 operator>(fvec4 other) const;
    ivec4 operator<(fvec4 other) const;
    ivec4 operator>=(fvec4 other) const;
    ivec4 operator<=(fvec4 other) const;
    operator ivec4() const;

    static ivec4 expandBitsToMask(int bitmask);

};


class ivec4 {
public:

    int32_t val[4];

    ivec4() {}
    ivec4(int v) : val {v,v,v,v} {}
    ivec4(int v1, int v2, int v3, int v4) : val {v1, v2, v3, v4} {}
    ivec4(const int* v) : val {v[0], v[1], v[2], v[3]} {}
    operator int() const
    {
        return *val;
    }

    int operator[](int i) const
    {
        switch (i) {
        case 0:
            return val[0];
        case 1:
            return val[1];
        case 2:
            return val[2];
        case 3:
            return val[3];
        }
        return 0;
    }

    void store(int* v) const
    {
        int i;
        for(i = 0; i<4; i++) {
            v[i] = val[i];
        }
    }

    ivec4 operator+(ivec4 other) const
    {
        int new_ivec4[4];
        int i;
        for(i = 0; i<4; i++) {
            new_ivec4[i] = val[i] + other[i];
        }
        return new_ivec4;
    }

    ivec4 operator-(ivec4 other) const
    {
        int new_ivec4[4];
        int i;
        for(i = 0; i<4; i++) {
            new_ivec4[i] = val[i] - other[i];
        }
        return ivec4(new_ivec4);
    }

    ivec4 operator*(ivec4 other) const
    {
        int new_ivec4[4];
        int i;
        for(i = 0; i<4; i++) {
            new_ivec4[i] = val[i] * other[i];
        }
        return ivec4(new_ivec4);
    }

    void operator+=(ivec4 other)
    {
        int i;
        for(i = 0; i<4; i++) {
            val[i] += other[i];
        }
    }

    void operator-=(ivec4 other)
    {
        int i;
        for(i = 0; i<4; i++) {
            val[i] -= other[i];
        }
    }

    void operator*=(ivec4 other)
    {
        int i;
        for(i = 0; i<4; i++) {
            val[i] *= other[i];
        }
    }

    ivec4 operator-() const
    {
        int new_ivec4[4];
        int i;
        for(i = 0; i<4; i++) {
            new_ivec4[i] = val[i] *(-1);
        }
        return ivec4(new_ivec4);
    }

    ivec4 operator&(ivec4 other) const
    {
        int new_ivec4[4];
        int i;
        for(i = 0; i<4; i++) {
            new_ivec4[i] = (int) (val[i] & other[i]);
        }
        return ivec4(new_ivec4);

    }

    ivec4 operator|(ivec4 other) const
    {
        int new_ivec4[4];
        int i;
        for(i = 0; i<4; i++) {
            new_ivec4[i] = (int) (val[i] | other[i]);
        }
        return ivec4(new_ivec4);
    }

    ivec4 operator==(ivec4 other) const
    {
        int new_ivec4[4];
        int i;
        for(i = 0; i<4; i++) {
            new_ivec4[i] = (int) (val[i] == other[i]);
        }
        return ivec4(new_ivec4);

    }

    ivec4 operator!=(ivec4 other) const
    {
        int new_ivec4[4];
        int i;
        for(i = 0; i<4; i++) {
            new_ivec4[i] = (int) (val[i] != other[i]);
        }
        return ivec4(new_ivec4);
    }

    ivec4 operator>(ivec4 other) const
    {
        int new_ivec4[4];
        int i;
        for(i = 0; i<4; i++) {
            new_ivec4[i] = (int) (val[i] > other[i]);
        }
        return ivec4(new_ivec4);
    }

    ivec4 operator<(ivec4 other) const
    {
        int new_ivec4[4];
        int i;
        for(i = 0; i<4; i++) {
            new_ivec4[i] = (int) (val[i] < other[i]);
        }
        return ivec4(new_ivec4);
    }

    ivec4 operator>=(ivec4 other) const
    {
        int new_ivec4[4];
        int i;
        for(i = 0; i<4; i++) {
            new_ivec4[i] = (int) (val[i] >= other[i]);
        }
        return ivec4(new_ivec4);
    }

    ivec4 operator<=(ivec4 other) const
    {
        int new_ivec4[4];
        int i;
        for(i = 0; i<4; i++) {
            new_ivec4[i] = (int) (val[i] <= other[i]);
        }
        return ivec4(new_ivec4);
    }

    operator fvec4() const;
};

inline fvec4::operator ivec4() const
{
    int new_ivec4[4];
    int i;
    for(i = 0; i<4; i++) {
        new_ivec4[i] = (int) val[i];
    }
    return ivec4(new_ivec4);
}

inline ivec4::operator fvec4() const
{
    float new_fvec4[4];
    int i;
    for(i = 0; i<4; i++) {
        new_fvec4[i] = (float) val[i];
    }
    return fvec4(new_fvec4);
}

inline ivec4 fvec4::expandBitsToMask(int bitmask)
{
    return ivec4(bitmask & 1 ? -1 : 0,
                 bitmask & 2 ? -1 : 0,
                 bitmask & 4 ? -1 : 0,
                 bitmask & 8 ? -1 : 0);
}

inline ivec4 fvec4::operator==(fvec4 other) const
{
    int new_ivec4[4];
    int i;
    for(i = 0; i<4; i++) {
        new_ivec4[i] = (int) (val[i] == other[i]);
    }
    return ivec4(new_ivec4);

}

inline ivec4 fvec4::operator!=(fvec4 other) const
{
    int new_ivec4[4];
    int i;
    for(i = 0; i<4; i++) {
        new_ivec4[i] = (int) (val[i] != other[i]);
    }
    return ivec4(new_ivec4);
}

inline ivec4 fvec4::operator>(fvec4 other) const
{
    int new_ivec4[4];
    int i;
    for(i = 0; i<4; i++) {
        new_ivec4[i] = (int) (val[i] > other[i]);
    }
    return ivec4(new_ivec4);

}

inline ivec4 fvec4::operator<(fvec4 other) const
{
    int new_ivec4[4];
    int i;
    for(i = 0; i<4; i++) {
        new_ivec4[i] = (int) (val[i] < (int)other[i]);
    }
    return ivec4(new_ivec4);
}

inline ivec4 fvec4::operator>=(fvec4 other) const
{
    int new_ivec4[4];
    int i;
    for(i = 0; i<4; i++) {
        new_ivec4[i] = (int) (val[i] >= (int)other[i]);
    }
    return ivec4(new_ivec4);
}

inline ivec4 fvec4::operator<=(fvec4 other) const
{
    int new_ivec4[4];
    int i;
    for(i = 0; i<4; i++) {
        new_ivec4[i] = (int) (val[i] <= (int)other[i]);
    }
    return ivec4(new_ivec4);
}

static inline fvec4 min(fvec4 v1, fvec4 v2)
{
    float new_fvec4[4];
    int i;
    for(i = 0; i<4; i++) {
        if(v1[i]<v2[i]) {
            new_fvec4[i] = v1[i];
        }
        else {
            new_fvec4[i] = v2[i];
        }
    }
    return fvec4(new_fvec4);
}

static inline fvec4 max(fvec4 v1, fvec4 v2)
{
    float new_fvec4[4];
    int i;
    for(i = 0; i<4; i++) {
        if(v1[i]>v2[i]) {
            new_fvec4[i] = v1[i];
        }
        else {
            new_fvec4[i] = v2[i];
        }
    }
    return fvec4(new_fvec4);
}

static inline fvec4 abs(fvec4 v)
{
    float new_fvec4[4];
    int i;
    for(i = 0; i<4; i++) {
        new_fvec4[i] = fabs(v[i]);
    }
    return fvec4(new_fvec4);
}

static inline fvec4 rsqrt(const fvec4& v)
{
    // Initial estimate of rsqrt().

    fvec4 y(1/sqrt(v));

    // Perform an iteration of Newton refinement.

    fvec4 x2 = 0.5f*v;
    y *= fvec4(1.5f)-x2*y*y;
    return y;
}


static inline fvec4 sqrt(fvec4 v)
{
    float new_fvec4[4];
    int i;
    for(i = 0; i<4; i++) {
        new_fvec4[i] = sqrt(v[i]);
    }
    return fvec4(new_fvec4);
}

static inline fvec4 exp(fvec4 v)
{
    float new_fvec4[4];
    int i;
    for(i = 0; i<4; i++) {
        new_fvec4[i] = exp(v[i]);
    }
    return fvec4(new_fvec4);
}

static inline fvec4 log(fvec4 v)
{
    float new_fvec4[4];
    int i;
    for(i = 0; i<4; i++) {
        new_fvec4[i] = log(v[i]);
    }
    return fvec4(new_fvec4);
}

static inline float dot3(fvec4 v1, fvec4 v2)
{
    float result[4];
    int i;
    for(i = 0; i<4; i++) {
        result[i] = v1[i]*v2[i];
    }
    return result[0] + result[1] + result[2];
}

static inline float dot4(fvec4 v1, fvec4 v2)
{
    float result[4];
    int i;
    for(i = 0; i<4; i++) {
        result[i] = v1[i]*v2[i];
    }
    return result[0] + result[1] + result[2] + result[3];
}

static inline fvec4 cross(fvec4 v1, fvec4 v2)
{
    return fvec4(v1[1]*v2[2] - v1[2]*v2[1],
                 v1[2]*v2[0] - v1[0]*v2[2],
                 v1[0]*v2[1] - v1[1]*v2[0], 0);
}

static inline float reduceAdd(fvec4 v)
{
    return dot4(v, fvec4(1.0f));
}

static inline void transpose(fvec4& v1, fvec4& v2, fvec4& v3, fvec4& v4)
{
    float t1[4] = {v1[0], v1[1], v1[2], v1[3]};
    float t2[4] = {v2[0], v2[1], v2[2], v2[3]};
    float t3[4] = {v3[0], v3[1], v3[2], v3[3]};
    float t4[4] = {v4[0], v4[1], v4[2], v4[3]};
    v1.val[0]=t1[0];
    v1.val[1]=t2[0];
    v1.val[2]=t3[0];
    v1.val[3]=t4[0];
    v2.val[0]=t1[1];
    v2.val[1]=t2[1];
    v2.val[2]=t3[1];
    v2.val[3]=t4[3];
    v3.val[0]=t1[2];
    v3.val[1]=t2[2];
    v2.val[2]=t3[2];
    v3.val[3]=t4[3];
    v4.val[0]=t1[3];
    v3.val[1]=t2[3];
    v2.val[2]=t3[3];
    v4.val[3]=t4[3];
}


static inline ivec4 min(ivec4 v1, ivec4 v2)
{
    int new_ivec4[4];
    int i;
    for(i = 0; i<4; i++) {
        if(v1[i]<v2[i]) new_ivec4[i] = v1[i];
        else new_ivec4[i] = v2[i];
    }
    return ivec4(new_ivec4);
}

static inline ivec4 max(ivec4 v1, ivec4 v2)
{
    int new_ivec4[4];
    int i;
    for(i = 0; i<4; i++) {
        if(v1[i]>v2[i]) new_ivec4[i] = v1[i];
        else new_ivec4[i] = v2[i];
    }
    return ivec4(new_ivec4);
}

static inline ivec4 abs(ivec4 v)
{
    int new_ivec4[4];
    int i;
    for(i = 0; i<4; i++) {
        new_ivec4[i] = abs(v[i]);
    }
    return ivec4(new_ivec4);
}

static inline bool any(ivec4 v)
{
    return (v[0]!= 0 || v[1] != 0 || v[2] != 0 || v[3] != 0);
}

static inline fvec4 operator+(float v1, fvec4 v2)
{
    float new_fvec4[4];
    int i;
    for(i = 0; i<4; i++) {
        new_fvec4[i] = v1+v2[i];
    }
    return fvec4(new_fvec4);
}

static inline fvec4 operator-(float v1, fvec4 v2)
{
    float new_fvec4[4];
    int i;
    for(i = 0; i<4; i++) {
        new_fvec4[i] = v1-v2[i];
    }
    return fvec4(new_fvec4);
}

static inline fvec4 operator-(fvec4 v2,float v1)
{
    float new_fvec4[4];
    int i;
    for(i = 0; i<4; i++) {
        new_fvec4[i] = v2[i]-v1;
    }
    return fvec4(new_fvec4);
}

static inline fvec4 operator*(fvec4 v2, float v1)
{
    float new_fvec4[4];
    int i;
    for(i = 0; i<4; i++) {
        new_fvec4[i] = v1*v2[i];
    }
    return fvec4(new_fvec4);
}

static inline fvec4 operator*(fvec4 v2, int v1)
{
    float new_fvec4[4];
    int i;
    for(i = 0; i<4; i++) {
        new_fvec4[i] = v1*v2[i];
    }
    return fvec4(new_fvec4);
}


static inline fvec4 operator/(float v1, fvec4 v2)
{
    float new_fvec4[4];
    int i;
    for(i = 0; i<4; i++) {
        new_fvec4[i] = v1/v2[i];
    }
    return fvec4(new_fvec4);
}


static inline fvec4 operator/(fvec4 v2, float v1)
{
    float new_fvec4[4];
    int i;
    for(i = 0; i<4; i++) {
        new_fvec4[i] = v2[i]/v1;
    }
    return fvec4(new_fvec4);
}

static inline fvec4 round(fvec4 v)
{
    float f[4];
    int i;
    for(i = 0; i<4; i++) {
        f[i] = round(v[i]);
    }
    return fvec4(f);
}

static inline fvec4 floor(fvec4 v)
{
    float f[4];
    int i;
    for(i = 0; i<4; i++) {
        f[i] = floor(v[i]);
    }
    return fvec4(f);
}

static inline void gatherVecPair(const float* table, ivec4 index, fvec4& out0, fvec4& out1)
{
    fvec4 t0(table + index[0]);
    fvec4 t1(table + index[1]);
    fvec4 t2(table + index[2]);
    fvec4 t3(table + index[3]);
    transpose(t0, t1, t2, t3);
    out0 = t0;
    out1 = t1;
}

static inline fvec4 reduceToVec3(fvec4 x, fvec4 y, fvec4 z)
{
    const auto nx = reduceAdd(x);
    const auto ny = reduceAdd(y);
    const auto nz = reduceAdd(z);
    return fvec4(nx, ny, nz, 0.0);
}

#endif