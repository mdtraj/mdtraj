// Math functions for old versions of visual studio
#ifndef MATH_PATCH_H__
#define MATH_PATCH_H__
#ifdef _MSC_VER
#if _MSC_VER < 1900 // prior to vs2015
//#pragma once 
inline float roundf(float x) {
   return x >= 0.0f ? floorf(x + 0.5f) : ceilf(x - 0.5f);
}

inline double cbrt(double x) {
    if (x<0)
        return -1.0 * pow(-1.0 * x, 1.0/3.0);
    else
        return pow(x, 1.0/3.0);
}
#endif // _MSC_VER < 1900

#ifndef isnan
#define isnan(x) ((x)!=(x))
#endif // isnan

#else // use C99 compliant header.
#include <math.h>
#endif // _MSC_VER
#endif // MATH_PATCH_H__
