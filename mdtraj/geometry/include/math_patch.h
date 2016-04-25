// Math functions for old versions of visual studio

#ifdef _MSC_VER
#if _MSC_VER < 1900 // prior to vs2015
float roundf(float x) {
   return x >= 0.0f ? floorf(x + 0.5f) : ceilf(x - 0.5f);
}
double cbrt(double x) {
   return pow(x, 1.0/3.0);
}
#endif // _MSC_VER < 1900

#ifndef isnan
#define isnan(x) ((x)!=(x))
#endif // isnan

#endif // _MSC_VER
