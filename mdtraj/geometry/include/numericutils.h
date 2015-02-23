#ifndef _NUMERICUTILS_H_
#define _NUMERICUTILS_H_
#ifdef __cplusplus
extern "C" {
#endif

int histogram(const float* data, const int n_data, const float* bins,
              const int n_bins, const float min_bin, const float max_bin,
              int* out);

#ifdef __cplusplus
}
#endif
#endif
