#ifndef __FANCY_INDEX_HPP__
#define __FANCY_INDEX_HPP__

void
fancy_index2d(const float* A, int nx, int ny,
              const int* indx, int nindx, const int* indy, int nindy,
              float* out);

#endif  //__FANCY_INDEX_H__
