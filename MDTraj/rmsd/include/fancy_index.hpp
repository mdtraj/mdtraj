#ifndef __FANCY_INDEX_HPP__
#define __FANCY_INDEX_HPP__

void
fancy_index2d(const float* A, long nx, long ny,
              const long* indx, long nindx, const long* indy, long nindy,
              float* out);

#endif  //__FANCY_INDEX_H__
