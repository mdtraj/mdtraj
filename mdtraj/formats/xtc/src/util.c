#include "util.h"

// defines PyArray_ENABLEFLAGS
#include "numpy/ndarraytypes.h"

// only if not defined (prior NumPy 1.7), provide impl
#if NPY_API_VERSION <= 0x00000006
#warning "numpy override enableflags"
void PyArray_ENABLEFLAGS(PyArrayObject *arr, int flags)
{
	arr->flags |= flags;
}

#endif
