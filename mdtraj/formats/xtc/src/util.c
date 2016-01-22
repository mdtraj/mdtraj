#include "util.h"
void PyArray_ENABLEFLAGS(PyArrayObject *arr, int flags)
{
    ((PyArrayObject_fields *)arr)->flags |= flags;
}
