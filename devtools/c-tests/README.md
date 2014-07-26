Memory Leakage Tests for c modules
-----------------------------------

test_geometry.c is a small c extecutable that calls the functions
in MDTraj/geometry/src/geometry.c without python. When run under
valgrind, this gives a good way to test for memory leakages and
invalid access.
