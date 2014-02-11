all: run

a.out: test_geometry.c
	gcc -g -O0 test_geometry.c ../../MDTraj/geometry/src/geometry.c ../../MDTraj/geometry/src/sasa.c -I../../MDTraj/geometry/include -msse3 -lm

run: a.out
	valgrind ./a.out
