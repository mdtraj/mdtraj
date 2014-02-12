all: valgrind clang-run

test-geom-gcc: test_geometry.c
	gcc -g -O0 test_geometry.c ../../MDTraj/geometry/src/geometry.c ../../MDTraj/geometry/src/sasa.c -I../../MDTraj/geometry/include -msse3 -lm -o test-geom-gcc

test-rmsd-gcc: test_rmsd.c
	gcc -g -O0 test_rmsd.c ../../MDTraj/rmsd/src/theobald_rmsd.c  -I../../MDTraj/rmsd/include/ -lm -o test-rmsd-gcc

valgrind: test-geom-gcc test-rmsd-gcc
	valgrind ./test-geom-gcc --error-exitcode=1
	valgrind ./test-rmsd-gcc --error-exitcode=1

test-clang: test_geometry.c
	clang -g -O0 -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls test_geometry.c ../../MDTraj/geometry/src/geometry.c ../../MDTraj/geometry/src/sasa.c -I../../MDTraj/geometry/include -msse3 -lm -o test-clang

clang-run: test-clang
	./test-clang
