all: valgrind clang-run

test-geom-gcc: test_geometry.c
	gcc -g -O0 test_geometry.c ../../MDTraj/geometry/src/geometry.c ../../mdtraj/geometry/src/sasa.c -I../../mdtraj/geometry/include -I../../mdtraj/geometry/src/kernels -mssse3 -lm -o test-geom-gcc

test-rmsd-gcc: test_rmsd.c
	gcc -g -O0 test_rmsd.c ../../mdtraj/rmsd/src/theobald_rmsd.c  -I../../mdtraj/rmsd/include/ -lm -mssse3 -lm -o test-rmsd-gcc

valgrind: test-geom-gcc test-rmsd-gcc
	valgrind ./test-geom-gcc --error-exitcode=1
	valgrind ./test-rmsd-gcc --error-exitcode=1

test-clang: test_geometry.c
	clang -g -O0 -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls test_geometry.c ../../mdtraj/geometry/src/geometry.c ../../mdtraj/geometry/src/sasa.c -I../../mdtraj/geometry/include -mssse3 -lm -o test-clang

clang-run: test-clang
	./test-clang
