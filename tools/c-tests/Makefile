all: valgrind clang-run

test-gcc: test_geometry.c
	gcc -g -O0 test_geometry.c ../../MDTraj/geometry/src/geometry.c ../../MDTraj/geometry/src/sasa.c -I../../MDTraj/geometry/include -msse3 -lm -o test-gcc

valgrind: test-gcc
	valgrind ./test-gcc

test-clang: test_geometry.c
	clang -g -O0 -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls test_geometry.c ../../MDTraj/geometry/src/geometry.c ../../MDTraj/geometry/src/sasa.c -I../../MDTraj/geometry/include -msse3 -lm -o test-clang

clang-run: test-clang
	./test-clang
