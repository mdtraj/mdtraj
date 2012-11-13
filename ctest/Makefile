all:	
	mkdir -p build
	gcc -I../MDTraj/xtc/include/ -I../MDTraj/xtc/ -c ../MDTraj/xtc/src/xdrfile.c -o build/xdrfile.o
	gcc -I../MDTraj/xtc/include/ -I../MDTraj/xtc/ -c ../MDTraj/xtc/src/xdrfile_xtc.c -o build/xdrfile_xtc.o
	gcc -I../MDTraj/xtc/include/ -I../MDTraj/xtc/ -c test.c -o build/test.o
	gcc build/xdrfile.o build/xdrfile_xtc.o build/test.o

local:
	gcc -DLOCAL -c test.c -o build/test.o
	gcc build/test.o
