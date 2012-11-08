# Minimal Makefile for building DCD and NAMDBIN file format shared
# libraries.

# Toni Giorgino at gmail.com - (c) Universitat Pompeu Fabra 2011


# The code distributed with this Makefile is bundled with VMD and
# licensed under the UIUC open source license (see LICENSE).  The
# plugins are being redistributed (along with their headers) to avoid
# depending on all of VMD's code base and build process. If you have
# trouble with this Makefile, download and build the latest VMD
# distribution.

# To integrate in MSMBUILDER, only dcdplugin_s.c and include files
# are required. They will be built by the setup.py script.


default: dcdplugin_s.so namdbinplugin_s.so

# Defaults for LINUXAMD64
CFLAGS=-I include -fPIC -Wall -O2 -m64
#CFLAGS=-I include -fPIC -Wall -g -m64


# xxx_s.c is an ugly hack to work-around the register-callback
# machinery: makes all static functions /directly/ callable from
# python. Normally they would be exposed through a the
# molfile_plugin_t structure containing function pointers. Using
# callbacks from python is cumbersome and requires to maintain a
# python version of the (non-trivial) struct.  THIS WILL BREAK ANY
# STATIC VARIABLE IN FUNCTION BODIES! On the other hand, there should
# be no function name conflict, since ctypes requires the caller to
# explicitly state the DLLs of each function to be called.


%_s.c: %.c
	sed s+static+/*static*/+ < $< > $@

%.o: %.c include/*.h *.h
	$(CC) -c $(CFLAGS)  $< -o $@

%.so: %.o
	$(CC) $(CFLAGS) -shared  -lm $< -o $@

clean:
	-rm *.so *.o


dcdplugin_test: dcdplugin.c
	gcc -lm -Wall -I include -D TEST_DCDPLUGIN -g $< -o $@

