


SHELL = /bin/sh

INSTALL = /usr/bin/install -c

BINDIR = /usr/local/bin

TOPDIR = /Users/mfgu/src/atomic/fac

CC = gcc
CFLAGS = -O2 -fPIC -fno-strength-reduce -no-cpp-precomp -fno-common
CPPFLAGS =  -I/Users/mfgu/src/atomic/fac -I/Users/mfgu/src/atomic/fac/faclib -I${PYINC}
LDFLAGS =  -L/Users/mfgu/src/atomic/fac
LIBS =  -lfac -lm -lcc_dynamic -L/sw/lib/gcc/powerpc-apple-darwin7.5.0/3.4.1 -L/sw/lib/gcc/powerpc-apple-darwin7.5.0/3.4.1/../../.. -lfrtbegin -lg2c -lSystem ${PYLIB}

all: mpy

mpy: 	mpy.c
	${CC} -o mpy ${CPPFLAGS} ${CFLAGS} mpy.c ${LDFLAGS} ${LIBS}

install:
	${INSTALL} -d ${BINDIR}
	${INSTALL} mpy ${BINDIR}

clean:
	rm -rf *.o *~ mpy



