


SHELL = /bin/sh

INSTALL = /usr/bin/install -c

BINDIR = /usr/local/bin

TOPDIR = /Users/mfgu/src/fac

CC = gcc
CFLAGS = -O2 -fPIC -fno-strength-reduce -no-cpp-precomp -fno-common
CPPFLAGS =  -I/Users/mfgu/src/fac -I/Users/mfgu/src/fac/faclib -I${PYINC}
LDFLAGS =  -L/Users/mfgu/src/fac
LIBS =  -lfac -lm -L/sw/lib/gcc4/lib/gcc/i386-apple-darwin8/4.2.0 -L/sw/lib/gcc4/lib/gcc/i386-apple-darwin8/4.2.0/../../.. -lgfortranbegin -lgfortran -lgcc_s.10.4 -lSystem ${PYLIB}

all: mpy

mpy: 	mpy.c
	${CC} -o mpy ${CPPFLAGS} ${CFLAGS} mpy.c ${LDFLAGS} ${LIBS}

install:
	${INSTALL} -d ${BINDIR}
	${INSTALL} mpy ${BINDIR}

clean:
	rm -rf *.o *~ mpy



