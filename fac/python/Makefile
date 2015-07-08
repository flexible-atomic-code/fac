


SHELL = /bin/sh

INSTALL = /sw/bin/ginstall -c

BINDIR = /usr/local/bin

TOPDIR = /disks/d4a/home/mfgu/src/fac/fac

CC = gcc
CFLAGS = -O2 -fPIC -fno-strength-reduce -fno-common
CPPFLAGS =  -I/disks/d4a/home/mfgu/src/fac/fac -I/disks/d4a/home/mfgu/src/fac/fac/faclib -I${PYINC}
LDFLAGS =  -L/disks/d4a/home/mfgu/src/fac/fac
LIBS =  -lfac -lm -L/sw/lib/gcc4.3/lib/gcc/i686-apple-darwin8/4.3.1 -L/sw/lib/gcc4.3/lib/gcc/i686-apple-darwin8/4.3.1/../../.. -lgfortranbegin -lgfortran ${PYLIB}

all: mpy

mpy: 	mpy.c
	${CC} -o mpy ${CPPFLAGS} ${CFLAGS} mpy.c ${LDFLAGS} ${LIBS}

install:
	${INSTALL} -d ${BINDIR}
	${INSTALL} mpy ${BINDIR}

clean:
	rm -rf *.o *~ mpy



