SRC = mpfun.f

OBJS = ${SRC:.f=.o}

.f.o:
	${FC} ${FFLAGS} $<

mpfun: ${OBJS}
	ar cr libfmpfun.a ${OBJS}
	ranlib libfmpfun.a

clean:
	rm -f *.o *.a *~
