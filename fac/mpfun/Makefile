SRC = mpfun.f

OBJS = ${SRC:.f=.o}

.f.o:
	${FC} -O -c $<

mpfun: ${OBJS}
	ar cr libmpfun.a ${OBJS}
	ranlib libmpfun.a

clean:
	rm -f *.o *.a
