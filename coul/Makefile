SRC = besljn.f clgam.f dcoul.f deltac.f fcoul.f scoul.f sum2f0.f

OBJS = ${SRC:.f=.o}

.f.o:
	${FC} -O -c $<

coul.a: ${OBJS}
	ar cr coul.a ${OBJS}
	ranlib coul.a

clean:
	rm -f *.o *.a
