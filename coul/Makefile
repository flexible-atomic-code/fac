SRC = clogam.f f11.f coulcc.f y5n.f cphamp.f dlogam.f argam.f besljn.f beslik.f coulrad.f

OBJS = ${SRC:.f=.o}

.f.o:
	${FC} ${FFLAGS} $<

coul: ${OBJS}
	ar cr libcoul.a ${OBJS}
	ranlib libcoul.a

clean:
	rm -f *.o *.a
