SRC = coulcc.f y5n.f dlogam.f besljn.f 

OBJS = ${SRC:.f=.o}

.f.o:
	${FC} -O -c $<

coul.a: ${OBJS}
	ar cr coul.a ${OBJS}
	ranlib coul.a

clean:
	rm -f *.o *.a
