SRC = uvip3p.f tripack.f rgbi3p.f sdbi3p.f

OBJS = ${SRC:.f=.o}

.f.o:
	${FC} -O -c $<

toms: ${OBJS}
	ar cr libtoms.a ${OBJS}
	ranlib libtoms.a

clean:
	rm -f *.o *.a
