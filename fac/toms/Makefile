SRC = uvip3p.f tripack.f rgbi3p.f sdbi3p.f tn.f

OBJS = ${SRC:.f=.o}

.f.o:
	${FC} ${FFLAGS} $<

toms: ${OBJS}
	ar cr libtoms.a ${OBJS}
	ranlib libtoms.a

clean:
	rm -f *.o *.a *~
