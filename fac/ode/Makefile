SRC = cfode.f dgbfa.f dgbsl.f dgefa.f dgesl.f ewset.f intdy.f lsode.f prepj.f solsy.f stode.f vnorm.f xerrwv.f

OBJS = ${SRC:.f=.o}

.f.o:
	${FC} ${FFLAGS} $<

ode: ${OBJS}
	ar cr libode.a ${OBJS}
	ranlib libode.a

clean:
	rm -f *.o *.a *~
