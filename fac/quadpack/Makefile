SRC = dqpsrt.f dqk15i.f dqelg.f dqagie.f dqagi.f

OBJS = ${SRC:.f=.o}

.f.o: 
	$(FC) ${FFLAGS} $<

quadpack: $(OBJS)
	ar cr libquadpack.a $(OBJS)
	ranlib libquadpack.a

clean:
	rm -f *.o *.a *~
