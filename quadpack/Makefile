SRC = dqpsrt.f dqk15i.f dqelg.f dqagie.f dqagi.f dqags.f dqagse.f dqk21.f dqng.f

OBJS = ${SRC:.f=.o}

.f.o: 
	$(FC) ${FFLAGS} $<

quadpack: $(OBJS)
	ar cr libfquadpack.a $(OBJS)
	ranlib libfquadpack.a

clean:
	rm -f *.o *.a *~
