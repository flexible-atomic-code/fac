SRC = dpmpar.f enorm.f lmder.f lmder1.f lmpar.f qrfac.f qrsolv.f chkder.f

OBJS =  ${SRC:.f=.o}

.f.o: 
	$(FC) ${FFLAGS} $<

minpack: $(OBJS)
	ar cr libminpack.a $(OBJS)
	ranlib libminpack.a

clean:
	rm -f *.o *.a *~
