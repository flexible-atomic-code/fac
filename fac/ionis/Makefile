SRC = intexp.f ionis.f recomb.f recomb_Fe.f

OBJS = ${SRC:.f=.o}

.f.o: 
	$(FC) ${FFLAGS} $<

ionis: $(OBJS)
	ar cr libionis.a $(OBJS)
	ranlib libionis.a

clean:
	rm -f *.o *.a *~
