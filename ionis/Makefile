SRC = intexp.f ionis.f recomb.f recomb_Fe.f cfit.f colfit.f phfit2.f rrfit.f drfit.f nrrfit.f ndrfit.f

OBJS = ${SRC:.f=.o}

.f.o: 
	$(FC) ${FFLAGS} $<

ionis: $(OBJS)
	ar cr libionis.a $(OBJS)
	ranlib libionis.a

clean:
	rm -f *.o *.a *~
