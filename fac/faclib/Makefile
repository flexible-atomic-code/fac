
SRC = pmalloc.c parser.c array.c interpolation.c grid.c coulomb.c config.c cfp.c angular.c rcfp.c recouple.c orbital.c radial.c dbase.c nucleus.c structure.c transition.c excitation.c recombination.c ionization.c init.c rates.c crm.c

OBJS = ${SRC:.c=.o}

.c.o: 
	$(CC) ${CFLAGS} $<

faclib: ${OBJS}
	ar cr libfac.a ${OBJS}
	ranlib libfac.a

clean: 
	rm -f *.o *.a *~



