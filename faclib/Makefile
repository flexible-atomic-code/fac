
SRC = array.c dbase.c config.c cfp.c angular.c rcfp.c recouple.c orbital.c radial.c nucleus.c structure.c transition.c excitation.c recombination.c ionization.c

OBJS = ${SRC:.c=.o}

.c.o: 
	$(CC) ${CFLAGS} $<

faclib.a: ${OBJS}
	ar cr faclib.a ${OBJS}
	ranlib faclib.a

clean: 
	rm -f *.o *.a



