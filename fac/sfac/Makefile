SRC = stoken.c sfac.c scrm.c

OBJS = ${SRC:.c=.o}

none:	 

sfac: 	none
	$(CC) ${CFLAGS} $(SRC)
	$(CC) -o sfac sfac.o stoken.o ${LDFLAGS} ${LIBS}
	$(CC) -o scrm scrm.o stoken.o ${LDFLAGS} ${LIBS}

clean:
	rm -rf *.o *~
	rm -rf sfac
	rm -rf scrm











