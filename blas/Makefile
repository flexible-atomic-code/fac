SRC = lsame.f dasum.f daxpy.f dcopy.f dscal.f ddot.f dgemm.f dgemv.f dger.f dnrm2.f drot.f dspmv.f dspr2.f dswap.f dtrmm.f dtrmv.f idamax.f xerbla.f

OBJS = ${SRC:.f=.o}

.f.o: 
	$(FC) ${FFLAGS} $<

blas: $(OBJS)
	ar cr libblas.a $(OBJS)
	ranlib libblas.a

clean:
	rm -f *.o *.a
