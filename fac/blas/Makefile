SRC = lsame.f dasum.f daxpy.f dcopy.f dscal.f ddot.f dgemm.f dgemv.f dger.f dnrm2.f drot.f dspmv.f dspr2.f dswap.f idamax.f xerbla.f

OBJS = ${SRC:.f=.o}

.f.o: 
	$(FC) -O -c $<

blas.a: $(OBJS)
	ar cr blas.a $(OBJS)
	ranlib blas.a

clean:
	rm -f *.o *.a