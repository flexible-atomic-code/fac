
SRC = dlacpy.f dlae2.f dlaed0.f dlaed1.f dlaed2.f dlaed3.f dlaed4.f dlaed5.f dlaed6.f dlaed7.f dlaed8.f dlaed9.f dlaeda.f dlaev2.f dlamch.f dlamrg.f dlansp.f dlanst.f dlapy2.f dlarf.f dlarfg.f dlartg.f dlascl.f dlaset.f dlasr.f dlasrt.f dlassq.f dopmtr.f dspevd.f dsptrd.f dstedc.f dsteqr.f dsterf.f ieeeck.f ilaenv.f

OBJS = ${SRC:.f=.o}

.f.o: 
	$(FC) -O -c $<

lapack.a: $(OBJS)
	ar cr lapack.a $(OBJS)
	ranlib lapack.a

clean:
	rm -f *.o *.a

