      subroutine njdbug(ibg)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      COMMON /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6

      ibug3 = ibg

      end

      subroutine cpydat(nd, idata, id)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      integer nd, idata(*), id
      
      PARAMETER (MZLR1=   5)
      PARAMETER (MZNR1=   5)
C
      PARAMETER (MXORB=MZNR1*MZNR1/2+MZLR1)
      PARAMETER (KFL1=100,KFL2=MXORB+4)
      PARAMETER (KFL2A=2*KFL2,KFL2B=4*KFL2,KFL2C=2*KFL2+2)
      PARAMETER (KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20,KFLS=12,
     A          KFLN=10,KFLV=40)

      LOGICAL FREE
      COMMON /COUPLE/M,N,J1(KFL1),J2(KFL2,3),J3(KFL2,3),FREE(KFL1)
C
      LOGICAL SUMVAR
      COMMON /ARGU/J6C,J7C,J8C,J9C,JWC,J6(KFL6),J7(KFL7),J8(KFL8),
     A       J9(KFL9),JW(6,KFLW),JDEL,LDEL(KFLW,2),SUMVAR(KFL1),MP
C
      COMMON /SUMARG/J6P(KFLV),J7P(KFLV),J8P(KFLV),J9P(KFLV),
     A       JWORD(6,KFLW),NLSUM,NBJ(KFLN),NB6J(KFLN),K6CP(KFLN),
     B       K7CP(KFLN),K8CP(KFLN),K9CP(KFLN),JSUM6(KFLS),
     C       JSUM4(KFLS,KFLW),JSUM5(KFLS,KFLW),INV6J(KFLW)
      
      if (id .eq. 1) then
         k = 1      
         j6c = idata(k)
         k = k + 1
         j7c = idata(k)
         k = k + 1
         j8c = idata(k)
         k = k + 1
         j9c = idata(k)
         k = k + 1
         jwc = idata(k)
         k = k + 1      
         do i = 1, kfl6
            j6(i) = idata(k)
            k = k + 1
         enddo
         do i = 1, kfl7
            j7(i) = idata(k)
            k = k + 1
         enddo
         do i = 1, kfl8
            j8(i) = idata(k)
            k = k + 1
         enddo
         do i = 1, kfl9
            j9(i) = idata(k)
            k = k + 1
         enddo
         do i = 1, kflw
            do j = 1, 6
               jw(j,i) = idata(k)
               k = k + 1
            enddo
         enddo
         jdel = idata(k)
         k = k + 1
         do i = 1, 2
            do j = 1, kflw
               ldel(j,i) = idata(k)
               k = k + 1
            enddo
         enddo
         do i = 1, kflv
            j6p(i) = idata(k)
            k = k + 1
         enddo
         do i = 1, kflv
            j7p(i) = idata(k)
            k = k + 1
         enddo
         do i = 1, kflv
            j8p(i) = idata(k)
            k = k + 1
         enddo
         do i = 1, kflv
            j9p(i) = idata(k)
            k = k + 1
         enddo
         do i = 1, kflw
            do j = 1, 6
               jword(j,i) = idata(k)
               k = k + 1
            enddo
         enddo
         nlsum = idata(k)
         k = k + 1
         do i = 1, kfln
            nbj(i) = idata(k)
            k = k + 1
         enddo
         do i = 1, kfln
            nb6j(i) = idata(k)
            k = k + 1
         enddo
         do i = 1, kfln
            k6cp(i) = idata(k)
            k = k + 1
         enddo
         do i = 1, kfln
            k7cp(i) = idata(k)
            k = k + 1
         enddo
         do i = 1, kfln
            k8cp(i) = idata(k)
            k = k + 1
         enddo
         do i = 1, kfln
            k9cp(i) = idata(k)
            k = k + 1
         enddo
         do i = 1, kfls
            jsum6(i) = idata(k)
            k = k + 1
         enddo
         do i = 1, kflw
            do j = 1, kfls
               jsum4(j,i) = idata(k)
               k = k + 1
            enddo
         enddo
         do i = 1, kflw
            do j = 1, kfls
               jsum5(j,i) = idata(k)
               k = k + 1
            enddo
         enddo
         do i = 1, kflw
            inv6j(i) = idata(k)
            k = k + 1
         enddo
      else
         k = 1      
         idata(k) = j6c
         k = k + 1
         idata(k) = j7c
         k = k + 1
         idata(k) = j8c
         k = k + 1
         idata(k) = j9c
         k = k + 1
         idata(k) = jwc
         k = k + 1      
         do i = 1, kfl6
            idata(k) = j6(i)
            k = k + 1
         enddo
         do i = 1, kfl7
            idata(k) = j7(i)
            k = k + 1
         enddo
         do i = 1, kfl8
            idata(k) = j8(i)
            k = k + 1
         enddo
         do i = 1, kfl9
            idata(k) = j9(i)
            k = k + 1
         enddo
         do i = 1, kflw
            do j = 1, 6
               idata(k) = jw(j,i)
               k = k + 1
            enddo
         enddo
         idata(k) = jdel
         k = k + 1
         do i = 1, 2
            do j = 1, kflw
               idata(k) = ldel(j,i)
               k = k + 1
            enddo
         enddo
         do i = 1, kflv
            idata(k) = j6p(i)
            k = k + 1
         enddo
         do i = 1, kflv
            idata(k) = j7p(i)
            k = k + 1
         enddo
         do i = 1, kflv
            idata(k) = j8p(i)
            k = k + 1
         enddo
         do i = 1, kflv
            idata(k) = j9p(i)
            k = k + 1
         enddo
         do i = 1, kflw
            do j = 1, 6
               idata(k) = jword(j,i)
               k = k + 1
            enddo
         enddo
         idata(k) = nlsum
         k = k + 1
         do i = 1, kfln
            idata(k) = nbj(i)
            k = k + 1
         enddo
         do i = 1, kfln
            idata(k) = nb6j(i)
            k = k + 1
         enddo
         do i = 1, kfln
            idata(k) = k6cp(i)
            k = k + 1
         enddo
         do i = 1, kfln
            idata(k) = k7cp(i)
            k = k + 1
         enddo
         do i = 1, kfln
            idata(k) = k8cp(i)
            k = k + 1
         enddo
         do i = 1, kfln
            idata(k) = k9cp(i)
            k = k + 1
         enddo
         do i = 1, kfls
            idata(k) = jsum6(i)
            k = k + 1
         enddo
         do i = 1, kflw
            do j = 1, kfls
               idata(k) = jsum4(j,i)
               k = k + 1
            enddo
         enddo
         do i = 1, kflw
            do j = 1, kfls
               idata(k) = jsum5(j,i)
               k = k + 1
            enddo
         enddo
         do i = 1, kflw
            idata(k) = inv6j(i)
            k = k + 1
         enddo         
      endif
      if (k .ge. nd) then
         write(*,*) "not enough space in idata array for cpydat"
      endif

      end

      subroutine njform(im, nb, ij2, ij3, ifree)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      integer ij2(*), ij3(*), ifree(*)

      logical fail
      
      PARAMETER (MZLR1=   5)
      PARAMETER (MZNR1=   5)
C
      PARAMETER (MXORB=MZNR1*MZNR1/2+MZLR1)
      PARAMETER (KFL1=100,KFL2=MXORB+4)
      PARAMETER (KFL2A=2*KFL2,KFL2B=4*KFL2,KFL2C=2*KFL2+2)
      PARAMETER (KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20,KFLS=12,
     A          KFLN=10,KFLV=40)

      LOGICAL FREE
      COMMON /COUPLE/M,N,J1(KFL1),J2(KFL2,3),J3(KFL2,3),FREE(KFL1)

      m = im
      n = nb
      k = 1
      do i = 1, im
         if (ifree(i) .eq. 0) then
            free(i) = .FALSE.
         else
            free(i) = .TRUE.
         endif
         j1(i) = 1
      enddo
      do i = 1, nb-1
         do j = 1, 3
            j2(i,j) = ij2(k)
            j3(i,j) = ij3(k)
            k = k + 1
         enddo
      enddo

      call njgraf(recup, fail)

      end

      subroutine njsum(ij1, recup)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      integer ij1(*)

      PARAMETER (MZLR1=   5)
      PARAMETER (MZNR1=   5)
C
      PARAMETER (MXORB=MZNR1*MZNR1/2+MZLR1)
      PARAMETER (KFL1=100,KFL2=MXORB+4)
      PARAMETER (KFL2A=2*KFL2,KFL2B=4*KFL2,KFL2C=2*KFL2+2)
      PARAMETER (KFL6=120,KFL7=150,KFL8=120,KFL9=40,KFLW=20,KFLS=12,
     A          KFLN=10,KFLV=40)

      LOGICAL FREE
      COMMON /COUPLE/M,N,J1(KFL1),J2(KFL2,3),J3(KFL2,3),FREE(KFL1)
C
      LOGICAL SUMVAR
      COMMON /ARGU/J6C,J7C,J8C,J9C,JWC,J6(KFL6),J7(KFL7),J8(KFL8),
     A       J9(KFL9),JW(6,KFLW),JDEL,LDEL(KFLW,2),SUMVAR(KFL1),MP
C
      COMMON /SUMARG/J6P(KFLV),J7P(KFLV),J8P(KFLV),J9P(KFLV),
     A       JWORD(6,KFLW),NLSUM,NBJ(KFLN),NB6J(KFLN),K6CP(KFLN),
     B       K7CP(KFLN),K8CP(KFLN),K9CP(KFLN),JSUM6(KFLS),
     C       JSUM4(KFLS,KFLW),JSUM5(KFLS,KFLW),INV6J(KFLW)

      k = 1
      do i = 1, m
         j1(i) = ij1(i) + 1
      enddo

      CALL GENSUM(J6C,J7C,J8C,J9C,JWC,J6,J7,J8,J9,JW,JDEL,LDEL,
     A     MP,J6P,J7P,J8P,J9P,JWORD,NLSUM,NBJ,NB6J,K6CP,K7CP,
     B     K8CP,K9CP,JSUM4,JSUM5,JSUM6,INV6J,RECUP)

      end
      
