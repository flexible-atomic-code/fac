      subroutine njdbug(ibg)
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6,IBUG7,IBUG8,IBUG9
!$OMP THREADPRIVATE(/DEBUG/)

      ibug3 = ibg

      end

      subroutine cpydat(nd, idata, id)

      write(*,*) 'not implemented'
      end

      subroutine njform(im, nb, ij2, ij3, ifree)  
      implicit double precision (a-h, o-z)    
      integer ij2(*), ij3(*), ifree(*)
      DIMENSION K6(40),K7(80),K8(40),KW(6,20),J2TEST(12),J3TEST(12)
      COMMON/COUPLE/M,N,J1(40),J2(12,3),J3(12,3)
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6,IBUG7,IBUG8,IBUG9
      COMMON/DIMEN/KFL1,KFL2,KFL3,KFL4,KFL5,KFL6,KFL7
      COMMON/INFORM/IREAD,IWRITE,IPUNCH
      COMMON/JCONST/J6C,J7C,J8C,JWC,ICOUNT,K6,K7,K8,KW,J2TEST,J3TEST
!$OMP THREADPRIVATE(/COUPLE/,/DEBUG/,/DIMEN/,/INFORM/,/JCONST/)      
      iread = 5
      iwrite = 6
      m = im
      n = nb
      k = 1
      do i = 1, nb-1
         do j = 1, 3
            j2(i,j) = ij2(k)
            j3(i,j) = ij3(k)
            k = k + 1
         enddo
      enddo

      call njsym(j6c,j7c,j8c,jwc,k6,k7,k8,kw,c,icount,j2test,j3test)

      end

      subroutine njsum(ij1, r)
      implicit double precision (a-h, o-z)
      integer ij1(*)

      DIMENSION K6(40),K7(80),K8(40),KW(6,20),J2TEST(12),J3TEST(12)
      COMMON/COUPLE/M,N,J1(40),J2(12,3),J3(12,3)
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6,IBUG7,IBUG8,IBUG9
      COMMON/DIMEN/KFL1,KFL2,KFL3,KFL4,KFL5,KFL6,KFL7
      COMMON/INFORM/IREAD,IWRITE,IPUNCH
      COMMON/JCONST/J6C,J7C,J8C,JWC,ICOUNT,K6,K7,K8,KW,J2TEST,J3TEST
!$OMP THREADPRIVATE(/COUPLE/,/DEBUG/,/DIMEN/,/INFORM/,/JCONST/)
      k = 1
      do i = 1, m
         j1(i) = ij1(i) + 1
      enddo
      
      call gensum(j6c,j7c,j8c,jwc,k6,k7,k8,kw,r,icount,j2test,j3test)
      end
      
