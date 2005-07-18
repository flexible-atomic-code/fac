      subroutine njdbug(ibg)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      COMMON /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6

      ibug3 = ibg

      end

      subroutine cpydat(nd, idata, id)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      write(*,*) 'not implemented'
      end

      subroutine njform(im, nb, ij2, ij3, ifree)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      integer ij2(*), ij3(*), ifree(*)

      PARAMETER(MANGM=60,MTRIAD=12,M2TRD=2*MTRIAD,M4TRD=4*MTRIAD)
      PARAMETER(M6J=20,MSUM=10,M3MNGM=3*MANGM,MANGMP=2*(MANGM/3))
      PARAMETER(MFACT=100,MTAB=30)
C
      LOGICAL FAIL,FREE
C
      COMMON/COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),FREE(MANGM)

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

      PARAMETER(MANGM=60,MTRIAD=12,M2TRD=2*MTRIAD,M4TRD=4*MTRIAD)
      PARAMETER(M6J=20,MSUM=10,M3MNGM=3*MANGM,MANGMP=2*(MANGM/3))
      PARAMETER(MFACT=100,MTAB=30)

      LOGICAL FREE
      COMMON/COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),FREE(MANGM)

      k = 1
      do i = 1, m
         j1(i) = ij1(i) + 1
      enddo
      
      call gensum(recup)

      end
      
