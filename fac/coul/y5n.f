C calculates the negtive energy coulomb wavefunction, 
C and the derivative.
C the asymptotic behavior is defined by Seaton, MNRAS 118, 504.

      subroutine y5n(lambda, eta0, x0, y5, y5p, ierr)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     lambda:    the orbital angular momentum. 
C     eta0:      2.0*z/k
C     x0:        k*r
C     y5:        wavefunction at r
C     y5p:       d(y5)/d(r) at r
C     ierr:      error code
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      implicit none

      real*8 HALFPI
      parameter (HALFPI = 1.5707963268D0)
      complex*16 CI
      parameter (CI = (0.0D0, 1.0D0)) 

      real*8 lambda, y5, y5p, x0, eta0
      complex*16 eta, x, zlmin, norm
      complex*16 fc(1), fcp(1), gc(1), gcp(1)
      complex*16 sig(1)
      integer kfn, ierr

      if (eta0 .gt. 0.0) then
         eta = dcmplx(0.0, eta0)
         kfn = 0
      else
         kfn = 1
      endif 
   
      x = dcmplx(0.0, x0)
      zlmin = dcmplx(lambda, 0.0)

      ierr = 0
      call coulcc(x, eta, zlmin, 1, fc, gc, fcp, gcp, sig, 
     +            11, kfn, ierr)
      norm = CI*((lambda-eta0)*HALFPI - sig(1))
      norm = exp(norm)
      y5 = dble(gc(1)*norm)
      y5p = -imag(gcp(1)*norm)
      end
