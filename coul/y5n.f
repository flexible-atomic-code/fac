C calculates the negtive energy coulomb wavefunction, 
C and the derivative.
C the asymptotic behavior is defined by Seaton, MNRAS 118, 504.

      subroutine y5n(lambda, eta0, x0, y5, y5p, norm, ierr)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     lambda:    the orbital angular momentum. 
C     eta0:      -z/k
C     x0:        k*r
C     y5:        wavefunction at r
C     y5p:       d(y5)/d(r) at r
C     norm:      normalization.
C     ierr:      error code
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      implicit none

      real*8 HALFPI
      parameter (HALFPI = 1.5707963268D0)
      integer MAXL, MAXN
      parameter (MAXL = 100, MAXN = 1)

      real*8 lambda, y5, y5p, x0, eta0, norm
      real*8 c, d, b, sb, cb
      complex*16 eta, x, zlmin
      complex*16 fc(MAXN), fcp(MAXN), gc(MAXN), gcp(MAXN)
      complex*16 sig(MAXN)
      integer kfn, ierr, k

      if (eta0 .gt. 0.0) then
         eta = dcmplx(0.0, eta0)
         kfn = 0
      else
         kfn = 1
      endif 

      x = dcmplx(0.0, x0)
      if (lambda .lt. MAXL+1) then
         zlmin = dcmplx(lambda, 0.0)
         k = 1
      else
         k = lambda - MAXL
         zlmin = dcmplx(lambda-k, 0.0)
         k = k+1
         if (k .gt. MAXN) then
            write(*,*) "Max recusion work array reached in y5n"
            ierr = -100
            return
         endif
      endif

      ierr = 1
      call coulcc(x, eta, zlmin, k, fc, gc, fcp, gcp, sig, 
     +            11, kfn, ierr)
      norm = imag(sig(k))
      b = (lambda-eta0)*HALFPI - dble(sig(k))
      sb = sin(b)
      cb = cos(b)
      c = dble(gc(k))
      d = imag(gc(k))
      y5 = c*cb - d*sb
      c = dble(gcp(k))
      d = imag(gcp(k))
      y5p = -(c*sb + d*cb)

      end
