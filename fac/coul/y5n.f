C calculates the negtive energy coulomb wavefunction, 
C and the derivative.
C the asymptotic behavior is defined by Seaton, MNRAS 118, 504.

      subroutine y5n(lambda, xi, eta0, x0, 
     +     y5, y5i, y5p, y5pi, ierr)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     lambda:    the orbital angular momentum. 
C     eta0:      -z/k, z negative for electron-ion system.
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

      real*8 lambda, xi, y5, y5i, y5p, y5pi, x0, eta0, c0, d0
      complex*16 b, c, d, a1, a2, clogam, logam, norm
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
         zlmin = dcmplx(lambda, xi)
         k = 1
      else
         k = lambda - MAXL
         zlmin = dcmplx(lambda-k, xi)
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
      c0 = imag(sig(k)) - HALFPI*xi
      d0 = (lambda-eta0)*HALFPI - dble(sig(k))
      norm = dcmplx(c0, d0)
      b = logam(2D-16)
      b = dcmplx(lambda, xi)
      a1 = clogam(eta0+b+1.0)
      a2 = clogam(eta0-b)
      b = norm - 0.5*(a1+a2)
      b = exp(b)/eta0
      c = b*gc(k)
      d = b*gcp(k)
      y5 = dble(c)
      y5i = imag(c)
      y5p = -imag(d)
      y5pi = dble(d)

      end
