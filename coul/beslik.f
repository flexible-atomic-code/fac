C calculates the modified bessel functions K with pure 
C imaginary order, and the derivatives, result is exponentially scaled

      subroutine beslik(lambda, x0, dk, dkp, ierr)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     lambda:    the order. 
C     x0:        argument.
C     di, dip:   I_lambda and derivative.
C     dk, dkp:   K_lambda and derivative.
C     ierr:      error code
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      implicit none

      integer MAXL, MAXN
      parameter (MAXL = 100, MAXN = 1)

      real*8 lambda, x0, dk, dkp
      complex*16 x, eta, zlmin
      complex*16 fc(MAXN), fcp(MAXN), gc(MAXN), gcp(MAXN)
      complex*16 sig(MAXN)
      integer kfn, ierr, k, mode

      kfn = 3
      x = dcmplx(x0, 0.0)
      zlmin = dcmplx(0.0, lambda)
      k = 1

      mode = -ierr
      ierr = 1
      call coulcc(x, eta, zlmin, k, fc, gc, fcp, gcp, sig, 
     +            mode, kfn, ierr)
      dk = dble(gc(k))
      if (abs(mode) .eq. 1) dkp = dble(gcp(k))

      end
