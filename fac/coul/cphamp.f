C     calculates the amplitude and phase of the positive energy 
C     Coulomb functions.

      subroutine cphamp(lambda, eta0, x0, a0, phase, a0p, ierr)
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     lambda: the orbital angular momentum
C     eta0:      -z/k, z negtive for electron ion collision.
C     x0:        k*r
C     a0:        amplitude at r
C     a0p:       derivative of amplitude at r
C     phase:     phase at r
C     ierr:      on output, error code
C                on input, 1--only amplitude
C                          2--only amplitude and phase
C                          3--amplitude, phase, and derivative.
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      implicit none
      integer MAXL, MAXN
      parameter (MAXL = 2500, MAXN = 1)
      
      real*8 lambda, eta0, x0, a0, a0p, phase
      real*8 dfc, dgc, dfcp, dgcp, PI
      parameter (PI = 3.14159265359D0)
      integer ierr, kfn, k, mode
      
      complex*16 eta, x, zlmin
      complex*16 fc(MAXN), fcp(MAXN), gc(MAXN), gcp(MAXN)
      complex*16 sig(MAXN)

      mode = ierr
      if (eta0 .gt. 0.0) then 
         eta = dcmplx(-eta0, 0.0)
         kfn = 0
      else
         kfn = 1
      endif

      x = dcmplx(x0, 0.0)
      if (lambda .lt. MAXL+1) then 
         zlmin = dcmplx(lambda, 0.0)
         k = 1
      else
         k = lambda - MAXL
         zlmin = dcmplx(lambda-k, 0.0)
         k = k+1
         if (k .gt. MAXN) then 
            write(*,*) "Max recusion work array reached in cphamp"
            ierr = -100
            return
         endif
      endif

      ierr = 1
      call coulcc(x, eta, zlmin, k, fc, gc, fcp, gcp, sig, 1, kfn, ierr)

      dfc = dble(fc(k))
      dgc = dble(gc(k))
      a0 = dfc*dfc + dgc*dgc
      if (mode .eq. 1) then
         a0 = sqrt(a0)
         return
      else 
         phase = atan2(dfc, dgc)
         if (mode .eq. 3) then
            dfcp = dble(fcp(k)) - dgc/a0
            dgcp = dble(gcp(k)) + dfc/a0
            a0 = sqrt(a0)
            a0p = sqrt(dfcp*dfcp + dgcp*dgcp)
            if ((dfcp .lt. 0.0 .and. dfc .gt. 0.0) .or. 
     +           (dfcp .gt. 0.0 .and. dfc .lt. 0.0)) then
               a0p = -a0p
            endif
         else
            a0 = sqrt(a0)
         endif
      endif

      end
      

