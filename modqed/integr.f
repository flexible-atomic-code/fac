c       =================================================
        function tint(k,p,q,a,b,r,v)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This function computes the matrix elements of r^k
c       tint= <p,q|r^k|a,b>
c       Numerical integration is performed by Simpson's rule
c
c       Input: 
c       p(i) - the upper radial component of the left wave function
c       q(i) - the lower radial component of the left wave function
c       a(i) - the upper radial component of the right wave function
c       b(i) - the lower radial component of the right wave function
c       r(i) - the radial grid
c       v(i) - the weights of the Simpson formula for nonuniform grid
c
c       Common parameters:
c       h    - grid step
c       ii   - Number of grid points
c
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /h/h/ii/ii
        real*8 p(maxii),q(maxii),a(maxii),b(maxii)
        real*8 r(maxii),v(maxii)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        i0=1
        r0=r(i0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        imax=ii
        do i=ii,1,-1
          d=dabs(p(i))+dabs(q(i))+dabs(a(i))+dabs(b(i))
          if (d.gt.1.d-20) goto 200
          imax=i
       enddo
       if (imax .lt. 3) then
          tint = 0d0
          return
       endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Integration from 0 to r1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
200     f1=(p(1)*a(1)+q(1)*b(1))*r(1)**k
        f5=(p(5)*a(5)+q(5)*b(5))*r(5)**k
c        if (f1 .gt. 1d-20) then
c           g=dlog(f5/f1)/dlog(r(5)/r(1))
c           t0=r(1)*f1/(g+1)
c        else
           t0=0d0
c        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        ds=0.d0
        i1=1
        i2=imax
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       If the number of points is even we use four-point
c       Simpson's rule for first 4 points. The number
c       of the rest points is odd
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (2*(imax/2).eq.imax) then
          i=i1
          f1=(p(i)*a(i)+q(i)*b(i))*v(i)*r(i)**k
          i=i1+1
          f2=(p(i)*a(i)+q(i)*b(i))*v(i)*r(i)**k
          i=i1+2
          f3=(p(i)*a(i)+q(i)*b(i))*v(i)*r(i)**k
          i=i1+3
          f4=(p(i)*a(i)+q(i)*b(i))*v(i)*r(i)**k
          ds=ds+3.d0/8.d0*h*(f1+3*f2+3*f3+f4)
          i1=i1+3
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Simpson's rule for the odd number of points.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dt1=0.d0
        do i=i1,i2,2
          fi=(p(i)*a(i)+q(i)*b(i))*v(i)*r(i)**k
          dt1=dt1+fi
        enddo
        f1=(p(i1)*a(i1)+q(i1)*b(i1))*v(i1)*r(i1)**k
        f2=(p(i2)*a(i2)+q(i2)*b(i2))*v(i2)*r(i2)**k
        dt1=2*dt1-f1-f2
c
        dt2=0.d0
        do i=i1+1,i2-1,2
          fi=(p(i)*a(i)+q(i)*b(i))*v(i)*r(i)**k
          dt2=dt2+fi
        enddo
        dt2=4*dt2
        dt=h*(dt1+dt2)/3.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        tint=ds+dt+t0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        function sint(c,p,q,a,b,r,v)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This function computes the matrix element
c       sint=<p,q|c(r)/r|a,b>.
c       Numerical integration is performed by Simpson's rule
c
c       Input: 
c       p(i) - the upper radial component of the left wave function
c       q(i) - the lower radial component of the left wave function
c       a(i) - the upper radial component of the right wave function
c       b(i) - the lower radial component of the right wave function
c       r(i) - the radial grid
c       v(i) - the weights of the Simpson formula for nonuniform grid
c
c       Common parameters:
c       h    - grid step
c       ii   - Number of grid points
c
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /h/h/ii/ii
        real*8 c(maxii),p(maxii),q(maxii),a(maxii),b(maxii)
        real*8 r(maxii),v(maxii)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        i0=1
        r0=r(i0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        imax=ii
        do i=ii,1,-1
          d=dabs(p(i))+dabs(q(i))+dabs(a(i))+dabs(b(i))+dabs(c(i))
          if (d.gt.1.d-20) goto 200
          imax=i-1
       enddo
       if (imax .lt. 3) then
          sint = 0d0
          return
       endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Integration from 0 to r1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
200     f1=c(1)*(p(1)*a(1)+q(1)*b(1))/r(1)
        f5=c(5)*(p(5)*a(5)+q(5)*b(5))/r(5)
c        if (f1 .gt. 1d-20) then
c           g=dlog(f5/f1)/dlog(r(5)/r(1))
c           t0=r(1)*f1/(g+1)
c        else
           t0=0d0
c        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        sint=0.d0
        i1=1
        i2=imax
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       If the number of points is even we use four-point
c       Simpson's rule for first 4 points. The number
c       of the rest points is odd
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (2*(imax/2).eq.imax) then
          i=i1
          f1=(p(i)*a(i)+q(i)*b(i))*c(i)*v(i)/r(i)
          i=i1+1
          f2=(p(i)*a(i)+q(i)*b(i))*c(i)*v(i)/r(i)
          i=i1+2
          f3=(p(i)*a(i)+q(i)*b(i))*c(i)*v(i)/r(i)
          i=i1+3
          f4=(p(i)*a(i)+q(i)*b(i))*c(i)*v(i)/r(i)
          sint=sint+3.d0/8.d0*h*(f1+3*f2+3*f3+f4)
          i1=i1+3
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Simpson's rule for the odd number of points.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dt1=0.d0
        do i=i1,i2,2
          fi=(p(i)*a(i)+q(i)*b(i))*c(i)*v(i)/r(i)
          dt1=dt1+fi
        enddo
        f1=(p(i1)*a(i1)+q(i1)*b(i1))*c(i1)*v(i1)/r(i1)
        f2=(p(i2)*a(i2)+q(i2)*b(i2))*c(i2)*v(i2)/r(i2)
        dt1=2*dt1-f1-f2
c
        dt2=0.d0
        do i=i1+1,i2-1,2
          fi=(p(i)*a(i)+q(i)*b(i))*c(i)*v(i)/r(i)
          dt2=dt2+fi
        enddo
        dt2=4*dt2
        dt=h*(dt1+dt2)/3.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        sint=sint+dt+t0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
