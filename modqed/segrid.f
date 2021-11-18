c       =================================================
        subroutine segrid
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine initializes the semi-logarithmic
c       radial grid r(i), which is used for the calculation
c       of the wave functions and the numerical integrations.
c
c       The nonuniform grid r is obtained from the uniform grid rho
c       by using the change of variables
c       rho(r)=al*r+bt*ln(r),
c       where al and bt are the grid parameters.
c
c       rho(r) is the uninform grid,
c       rho_i=rho_1+h*(i-1), where h is the grid spacing. 
c
c       ii   .... number of grid points (ii=maxii-20), where
c                 maxii is given in the include file
c                 'qedmod.inc'
c       r(i) .... non-unform radial grid (i=1,ii)
c       r1   .... the first grid point closest to the origin
c       r2   .... the last grid point (the size of the box), defined
c                 in the atom_data subroutine r2=r(ii).
c       h    .... Semi-logarithmic grid step                
c       z    .... Nuclear charge
c       rnucl ... Radius of the nuclear sphere 
c       knucl  ..... Nuclear model (0 - point, 1 - Fermi model,
c                    2 - uniform sphere)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /z/z
        common /ii/ii/r1/r1/r2/r2/h/h/al/al/bt/bt
        common /r/r(maxii)/v/v(maxii)/rho/rho(maxii)
        common /knucl/knucl/inucl/inucl/rnucl/rnucl
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        r1=dexp(-6.d0)*(256.d0/maxii)**2/5.d0
        if (z.gt.0.99) then
          r1=r1/z
        else
          r1=r1/100.d0
        endif
        if (rnucl.gt.0.and.r1.gt.0.10d0*rnucl) r1=0.10d0*rnucl
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        rmax=dabs(r2)
        imax=maxii-20
        imax=((imax+1)/2)*2-1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Semi-logarithmic grid: ro=al*r+bt*ln(r)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        ii=imax
        bt=1.d0
        h=bt*dlog(r2/r1)/((ii-1)*0.67)
        al=(h*(imax-1)-bt*dlog(r2/r1))/(r2-r1)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (al.ge.0.d0) goto 210
        ii=(bt*dlog(r2/r1))/h+1
c        write( *,99) imax,ii,h,al,bt,r1,r2
c99      format(/'  ii > Imax'/'  Imax=',i6/'  ii  =',i6/
c     1  '  h   =',f9.4,/'  al  =',f9.3/,'  bt  =',f9.3/
c     2  '  r1  =',e10.3,'  r2  =',f9.3)
        call exit1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Correction for the Sphere nuclear model
c       r(inucl)=rnucl
c       - - - - - - - - - - - - - - - - - - - - - - - - -
210     inucl=0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (knucl.eq.2) then
          ro1=al*r1+bt*dlog(r1)
          ron=al*rnucl+bt*dlog(rnucl)
          d=ron-ro1
          inucl=(d/h+0.0001d0)+1
          if (inucl.lt.11) then
            inucl=11
            dr=d-(inucl-1)*h
            ro1=ro1+dr
            r0=dexp(ro1)/bt
            r1=rr(ro1,r0)
          endif
          dr=(rnucl-r1)/(r2-r1)
          h=(bt*dlog(rnucl/r1)-bt*dlog(r2/r1)*dr)/
     &    ((inucl-1)-(ii-1)*dr)
          al=((inucl-1)*h-bt*dlog(rnucl/r1))/(rnucl-r1)
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,maxii
          r(i)=0.d0
          v(i)=0.d0
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call tbr(ii,r1,h,al,bt,r,v)
        do i=1,ii
           rho(i) = al*r(i) + bt*dlog(r(i))
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c        write( *,5) ii,r1,r2
c        write(11,5) ii,r1,r2
c5       format (/2x,'Semi-logarithmic radial grid:',
c     1  /2x,'npoints    =',i6,/2x,'rmin    =',e13.4,
c     2  /2x,'rmax    =',e13.4)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
1000    return
        end

        subroutine rgmqed(a, b)
        implicit real*8 (a-h,o-z)
        common /ii/ii/r1/r1/r2/r2/h/h/al/al/bt/bt

        a = al
        b = bt
        return
        end
      
c       =================================================
        subroutine tbr(ii,r1,h,al,bt,r,v)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Initialization of the Semi-logarithmic radial grid
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dimension r(*),v(*)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        r(1)=r1
        d=r1
        t=al*d+bt
        v(1)=d/t
        p=al*d+bt*dlog(d)
        do 10 i=2,ii
        p=p+h
200     t=al*d+bt*dlog(d)
        t=(p-t)/(al*d+bt)
        d=d*(1.d0+t)
        if (dabs(t).gt.0.5d-11) goto 200
        t=al*d+bt
        v(i)=d/t
        r(i)=d
10      continue
        return
        end
c       =================================================
        function rr(ro,r0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine solve the equation ro=al*r+bt*ln(r)
c       and define the value of r for the given ro.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /al/al/bt/bt
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        r=r0
200     t=al*r+bt*dlog(r)
        t=(ro-t)/(al*r+bt)
        if (t.le.-1.d0) t=-0.5d0
        r=r*(1.d0+t)
        if (dabs(t).gt.0.1d-8) goto 200
        rr=r
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
