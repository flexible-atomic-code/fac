c       =================================================
        subroutine uehling(maxii,r,unuc,uehl)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine computes Uehling vacuum-polarization
c       potential by using the accurate approximation method,
c       from the paper:
c       L.W.Fullerton and G.A. Ringer, PRA 13, 1283 (1976)
c       This approximation may be applied to point-charge
c       or to finite charge distribution nuclear models.
c
c       Input parameters:
c       r(i)    .... Radial grid
c       maxii   .... the size of the r array
c       unuc(i) .... Nonpoint part of the nuclear potential
c                    multiplied by r(i).
c
c       Output parameters:
c       uehl(i) .... Uehling potential multiplied by r(i).
c
c       Common block parameters:
c       cl     ..... Speed of light
c       z      ..... Nuclear charge
c       knucl  ..... Nuclear model (0 - point, 1 - Fermi model,
c                    2 - uniform sphere)
c       ii     ..... Number of grid points
c       h, al, bt ..... Parameter of the radial grid
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /cl/cl/z/z/ii/ii/h/h/al/al/bt/bt
        common /inucl/inucl/knucl/knucl
        real*8 r(maxii),unuc(maxii),uehl(maxii)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (pi=3.1415926535897932385d0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation Uehling potential times on r
c       - - - - - - - - - - - - - - - - - - - - - - - - -
	Do i=1,maxii
          uehl(i)=0.D0
	End do
c       - - - - - - - - - - - - - - - - - - - - - - - - -
	coef1=-2.d0/3.d0/cl**2/(4*pi)
	coef2=-z*8.d0/3.d0/cl/(4*pi)
	ir=0
 500    ir=ir+1
        rr=r(ir)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (knucl.ne.0) then
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Finite charge nuclear distribution
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        npoints=4*maxii
        r0=r(1)/10
        rmax=r(inucl)
        tmax=dlog(rmax)
        tmin=dlog(r0)
        hi=(tmax-tmin)/(npoints-1)
        fint=0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,npoints
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          t=hi*(i-1)+tmin
          ri=dexp(t)
          x1=2*cl*dabs(rr-ri)
          x2=2*cl*dabs(rr+ri)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          d1=dk0(x1)
          d2=dk0(x2)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          fu=d1-d2
          fr=fu*ro_nucl(ri)*ri
          fr=fr*ri
          fint=fint+fr*hi
          if (i.eq.npoints-2) then
            fr1=fr
          endif
          if (i.eq.npoints-1) then
            fr2=fr
          endif
          if (i.eq.npoints) then
            fr3=fr
          endif
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dfr=(3*fr3-4*fr2+fr1)/(2*hi)
        fint=fint-0.5d0*hi*fr3-hi**2/12.d0*dfr
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        Uehl(ir)=fint*coef1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        else
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Point nuclear distribution
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          x=2*cl*rr
          Uehl(ir)=dk1(x)*coef2
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
	If (dabs(uehl(ir)/r(ir)).gt.1.D-40.and.ir.lt.ii) go to 500
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        imax=ir
        uehl(ii+3)=imax
        uehl(ii+4)=0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        real*8 function dk0(x)
        implicit real*8(a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Fullerton'a and Ringer'a, PRA 13, 1283 (1976)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dk0=0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        a0= 0.88357293375d0
        a1=-0.28259817381d0
        a2=-0.58904879578d0
        a3= 0.12500133434d0
        a4=-0.032729913852d0
        a5= 0.0082888574511d0
        a6=-0.0000103277658d0
        a7= 0.0000636436689d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        b0=-319.999594323d0
        b1= 2.53900995981d0
        b2= 1.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        c0=-319.999594333d0
        c1= 2.53901020662d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        d0=5.018065179d0
        d1=71.51891262d0
        d2=211.6209929d0
        d3=31.40327478d0
        d4=-1.d0 
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        e0=2.669207401d0
        e1=51.72549669d0
        e2=296.9809720d0
        e3=536.4324164d0
        e4=153.5335924d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if(x.ge.0.d0.and.x.le.1.d0) then       
          dp7=a0+a1*x+a2*x**2+a3*x**3+a4*x**4+a5*x**5
     #        +a6*x**6+a7*x**7
          dp21=(b0+b1*x**2+b2*x**4)/(c0+c1*x**2)
          if (x.gt.0) then
            dk0=dp7+x*dp21*dlog(x)
          else
            dk0=dp7
          endif
        else
          dk0=dexp(-x)*(d0+d1/x+d2/x**2+d3/x**3+d4/x**4)
     #        /(dsqrt(x))**3/(e0+e1/x+e2/x**2+e3/x**3+e4/x**4)
        end if
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        real*8 function dk1(x)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8(a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Fullerton'a and Ringer'a, PRA 13, 1283 (1976)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dk1=0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        a0=-0.71740181754d0
        a1= 1.1780972274d0
        a2=-0.37499963087d0
        a3= 0.13089675530d0
        a4=-0.038258286439d0
        a5=-0.0000242972873d0
        a6=-0.0003592014867d0
        a7=-0.0000171700907d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        b0=-64.0514843293d0
        b1= 0.711722714285d0
        b2= 1.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        c0= 64.0514843293d0
        c1= -0.711722686403d0
        c2=0.0008042207748
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        d0=217.2386409d0
        d1=1643.364528d0
        d2=2122.244512d0
        d3=-45.12004044d0
        d4=1.d0 
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        e0=115.5589983d0
        e1=1292.191441d0
        e2=3831.198012d0
        e3=2904.410075d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if(x.ge.0.d0.and.x.le.1.d0) then       
          dp7=a0+a1*x+a2*x**2+a3*x**3+a4*x**4+a5*x**5
     #        +a6*x**6+a7*x**7
          dp21=(b0+b1*x**2+b2*x**4)/(c0+c1*x**2+c2*x**4)
          dk1=dp7+dp21*dlog(x)
        else
          dk1=dexp(-x)*(d0+d1/x+d2/x**2+d3/x**3+d4/x**4)
     #        /(dsqrt(x))**3/(e0+e1/x+e2/x**2+e3/x**3)
        end if
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
