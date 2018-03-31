c       =================================================
        subroutine nucl(maxii,r,unuc)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine computes "nonpoint" part of the
c       nuclear potential, i.e., the deviation of the nuclear
c       potential from the point-nucleus Coulomb potential
c
c       Input parameters:
c       r(i)    .... the radial grid
c       maxii   .... the size of the r array
c
c       Output parameters:
c       unuc(i) .... Nonpoint part of the nuclear potential
c                    multiplied by r(i).
c
c       Common-block parameters:
c       ii     ..... Number of grid points
c       cl     ..... Speed of light 
c       z      ..... Nuclear charge
c       knucl  ..... Nuclear model (0 - point, 1 - Fermi model,
c                    2 - uniform sphere)
c       rnucl  ..... Radius of the nuclear sphere 
c       vnuc(n) .... The Taylor expansion coefficients of the
c                    nuclear potential around 0 (n=0,..,nmax).
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
        common /z/z/cl/cl/ii/ii
        real*8 r(maxii),unuc(maxii)
        common /inucl/inucl/knucl/knucl/rnucl/rnucl
     1  /vnuc/vnuc(20)
        common /nmax/nmax
        common /fermi/aa,tt,cc,ro0,rnucl0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (pi=3.1415926535897932385d0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        rnucl0=rnucl/dsqrt(5.d0/3.d0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do n=0,nmax
          i=n+1
          vnuc(i)=0.d0
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Point nuclear
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (knucl.eq.0) then
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c          write( *,5)
c          write(11,5)
c5         format(/2x,'Point Nuclear Charge'/2x,20('*'))
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          vnuc(1)=-z
          goto 1000
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        fermi2at=1.d-5/0.529177249d0
        rfermi=rnucl/fermi2at
        rms_au=rnucl/dsqrt(5.d0/3.d0)
        rms_fermi=rfermi/dsqrt(5.d0/3.d0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Volume distribution
c       - - - - - - - - - - - - - - - - - - - - - - - - -
200     if (knucl.eq.2) then
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c          write( *,15)
c          write(11,15)
c15        format(/2x,'Volume Distribution of the Nuclear Charge'
c     &    /2x,41('*'))
c          write( *,25) rnucl,rfermi,rms_au,rms_fermi
c          write(11,25) rnucl,rfermi,rms_au,rms_fermi
c25        format 
c     &    ( 2x,'Rnucl     =',e14.6,1x,'[a.u.] =',f8.4,1x,'[Fermi]',
c     &     /2x,'Sqrt<R^2> =',e14.6,1x,'[a.u.] =',f8.4,1x,'[Fermi]')
c
          vnuc(2)=-1.5d0*z/rnucl
          vnuc(4)= 0.5d0*z/rnucl**3
          goto 1000
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (knucl.eq.1) then
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Fermi distribution
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       1 fermi = 10**(-13) cm
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        fermi2at=1.d-5/0.529177249d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c        write( *,35)
c        write(11,35)
c35      format(/2x,'Fermi distribution of the nuclear charge:')
c
c        write( *,25) rnucl,rfermi,rms_au,rms_fermi
c        write(11,25) rnucl,rfermi,rms_au,rms_fermi
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Setting nuclear parameters
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        tt=2.30d0
        aa=tt/(4.d0*dlog(3.d0))*fermi2at
        cc=5.d0/3.d0*rnucl0**2-7.d0/3.d0*pi**2*aa**2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (cc.le.0.d0) then
          write( *,'(/2x,a,1x,a/2x,a)') '*** Warning!!!',
     &    'Z too small. Fermi Distribution does not work ***.',
     &    'The program use Sphere Ditribution.' 
          write(11,'(/2x,a,1x,a/2x,a)') '*** Warning!!!',
     &    'Z too small. Fermi Distribution does not work ***.',
     &    'The program use Sphere Ditribution.'
          knucl=2
          goto 200
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        cc=dsqrt(cc)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c        write( *,45) cc,cc/fermi2at,tt*fermi2at,tt
c        write(11,45) cc,cc/fermi2at,tt*fermi2at,tt
c45      format (2x,'Parameter c=',e13.6,1x,'[a.u.] =',f8.4,
c     &  1x,'[Fermi]',/2x,'Parameter t=',e13.6,1x,
c     &  '[a.u.] =',f8.4,1x,'[Fermi]')
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dnorm=cc**3*(1+(pi*aa/cc)**2-6*(aa/cc)**3*sk(3,-cc/aa))/3
        ro0=z/dnorm
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation of the Taylor coefficients
c       - - - - - - - - - - - - - - - - - - - - - - - - -
	B0=ro0/(1+dexp(-cc/aa))
	B1=-ro0/((1+dexp(-cc/aa))**2)*dexp(-cc/aa)/aa
        B2=B1/(1+dexp(-cc/aa))*(1.d0-dexp(-cc/aa))/aa
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       RM(n)=<ro(r)*r^n>/ro0;  ro0=2*ro(c)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
	rm1=rm(1,cc,aa)
        vnuc(2)=-rm1*ro0
        vnuc(4)=B0/6.d0
        vnuc(5)=B1/12.d0
        vnuc(6)=B2/40.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (knucl.eq.3) then
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Gaussian distribution
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       V(r)=-Z*erf(dsqrt(1.5)*r/Rms)/r
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c        write( *,65)
c        write(11,65)
c65      format(/2x,'Gaussian Distribution of the Nuclear Charge'
c     &  /2x,43('*'))
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dn0=-2*z/dsqrt(pi)
        a0=dsqrt(1.5d0)/rnucl0
        vnuc(2)= dn0*a0
        vnuc(4)=-dn0*a0**3/3.d0
        vnuc(6)= dn0*a0**5/10.d0
        vnuc(8)=-dn0*a0**7/42.d0
        vnuc(10)=dn0*a0**9/216.d0
c       vnuc(12)=-dn0*a0**11/1320.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (dabs(vnuc(1)).gt.dabs(cl)) then
c          write( *,'(/2x,a)') '*** Error: Vnuc(1) > |Cl| ***'
c          write(11,'(/2x,a)') '*** Error: Vnuc(1) > |Cl| ***'
        call exit1
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
1000    do i=1,maxii
          unuc(i)=0.d0
        enddo
        if (knucl.eq.0) return
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (knucl.eq.2) then
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          if (inucl.lt.1) return
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          do i=1,inucl
            vi=0.d0
            do n=0,nmax
              vi=vi+vnuc(n+1)*r(i)**n
            enddo
            unuc(i)=vi+z
          enddo
          goto 1100
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Fermi and Gaussian distributions        
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,ii
          if (knucl.eq.1) vi=vpot_fermi(z,cc,aa,r(i))
          if (knucl.eq.3) vi=vpot_gauss(z,rnucl0,r(i))
          unuc(i)=z+vi
        enddo
        inucl=ii
        do i=ii,1,-1
          if (dabs(unuc(i)).gt.1.d-13) goto 1100
          unuc(i)=0.d0
          inucl=i
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
 1100   unuc(ii+3)=inucl
        unuc(ii+4)=0.d0
        unuc(ii+5)=z+vnuc(1)
        do m=1,nmax-1
          i=ii+5+m
          unuc(i)=vnuc(m+1)*r(1)**m
        enddo
        unuc(ii+5+nmax)=0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        real*8 function ro_nucl(r)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8(a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /knucl/knucl/rnucl/rnucl/z/z
        common /fermi/aa,tt,cc,ro0,rnucl0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine computes the charge density
c       ditribution of the nucleus
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (pi=3.1415926535897932385d0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Volume distribution
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (knucl.eq.2) then
          if(r.lt.rnucl+1.d-10)then
            ro_nucl=3.d0*z/rnucl**3
          else
            ro_nucl=0.d0
          endif
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Fermi distribution
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (knucl.eq.1) then
          x=(r-cc)/aa
          if (x.le.200) then 
            ro_nucl=ro0/(1.d0+dexp(x))
          else
            ro_nucl=0.d0
          endif
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Gaussian distribution
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (knucl.eq.3) then
          ro_nucl=z*4*pi*(1.5d0/(pi*rnucl0**2))**1.5d0*
     &             dexp(-1.5d0*(r/rnucl0)**2)
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================              
        function vpot_fermi(z,c,a,r)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This routine computes the nuclear potential for
c       the Fermi distrubution of the nclear density.
c       Details see in: F.A.Parpia, A.K.Mohanty, PRA,
c       v.46, 3735 (1992).
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (A-H,O-Z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (pi=3.1415926535897932385d0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dn=1+(pi*a/c)**2-6d0*((a/c)**3)*SK(3,-c/a)
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
        if (r.lt.c) then
        B1=-6d0*((a/c)**3)*SK(3,-c/a)+6d0*((a/c)**3)*SK(3,(r-c)/a)
        B2=r/c*(3d0/2d0+((pi*a/c)**2)/2d0-3d0*((a/c)**2)*SK(2,(r-c)/a))
        B3=-((r/c)**3)/2d0
        V=Z/dn*(B1+B2+B3)
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
        else
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
        B1=dn
        B2=3d0*((a/c)**2)*(r/c*SK(2,(c-r)/a)+2d0*a/c*SK(3,(c-r)/a))
        V=Z/dn*(B1+B2)
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -      	
        vpot_fermi=-V
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
	return
        end
c       =================================================              
        function vpot_gauss(z,rms,r)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Evaluation of the Gauss potential
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (A-H,O-Z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        external erf
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        vpot_gauss=-z*erf(dsqrt(1.5d0)*r/rms)
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
	return
        end
c       =================================================
        function sk(k,x)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Auxiliary function for the evaluationb of the Fermi 
c       potential
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (error=1d-20)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dk=K
        N=0
        sk=0d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
10      N=N+1
        dn=N
        add=DEXP(dN*X)/(dN**dK)*(-1)**N
	if(N.eq.1) add1=add
        sk=sk+add
        if (add.eq.0d0) goto 1000
        if (dabs(add/add1).lt.error) goto 1000
        goto 10
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
1000    return
        end
c       =================================================
        function rm(M,C,A)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (A-H,O-Z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (error=1.d-4)
        parameter (pi=3.1415926535897932385d0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (m.eq.1) then
          rm=c*c/2+a**2*sk(2,-c/a)+2*a**2*(pi**2/12)
        return
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        rM=c**(M+1)/(M+1)
        rm1=c**(M+1)/(M+1)
        N=0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
10      N=N+1
        dn=N
        add=((-1)**N)*(fxneax(M,dN/A,C)-DEXP(-dN*C/A)*
     &  fxneax(M,dN/A,0d0))
        rm=rm+add
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if(dabs(add/rm1).lt.error) goto 20
        goto 10
c       - - - - - - - - - - - - - - - - - - - - - - - - -       
20      N=0
        addd=add
c       - - - - - - - - - - - - - - - - - - - - - - - - -             
30      N=N+1
        dN=N
        add=((-1)**(N+1))*(-fxneax(M,-dN/A,C))
        IF(N.eq.1) RM1=add
        RM=RM+add
c       - - - - - - - - - - - - - - - - - - - - - - - - -             
        if(dabs(add).lt.dabs(addd)) goto 1000
        goto 30
c       - - - - - - - - - - - - - - - - - - - - - - - - -
1000    return
        end
c       =================================================
        FUNCTION fxneax(N,A,X)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (A-H,O-Z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        SUM=0d0
        DO 10 l=0,N
c       SUM=SUM+dfac(n)/dfac(n-l)*((A*X)**(n-l))*(-1)**l
        SUM=SUM+dfac(n)/dfac(n-l)*(stpdi(A*X,n-l))*(-1)**l
10      CONTINUE
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        fxneax=SUM/(A**(N+1))
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        RETURN
        END
c       =================================================
        function dfac(n)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Evaluation of factorial
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit double precision (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (one=1.0D+00)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if(n.eq.0) then
          dfac=one
	  return
        end if
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        df= one
        ns=n
        do i=ns,1,-1
          df=df*i
        enddo
        dfac=df
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================       
        function stpdi(X,I)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        IMPLICIT REAL*8 (A-H,O-Z)      
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        IF(I.EQ.0) THEN
          stpdi=1.d0
          goto 1000
        ELSE
          stpdi=x**i
          goto 1000
        END IF 
c       - - - - - - - - - - - - - - - - - - - - - - - - -
1000    return
        END
c       =================================================              
        DOUBLE PRECISION FUNCTION erf(x)
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
C       Evaluation of the real error function.
c       This subroutine was taken from the Naval Surface
C       Warfare Center (NSWC) Mathematics Library (1993 version).
c       http://www.ualberta.ca/CNS/RESEARCH/Software/
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        DOUBLE PRECISION x
C       ..
C       .. Local Scalars ..
        DOUBLE PRECISION ax,bot,c,t,top,x2
C       ..
C       .. Local Arrays ..
        DOUBLE PRECISION a(5),b(3),p(8),q(8),r(5),s(4)
C       ..
C       .. Intrinsic Functions ..
        INTRINSIC abs,exp,sign
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
C       Data statements ..
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
        DATA c/.564189583547756D0/
        DATA a(1)/.771058495001320D-04/,a(2)/-.133733772997339D-02/,
     +       a(3)/.323076579225834D-01/,a(4)/.479137145607681D-01/,
     +       a(5)/.128379167095513D+00/
        DATA b(1)/.301048631703895D-02/,b(2)/.538971687740286D-01/,
     +       b(3)/.375795757275549D+00/
        DATA p(1)/-1.36864857382717D-07/,p(2)/5.64195517478974D-01/,
     +       p(3)/7.21175825088309D+00/,p(4)/4.31622272220567D+01/,
     +       p(5)/1.52989285046940D+02/,p(6)/3.39320816734344D+02/,
     +       p(7)/4.51918953711873D+02/,p(8)/3.00459261020162D+02/
        DATA q(1)/1.00000000000000D+00/,q(2)/1.27827273196294D+01/,
     +       q(3)/7.70001529352295D+01/,q(4)/2.77585444743988D+02/,
     +       q(5)/6.38980264465631D+02/,q(6)/9.31354094850610D+02/,
     +       q(7)/7.90950925327898D+02/,q(8)/3.00459260956983D+02/
        DATA r(1)/2.10144126479064D+00/,r(2)/2.62370141675169D+01/,
     +       r(3)/2.13688200555087D+01/,r(4)/4.65807828718470D+00/,
     +       r(5)/2.82094791773523D-01/
        DATA s(1)/9.41537750555460D+01/,s(2)/1.87114811799590D+02/,
     +       s(3)/9.90191814623914D+01/,s(4)/1.80124575948747D+01/
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
C       Executable Statements ..
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
        ax = abs(x)
        IF (ax.GT.0.5D0) GO TO 10
        t = x*x
        top = ((((a(1)*t+a(2))*t+a(3))*t+a(4))*t+a(5)) + 1.0D0
        bot = ((b(1)*t+b(2))*t+b(3))*t + 1.0D0
        erf = x* (top/bot)
        RETURN
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
   10   IF (ax.GT.4.0D0) GO TO 20
        top = ((((((p(1)*ax+p(2))*ax+p(3))*ax+p(4))*ax+p(5))*ax+
     +        p(6))*ax+p(7))*ax + p(8)
        bot = ((((((q(1)*ax+q(2))*ax+q(3))*ax+q(4))*ax+q(5))*ax+
     +        q(6))*ax+q(7))*ax + q(8)
        erf = 0.5D0 + (0.5D0-exp(-x*x)*top/bot)
        IF (x.LT.0.0D0) erf = -erf
        RETURN
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
   20   IF (ax.GE.5.8D0) GO TO 30
        x2 = x*x
        t = 1.0D0/x2
        top = (((r(1)*t+r(2))*t+r(3))*t+r(4))*t + r(5)
        bot = (((s(1)*t+s(2))*t+s(3))*t+s(4))*t + 1.0D0
        erf = (c-top/ (x2*bot))/ax
        erf = 0.5D0 + (0.5D0-exp(-x2)*erf)
        IF (x.LT.0.0D0) erf = -erf
        RETURN
c       - - - - - - - - - - - - - - - - - - - - - - - - -      
   30   erf = sign(1.0D0,x)
        RETURN
        END
