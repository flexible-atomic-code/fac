c       =================================================
        subroutine dirac(n,kappa,p,q)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine computes radial Dirac wavefunctions
c       for an extended nuclear potential by solving 
c       numerically the Dirac equation
c
c       Input parameters:
c       n ......... Principal quantum number
c       kappa ..... Relativistic angular quantum number
c
c       Output parameters:
c       p(i)  ..... Large component P(r) of the wavefunction 
c       q(i)  ..... Small component Q(r) of the wavefunction
c       i=1,...,ii
c
c       Common-block parameters:
c       z     ..... Nuclear charge
c       cl    ..... Speed of light
c       e     ..... one-electron energy
c       ii    ..... Number of grid points
c       h     ..... grid step
c       r(i)  ..... radial grid
c       y(i)  ..... Coulomb (or external) potential multiplied
c                   by r(i)
c       unuc(i) ... Nonpoint part of the Nuclear potential
c                   multiplied by r(i)
c       cp(i),cq(i)..corrections for the finite difference scheme
c       knucl  ..... Nuclear model (0 - point, 1 - Fermi model,
c                    2 - uniform sphere)
c       vnuc(n) .... The Taylor expansion coefficients of the
c                    nuclear potential around 0 (n=0,..,nmax).
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /cl/cl/z/z
        common /ii/ii/h/h/r1/r1/al/al/bt/bt
        common /a/a(maxii)/b/b(maxii)/c/c(maxii)
     4  /cp/cp(maxii)/cq/cq(maxii)/y/y(maxii)
     5  /r/r(maxii)/v/v(maxii)/unuc/unuc(maxii)
        common /nmax/nmax
        common /knucl/knucl/inucl/inucl/rnucl/rnucl
     &  /vnuc/vnuc(20)
        common /niter/niter/nit/nit/mi/m1,m2,m3
        real*8 p(maxii),q(maxii)
        real*8 w(maxii)
        real*8 fac(30)
        character*1 let(11)
        data
     1  let /'s','p','d','f','g','h','i','k','l','m','n'/
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        nmax=9
        l=(iabs(2*kappa+1)-1)/2
        k=kappa
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,ii
          y(i)=-z
        enddo
        do m=0,nmax
          i=ii+5+m
          y(i)=0.d0
        enddo
        y(ii+5)=0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        eps0=1.d-10*(256.d0/maxii)**2
        if (eps0.lt.1.d-14) eps0=1.d-14
        tol=1.d-50
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        fac(1)=0.d0
        do i=2,30
          d=i-1.d0
          fac(i)=fac(i-1)+dlog(d)
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,maxii
          p(i)=0.d0
          q(i)=0.d0
        enddo
        e=(z/n)**2
        p(ii+1)=e
        p(ii+2)=1.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       m1 - First turning point
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dm=z*z-e*(l+0.5d0)**2
        if (dm.gt.0.1) then
          rm1=(l+0.5d0)**2/(z+dsqrt(dm))
          m1=(al*(rm1-r(1))+bt*dlog(rm1/r(1)))/h+1.d0
        else
          m1=11
        endif
        if (m1.lt.3) m1=3
        m1=(m1/2)*2+1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        p(ii+15)=m1
        p(ii+3)=ii
        q(ii+3)=ii
        gam=dsqrt(1.d0-(vnuc(1)/cl)**2)
        p(ii+4)=gam
        q(ii+4)=gam
        p(ii+5)=0.d0
        q(ii+5)=0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dz=z
        if (z.lt.0.1d0) dz=1.d0
        a0=dexp(dlog(2.d0*dz/n)+(l+0.5d0)*
     1  dlog(2.d0*dz/n)-
     1  fac(l+l+2)+0.5d0*(fac(n+l+1)-fac(n-l)))/dsqrt(2.d0*n)
        if (kappa.gt.0)
     1  q(ii+5)=-a0*(l+kappa+1)/(2*cl)
        if (kappa.lt.0) p(ii+5)=a0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        c1=0.01d0
        eps=eps0*0.1d0
        hh=0.5d0*h
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        ee1=0.d0
        dd1=0.d0
        niter=0
800     niter=niter+1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (niter.gt.200) then
          write( *,'(/2x,a)') 'Niter > 200 in dirac'
          write(11,'(/2x,a)') 'Niter > 200 in dirac'
          call exit1
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,ii
          cp(i)=0.d0
          cq(i)=0.d0
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Boundary conditions at the origin
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call origin(kappa,p,q)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        gam=p(ii+4)
        e=p(ii+1)
        e0=e
        m1=p(ii+15)+c1
        imax=ii
        iemax=0
        emax=0.d0
        iemin=0
        emin=0.d0
        iflag=0
        ifail=0
        istep=0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation of the homogeneous part
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        tn=dsqrt(1.d0/p(ii+2))
        t=hh/cl
        t1=t*tn
        s=0.5d0*e
        hk=hh*kappa
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 20 i=1,ii
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        c(i)=hk
        yi=y(i)+unuc(i)
        a(i)=t*(s*r(i)+yi)
        d=v(i)/r(i)
        c(i)=c(i)*d
        a(i)=a(i)*d
c       - - - - - - - - - - - - - - - - - - - - - - - - -
20      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation of the point m3
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        m3=ii
        s=1.d0/(h*h)
        do 30 i=1,ii
        j=ii+1-i
        if (e*v(j)*v(j).ge.s) goto 30
        m3=j
        goto 320
30      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
320     if (niter.gt.1) then
          if (knucl.ne.1.or.inucl.le.1) then
            call correct(1,ii-2,tn,p,q,cp,cq)
          else
            call correct(1,inucl,tn,p,q,cp,cq)
            call correct(inucl,ii-2,tn,p,q,cp,cq)
          endif
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        n1=0
        n2=0
        n3=0
        nt1=0
        nit=0
        nt3=0
        nj=0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       The begining of the iteration procedure
c       - - - - - - - - - - - - - - - - - - - - - - - - -
600     if (nit.le.380) goto 500
        if(nt1.eq.0.and.nt3.eq.0)
     1  write( *,5) iflag,m1,m2,m3,nj,e0,e
        if (nt1.eq.0.and.nt3.eq.0)
     1  write(11,5) iflag,m1,m2,m3,nj,e0,e
5       format(' Flag=',i1,2x,'m1=',i4,2x,'m2=',i4,2x,
     1  'm3=',i4,2x,'nj =',i2,2x,'e =',2e15.8)
        if (nt1.ne.0.or.nt3.ne.0)
     1  write( *,15) iflag,m1,m2,m3,nj,e0,dp
        if (nt1.ne.0.or.nt3.ne.0)
     1  write(11,15) iflag,m1,m2,m3,nj,e0,dp
15      format(' Flag=',i1,2x,'m1=',i4,2x,'m2=',i4,2x,
     1  'm3=',i4,2x,'nj=',i2,2x,'e =',e15.8,2x,'df=',e12.4)
        if (nit.le.450) goto 500
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,m2
          write(11,'(6e14.6)') r(i),p(i),q(i),y(i),cp(i),cq(i)
        enddo
        j=2*iabs(kappa)-1
        write( *,25) n,let(l+1),j
        write(11,25) n,let(l+1),j
25      format(/2x,'Divergancy for the shell:',i3,a1,i2,'/2')
        write( *,35) niter
        write(11,35) niter
35      format(2x,'nit>250',2x,'Niter=',i5)
        close (unit=11)
        call exit1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
500     nit=nit+1
        nj=0
        d=0.5d0*(e-e0)*hh/cl
        s=cl*h
        st=0.d0
        do i=1,ii
          t=a(i)+d*v(i)
          a(i)=t
          b(i)=s*v(i)-a(i)
          w(i)=c(i)*c(i)+a(i)*b(i)
          if (st.gt.t) st=t
        enddo
        if (st.lt.0.d0) goto 350
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       iflag=1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        iflag=1
        n1=n1+1
        if (iemax.eq.0) emax=e
        if (iemax.eq.1.and.e.lt.emax) emax=e
        iemax=1
        e0=e
        e=e*(1.d0-0.5d0*n1/(n1+1.d0))
        if (iemin.eq.1.and.e.le.emin) e=0.5d0*(e0+emin)
        goto 600
350     n1=0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       m2 - Second turning point
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (nt1.gt.0.and.nt3.gt.0) goto 380
        im=m1+m3
        s=0.5d0*kappa*(kappa+1)*hh/cl
        do 70 i=m1,m3
        j=im-i
        rj=r(j)
        df=a(j)+s*v(j)/(rj*rj)
        if (df.gt.0.d0.and.w(j).gt.0.d0) goto 70
        m2=j
        goto 370
70      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       iflag=2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        iflag=2
        if (iemax.eq.0) emax=e
        if (iemax.eq.1.and.e.lt.emax) emax=e
        iemax=1
        e0=e
        e=e*0.8d0
        if (iemin.eq.1.and.e.le.emin) e=0.5d0*(e0+emin)
        goto 600
370     if (m2.lt.m3-2) goto 380
        if (iemin.eq.0) emin=e
        if (iemin.eq.1.and.e.gt.emin) emin=e
        iemin=1
        e0=e
        e=e*1.2d0
        if (iemax.eq.1.and.e.ge.emax) e=0.5d0*(e0+emax)
        goto 600
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Boundary conditions at the origin
c       - - - - - - - - - - - - - - - - - - - - - - - - -
380     i0=1
        pi=0.d0
        qi=0.d0
        do m=0,nmax
          j=ii+5+m
          pi=pi+p(j)
          qi=qi+q(j)
        enddo
        st=r(i0)**gam
        pi=pi*st
        qi=qi*st
        p(i0)=(1.d0+c(i0))*pi+b(i0)*qi
        q(i0)=(1.d0-c(i0))*qi+a(i0)*pi
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation functions p(i),q(i) at points (i=1,m2)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        i1=i0+1
        pj=p(i0)
        qj=q(i0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 120 i=i1,m2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        j=i-1
        dw=w(j)
        dc=c(j)
        da=a(j)
        db=b(j)
        dt=0.5d0*(1.d0-dw)
        dq=((dw+dc)*qj-da*pj)/dt+cq(j)
        dp=((dw-dc)*pj-db*qj)/dt+cp(j)
        pj=pj+dp
        qj=qj+dq
        p(i)=pj
        q(i)=qj
c       - - - - - - - - - - - - - - - - - - - - - - - - -
120     continue
        pm2=pj
        qm2=qj
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Testing number of nodes.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        nj=n-l-1
        do i=2,m2
          if (p(i)*p(i-1).lt.0) nj=nj-1
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (nj.eq.0) goto 390
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       iflag=3
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        iflag=3
        i=1
        if (nj.lt.0) i=-1
        n2=n2+1
        e0=e
        j=iabs(nj)
        if (j.lt.n2) j=n2
        e=e*(1.d0-0.04d0*(i*j)/(j+2.d0))
        if (e.gt.2*cl**2) e=0.5d0*(e0+2*cl**2)
c        if (e.gt.2*cl**2) e=2*cl**2
        if (iemin.eq.1.and.e.le.emin) e=0.5d0*(emin+e0)
        if (iemax.eq.1.and.e.ge.emax) e=0.5d0*(emax+e0)
        if (nj.gt.0.and.iemax.eq.0) emax=e0
        if (nj.gt.0.and.iemax.eq.1.and.e.lt.emax) emax=e0
        if (nj.gt.0) iemax=1
        if (nj.lt.0.and.iemin.eq.0) emin=e0
        if (nj.lt.0.and.iemin.eq.1.and.e.gt.emin) emin=e0
        if (nj.lt.0) iemin=1
        if (iemin.eq.0.or.iemax.eq.0) goto 600
        if (emax-emin.gt.eps*emax*0.00001) goto 600
        ifail=1
        goto 390
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Boundary conditions at the infinity
c       - - - - - - - - - - - - - - - - - - - - - - - - -
390     im=ii-2
        i1=im+imax
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 130 i=im,imax
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        j=i1-i
        dw=w(j)
        if (dw.gt.0.d0) then
          dw=dsqrt(dw)
        else
          write( *,45)
          write(11,45)
45        format(2x,'Warning!!  R2 too small or',
     1    ' Bound state does not exist.')
          call exit1
        endif
        db=b(j)
        dc=c(j)
        pj=(dw-dc)/db
        qj=cq(j)/(dw+cp(j))
        qj=0.d0
        p(j)=pj
        q(j)=qj
c       - - - - - - - - - - - - - - - - - - - - - - - - -
130     continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation the ratio q(i)/p(i) at points (i=m2,ii)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        qj=qj*(1.d0-dc-pj*db)
        p(im)=pj
        q(im)=qj
        j=im-1
400     dc=c(j)
        da=a(j)
        db=b(j)
        dw=w(j)
        dq=dc+dc
        dt=1.d0+dw+dq+2*db*pj
        if (j.le.m3) dp=(da+da+(1.d0+dw-dq)*pj)/dt
        if (j.gt.m3) dp=(dsqrt(dw)-dc)/db
        qj=(qj+pj*cp(j)-cq(j))/dt*(1.d0-dw)
        pj=dp
        if (j.eq.m2) goto 410
        p(j)=pj
        q(j)=qj
        j=j-1
        goto 400
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Sewing the wave functions at the point m2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
410     dq=qm2
        qm2=pj*pm2+qj
        dp=qm2-dq
        dq=dq+kappa/r(m2)*pm2/(2*cl)
        t=dabs(dp)/(dabs(dq)+eps/cl)
c       if (t.lt.5*eps.and.nit.gt.2) goto 450
        if (t.lt.5*eps) goto 450
        if (dp*b(m2).le.0.d0) goto 420
        j=1
        e3=e
        dt3=dp
        nt3=1
        goto 430
420     e1=e
        dt1=dp
        nt1=1
        j=-1
430     if (nt1.eq.nt3) goto 440
        istep=0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       iflag=4
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        iflag=4
        i=n-l-1
        n3=n3+1
        d=i-2*(i/2)-0.5d0
        t=0.01d0
        e0=e
        d=t*(j*n3)*d/(n3+3.d0)
        e=e*(1.d0+d)
c        if (e.gt.2*cl**2) e=0.5d0*(e0+2*cl**2)
         goto 600
        if (iemin.eq.1.and.e.le.emin) e=0.5d0*(e0+emin)
        if (iemax.eq.1.and.e.ge.emax) e=0.5d0*(e0+emax)
        if (d.le.0.d0.and.iemax.eq.0) emax=e0
        if (d.le.0.d0.and.iemax.eq.1.and.e.le.emax) emax=e0
        if (d.le.0.d0) iemax=1
        if (d.ge.0.d0.and.iemin.eq.0) emin=e0
        if (d.ge.0.d0.and.iemin.eq.1.and.e.ge.emin) emin=e0
        if (d.ge.0.d0) iemin=1
        if (iemin.eq.0.or.iemax.eq.0) goto 600
        if (emax-emin.gt.eps*emax*0.00001) goto 600
        ifail=1
        goto 450
440     n3=0
        e0=e
        e=(dt3*e1-dt1*e3)/(dt3-dt1)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        istep=istep+1
        if (istep.gt.3) then
           e=0.5d0*(e1+e3)
           istep=0
        endif
        if (nt3.gt.5.or.nt1.gt.5) then
          e=0.5d0*(e1+e3)
          nt3=1
          nt1=1
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       iflag=5
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        iflag=5
        goto 600
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation functions p(i),q(i) at points (i=m2,ii)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
450     i2=m2
        im=ii-2
        do 140 i=m2+1,im
        j=i-1
        d=c(j)+c(j)
        s=b(j)+b(j)
        t=((1.d0+w(j)-d)*p(j)-s*q(j))/(1.d0-w(j))
        pj=t+cp(j)
        qj=p(i)*pj+q(i)
        if (dabs(pj).lt.1.d-20) goto 460
        i2=i
        p(i)=pj
        q(i)=qj
140     continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Reconstruction functions p(i),q(i) at points (i=1,m2)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
460     do j=1,i2
          pj=p(j)
          qj=q(j)
          t=1.d0-w(j)
          p(j)=((1.d0-c(j))*pj-b(j)*qj)/t
          q(j)=((1.d0+c(j))*qj-a(j)*pj)/t
        enddo
        if (i2.lt.im) goto 470
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        i2=im
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 160 i=im+1,ii
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        d=(r(i)-r(i-1))/(hh*v(i))
        s=dsqrt(w(i))
        t=cp(i)-s
        t1=t*d
        if (dabs(t1).lt.0.01d0) t2=-d
        if (dabs(t1).ge.0.01d0) t2=(1.d0-dexp(t1))/t
        t1=p(i-1)*dexp(-s*d)
        t2=t2*b(i)*q(i)
        pj=t1+t2
        if (dabs(pj).lt.tol) goto 470
        i2=i
        t2=t1
        q(i)=pj*p(i)+q(i)
        p(i)=pj
c       - - - - - - - - - - - - - - - - - - - - - - - - -
160     continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        goto 1000
c       - - - - - - - - - - - - - - - - - - - - - - - - -
470     imax=i2
        if (imax.eq.ii) goto 1000
        do i=imax+1,ii
          p(i)=0.d0
          q(i)=0.d0
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
1000    p(ii+1)=e
        p(ii+3)=imax
        p(ii+18)=m2
        q(ii+1)=e
        q(ii+3)=imax
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Normalization
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        d=tint(0,p,q,p,q,r,v)
        d=1.d0/dsqrt(d)
        do i=1,imax
          p(i)=p(i)*d
          q(i)=q(i)*d
        enddo
        do m=0,nmax
          i=m+ii+5
          p(i)=p(i)*d
          q(i)=q(i)*d
        enddo
        p(ii+2)=1.d0
        q(ii+2)=1.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (dabs(p(imax)).gt.1.d-3) then
        write( *,65) kappa,imax,p(imax)
        write(11,65) kappa,imax,p(imax)
65      format(2x,'Warning!!  R2 too small. Kappa =',i2,
     1  '  imax=',i5,'  p(imax) =',e9.2)
        ifail=1
        if (ifail.ne.0) then
          write( *,75) kappa
          write(11,75) kappa
75        format(2x,'Warning!!  Failed to find eigenvalue.',
     1    '  Kappa =',i2)
        write(*,*) m1,m2,m3,imax
          stop
        endif
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        ee0=0.5d0*p(ii+1)
        dd0=(p(ii+5)**2+q(ii+5)**2)/p(ii+2)
        de=dabs((ee0-ee1)/(ee0+ee1+eps))
        dn=dabs((dd0-dd1)/(dd0+dd1+eps))
        del=de
        if (dn.gt.de) del=dn
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        ee1=ee0
        dd1=dd0
        if (niter.ge.3) then
          if (del.gt.1.d-10.or.niter.le.2) goto 800
        else
          goto 800
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine correct(imin,imax,tn,p,q,cp,cq)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine computes the correction term to
c       the finite-difference scheme for solving the
c       radial equations. The procedure uses the radial wave
c       function P,Q obtained in the previous iteration
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /h/h
        real*8 p(maxii),q(maxii),cp(maxii),cq(maxii)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        s=tn/3.d0
        t=tn/8.d0
        d=tn/120.d0
        i0=imin+2
        im=imax-3
        do 10 i=i0,im
        j =i+1
        j1=j+1
        j2=j1+1
        i1=i-1
        i2=i1-1
        tp=s*(p(j)-p(i))-t*(p(j1)-p(i1))+d*(p(j2)-p(i2))
        tq=s*(q(j)-q(i))-t*(q(j1)-q(i1))+d*(q(j2)-q(i2))
c        tp=tn/12.d0*(-p(j1)+3.d0*p(j)-3.d0*p(i)+p(i1))
c        tq=tn/12.d0*(-q(j1)+3.d0*q(j)-3.d0*q(i)+q(i1))
        cp(i)=cp(i)+tp
        cq(i)=cq(i)+tq
c       - - - - - - - - - - - - - - - - - - - - - - - - -
10      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        i=imin+1
        j =i+1
        j1=j+1
        j2=j1+1
        j3=j2+1
        i1=i-1
        tp=d*(9.d0*p(i1)-25.d0*p(i)+20.d0*p(j)-
     1  5.d0*p(j2)+p(j3))
        tq=d*(9.d0*q(i1)-25.d0*q(i)+20.d0*q(j)-
     1  5.d0*q(j2)+q(j3))
        cp(i)=cp(i)+tp
        cq(i)=cq(i)+tq
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        i=imin
        j =i+1
        j1=j+1
        j2=j1+1
        j3=j2+1
        j4=j3+1
        tp=d*(29.d0*p(i)-115.d0*p(j)+180.d0*p(j1)-
     1  140.d0*p(j2)+55.d0*p(j3)-9.d0*p(j4))
        tq=d*(29.d0*q(i)-115.d0*q(j)+180.d0*q(j1)-
     1  140.d0*q(j2)+55.d0*q(j3)-9.d0*q(j4))
        cp(i)=cp(i)+tp
        cq(i)=cq(i)+tq
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        i=imax-2
        i1=i+1
        i2=i1+1
        j =i-1
        j1=j-1
        j2=j1-1
        tp=d*(9.d0*p(i2)-25.d0*p(i1)+20.d0*p(i)-
     1  5.d0*p(j1)+p(j2))
        tq=d*(9.d0*q(i2)-25.d0*q(i1)+20.d0*q(i)-
     1  5.d0*q(j1)+q(j2))
        cp(i)=cp(i)-tp
        cq(i)=cq(i)-tq
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        i=imax-1
        i1=i+1
        j =i-1
        j1=j-1
        j2=j1-1
        j3=j2-1
        tp=d*(29.d0*p(i1)-115.d0*p(i)+180.d0*p(j)-
     1  140.d0*p(j1)+55.d0*p(j2)-9.d0*p(j3))
        tq=d*(29.d0*q(i1)-115.d0*q(i)+180.d0*q(j)-
     1  140.d0*q(j1)+55.d0*q(j2)-9.d0*q(j3))
        cp(i)=cp(i)-tp
        cq(i)=cq(i)-tq
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine origin(kappa,p,q)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /ii/ii/cl/cl/z/z
        common /unuc/unuc(maxii)/y/y(maxii)/r/r(maxii)
        common /nmax/nmax
        real*8 p(maxii),q(maxii)
        dimension yy(20),p1(20),q1(20),vp(20),vq(20)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine calculates the coefficients of the
c       Taylor series expansion of the wavefunctions P,Q,
c       around the origin (r=0) and stores them in the
c       arrays p(ii+5+n), q(ii+5+n) where n=0,...nmax.
c       and ii - number of the grid points.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        c1=0.01d0
        e=p(ii+1)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        yy(1)=0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        ny=y(ii+4)+c1
        do n=0,nmax
          i=n+1
          yy(i+ny)= y(ii+n+5)/r(1)**n
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        vp(1)=unuc(ii+5)-z
        vq(1)=unuc(ii+5)-z
        do m=1,nmax
          i=ii+5+m
          ui=unuc(i)/r(1)**m
          vp(m+1)=yy(m+1)+ui
          vq(m+1)=yy(m+1)+ui
        enddo
c
        vp(2)=vp(2)+(0.5d0*e-2*cl*cl)
        vq(2)=vq(2)+0.5d0*e
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        gam=dsqrt(kappa**2-(vp(1)/cl)**2)
        t=iabs(kappa)
        s=gam+t
        d=dsqrt((p(ii+5)**2+q(ii+5)**2)*0.5d0*s/
     1  (t*p(ii+2)))
        t=d*dabs(vp(1))/(s*cl)
        if (kappa.gt.0) goto 350
        p1(1)=d
        q1(1)=t
        goto 360
350     p1(1)=t
        q1(1)=-d
360     do n=1,nmax
          i=n+1
          d=n*(n+2*gam)*cl*cl
          t1=(n+gam-kappa)*cl
          t2=(n+gam+kappa)*cl
          s1=0.d0
          s2=0.d0
          do m=1,n
            j=m+1
            s1=s1+p1(i-m)*vq(j)
            s2=s2+q1(i-m)*vp(j)
          enddo
          p1(i)=(-vp(1)*s1+t1*s2)/d
          q1(i)=(-vq(1)*s2-t2*s1)/d
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        p(ii+4)=gam
        q(ii+4)=gam
        do n=0,nmax
          i=n+1
          p(ii+5+n)=p1(i)*r(1)**n
          q(ii+5+n)=q1(i)*r(1)**n
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
1000    return
        end
