c       =================================================
        subroutine init_se
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine initializes some data needed to
c       construct self-energy QED potential. In partuclar,
c       it computes the diagonal and non-diagonal SE matrix
c       elements by means of subroutine FSE_dat and stores them
c       in the array qed_matrix.
c
c       The quantum numbers listed below define the quantum
c       numbers of the projected wave functions used by the
c       nonlocal part of SE model potential.
c
c       Parameters of the projected wave functions:
c       ns_proj     - number of atomic shells
c       maxns       - maximal number of atomic shells
c                     (see file 'qedmod.inc')
c       ni          - shell number (ni=1,,,ns)
c       nn_proj(ni) - principal quantum number
c       ll_proj(ni) - orbital quantum number
c       jj_proj(ni) - angular quantum number
c       kk_proj(ni) - relativistic quantum numbers
c       num_kappa   - total number of different relativistic
c                     angular quantum numbers
c
c       z    .... Nuclear charge
c       cl   ..... Speed of light
c       knucl  ..... Nuclear model (0 - point, 1 - Fermi model,
c                    2 - uniform sphere)
c
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /cl/cl/z/z/ns_proj/ns_proj/knucl/knucl
        common /num_kappa/num_kappa
        common /nn_proj/nn_proj(maxns)/ll_proj/ll_proj(maxns)
     &  /jj_proj/jj_proj(maxns)/kk_proj/kk_proj(maxns)
        common /qed_matrix/qed_matrix(maxns,maxns)
        common /kk_se/kk_se(6)/ll_se/ll_se(6)/nn_se/nn_se(6)
        common /se_loc/se_loc(6)/num_loc/num_loc(-3:3)
        common /nshell/nshell(-5:5,10)
        common /nn_max/nn_max(0:3)
        parameter (np=5)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (pi=3.1415926535897932385d0)
        alpha=1.d0/cl
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Default value of the max principal quantum number
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        nma=3
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do l=0,3
          if (nma.gt.0) then
            nn_max(l)=nma
          else
            nn_max(l)=iabs(nma)+l
          endif
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        kk_se(1)=-1
        kk_se(2)=1
        kk_se(3)=-2
        kk_se(4)=2
        kk_se(5)=-3
        kk_se(6)=3
        num_kappa=5
        ni=0
        do il=1,num_kappa
          kappa=kk_se(il)
          l=(iabs(2*kappa+1)-1)/2
          ll_se(il)=l
          nn_se(il)=l+1
          do n=l+1,nn_max(l)
            ni=ni+1
            kk_proj(ni)=kappa
            ll_proj(ni)=l
            jj_proj(ni)=2*iabs(kappa)-1
            nn_proj(ni)=n
          enddo
        enddo
        ns_proj=ni
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        num_loc(-3)=5
        num_loc(-2)=3
        num_loc(-1)=1
        num_loc( 0)=0
        num_loc( 1)=2
        num_loc( 2)=4
        num_loc( 3)=6
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do kappa=-5,5
        do n=1,10
          nshell(kappa,n)=0
        enddo
        enddo
        do ni=1,ns_proj
          kappa=kk_proj(ni)
          n=nn_proj(ni)
          nshell(kappa,n)=ni
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation of the qed_matrix
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call se_pnt_stor()
        call se_fn_stor()
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do ni=1,ns_proj
        do nj=ni,ns_proj
          qed_matrix(ni,nj)=0.d0
          kappa=kk_proj(ni)
          nn1=nn_proj(ni)
          nn2=nn_proj(nj)
          if (kk_proj(nj).eq.kappa) then
            iz=z+0.5d0
            call FSEdat(kappa,nn1,nn2,iz,FSE_pnt,FSE_ext)
            dn1=nn1
            dn2=nn2
            dn=dsqrt(dn1*dn2)
            coef=alpha/pi*(z*alpha)**4/dn**3*cl**2
            if (knucl.eq.0) then
              qed_matrix(ni,nj)=fse_pnt*coef
              qed_matrix(nj,ni)=fse_pnt*coef
            else
              qed_matrix(ni,nj)=fse_ext*coef
              qed_matrix(nj,ni)=fse_ext*coef
            endif
            do il=1,num_kappa
              n=nn_se(il)
              if (kappa.eq.kk_se(il).and.nn1.eq.n.and.nn2.eq.n) then
                se_loc(il)=qed_matrix(ni,nj)
              endif
            enddo
          endif
        enddo
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine local_se_pot
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine computes the local part Vloc(r)
c       of the Self-Energy (SE) QED potential.
c       Local potential multiplied on r is store in the
c       array 'vsemi_loc(i,kappa)' for each value of the
c       relativistic quantum number kappa, where 'i' is
c       a number of grid point. 
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /cl/cl/z/z/r1/r1/ii/ii/h/h
        common /p/p(maxii)/q/q(maxii)/y/y(maxii)
        common /r/r(maxii)/v/v(maxii)
        common /ee_loc/ee_loc(maxns)
        common /ns_proj/ns_proj/num_kappa/num_kappa        
        common /nn_proj/nn_proj(maxns)/ll_proj/ll_proj(maxns)
     &  /jj_proj/jj_proj(maxns)/kk_proj/kk_proj(maxns)
        common /vloc/vloc(maxii)/ww/ww(maxii)
        common /vsemi_loc/vsemi_loc(maxii,6)
        common /wsemi_loc/wsemi_loc(maxii,6)
        common /kk_se/kk_se(6)/ll_se/ll_se(6)/nn_se/nn_se(6)
        common /se_loc/se_loc(6)
        common /nshell/nshell(-5:5,10)        
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (pi=3.1415926535897932385d0)
        alpha=1.d0/cl
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do il=1,num_kappa
          do i=1,maxii
            vsemi_loc(i,il)=0.d0
            wsemi_loc(i,il)=0.d0
          enddo
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 10 il=1,num_kappa
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        kappa=kk_se(il)
        l=ll_se(il)
        n=nn_se(il)
        ni=nshell(kappa,n)
        if (ni.eq.0) goto 10
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,ii
          ri=r(i)
          vloc(i)=dexp(-ri*cl)*ri
          ww(i)=  dexp(-ri*cl)*ri
        enddo
        do i=ii,1,-1
          if (dabs(vloc(i)).gt.1.d-40) goto 200
          imax=i-1
          vloc(i)=0.d0
        enddo
c
 200    vloc(ii+3)=imax
        vloc(ii+4)=1.d0
        ww(ii+3)=ii
        ww(ii+4)=1.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call read_func(ni,p,q,2)
        de=sint(vloc,p,q,p,q,r,v)
        vloc0=se_loc(il)/de
        de=sint(ww,p,q,p,q,r,v)
        w0=1.d0/de
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,ii
          vloc(i)=vloc0*vloc(i)
          ww(i)=w0*ww(i)
        enddo
        do i=1,maxii
          vsemi_loc(i,il)=vloc(i)
          wsemi_loc(i,il)=ww(i)
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
10      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine nonlocal_se_pot
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine computes the nonlocal (separable)
c       part of the Self-Energy (SE) potential. This
c       potential has the form:
c       V_nonloc = \sum_{a,b} |a> D_ab <b|, where a,b
c       are projected wave functions and
c       Dab=Sab^{-1} * QED_matr * Sab^{-1}
c       The projected wave functions and matrix D_ab are
c       stored in the external file
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /cl/cl/z/z/r1/r1/ns_proj/ns_proj
        common /p/p(maxii)/q/q(maxii)/a/a(maxii)/b/b(maxii)
     &  /cp/cp(maxii)/cq/cq(maxii)/r/r(maxii)/v/v(maxii)
        common /nn_proj/nn_proj(maxns)/ll_proj/ll_proj(maxns)
     &  /jj_proj/jj_proj(maxns)/kk_proj/kk_proj(maxns)
        common /ii/ii/h/h/al/al/bt/bt
        common /uq/uq(maxii,maxns)/up/up(maxii,maxns)
        common /Dab/Dab(maxns,maxns)/Gab/Gab(maxns,maxns)
     &  /sab/sab(maxns,maxns)/vab/vab(maxns,maxns)
        common /qed_matrix/qed_matrix(maxns,maxns)
        common /diag/diag(maxns)/diag_qed/diag_qed(maxns)
     &  /diag_inv/diag_inv(maxns)
        common /de_qed_loc0/de_qed_loc0(maxns)/de_qed0/de_qed0(maxns)
        common /vloc/vloc(maxii)/ww/ww(maxii)
        common /vsemi_loc/vsemi_loc(maxii,6)
        common /wsemi_loc/wsemi_loc(maxii,6)
        common /num_loc/num_loc(-3:3)
        common /nshell/nshell(-5:5,10)
        real*8 p1(maxii),q1(maxii)
        real*8 p2(maxii),q2(maxii)
        real*8 p3(maxii),q3(maxii)
        integer iwrk1(maxns),iwrk2(maxns)
        parameter (np=5)
        character*1 let(11)
        data
     1  let /'s','p','d','f','g','h','i','k','l','m','n'/
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (pi=3.1415926535897932385d0)
        alpha=1.d0/cl
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Functions up, uq
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do ni=1,ns_proj
          kappa=kk_proj(ni)
          n=nn_proj(ni)
          l=ll_proj(ni)
          call read_func(ni,p,q,2)
          il=num_loc(kappa)
          do i=1,maxii
            vloc(i)=vsemi_loc(i,il)       
            ww(i)=wsemi_loc(i,il)
          enddo
          do i=1,maxii
            a(i)=p(i)
            b(i)=q(i)
          enddo
          do i=1,ii
            a(i)=ww(i)*p(i)
            b(i)=ww(i)*q(i)
          enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          k=nn_proj(ni)-ll_proj(ni)
            do i=1,ii
              ppi=p(i)
              qqi=q(i)
              dw=dexp(-r(i)*z*2)*r(i)
              factor_p=(1.d0-(-1.d0)**k)/2.d0
              factor_q=(1.d0+(-1.d0)**k)/2.d0
              a(i)= factor_p*dw*ppi
              b(i)= factor_q*dw*qqi
            enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          gam=p(ii+4)+ww(ii+4)
          a(ii+4)=gam
          b(ii+4)=gam
          a(ii+3)=ii
          b(ii+3)=ii
c
          dn=tint(-1,p,q,a,b,r,v)
          do i=1,ii
            up(i,ni)=a(i)/dn
            uq(i,ni)=b(i)/dn
          enddo
          imax=ii
          do i=ii,1,-1
            if (dabs(up(i,ni))+dabs(uq(i,ni)).gt.1.d-40) goto 200
            imax=i-1
            up(i,ni)=0.d0
            uq(i,ni)=0.d0
          enddo
 200      up(ii+3,ni)=imax
          uq(ii+3,ni)=imax
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Matrix Sab and Vab
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 10 ni=1,ns_proj
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        kappa=kk_proj(ni)
        il=num_loc(kappa)
        do i=1,maxii
          vloc(i)=vsemi_loc(i,il)
          ww(i)=wsemi_loc(i,il)
        enddo
        call read_func(ni,p,q,2)
        do nj=1,ns_proj
          sab(ni,nj)=0.d0
          vab(ni,nj)=0.d0
          if (kk_proj(ni).eq.kk_proj(nj)) then
            call read_func(nj,cp,cq,2)
            de=sint(vloc,p,q,cp,cq,r,v)
            vab(ni,nj)=de
            do i=1,maxii
              a(i)=up(i,nj)
              b(i)=uq(i,nj)
            enddo
            sab(ni,nj)=tint(-1,p,q,a,b,r,v)
          endif
        enddo
        de_qed_loc0(ni)=vab(ni,ni)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
10      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Matrix Inversion
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do ni=1,ns_proj
        do nj=1,ns_proj
          Gab(ni,nj)=Sab(ni,nj)
        enddo
        enddo
        call dminv(Gab,ns_proj,maxns,det,iwrk1,iwrk2)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Inversion test
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        del=0.d0
        do ni=1,ns_proj
        do nj=1,ns_proj
          aij=0.d0
          do nk=1,ns_proj
            aij=aij+Gab(ni,nk)*Sab(nk,nj)
          enddo
          if (ni.eq.nj) aij=aij-1.d0
          del=del+dabs(aij)/ns_proj
        enddo
        enddo
        if (abs(del).lt.1.d-10) then
c           write( *,'(/2x,a,e12.2,a)') 'Inversion test:',del,"  OK"
c           write(11,'(/2x,a,e12.2,a)') 'Inversion test:',del,"  OK"
        else
c           write( *,'(/2x,a,e12.2,a)') 'Inversion test:',del
c           write(11,'(/2x,a,e12.2,a)') 'Inversion test:',del
        end if
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       D-matrix
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call dab_matrix(qed_matrix)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do ni=1, ns_proj
          delt=0.d0
          do nj=1, ns_proj
          do nk=1, ns_proj
            delt=delt+Sab(ni,nj)*Dab(nj,nk)*Sab(ni,nk)
          enddo
          enddo
          de_qed0(ni)=Vab(ni,ni)+delt
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
 1000   return
        end
c       =================================================
        subroutine dab_matrix(qed_matrix)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       The subroutine computes Dab matrix 
c       Dab=Sab^{-1} * QED_matr * Sab^{-1}
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /ns_proj/ns_proj
        common /nn_proj/nn_proj(maxns)/ll_proj/ll_proj(maxns)
     &  /jj_proj/jj_proj(maxns)/kk_proj/kk_proj(maxns)
        common /Dab/Dab(maxns,maxns)/Gab/Gab(maxns,maxns)
     &  /sab/sab(maxns,maxns)/vab/vab(maxns,maxns)
        real*8 qed_matrix(maxns,maxns)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       D-matrix
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do ni=1,ns_proj
        do nj=1,ns_proj
          dab(ni,nj)=0.d0
        enddo
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 10 ni=1,ns_proj
        do 20 nj=1,ns_proj
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        Dab(ni,nj)=0.d0
        if (kk_proj(nj).ne.kk_proj(ni)) goto 20
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 30 nk=1,ns_proj
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (kk_proj(nk).ne.kk_proj(ni)) goto 30
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 40 nl=1,ns_proj
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (kk_proj(nl).ne.kk_proj(ni)) goto 40
        qkl=qed_matrix(nk,nl)
        qlk=qed_matrix(nl,nk)
        dkl=0.5d0*(qkl+qlk)-vab(nk,nl)
        Dab(ni,nj)=Dab(ni,nj)+Gab(ni,nk)*dkl*Gab(nj,nl)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
40      continue
30      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
20      continue
10      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine write_pot(output_filename)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine store the radial grid, Uehling and
c       Wichmann-Kroll potentials, local potential and
c       unlocal part of the model potential in the
c       output_filename.
c       The name of this file is input parameter.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine store Dab matrix 
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /cl/cl/z/z/ns_proj/ns_proj/num_kappa/num_kappa
        common /ii/ii/h/h/al/al/bt/bt
        common /kk_se/kk_se(6)
        common /p/p(maxii)/q/q(maxii)/r/r(maxii)/v/v(maxii)
        common /nn_proj/nn_proj(maxns)/ll_proj/ll_proj(maxns)
     &  /jj_proj/jj_proj(maxns)/kk_proj/kk_proj(maxns)
        common /vloc/vloc(maxii)/unuc/unuc(maxii)
        common /uehl/uehl(maxii)/v_wk/v_wk(maxii)
        common /knucl/knucl/inucl/inucl/rnucl/rnucl
        common /num_loc/num_loc(-3:3)
        common /vsemi_loc/vsemi_loc(maxii,6)
        common /uq/uq(maxii,maxns)/up/up(maxii,maxns)
        common /Dab/Dab(maxns,maxns)
        character*16 output_filename
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c        open (unit=13,file=output_filename,status='unknown')
c
c        write(13,'(i22,16x,a)') ii,'ii'
c        write(13,'(f22.14,16x,a)') z,'z'
c
c        write(13,'(3x,5(i5))') (kk_se(il),il=1,num_kappa)
c
c        write(13,'(16x,a,7x,5(8x,a,i2))') 'r',
c     &  ('vloc  kappa=',kk_se(il),il=1,num_kappa)
c        do i=1,ii
c          write (13,15) i,r(i), (vsemi_loc(i,il),il=1,num_kappa)
c        enddo
c15      format(i5,7d22.12)
c        write(13,*)
c        write (13,25) (vsemi_loc(ii+3,il),il=1,num_kappa)
c        write (13,25) (vsemi_loc(ii+4,il),il=1,num_kappa)
c25      format(27x,6d22.12)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c        write(13,'(15x,a,19x,a,15x,a,18x,a)') 'r','v_uehling',
c     &  'v_wk','unuc'
c        do i=1,ii
c          write (13,35) i,r(i),uehl(i),v_wk(i),unuc(i)
c        enddo
c35      format(i5,4d22.12)
c        write (13,*)
c        write (13,45) uehl(ii+3),v_wk(ii+3),unuc(ii+3)
c        write (13,45) uehl(ii+4),v_wk(ii+4),unuc(ii+4)
c45      format(27x,3d22.12)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c        write(13,'(/2x,i6,16x,a)') ns_proj,'ns_proj'
c        do ni=1, ns_proj
c          write(13,'(5i4,16x,a)') ni,nn_proj(ni),ll_proj(ni),
c     &    jj_proj(ni),kk_proj(ni),'ni, n, l, j, k'
c        enddo
c        write(13,'(/22x,a)') 'Matrix Dab'
c        do ni=1, ns_proj
c        do nj=1, ns_proj
c         write(13,85) ni,nj, Dab(ni,nj)
c        enddo
c        enddo
c85      format(2i6, e23.13)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c        do ni=1,ns_proj
c         write(13,'(2x,a,i2,11x,a,18x,a)') 'ni=',ni,'u_p','u_q'
c          do i=1, ii
c            write(13,95) i,up(i,ni),uq(i,ni)
c95        format(i5,2x,2e23.13)
c          enddo
c          write(13,*)
c          write(13,105) up(ii+3,ni), uq(ii+3,ni)
c          write(13,105) up(ii+4,ni), uq(ii+4,ni)
c105       format(7x,2e23.13)
c        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c        write(13,*)
c       write(13,'(i6,20x,a)') knucl,'knucl'
c        write(13,'(i6,20x,a)') inucl,'inucl'
c        write(13,'(d22.14,4x,a)') rnucl,'rnucl'
c        write(13,'(d22.14,4x,a)') al,'al'
c        write(13,'(d22.14,4x,a)') bt,'bt'
c        write(13,'(d22.14,4x,a)') h,'h'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c        close(unit=13)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine dminv(a,n,nmax,d,l,m)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Invertion of arbitrary matrix:  C=A**(-1),
c       using Gauss-Gordan method.
c       This subroutine is a sligtly modified version of the
c       minv subroutine from the IBM Application Program:
c       Scientific Subroutine Package (SSP)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       l,m - work arrays (dimension - n)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dimension a(nmax,nmax),l(nmax),m(nmax)
c       double precision a,d,biga,hold
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Searching maximal element.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        d=1.d0
        do 80 k=1,n
        l(k)=k
        m(k)=k
        biga=a(k,k)
        do 20 j=k,n
        do 20 i=k,n
10      if (dabs(biga)-dabs(a(i,j))) 15,20,20
15      biga=a(i,j)
        l(k)=i
        m(k)=j
20      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Changing the strings
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        j=l(k)
        if (j-k) 35,35,25
25      do 30 i=1,n
        hold=-a(k,i)
        a(k,i)=a(j,i)
30      a(j,i)=hold
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Changing the columns
c       - - - - - - - - - - - - - - - - - - - - - - - - -
35      i=m(k)
        if (i-k) 45,45,38
38      jp=n*(i-1)
        do 40 j=1,n
        hold=-a(j,k)
        a(j,k)=a(j,i)
40      a(j,i)=hold
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Deviding column on the leader element.
c       The leader element is in 'biga'.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
45      if (biga) 48,46,48
46      d=0.d0
        return
48      biga=-biga
        do 55 i=1,n
        if (i-k) 50,55,50
50      a(i,k)=a(i,k)/biga
55      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        biga=-biga
        do 65 i=1,n
        hold=a(i,k)
        ij=i-n
        do 65 j=1,n
        ij=ij+n
        if (i-k) 60,65,60
60      if (j-k) 62,65,62
62      kj=ij-i+k
        a(i,j)=hold*a(k,j)+a(i,j)
65      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Deviding the string on the leader element.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        kj=k-n
        do 75 j=1,n
        kj=kj+n
        if (j-k) 70,75,70
70      a(k,j)=a(k,j)/biga
75      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Composition of the leader elements.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        d=d*biga
        a(k,k)=1.d0/biga
80      continue
        k=n
100     k=k-1
        if (k) 150,150,105
105     i=l(k)
        if (i-k) 120,120,108
108     do 110 j=1,n
        hold=a(j,k)
        a(j,k)=-a(j,i)
110     a(j,i)=hold
120     j=m(k)
        if (j-k) 100,100,125
125     do 130 i=1,n
        hold=a(k,i)
        a(k,i)=-a(j,i)
130     a(j,i)=hold
        goto 100
c       - - - - - - - - - - - - - - - - - - - - - - - - -
150     return
        end
