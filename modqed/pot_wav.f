      subroutine locsep(i,nr,ur,vr)
      implicit real*8 (a-h,o-z)
      include 'qedmod.inc'
      common /ii/ii
      real*8 ur(nr), vr(nr), vt(maxii)
      common /kk_se/kk_se(6)
      common /r/r(maxii)/v/v(maxii)/rho/rho(maxii)
      common /ns_proj/ns_proj/num_kappa/num_kappa
      common /vsemi_loc/vsemi_loc(maxii,5)

      il = i+1
      if (il .gt. 0 .and. il .le. num_kappa) then
         do i=1,ii
            vt(i)=vsemi_loc(i,il)
         enddo
         call uvip3p(3, ii, rho, vt, nr, ur, vr)
      else
         do i=1,nr
            vr(i)=0d0
         enddo
      endif
      end
      
c     =================================================
        subroutine se_pot_wav(n,kappa,iw,r,pp,qq,cp,cq)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine computes the functions cp(i),cq(i),
c       which are the result of model potentail acting on
c       the wave functions p(i) and q(i) respectively.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /z/z/cl/cl
        common /ns_proj/ns_proj/num_kappa/num_kappa
        common /kk_proj/kk_proj(maxns)
        common /uq/uq(maxii,maxns)/up/up(maxii,maxns)
        common /kk_se/kk_se(6)
        common /vsemi_loc/vsemi_loc(maxii,5)
        common /Dab/Dab(maxns,maxns)
        common /v/v(maxii)
        common /ii/ii
        real*8 r(maxii)
        real*8 pp(maxii),qq(maxii)
        real*8 cp(maxii),cq(maxii)
        real*8 dproj(maxns)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        gam=dsqrt(kappa**2-(z/cl)**2)
        pp(ii+3)=ii
        qq(ii+3)=ii
        pp(ii+4)=gam
        qq(ii+4)=gam
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Local potential
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        il=0
        do i=1,num_kappa
           if (kappa.eq.kk_se(i)) then
              il=i
           endif
        enddo
       
        if (il.gt.0 .and. iw .eq. 0) then
           do i=1,ii
              cp(i)=vsemi_loc(i,il)*pp(i)/r(i)
              cq(i)=vsemi_loc(i,il)*qq(i)/r(i)
           enddo
        else
           do i=1,maxii
              cp(i)=0.d0
              cq(i)=0.d0
           enddo
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       non-local
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do nj=1, ns_proj
          dproj(nj)=0.d0
          if (kappa.eq.kk_proj(nj)) then
            dproj(nj)=tint(-1,pp,qq,up(1,nj),uq(1,nj),r,v)
          endif
        enddo
c
        do nj=1, ns_proj
          dnonloc=0.d0
          do nk=1, ns_proj
            dnonloc=dnonloc+Dab(nj,nk)*dproj(nk)
          enddo
          do i=1,ii
            cp(i)=cp(i)+dnonloc*up(i,nj)/r(i)
            cq(i)=cq(i)+dnonloc*uq(i,nj)/r(i)
          enddo
        enddo
c
        cp(ii+3)=ii
        cq(ii+3)=ii
        if (il .gt. 0) then
           cp(ii+4)=vsemi_loc(ii+4,il)+gam
           cq(ii+4)=vsemi_loc(ii+4,il)+gam
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
