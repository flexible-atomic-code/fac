c       =================================================
        subroutine read_pot(data_filename)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine read the radial grid, Uehling and
c       Wichmann-Kroll potentials, local potential and
c       unlocal parts of the model potential from data_filename.
c       The name of this file is input parameter.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
        character*(*) data_filename
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /z/z
        common /ii/ii/h/h/r1/r1/r2/r2/al/al/bt/bt
        common /r/r(maxii)/v/v(maxii)
        common /knucl/knucl/inucl/inucl/rnucl/rnucl
        common /ns_proj/ns_proj/num_kappa/num_kappa
        common /nn_proj/nn_proj(maxns)/ll_proj/ll_proj(maxns)
     &  /jj_proj/jj_proj(maxns)/kk_proj/kk_proj(maxns)
        common /uq/uq(maxii,maxns)/up/up(maxii,maxns)
        common /vloc/vloc(maxii)/unuc/unuc(maxii)
        common /uehl/uehl(maxii)/v_wk/v_wk(maxii)
        common /kk_se/kk_se(6)
        common /vsemi_loc/vsemi_loc(maxii,6)
        common /Dab/Dab(maxns,maxns)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        num_kappa=5
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Read QED potential
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        write(*,*) " Reading model QED potential from ",data_filename
        open(unit=13,file=data_filename,status='old',ERR=333)
        go to 334
 333    write(*,'(2x,a,2x,a)') "Error in opening file",data_filename
        stop
 334    continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        read(13,*) ii
        read(13,*) z
          read (13,*) (kk_se(il),il=1,num_kappa)
        read(13,*)
        do i=1,ii
          read (13,*) i1,r(i),(vsemi_loc(i,il),il=1,num_kappa)
        enddo
        read (13,*)
        read (13,*) (vsemi_loc(ii+3,il),il=1,num_kappa)
        read (13,*) (vsemi_loc(ii+4,il),il=1,num_kappa)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        read(13,*)
        do i=1,ii
          read (13,*) i1,ri,uehl(i),v_wk(i),unuc(i)
        enddo
        read (13,*)
        read (13,*) uehl(ii+3),v_wk(ii+3),unuc(ii+3)
        read (13,*) uehl(ii+4),v_wk(ii+4),unuc(ii+4)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        read(13,*)
        read(13,*) ns_proj
        do ni=1, ns_proj
          read(13,*) ni1,nn_proj(ni),ll_proj(ni),jj_proj(ni),
     &    kk_proj(ni)
        enddo
c
        read(13,*)
        read(13,*)
        do ni=1, ns_proj
        do nj=1, ns_proj
          read(13,*) ni1,nj1, Dab(ni,nj)
        enddo
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do ni=1, ns_proj
          read(13,*)
          do i=1,ii
            read(13,*) i1,up(i,ni),uq(i,ni)
          enddo
          read (13,*)
          read (13,*) up(ii+3,ni), uq(ii+3,ni)
          read (13,*) up(ii+4,ni), uq(ii+4,ni)
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        read(13,*)
        read(13,*) knucl
        read(13,*) inucl
        read(13,*) rnucl
        read(13,*) al
        read(13,*) bt
        read(13,*) h
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        r1=r(1)
        r2=r(ii)
        do i=1,ii
          v(i)=r(i)/(al*r(i)+bt)
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        close(unit=13)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
