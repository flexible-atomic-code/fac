c =================================================
        subroutine iniqed(z0, nmax0, knucl0, rms)
c - - - - - - - - - - - - - - - - - - - - - - - - -
c
c  This is the main part of the potgen routine. potgen
c  calculates the auxiliar data necessary for the evaluation 
c  of the matrix elements of the model QED operator and 
c  stores them into an external file on the disk. 
c
c  The model QED operator has the form: 
c  Vmod= Vloc+ \sum_{a,b} |a> Dab_{a,b} <b|,
c  where Vloc is the local part of the operator (some potential
c  which is stored on the radial grid), |a> and |b> are the
c  projected wave function and D_ab is the QED matrix. 
c  Vloc, |a,b>, and Dab are calculated and stored in the 
c  external file.
c
c  Input parameters:
c  z -- nuclear charge
c  knucl -- nuclear model 
c        knucl=0 - the point nucleus;
c        knucl=1 - the extended nucleus (Fermi model);
c  output_filename -- the name of the output file in which
c                     QED potential is stored
c
c - - - - - - - - - - - - - - - - - - - - - - - - -
c  The code uses the following subroutines:
c
c  grid        - Subroutine initializes nonuniform semi-logarithmic
c                radial grid.
c
c  init_se -     This subroutine initializes auxiliar data needed to
c                construct self-energy potential. It also
c                computes the diagonal and non-diagonal SE matrix
c                elements using subroutine FSE_dat.
c
c  FSE_dat  -    This subroutine obtaines diagonal and non-diagonal
c                self-energy matrix elements, by interpolating the
c                prestored tabulated data from the paper
c                V.M. Shabaev, I.I. Tupitsyn and V. A. Yerokhin,
c                PRA, 88, 012513 (2013).
c
c  local_se_pot - This subroutine computes the local part
c                of the self-energy potential.
c                Local potential multiplied on r is stored in the
c                array 'vsemi_loc(i,kappa)' for each value of the
c                relativistic quantum number kappa, where 'i' numerates
c                the grid points. 
c
c  nonlocal_se_pot -This subroutine computes the nonlocal
c                (separable) part of the self-energy potential.
c                This potential has the form:
c                V_nonloc = \sum_{a,b} |a> D_ab <b|, where a,b
c                are projected wave functions.
c                The projected wave functions and matrix D_ab are
c                stored in the external file
c
c  uehling  -    This subroutine computes Uehling vacuum-polarization
c                potential on a grid.
C                The subroutine uses accurate approximation method,
c                described in the paper:
c                L.W.Fullerton and G.A. Ringer, PRA 13, 1283 (1976)
c
c  wk  -         This subroutine computes Wichmann-Kroll point nuclear
c                vacuum-polarization potential on a grid and.
C                The subroutine use the analytical approximation formulas,
c                obtained in the paper:
c                A G Fainshtein, N L Manakov and A A Nekipelov,
c                J.Phys.E: At.Mol.Opt.Phys. v.23 (1990) 559-569. 
c
c  nucl  -       This subroutine computes the nuclear potential.
c
c  Dirac -       This subroutine computes radial Dirac
c                wavefunctions for an arbitrary nuclear charge distribution.
c
c  tint  -       This routine performs the numerical integration in 
c                the matrix elements <p,q|r^k|a,b>
c
c  sint  -       This routine performs the numerical integration in 
c                the matrix elements of the potential v(r)/r, 
c                <p,q|v(r)/r|a,b>.
c
c- - - - - - - - - - - - - - - - - - - - - - - - -
c
c       Common-block variables:
c
c       z    .... Nuclear charge
c       cl    ..... Speed of light 
c       ii    ..... Number of radial grid points (ii=maxii-20), where
c                   maxii is given in the include file
c       r(i)  ..... non-unform radial grid (i=1,ii)
c       r1    ..... the first grid point 
c       r2    ..... the last grid point (the size of the box)
c       h, al, bt  ..... grid parameters
c       p(i)  ..... Large component of the wavefunction (p(r)= G(r))
c       q(i)  ..... Small component of the wavefunction (q(r)=-F(r))
c       y(i)  ..... Coulomb (or external) potential multiplied
c       unuc(i) ... Nonpoint part of the Nuclear potential
c                   multiplied by r(i)
c       uehl(i) ... Uehling potential multiplied by r(i).
c       wk(i)   ... Wichmann-Kroll potential multiplied by r(i).
c
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
c     - - - - - - - - - - - - - - - - - - - - - - - - -
        real*8 cp(maxii),cq(maxii)
        common /cl/cl/z/z
        common /ii/ii/h/h/r1/r1/r2/r2/al/al/bt/bt
        common /p/p(maxii)/q/q(maxii)/r/r(maxii)/v/v(maxii)
        common /uehl/uehl(maxii)/v_wk/v_wk(maxii)
        common /nmax/nmax/knucl/knucl/inucl/inucl
     1  /rnucl/rnucl/unuc/unuc(maxii)
        common /niter/niter/mi/m1,m2,m3/nit/nit
        common /ns_proj/ns_proj
        common /nn_proj/nn_proj(maxns)/ll_proj/ll_proj(maxns)
     &  /jj_proj/jj_proj(maxns)/kk_proj/kk_proj(maxns)
        common /de_qed0/de_qed0(maxns)
        common /qed_matrix/qed_matrix(maxns,maxns)
        common /num_kappa/num_kappa
        common /kk_se/kk_se(6)
        common /rfn/rfn(6)
        common /Rrms/Rrms(120)/name_at/name_at(120)/nprms/rms0
        character*16 output_filename
        character*1 let(11)
        data
     1  let /'s','p','d','f','g','h','i','k','l','m','n'/
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       file xx:potgen.res
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c        open(unit=11,file='potgen.res',status='unknown')
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        nmax=nmax0
        z = z0
        knucl = knucl0
        rms0 = rms
        nnm = 5
c        write(*,*) nmax, z, knucl, nmax0, z0, knucl0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c        print *,"input the nuclear charge Z"
c        read(*,*) z
c        if (z.lt.10.or.z.gt.120) stop "Z should be from range [10,120]"
c        print *,"input the nuclear model (0: point nucleus;",
c     >              " 1: extended nucleus)"
c        read(*,*) knucl
c        if ((knucl.lt.0).or.(knucl.gt.1)) then
c           print *,
c     > "knucl is not recognized. Extended nucleus (knucl=1) is used"
c           knucl = 1
c        end if
c        print *,"input the name of the auxiliary data file"
c        read *,output_filename
c        write( *,1)
c        write(11,1)
c 1      format(/2x,"Potgen: calculation of the model QED potential",/2x)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call atom_data
        call segrid
        call nucl(maxii,r,unuc)
        call init_se
        do ik = 1,num_kappa
           kappa = kk_se(ik)
           if (kappa .ne. 0) then
              call dirac(nnm, kappa, p, q)
              do ir=1,ii
                 if (p(ir) .gt. 0 .and. p(ir+1) .le. 0) goto 10
                 if (p(ir) .lt. 0 .and. p(ir+1) .ge. 0) goto 10
              enddo
 10           rfn(ik) = r(ir-1)
c              write(*,*) kappa, ik, ir, rfn(ik)
           endif
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Hydrogen like wave functions
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c        write( *,5)
c        write(11,5)
c5       format(/2x,'Calculating hydrogenic wavefunctions:')
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c        write( *,15)
c        write(11,15)
c15      format(2x,47('=')/10x,'jj',6x,'energy',7x,'Iterations')
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do ni=1,ns_proj
          n=nn_proj(ni)
          l=ll_proj(ni)+1
          kappa=kk_proj(ni)
          call dirac(n,kappa,p,q)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          e=0.5d0*p(ii+1)
          d=p(ii+5)**2+q(ii+5)**2
c          write( *,25) ni,n,let(l),jj_proj(ni),e,niter
c          write(11,25) ni,n,let(l),jj_proj(ni),e,niter
c25        format (i3,i4,a1,i2,'/2',f16.8,i9)
          call write_func(ni,p,q,2)
c          write(*,*) r(1),p(1),q(1),r(50),p(50),q(50)
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c        write( *,35)
c        write(11,35)
c35      format(2x,47('='))
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c	write( *,'(/2x,a)') "Calculating Uehling potential..." 
c        call uehling(maxii,r,unuc,uehl)
c	write( *,'(2x,a)') "Calculating Wichmann-Kroll potential..." 
c        call wk(maxii,r,v_wk)
c	write( *,'(2x,a)') 
c     &  "Calculating local part of self-energy potential..."     
        call local_se_pot
c	write( *,'(2x,a)') 
c     & "Calculating nonlocal part of self-energy potential..."     
        call nonlocal_se_pot
        idebug = 0
        if (idebug .gt. 0) then
c       - - - - - - - - - - - - - - - - - - - - - - - - -
           write( *,45)
c     write(11,45)
 45        format(/2x,'Control check of matrix elements with the ',
     &          'hydrogenic wave functions: ')
           write( *,55)
c     write(11,55)
 55        format(9x,'State',10x,'Uehl',12x,'WK',12x,'SE',
     &          8x,'SE(table data)')
c     - - - - - - - - - - - - - - - - - - - - - - - - -
           pi=3.1415926535897932385d0
           alpha=1.d0/cl
           alz=z/cl
           do ni=1,ns_proj
              n=nn_proj(ni)
              kappa=kk_proj(ni)
              j=2*iabs(kappa)-1
              l=(iabs(2*kappa+1)-1)/2
              coef=alpha/pi*alz**4/n**3*cl**2
              call read_func(ni,p,q,2)
c     de=sint(uehl,p,q,p,q,r,v)
c     dw=sint(v_wk,p,q,p,q,r,v)
c     falz_ue=de/coef
c     falz_wk=dw*cl/(1.d0/pi/n**3 *alz**4)
              call se_pot_wav(n,kappa,0,r,p,q,cp,cq)
              dse0 = tint(0,p,q,cp,cq,r,v)/coef
              falz_se=de_qed0(ni)/coef
              fse_int=qed_matrix(ni,ni)/coef
              write( *,65) ni,n,let(l+1),j,falz_ue,falz_wk,
     &             falz_se,fse_int,dse0,dse0*coef
c     write(11,65) ni,n,let(l+1),j,falz_ue,falz_wk,
c     &    falz_se,fse_int,dse0
 65           format(2x,i3,i5,a1,i2,'/2',6f15.6)
           enddo
c     - - - - - - - - - - - - - - - - - - - - - - - - -
c     write( *,75) output_filename
c     write(11,75) output_filename
c     75     format(/2x, 'Writing model QED potential into the file: ',a)
c     call write_pot(output_filename)
c     - - - - - - - - - - - - - - - - - - - - - - - - -
c     stop
        endif
        end
c       =================================================
        subroutine read_func(ni,v1,v2,nrec)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       reading a vector from the common block
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter(maxff=2*(2*maxns+1))
        common /ff/ff(maxii,maxff)
        integer ni,nrec
        real*8 v1(maxii),v2(maxii)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        nr1=2*ni-1
        nr2=nr1+1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,maxii
          v1(i)=ff(i,nr1)
        enddo
        if (nrec.eq.1) goto 1000
        do i=1,maxii
          v2(i)=ff(i,nr2)
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
1000    return
        end
      
      subroutine radfnd(ik, rn)
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
        common /num_kappa/num_kappa
        common /kk_se/kk_se(6)
        common /rfn/rfn(6)
        rn = 0d0
        il = ik+1
        if (il .gt. 0 .and. il .le. num_kappa) then
           rn = rfn(il)
           return
        endif        
        return
        end
      
c       =================================================
        subroutine write_func(ni,v1,v2,nrec)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       writing a vector to the common block
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter(maxff=2*(2*maxns+1))
        common /ff/ff(maxii,maxff)
        integer ni,nrec
        real*8 v1(maxii),v2(maxii)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        nr1=2*ni-1
        nr2=nr1+1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,maxii
          ff(i,nr1)=v1(i)
        enddo
        if (nrec.eq.1) goto 1000
        do i=1,maxii
          ff(i,nr2)=v2(i)
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
1000    return
        end
c       =================================================
        subroutine exit1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call exit(1)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
