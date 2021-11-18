c=======================================================================
c This is an example #1 of how to use the model QED potential in actual
c calculations. 
c
c The code reads the model QED potential from the disk and calculates 
c the matrix elements of the model QED operator with the hydrogenic (Dirac) 
c wave functions for the point or the extended nucleus. 
c
c=======================================================================
        subroutine genqed(n0,n1,ka,iw,nr,r0,p0,q0,p1,q1,dse)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
c       - - - - - - - - - - - - - - - - - - - - - - - - -

        character*7 GetState
c        character*72 data_filename
        real*8 p(maxii),q(maxii) 
        real*8 pd(maxii),qd(maxii) 
        real*8 cp(maxii),cq(maxii) 

        common /z/z     ! nuclear charge
        common /knucl/knucl  ! nuclear model (0,1)
        common /r/r(maxii)/rho/rho(maxii)  ! radial grid
        common /ii/ii/h/h/r1/r1/r2/r2/al/al/bt/bt
        common /v/v(maxii) ! weights of Simpson quadrature formula on radial grid
        common /uehl/uehl(maxii) ! Uehling potential on radial grid
        common /v_wk/v_wk(maxii) ! WK potential on radial grid
        common /cl/cl            ! speed of light
        common /ns_proj/ns_proj
        common /nn_proj/nn_proj(maxns)/ll_proj/ll_proj(maxns)
     &  /jj_proj/jj_proj(maxns)/kk_proj/kk_proj(maxns)
c
c-
c     Some constants
c--
        cl=137.03599911d0
        pi=3.1415926535897932385d0
        alpha=1.d0/cl
c
c-
c     Input
c--
c        print *,
c     & "input the name of the data file with model QED potential"
c        read (*,'(a)') data_filename
c-
c     Reading model QED potential from disk
c--
c        call read_pot(data_filename)
c
c-
c     Output
c--
c        write( *,1)
c 1      format(/2x,"Potuse: An exemplary calculation",/2x,
c     & "Matrix elements of the model QED ",
c     & "potential with hydrogenic wave functions")
c        write( *,2) z
c 2      format (/2x,'Z=',f10.4)
c        if (knucl.eq.0) then
c           write (*,*) " Point nuclear charge"
c        else if (knucl.eq.1) then
c           write (*,*) " Extended nuclear charge"
c        end if
c        write (*,*) " Units are F(Z\alpha)"
c        write( *,55)
c55      format(/6x,'State',9x,'Ueling',11x,'WK',13x,'SE',11x,'QEDMOD')

c
c-
c     Calculation of matrix elements of the model QED operator with
c      Dirac wave functions
c--
        alz=z/cl

        l=(iabs(2*ka+1)-1)/2
        if (l .gt. 2) then
           dse = 0d0
           return
        endif
        call uvip3p(3, nr, r0, p0, ii, rho, p)
        call uvip3p(3, nr, r0, q0, ii, rho, q)
        do i=1,ii
           q(i)=-q(i)
        enddo
! the result of the SE operator acting on the Dirac wave function
! is calculated and stored on a grid
        call se_pot_wav(n0,ka,iw,r,p,q,cp,cq)
        
!     matrix element of the SE operator
        if (n0 .ne. n1) then
           call uvip3p(3, nr, r0, p1, ii, rho, p)
           call uvip3p(3, nr, r0, q1, ii, rho, q)
           do i=1,ii
              q(i) = -q(i)
           enddo
c           do ni=1,ns_proj
c              ki=kk_proj(ni)
c              if (ki .eq. ka) then
c                 call read_func(ni,pd,qd,2)
c                 wi=tint(0,p,q,pd,qd,r,v)
c                 dse=dse+wi*tint(0,pd,qd,cp,cq,r,v)
c              endif
c           enddo
           dse=tint(0,p,q,cp,cq,r,v)
        else     
           dse=tint(0,p,q,cp,cq,r,v)
           coef = cl/pi*alz**4/n0**3
           dse = dse/coef
        endif
! matrix element of the Ueling potential
c        due=sint(uehl,p,q,p,q,r,v)
        
! matrix element of the Wichmann-Kroll potential
c        dwk=sint(v_wk,p,q,p,q,r,v)*cl**2
c        if (n0 .eq. n1) then
c           coef = cl/pi*alz**4/n0**3
c           dse = dse/coef
c        endif
c        write( *,65) GetState(n,kappa)
c     &       ,due/coef,dwk/coef,dse/coef
c     &       ,(due+dwk+dse)/coef
c 65     format(5x,a7,4f15.6)
        end
c=======================================================================
c
c=======================================================================
        character*7 function GetState(n,kappa)
c
c Generates a spectroscopical notation of the state
c
        implicit real*8 (a-h,o-z)
        character*1 ch
*
        j=2*iabs(kappa)-1
        l=(iabs(2*kappa+1)-1)/2

        if (l.eq.0) then
           ch = "s"
        else if (l.eq.1) then
           ch = "p"
        else if (l.eq.2) then
           ch = "d"
        else if (l.eq.3) then
           ch = "f"
        else
           ch = "X"
        end if

        write (GetState,"(i2,a1,i2,a2)") n,ch,j,"/2"
*
        return
        end
c       =================================================
