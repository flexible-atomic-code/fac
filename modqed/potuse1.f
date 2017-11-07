c=======================================================================
c This is an example #1 of how to use the model QED potential in actual
c calculations. 
c
c The code reads the model QED potential from the disk and calculates 
c the matrix elements of the model QED operator with the hydrogenic (Dirac) 
c wave functions for the point or the extended nucleus. 
c
c=======================================================================
        program main_potuse
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
c       - - - - - - - - - - - - - - - - - - - - - - - - -

        character*7 GetState
        character*72 data_filename
        real*8 p(maxii),q(maxii) 
        real*8 cp(maxii),cq(maxii) 

        common /z/z     ! nuclear charge
        common /knucl/knucl  ! nuclear model (0,1)
        common /r/r(maxii)  ! radial grid
        common /v/v(maxii) ! weights of Simpson quadrature formula on radial grid
        common /uehl/uehl(maxii) ! Uehling potential on radial grid
        common /v_wk/v_wk(maxii) ! WK potential on radial grid
        common /cl/cl  ! speed of light
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
        print *,
     & "input the name of the data file with model QED potential"
        read (*,'(a)') data_filename
c-
c     Reading model QED potential from disk
c--
        call read_pot(data_filename)
c
c-
c     Output
c--
        write( *,1)
 1      format(/2x,"Potuse: An exemplary calculation",/2x,
     & "Matrix elements of the model QED ",
     & "potential with hydrogenic wave functions")
        write( *,2) z
 2      format (/2x,'Z=',f10.4)
        if (knucl.eq.0) then
           write (*,*) " Point nuclear charge"
        else if (knucl.eq.1) then
           write (*,*) " Extended nuclear charge"
        end if
        write (*,*) " Units are F(Z\alpha)"
        write( *,55)
55      format(/6x,'State',9x,'Ueling',11x,'WK',13x,'SE',11x,'QEDMOD')

c
c-
c     Calculation of matrix elements of the model QED operator with
c      Dirac wave functions
c--
        alz=z/cl

        do kapabs = 1,3
           do isig = -1,1,2
              kappa = isig*kapabs
              l=(iabs(2*kappa+1)-1)/2

              do n = 1,3
                 if (n.lt.l+1) cycle

                   ! the Dirac wave function for point and extended nucleus
                   ! is calculated and stored on a grid
                 call dirac(n,kappa,p,q)

                   ! the result of the SE operator acting on the Dirac wave function
                   ! is calculated and stored on a grid
                 call se_pot_wav(n,kappa,r,p,q,cp,cq)

                   ! matrix element of the SE operator                    
                 dse=tint(0,p,q,cp,cq,r,v)

                   ! matrix element of the Ueling potential
                 due=sint(uehl,p,q,p,q,r,v)

                   ! matrix element of the Wichmann-Kroll potential
                 dwk=sint(v_wk,p,q,p,q,r,v)*cl**2
                 
                 coef = cl/pi*alz**4/n**3
                 
                 write( *,65) GetState(n,kappa)
     &                ,due/coef,dwk/coef,dse/coef
     &                ,(due+dwk+dse)/coef
 65              format(5x,a7,4f15.6)
              enddo
           enddo
        enddo
c--
        stop
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
        subroutine exit1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call exit(1)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
