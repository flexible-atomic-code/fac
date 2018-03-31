c=======================================================================
c
c This is an example #2 of how to use the model QED potential in actual
c calculations. 
c
c The code reads the model QED potential and the one-electron wave 
c functions from the disk and calculates the expectation values of the
c model QED operator with these wave functions.
c
c=======================================================================
        program main_potuse
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
c       - - - - - - - - - - - - - - - - - - - - - - - - -

        character*1 ch
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
        common /ii/ii ! number of grid points
        common /cl/cl  ! speed of light

        logical OK
        character*32 dir_name
        parameter (numwfMax = 10)
        character*7 states(numwfMax)
        real*8 due(numwfMax),dwk(numwfMax),dse(numwfMax)
        common /conf/nc,KC(numwfMax),NnC(numwfMax)
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
        read *,data_filename

        print *,
     & "input the directory name where the wave functions are stored"
        read *,dir_name

        OK = .TRUE.
        nc = 0
        do while (OK)
           nc = nc + 1
           write (*,"(a,i1)") 
     & "input kappa and n of the reference state #",nc
           read *,KC(nc),NnC(nc)
           write (*,*) "next state? (y/n)"
           read *,ch
           if (ch.ne."y".and.ch.ne."Y") OK = .FALSE.
        enddo
c
c-
c     Reading model QED potential from disk
c--
        call read_pot(data_filename)
c
c-
c     Reading the wave functions from disk
c---
        do i = 1,nc
           call wfunc_read(i,dir_name)
        enddo
c
c-
c     Output
c--
        write( *,1)
 1      format(/2x,"Potuse: An exemplary calculation",/2x,
     & "Matrix elements of the model QED ",
     & "potential with wave functions prestored on the disk")
        write( *,2) z
 2      format (/2x,'Z=',f10.4)
        if (knucl.eq.0) then
           write (*,*) " Point nuclear charge"
        else if (knucl.eq.1) then
           write (*,*) " Extended nuclear charge"
        end if
        write (*,*) " Units are F(Z\alpha)"
        write( *,55)
55      format(/6x,'State',13x,'Ueling',11x,'WK',13x,'SE',11x,'QEDMOD')
c
c-
c     Calculation of matrix elements of the model QED operator with
c      Dirac wave functions
c--
        alz=z/cl

        do jj = 1,nc
           kappa = KC(jj)
           n = NnC(jj)
           states(jj) = GetState(n,kappa)

             ! wave function is interpolated and stored on a grid
           do i=1,ii
              call wfunc_interpolate(kappa,n,r(i),p(i),q(i))
           enddo

             ! the result of the SE operator acting on the Dirac wave function
             ! is calculated and stored on a grid
           call se_pot_wav(n,kappa,r,p,q,cp,cq)

             ! matrix element of the SE operator                    
           dse(jj)=tint(0,p,q,cp,cq,r,v)

             ! matrix element of the Ueling potential
           due(jj)=sint(uehl,p,q,p,q,r,v)

             ! matrix element of the Wichmann-Kroll potential
           dwk(jj)=sint(v_wk,p,q,p,q,r,v)*cl**2
                 
           coef = cl/pi*alz**4/n**3
                 
           write( *,65) states(jj)
     &            ,due(jj)/coef,dwk(jj)/coef,dse(jj)/coef
     &            ,(due(jj)+dwk(jj)+dse(jj))/coef
 65        format(5x,a7,4x,4f15.6)
        enddo
*
        do jj = 2,nc
           write( *,66) states(1)//"--"//states(jj)
     &            ,(due(1)-due(jj))/coef
     &            ,(dwk(1)-dwk(jj))/coef
     &            ,(dse(1)-dse(jj))/coef
     &            ,(due(1)+dwk(1)+dse(1)-(due(jj)+dwk(jj)+dse(jj)))/coef
 66        format(a,4f15.6)
        enddo
c--
        stop
        end
c=======================================================================
c
c=======================================================================
        subroutine wfunc_interpolate(KA,NnA,x,GA,FA)
c
c 5-point interpolation of the stored wave function
c
        implicit real*8 (a-h,o-z)
        parameter (numRadPo = 10000)
        real*8 Fg(0:10),Ff(0:10)

        parameter (numwfMax = 10)
        common /conf/nc,KC(numwfMax),NnC(numwfMax)
        common /DF_/xr(numRadPo,numwfMax)
     &             ,Gc(numRadPo,numwfMax)
     &             ,Fc(numRadPo,numwfMax)
        common /DF /npoints(numwfMax)

        save left
        data left/0/
c--
        inum = 0
        do j = 1,nc
           if ((KA.eq.KC(j)).and.(NnA.eq.NnC(j))) inum = j
        enddo
        if (inum.eq.0) then
           print *,"kappa = ",KA," n = ",NnA
           stop 'wf_DF_int: 1'
        end if
c-
        GA = 0.d0
        FA = 0.d0
        nmax = npoints(inum)

        if (x.gt.xr(nmax,inum)) return
c-
        if ((left.ge.nmax).or.(left.lt.1)) left=1
        if (x.lt.xr(left,inum)) left=1

        do while ((left.lt.nmax).and.(x.gt.xr(left+1,inum)))
           left = left+1
        enddo
        
        imin = max(1,   left-2)
        imax = min(nmax,left+2)
        if (imin.eq.1) imax = 3
        if (imax.eq.nmax) imin = nmax-4
        
        do j = imin,imax
           Fg(j-imin) = Gc(j,inum)
           Ff(j-imin) = Fc(j,inum)
        enddo
        
        do j = imin+1,imax
           do i = j,imax
              Fg(i-imin) = (Fg(j-imin-1)-Fg(i-imin))
     >                                /(xr(j-1,inum)-xr(i,inum))
              Ff(i-imin) = (Ff(j-imin-1)-Ff(i-imin))
     >                                /(xr(j-1,inum)-xr(i,inum))
           enddo
        enddo

        GA = Fg(imax-imin)
        FA = Ff(imax-imin)
        do i = imax-1,imin,-1
           GA = Fg(i-imin)+ (x-xr(i,inum))* GA
           FA = Ff(i-imin)+ (x-xr(i,inum))* FA
        enddo
c--   
        return
        end
c=======================================================================
c
c=======================================================================
        subroutine wfunc_read(nwf,dir_name)
c
c Reading the wave function from the disk and storing it into common blocks
c
        implicit real*8 (a-h,o-z)
        parameter (nwfMax = 10)
        parameter (numRadPo = 10000)

        character*32 dir_name
        character*32 file_name
        common /conf/nc,KC(nwfMax),NnC(nwfMax)

        logical found
        
        common /DF_/xr(numRadPo,nwfMax)
     &             ,Gc(numRadPo,nwfMax)
     &             ,Fc(numRadPo,nwfMax)
        common /DF /npoints(nwfMax)
c--   
        kappa = KC(nwf)
        n = NnC(nwf)
        j=2*iabs(kappa)-1
        l=(iabs(2*kappa+1)-1)/2

         ! generate file name where wave function is stored
        file_name = ""
        write (file_name(1:1),"(i1)") n
        if (l.eq.0) then
           file_name(2:2) = "s"
        else if (l.eq.1) then
           file_name(2:2) = "p"
        else if (l.eq.2) then
           file_name(2:2) = "d"
        else
           file_name(2:2) = "X"
        end if
        if (j > 2*l) then
           file_name(3:8) = "_p.dat"
        else
           file_name(3:8) = "_m.dat"
        end if
        file_name = dir_name(1:length(dir_name))//"/"//file_name
        
        inquire (file=file_name,exist=found)
        if (.not.found)  then
           write (*,*) "wfunc_read: file ",file_name," is not found"
           stop 
        end if
        open (23,file=file_name,status='old')
*     
        write (*,"(2x,a,i2,a,i1)")
     & "Reading wave function from "//file_name(1:length(file_name))//
     & " with kappa = ",kappa," n = ",n
        read (23,*) KC_,NnC_,ea_
        i = 0
        do while (.TRUE.)
           read (23,FMT=*,END=111,ERR=111) r_,Gn_,Fn_
           i = i+1
           if (i.gt.numRadPo) stop 'wf_DF_stor: numRadPo is too small'
           xr(i,nwf) = r_
           Gc(i,nwf) = Gn_ 
           Fc(i,nwf) = Fn_ 
        enddo
 111    npoints(nwf) = i
        close (23)
c--   
        return
        end
c=======================================================================
c
c=======================================================================
      FUNCTION LENGTH (STRING)
      CHARACTER*(*) STRING
*
      LRIGHT = LEN (STRING)
      DO I = LRIGHT,1,-1
         IF (STRING(I:I) .NE. ' ') THEN
            LENGTH = I
            return
         ENDIF
      enddo
      LENGTH = 0
*
      return
      END
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
