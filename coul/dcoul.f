C     calculate dirac coulomb function for positive and negtive energies
C     INPUT:
C     z Nuclear charge. positive for electron-ion system.
C     e energy.
C     k kappa
C     r radius
C     OUTPUT:
C     p large component of regular solution
C     q small component of regular solution
C     p1 iregular large for e > 0. ignore for e < 0
C     q1 iregular small for e > 0, ignore for e < 0
C     ierr error code returned by coulcc
      subroutine dcoul(z, e, k, r, p, q, p1, q1, ierr)
      implicit none     
      integer k, ierr, kfn, inorm
      double precision z, e, r, p, q, p1, q1, c, ki, zp, gam
      double precision lambda, y, qi, x0, b1, b2, np
      complex*16 x, eta, zlmin, omega, a, pp, qq, mu, nu, IONE
      complex*16 fc(1), gc(1), fcp(1), gcp(1), sig(1), clgam, lam0
      double precision SL, SL2, TSL2, ALPHA
      PARAMETER (SL=137.036D0,SL2=SL*SL,TSL2=SL2+SL2,ALPHA=1.0D0/SL)
      real*8 HALFPI
      parameter (HALFPI = 1.5707963268D0)
      
      inorm = ierr
      IONE = dcmplx(0.0, 1.0)
      c = 1.0+0.5*e/SL2
      ki = sqrt(2.0*abs(e)*c)
      zp = z*ALPHA
      gam = sqrt(k*k - zp*zp)
      lambda = gam - 0.5
      qi = sqrt(c/ki)
      y = (1.0+e/SL2)*z/ki
      
      x0 = ki*r      
      if (e .lt. 0) then
         x = dcmplx(0.0, x0)
         if (inorm .eq. 1) then
            eta = dcmplx(1D-5*(0.5+y), 0.5+y)
         else
            eta = dcmplx(0.0, 0.5+y)
         endif
         mu = dcmplx(k - z/ki, 0.0)
         nu = dcmplx(0.5+y-x0, 0.0)
      else
         x = dcmplx(x0, 0.0)
         eta = dcmplx(-y, 0.5)
         mu = dcmplx(k, -z/ki)
         nu = IONE*(x - eta)
      endif

      zlmin = dcmplx(lambda, 0.0)
      ierr = 1
      kfn = 0
      
      call coulcc(x, eta, zlmin, 1, fc, gc, fcp, gcp, sig, 
     +     11, kfn, ierr)
      if (e .lt. 0) then
         omega = IONE*(HALFPI*(lambda - y - 0.5) - sig(1))
         if (inorm .gt. 0 .and. z > 0) then
            np = y - gam
            b1 = sqrt(z*c*(z/ki-k))*ki/z
            lam0 = np+1.0
            b2 = dble(clgam(lam0))
            lam0 = np + 1.0 + 2.0*gam
            b2 = b2 + dble(clgam(lam0))
            b2 = -0.5*b2
            omega = omega + b2
            a = exp(omega)*b1
            a = a/(mu*sqrt(2.0*x0))
         else
            a = exp(IONE*dimag(omega))
            a = a/mu
            a = a*qi/sqrt(2.0*x0)
         endif
         pp = a*((mu + nu)*gc(1) - x*gcp(1))
         qq = (ALPHA*e/ki)*a*((mu - nu)*gc(1) + x*gcp(1))
         p = dble(pp)
         q = dble(qq)
         p1 = dble(omega)
         q1 = dimag(omega)
      else
         a = qi/sqrt(2.0*IONE*mu*x)
         pp = a*((mu + nu)*gc(1) - x*gcp(1))
         qq = (IONE*e*ALPHA/ki)*a*((mu - nu)*gc(1) + x*gcp(1))
         p1 = dble(pp)
         q1 = dble(qq)
         p = dimag(pp)
         q = dimag(qq)
      endif

      end
      
C  **************************************************************       
C                       SUBROUTINE DCOUL1                                
C  **************************************************************       
      SUBROUTINE DCOUL1(Z1,E,K,R,FU,FL,GU,GL,IERR)                         
C                                                                       
C     THIS SUBROUTINE COMPUTES RADIAL DIRAC-COULOMB WAVE FUNC-          
C  TIONS FOR FREE STATES.                                               
C                                                                       
C  **** ALL QUANTITIES IN ATOMIC UNITS.                                 
C                                                                       
C  INPUT ARGUMENTS:                                                     
C     Z ........ FIELD STRENGTH, I.E. VALUE OF R*V(R) (ASSUMED          
C                CONSTANT).                                             
C     E ........ PARTICLE KINETIC ENERGY (POSITIVE).                    
C     K ........ ANGULAR MOMENTUM QUANTUM NUMBER KAPPA (.NE.0).         
C     R ........ RADIAL DISTANCE (POSITIVE).                            
C                                                                       
C  OUTPUT ARGUMENTS:                                                    
C     FU, FL ... UPPER AND LOWER COMPONENTS OF THE REGULAR RADIAL       
C                COULOMB FUNCTION.                                      
C     GU, GL ... UPPER AND LOWER COMPONENTS OF THE IRREGULAR RA-        
C                DIAL COULOMB FUNCTION.                                 
C     ERR ...... ACCURACY OF THE COMPUTED FUNCTIONS (RELATIVE UN-       
C                CERTAINTY).                                            
C  OUTPUT THROUGH COMMON/OCOUL/:                                        
C     WAVNUM ... WAVE NUMBER.                                           
C     ETA ...... SOMMERFELD'S PARAMETER.                                
C     DELTA .... COULOMB PHASE SHIFT (MODULUS 2*PI).                    
C                                                                       
C     RADIAL FUNCTIONS ARE NORMALIZED SO THAT, FOR LARGE R, THE         
C  UPPER COMPONENT OSCILLATES WITH UNIT AMPLITUDE.                      
C                                                                       
C     OTHER SUBPROGRAMS REQUIRED: SUBROUTINES FCOUL AND SUM2F0,         
C                                 AND FUNCTION CLGAM.                   
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)           
      PARAMETER (PI=3.1415926535897933D0,PIH=0.5D0*PI,TPI=PI+PI,        
     1  SL=137.036D0,SL2=SL*SL,TSL2=SL2+SL2,ALPHA=1.0D0/SL)             
      COMMON/OFCOUL/DELTA0                                              
      COMMON/OCOUL/WAVNUM,ETA,DELTA      
!$OMP THREADPRIVATE(/OFCOUL/,/OCOUL/)
C               
      Z = -Z1
      IERR = 0
      IF(E.LT.0.0001D0.OR.K.EQ.0) THEN                                  
        IERR = 2
        FU=0.0D0                                                        
        FL=0.0D0                                                        
        GU=0.0D0                                                        
        GL=0.0D0                                                        
        ERR=1.0D0                                                       
        IF(E.LT.0.0001D0) WRITE(6,2101)                                 
 2101   FORMAT(1X,'*** ERROR IN DCOUL: E IS TOO SMALL.')                
        IF(K.EQ.0) WRITE(6,2102)                                        
 2102   FORMAT(1X,'*** ERROR IN DCOUL: K.EQ.0.')                        
        RETURN                                                          
      ENDIF                                                             
C                                                                       
C  ****  PARAMETERS.                                                    
C                                                                       
      IF(DABS(Z).GT.0.00001D0) THEN                                     
        ZETA=Z*ALPHA                                                    
        ICAL=0                                                          
      ELSE                                                              
        ZETA=0.0D0                                                      
        ICAL=1                                                          
      ENDIF                                                             
      RLAMBS=K*K-ZETA*ZETA                                              
      RLAMB=DSQRT(RLAMBS)                                               
      W=E+SL2                                                           
      PC=DSQRT(E*(E+TSL2))                                              
      WAVNUM=PC/SL                                                      
      ETA=ZETA*W/PC                                                     
      X=WAVNUM*R                                                        
      RLA=DSQRT(RLAMBS+ETA*ETA)                                         
      IF(ICAL.EQ.1) GO TO 1                                             
C                                                                       
C  ************  COULOMB FUNCTIONS.                                     
C                                                                       
      DELTA0=1.0D30                                                     
      RLAMB1=RLAMB-1.0D0                                                
      CALL FCOUL(ETA,RLAMB1,X,FM1,FPM1,GM1,GPM1,ERR)                    
      IF(ERR.GE.1.0D-6) RETURN                                          
      SLA=(RLAMB/X)+(ETA/RLAMB)                                         
      F=RLAMB*(SLA*FM1-FPM1)/RLA                                        
      G=RLAMB*(SLA*GM1-GPM1)/RLA                                        
C                                                                       
C  ****  DIRAC-COULOMB WAVE FUNCTIONS AND PHASE SHIFT.                  
C                                                                       
      P1=K+RLAMB                                                        
      P2=RLAMB*SL2-K*W                                                  
      RNUR=ZETA*(W+SL2)                                                 
      RNUI=-P1*PC                                                       
      RNORM=1.0D0/(DSQRT(RNUR*RNUR+RNUI*RNUI)*RLAMB)                    
      IF(K.GT.0) THEN                                                   
        L=K                                                             
      ELSE                                                              
        L=-K-1                                                          
      ENDIF                                                             
C                                                                       
      IF(DELTA0.GT.1.0D29) DELTA0=DELTAC(ETA,RLAMB)                     
      RNU=DATAN2(RNUI,RNUR)                                             
      DELTA=RNU-(RLAMB-L-1)*PIH+DELTA0                                  
      IF(Z.LT.0.0D0.AND.K.LT.0) THEN                                    
        RNORM=-RNORM                                                    
        DELTA=DELTA-PI                                                  
      ENDIF                                                             
      IF(DELTA.GE.0.0D0) THEN                                           
        DELTA=DMOD(DELTA,TPI)                                           
      ELSE                                                              
        DELTA=-DMOD(-DELTA,TPI)                                         
      ENDIF                                                             
      Q2=P1*P2*RNORM                                                    
      Q1=RLA*PC*RNORM                                                   
      P1=P1*Q1                                                          
      Q1=ZETA*Q1                                                        
      P2=ZETA*P2*RNORM                                                  
C                                                                       
      FU=P1*F+P2*FM1                                                    
      GU=P1*G+P2*GM1                                                    
      FL=Q1*F+Q2*FM1                                                    
      GL=Q1*G+Q2*GM1    
      FU = -FU
      GU = -GU
      RETURN                                                            
C                                                                       
C  ************  Z=0. SPHERICAL BESSEL FUNCTIONS.                       
C                                                                       
    1 CONTINUE                                                          
      RLAMB=IABS(K)                                                     
      CALL FCOUL(0.0D0,RLAMB,X,F,FP,G,GP,ERR)                           
      IF(ERR.GE.1.0D-6) THEN
         IERR = 1
         RETURN  
      ENDIF                                        
      FM1=(RLAMB*F/X)+FP                                                
      GM1=(RLAMB*G/X)+GP                                                
      FACT=DSQRT(E/(E+TSL2))                                            
      IF(K.LT.0) THEN                                                   
        FU=FM1                                                          
        GU=GM1                                                          
        FL=FACT*F                                                       
        GL=FACT*G                                                       
      ELSE                                                              
        FU=F                                                            
        GU=G                                                            
        FL=-FACT*FM1                                                    
        GL=-FACT*GM1                                                    
      ENDIF                                                             
      DELTA=0.0D0                                                       
      FU = -FU
      GU = -GU
      RETURN                                                            
      END                                                               
C  **************************************************************       
C                       SUBROUTINE FCOUL                                
C  **************************************************************       
      SUBROUTINE FCOUL(ETA,RLAMB,X,F,FP,G,GP,ERR)                       
C                                                                       
C     CALCULATION OF COULOMB FUNCTIONS FOR REAL ETA, RLAMB.GT.-1        
C  AND X LARGER THAN, OR OF THE ORDER OF XTP0 (THE TURNING POINT        
C  FOR RLAMB=0). STEED'S CONTINUED FRACTION METHOD IS COMBINED          
C  WITH RECURSION RELATIONS AND AN ASYMPTOTIC EXPANSION. THE            
C  OUTPUT VALUE ERR=1.0D0 INDICATES THAT THE ADOPTED EVALUATION         
C  ALGORITHM IS NOT APPLICABLE (X IS TOO SMALL).                        
C                                                                       
C  INPUT ARGUMENTS:                                                     
C     ETA ...... SOMMERFELD'S PARAMETER.                                
C     RLAMB .... ANGULAR MOMENTUM.                                      
C     X ........ VARIABLE (=WAVE NUMBER TIMES RADIAL DISTANCE).         
C                                                                       
C  OUTPUT ARGUMENTS:                                                    
C     F, FP .... REGULAR FUNCTION AND ITS DERIVATIVE.                   
C     G, GP .... IRREGULAR FUNCTION AND ITS DERIVATIVE.                 
C     ERR ...... RELATIVE NUMERICAL UNCERTAINTY. A VALUE OF THE         
C                ORDER OF 10**(-N) MEANS THAT THE CALCULATED            
C                FUNCTIONS ARE ACCURATE TO N DECIMAL FIGURES.           
C                THE MAXIMUM ACCURACY ATTAINABLE WITH DOUBLE            
C                PRECISION ARITHMETIC IS ABOUT 1.0D-15.                 
C                                                                       
C     OTHER SUBPROGRAMS REQUIRED: SUBROUTINE SUM2F0 AND                 
C                                 FUNCTIONS DELTAC AND CLGAM.           
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)           
      PARAMETER (PI=3.1415926535897932D0,PIH=0.5D0*PI,TPI=PI+PI,        
     1  EPS=1.0D-16,TOP=1.0D5,NTERM=1000)                               
      COMMON/OFCOUL/DELTA       
!$OMP THREADPRIVATE(/OFCOUL/)                                              
C                                                                       
      IF(RLAMB.LT.-0.999D0) THEN                                        
        WRITE(6,'(1X,''*** ERROR IN RCOUL: RLAMB.LT.-0.999'')')         
        STOP                                                            
      ENDIF                                                             
      IF(X.LT.EPS) GO TO 10                                             
C                                                                       
C  ****  NUMERICAL CONSTANTS.                                           
C                                                                       
      CI=DCMPLX(0.0D0,1.0D0)                                            
      CI2=2.0D0*CI                                                      
      CIETA=CI*ETA                                                      
      X2=X*X                                                            
      ETA2=ETA*ETA                                                      
C                                                                       
C  ****  TURNING POINT (XTP). (44)                                      
C                                                                       
      IF(RLAMB.GE.0.0D0) THEN                                           
        XTP=ETA+DSQRT(ETA2+RLAMB*(RLAMB+1.0D0))                         
      ELSE                                                              
        XTP=EPS                                                         
      ENDIF                                                             
      ERRS=10.0D0                                                       
      IF(X.LT.XTP) GO TO 1                                              
C                                                                       
C  ************  ASYMPTOTIC EXPANSION. (71-75)                          
C                                                                       
C  ****  COULOMB PHASE-SHIFT.                                           
      DELTA=DELTAC(ETA,RLAMB)                                           
C                                                                       
      CPA=CIETA-RLAMB                                                   
      CPB=CIETA+RLAMB+1.0D0                                             
      CPZ=CI2*X                                                         
      CALL SUM2F0(CPA,CPB,CPZ,C2F0,ERR1)                                
      CQA=CPA+1.0D0                                                     
      CQB=CPB+1.0D0                                                     
      CALL SUM2F0(CQA,CQB,CPZ,C2F0P,ERR2)                               
      C2F0P=CI*C2F0P*CPA*CPB/(2.0D0*X2)                                 
C  ****  FUNCTIONS.                                                     
      THETA=X-ETA*DLOG(2.0D0*X)-RLAMB*PIH+DELTA                         
      IF(THETA.GT.1.0D4) THETA=DMOD(THETA,TPI)                          
      CEITH=CDEXP(CI*THETA)                                             
      CGIF=C2F0*CEITH                                                   
      G=CGIF                                                            
      F=-CI*CGIF                                                        
C  ****  DERIVATIVES.                                                   
      CGIFP=(C2F0P+CI*(1.0D0-ETA/X)*C2F0)*CEITH                         
      GP=CGIFP                                                          
      FP=-CI*CGIFP                                                      
C  ****  GLOBAL UNCERTAINTY. THE WRONSKIAN MAY DIFFER FROM 1 DUE        
C        TO TRUNCATION AND ROUNDOFF ERRORS.                             
      ERR=DMAX1(ERR1,ERR2,DABS(G*FP-F*GP-1.0D0))                        
      IF(ERR.LE.EPS) RETURN                                             
      ERRS=ERR                                                          
C                                                                       
C  ************  STEED'S CONTINUED FRACTION METHOD.                     
C                                                                       
    1 CONTINUE                                                          
      CIETA2=CIETA+CIETA                                                
      ETAX=ETA*X                                                        
C                                                                       
C  ****  CONTINUED FRACTION FOR F. (60-70)                              
C                                                                       
      INULL=0                                                           
      RLAMBN=RLAMB+1.0D0                                                
      A1=-(RLAMBN+1.0D0)*(RLAMBN**2+ETA2)*X/RLAMBN                      
      B0=(RLAMBN/X)+(ETA/RLAMBN)                                        
      B1=(2.0D0*RLAMBN+1.0D0)*(RLAMBN*(RLAMBN+1.0D0)+ETAX)              
      FA3=B0                                                            
      FA2=B0*B1+A1                                                      
      FB3=1.0D0                                                         
      FB2=B1                                                            
      RF=FA3                                                            
C                                                                       
      DO 2 N=2,NTERM                                                    
      RFO=RF                                                            
      DAF=DABS(RF)                                                      
      RLAMBN=RLAMB+N                                                    
      AN=-(RLAMBN**2-1.0D0)*(RLAMBN**2+ETA2)*X2                         
      BN=(2.0D0*RLAMBN+1.0D0)*(RLAMBN*(RLAMBN+1.0D0)+ETAX)              
      FA1=FA2*BN+FA3*AN                                                 
      FB1=FB2*BN+FB3*AN                                                 
      TST=DABS(FB1)                                                     
C                                                                       
      IF(TST.LT.1.0D-25) THEN                                           
        IF(INULL.GT.0) STOP                                             
        INULL=1                                                         
        FA3=FA2                                                         
        FA2=FA1                                                         
        FB3=FB2                                                         
        FB2=FB1                                                         
        RF=RFO                                                          
      ELSE                                                              
        FA3=FA2/TST                                                     
        FA2=FA1/TST                                                     
        FB3=FB2/TST                                                     
        FB2=FB1/TST                                                     
        RF=FA2/FB2                                                      
        IF(DABS(RF-RFO).LT.EPS*DAF) GO TO 3                             
      ENDIF                                                             
    2 CONTINUE                                                          
    3 CONTINUE                                                          
      IF(DAF.GT.1.0D-25) THEN                                           
        ERRF=DABS(RF-RFO)/DAF                                           
      ELSE                                                              
        ERRF=EPS                                                        
      ENDIF                                                             
      IF(ERRF.GT.ERRS) THEN                                             
        ERR=ERRS                                                        
        RETURN                                                          
      ENDIF                                                             
C                                                                       
C  ****  DOWNWARD RECURSION FOR F AND FP. ONLY IF RLAMB.GT.1 AND        
C        X.LT.XTP. (48,49)                                              
C                                                                       
      RLAMB0=RLAMB                                                      
      IF(X.GE.XTP.OR.RLAMB0.LT.1.0D0) THEN                              
        ISHIFT=0                                                        
      ELSE                                                              
        FT=1.0D0                                                        
        FTP=RF                                                          
        IS0=RLAMB0+1.0D-6                                               
        TST=X*(X-2.0D0*ETA)                                             
        DO 4 I=1,IS0                                                    
        ETARL0=ETA/RLAMB0                                               
        RL=DSQRT(1.0D0+ETARL0**2)                                       
        SL=(RLAMB0/X)+ETARL0                                            
        RLAMB0=RLAMB0-1.0D0                                             
        FTO=FT                                                          
        FT=(SL*FT+FTP)/RL                                               
        FTP=SL*FT-RL*FTO                                                
        IF(FT.GT.1.0D10) THEN                                           
          FTP=FTP/FT                                                    
          FT=1.0D0                                                      
        ENDIF                                                           
        RL1T=RLAMB0*(RLAMB0+1.0D0)                                      
        IF(TST.GT.RL1T) THEN                                            
          ISHIFT=I                                                      
          GO TO 5                                                       
        ENDIF                                                           
    4   CONTINUE                                                        
        ISHIFT=IS0                                                      
    5   CONTINUE                                                        
        XTPC=ETA+DSQRT(ETA2+RL1T)                                       
        RFM=FTP/FT                                                      
      ENDIF                                                             
C                                                                       
C  ****  CONTINUED FRACTION FOR P+CI*Q WITH RLAMB0. (76-79)             
C                                                                       
      INULL=0                                                           
      CAN=CIETA-ETA2-RLAMB0*(RLAMB0+1.0D0)                              
      CB0=X-ETA                                                         
      CBN=2.0D0*(X-ETA+CI)                                              
      CFA3=CB0                                                          
      CFA2=CB0*CBN+CAN                                                  
      CFB3=1.0D0                                                        
      CFB2=CBN                                                          
      CPIQ=CFA3                                                         
C                                                                       
      DO 6 N=2,NTERM                                                    
      CPIQO=CPIQ                                                        
      DAPIQ=CDABS(CPIQ)                                                 
      CAN=CAN+CIETA2+(N+N-2)                                            
      CBN=CBN+CI2                                                       
      CFA1=CFA2*CBN+CFA3*CAN                                            
      CFB1=CFB2*CBN+CFB3*CAN                                            
      TST=CDABS(CFB1)                                                   
C                                                                       
      IF(TST.LT.1.0D-25) THEN                                           
        IF(INULL.GT.0) STOP                                             
        INULL=1                                                         
        CFA3=CFA2                                                       
        CFA2=CFA1                                                       
        CFB3=CFB2                                                       
        CFB2=CFB1                                                       
        CPIQ=CPIQO                                                      
      ELSE                                                              
        CFA3=CFA2/TST                                                   
        CFA2=CFA1/TST                                                   
        CFB3=CFB2/TST                                                   
        CFB2=CFB1/TST                                                   
        CPIQ=CFA2/CFB2                                                  
        IF(CDABS(CPIQ-CPIQO).LT.EPS*DAPIQ) GO TO 7                      
      ENDIF                                                             
    6 CONTINUE                                                          
    7 CONTINUE                                                          
      IF(DAPIQ.GT.1.0D-25) THEN                                         
        ERRPIQ=CDABS(CPIQ-CPIQO)/DAPIQ                                  
      ELSE                                                              
        ERRPIQ=EPS                                                      
      ENDIF                                                             
      IF(ERRPIQ.GT.ERRS) THEN                                           
        ERR=ERRS                                                        
        RETURN                                                          
      ENDIF                                                             
      CPIQ=CI*CPIQ/X                                                    
C                                                                       
      RP=CPIQ                                                           
      RQ=-CI*CPIQ                                                       
      IF(RQ.LE.1.0D-25) GO TO 10                                        
      ERR=DMAX1(ERRF,ERRPIQ)                                            
C                                                                       
C  ****  INVERTING STEED'S TRANSFORMATION. (57,58)                      
C                                                                       
      IF(ISHIFT.LT.1) THEN                                              
        RFP=RF-RP                                                       
        F=DSQRT(RQ/(RFP**2+RQ**2))                                      
        IF(FB2.LT.0.0D0) F=-F                                           
        FP=RF*F                                                         
        G=RFP*F/RQ                                                      
        GP=(RP*RFP-RQ**2)*F/RQ                                          
        IF(X.LT.XTP.AND.G.GT.TOP*F) GO TO 10                            
      ELSE                                                              
        RFP=RFM-RP                                                      
        FM=DSQRT(RQ/(RFP**2+RQ**2))                                     
        G=RFP*FM/RQ                                                     
        GP=(RP*RFP-RQ**2)*FM/RQ                                         
        IF(X.LT.XTPC.AND.G.GT.TOP*FM) GO TO 10                          
C  ****  UPWARD RECURSION FOR G AND GP (IF ISHIFT.GT.0). (50,51)        
        DO 8 I=1,ISHIFT                                                 
        RLAMB0=RLAMB0+1.0D0                                             
        ETARL0=ETA/RLAMB0                                               
        RL=DSQRT(1.0D0+ETARL0**2)                                       
        SL=(RLAMB0/X)+ETARL0                                            
        GO=G                                                            
        G=(SL*GO-GP)/RL                                                 
        GP=RL*GO-SL*G                                                   
        IF(G.GT.1.0D35) GO TO 10                                        
    8   CONTINUE                                                        
    9   W=RF*G-GP                                                       
        F=1.0D0/W                                                       
        FP=RF/W                                                         
      ENDIF                                                             
C  ****  THE WRONSKIAN MAY DIFFER FROM 1 DUE TO ROUNDOFF ERRORS.        
      ERR=DMAX1(ERR,DABS(FP*G-F*GP-1.0D0))                              
      RETURN                                                            
C                                                                       
   10 F=0.0D0                                                           
      FP=0.0D0                                                          
      G=1.0D35                                                          
      GP=-1.0D35                                                        
      ERR=1.0D0                                                         
      RETURN                                                            
      END                                                               
C  **************************************************************       
C                       SUBROUTINE SUM2F0                               
C  **************************************************************       
      SUBROUTINE SUM2F0(CA,CB,CZ,CF,ERR)                                
C                                                                       
C     SUMMATION OF THE 2F0(CA,CB;CS) HYPERGEOMETRIC ASYMPTOTIC          
C  SERIES. THE POSITIVE AND NEGATIVE CONTRIBUTIONS TO THE REAL          
C  AND IMAGINARY PARTS ARE ADDED SEPARATELY TO OBTAIN AN ESTIMATE       
C  OF ROUNDING ERRORS.                                                  
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)           
      PARAMETER (EPS=1.0D-16,ACCUR=0.5D-15,NTERM=75)                    
      RRP=1.0D0                                                         
      RRN=0.0D0                                                         
      RIP=0.0D0                                                         
      RIN=0.0D0                                                         
      CDF=1.0D0                                                         
      ERR2=0.0D0                                                        
      ERR3=1.0D0                                                        
      DO 1 I=1,NTERM                                                    
      J=I-1                                                             
      CDF=CDF*(CA+J)*(CB+J)/(I*CZ)                                      
      ERR1=ERR2                                                         
      ERR2=ERR3                                                         
      ERR3=CDABS(CDF)                                                   
      IF(ERR1.GT.ERR2.AND.ERR2.LT.ERR3) GO TO 2                         
      AR=CDF                                                            
      IF(AR.GT.0.0D0) THEN                                              
        RRP=RRP+AR                                                      
      ELSE                                                              
        RRN=RRN+AR                                                      
      ENDIF                                                             
      AI=DCMPLX(0.0D0,-1.0D0)*CDF                                       
      IF(AI.GT.0.0D0) THEN                                              
        RIP=RIP+AI                                                      
      ELSE                                                              
        RIN=RIN+AI                                                      
      ENDIF                                                             
      CF=DCMPLX(RRP+RRN,RIP+RIN)                                        
      AF=CDABS(CF)                                                      
      IF(AF.GT.1.0D25) THEN                                             
        CF=0.0D0                                                        
        ERR=1.0D0                                                       
        RETURN                                                          
      ENDIF                                                             
      IF(ERR3.LT.1.0D-25*AF.OR.ERR3.LT.EPS) THEN                        
         ERR=EPS                                                        
         RETURN                                                         
      ENDIF                                                             
    1 CONTINUE                                                          
C  ****  ROUNDOFF ERROR.                                                
    2 CONTINUE                                                          
      TR=DABS(RRP+RRN)                                                  
      IF(TR.GT.1.0D-25) THEN                                            
        ERRR=(RRP-RRN)*ACCUR/TR                                         
      ELSE                                                              
        ERRR=1.0D0                                                      
      ENDIF                                                             
      TI=DABS(RIP+RIN)                                                  
      IF(TI.GT.1.0D-25) THEN                                            
        ERRI=(RIP-RIN)*ACCUR/TI                                         
      ELSE                                                              
        ERRI=1.0D0                                                      
      ENDIF                                                             
C  ****  ... AND TRUNCATION ERROR.                                      
      IF(AR.GT.1.0D-25) THEN                                            
      ERR=DMAX1(ERRR,ERRI)+ERR2/AF                                      
      ELSE                                                              
      ERR=DMAX1(ERRR,ERRI)                                              
      ENDIF                                                             
      RETURN                                                            
      END                                                               
C  **************************************************************       
C                         FUNCTION DELTAC                               
C  **************************************************************       
      FUNCTION DELTAC(ETA,RLAMB)                                        
C                                                                       
C     CALCULATION OF COULOMB PHASE SHIFT (MODULUS 2*PI). (47)           
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)           
      PARAMETER (PI=3.1415926535897932D0,TPI=PI+PI)                     
      CI=DCMPLX(0.0D0,1.0D0)                                            
C  ****  COULOMB PHASE-SHIFT.                                           
      DELTAC=-CI*CLGAM(RLAMB+1.0D0+CI*ETA)                              
      IF(DELTAC.GE.0.0D0) THEN                                          
        DELTAC=DMOD(DELTAC,TPI)                                         
      ELSE                                                              
        DELTAC=-DMOD(-DELTAC,TPI)                                       
      ENDIF                                                             
      RETURN                                                            
      END                                                               
C  **************************************************************       
C                       FUNCTION CLGAM                                  
C  **************************************************************       
      FUNCTION CLGAM(CZ)                                                
C                                                                       
C     THIS FUNCTION GIVES LOG(GAMMA(CZ)) FOR COMPLEX ARGUMENTS.         
C                                                                       
C   REF.: M. ABRAMOWITZ AND I.A. STEGUN, 'HANDBOOK OF MATHEMATI-        
C         CAL FUNCTIONS'. DOVER, NEW YORK (1974). PP 255-257.           
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)           
      PARAMETER (PI=3.1415926535897932D0)                               
      CZA=CZ                                                            
      ICONJ=0                                                           
      AR=CZA                                                            
      CLGAM=36.84136149D0                                               
      IF(CDABS(CZA).LT.1.0D-16) RETURN                                  
C                                                                       
      AI=CZA*DCMPLX(0.0D0,-1.0D0)                                       
      IF(AI.GT.0.0D0) THEN                                              
        ICONJ=0                                                         
      ELSE                                                              
        ICONJ=1                                                         
        CZA=DCONJG(CZA)                                                 
      ENDIF                                                             
C                                                                       
      CZFAC=1.0D0                                                       
      CZFL=0.0D0                                                        
    1 CZFAC=CZFAC/CZA                                                   
      IF(CDABS(CZFAC).GT.1.0D8) THEN                                    
        CZFL=CZFL+CDLOG(CZFAC)                                          
        CZFAC=1.0D0                                                     
      ENDIF                                                             
      CZA=CZA+1.0D0                                                     
      AR=CZA                                                            
      IF(CDABS(CZA).LT.1.0D-16) RETURN                                  
      IF(CDABS(CZA).GT.15.0D0.AND.AR.GT.0.0D0) GO TO 2                  
      GO TO 1                                                           
C  ****  STIRLING'S EXPANSION OF CDLOG(GAMMA(CZA)).                     
    2 CZI2=1.0D0/(CZA*CZA)                                              
      CZS=(43867.0D0/244188.0D0)*CZI2                                   
      CZS=(CZS-3617.0D0/122400.0D0)*CZI2                                
      CZS=(CZS+1.0D0/156.0D0)*CZI2                                      
      CZS=(CZS-691.0D0/360360.0D0)*CZI2                                 
      CZS=(CZS+1.0D0/1188.0D0)*CZI2                                     
      CZS=(CZS-1.0D0/1680.0D0)*CZI2                                     
      CZS=(CZS+1.0D0/1260.0D0)*CZI2                                     
      CZS=(CZS-1.0D0/360.0D0)*CZI2                                      
      CZS=(CZS+1.0D0/12.0D0)/CZA                                        
      CLGAM=(CZA-0.5D0)*CDLOG(CZA)-CZA+9.1893853320467274D-1+CZS        
     1     +CZFL+CDLOG(CZFAC)                                           
      IF(ICONJ.EQ.1) CLGAM=DCONJG(CLGAM)                                
      RETURN                                                            
      END                       
