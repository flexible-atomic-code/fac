C  **************************************************************      
C                       SUBROUTINE DCOUL                               
C  **************************************************************      
      SUBROUTINE DCOUL(Z,E,K,R,FU,FL,GU,GL,ERR)                        
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
C                                                                      
      IF(E.LT.0.0001D0.OR.K.EQ.0) THEN                                 
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
      RETURN                                                           
C                                                                      
C  ************  Z=0. SPHERICAL BESSEL FUNCTIONS.                      
C                                                                      
    1 CONTINUE                                                         
      RLAMB=IABS(K)                                                    
      CALL FCOUL(0.0D0,RLAMB,X,F,FP,G,GP,ERR)                          
      IF(ERR.GE.1.0D-6) RETURN                                         
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
      RETURN                                                           
      END                               
