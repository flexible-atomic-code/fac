CCCCCCCCCC         COULOMB AND BESSEL FUNCTIONS        CCCCCCCCCC      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C  **************************************************************      
C                       SUBROUTINE SCOUL                               
C  **************************************************************      
      SUBROUTINE SCOUL(Z,E,L,R,F,FP,G,GP,ERR)                          
C                                                                      
C     THIS SUBROUTINE COMPUTES RADIAL SCHRODINGER-COULOMB WAVE         
C  FUNCTIONS FOR FREE STATES.                                          
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
C     F, FP .... REGULAR RADIAL SCHRODINGER-COULOMB FUNCTION AND       
C                ITS DERIVATIVE.                                       
C     G, GP .... IRREGULAR RADIAL SCHRODINGER-COULOMB FUNCTION         
C                AND ITS DERIVATIVE.                                   
C     ERR ...... ACCURACY OF THE COMPUTED FUNCTIONS (RELATIVE UN-      
C                CERTAINTY).                                           
C  OUTPUT THROUGH COMMON/OCOUL/:                                       
C     WAVNUM ... WAVE NUMBER.                                          
C     ETA ...... SOMMERFELD'S PARAMETER.                               
C     DELTA .... COULOMB PHASE SHIFT (MODULUS 2*PI).                   
C                                                                      
C     RADIAL FUNCTIONS ARE NORMALIZED SO THAT, FOR LARGE R, THEY       
C  OSCILLATE WITH UNIT AMPLITUDE.                                      
C                                                                      
C     OTHER SUBPROGRAMS REQUIRED: SUBROUTINES FCOUL AND SUM2F0,        
C                                 AND FUNCTION CLGAM.                  
C                                                                      
      IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)          
      COMMON/OFCOUL/DELTA0                                             
      COMMON/OCOUL/WAVNUM,ETA,DELTA                                    
C                                                                      
      IF(E.LT.0.0001D0.OR.L.LT.0) THEN                                 
        F=0.0D0                                                        
        FP=0.0D0                                                       
        G=1.0D35                                                       
        GP=-1.0D35                                                     
        ERR=1.0D0                                                      
        IF(E.LT.0.0001D0) WRITE(6,2101)                                
 2101   FORMAT(1X,'*** ERROR IN SCOUL: E IS TOO SMALL.')               
        IF(L.LT.0) WRITE(6,2102)                                       
 2102   FORMAT(1X,'*** ERROR IN SCOUL: L.LT.0.')                       
        RETURN                                                         
      ENDIF                                                            
C                                                                      
C  ****  PARAMETERS.                                                   
C                                                                      
      WAVNUM=DSQRT(E+E)                                                
      IF(DABS(Z).GT.0.00001D0) THEN                                    
        ETA=Z/WAVNUM                                                   
        ICAL=0                                                         
      ELSE                                                             
        ETA=0.0D0                                                      
        ICAL=1                                                         
      ENDIF                                                            
      RLAMB=L                                                          
      X=WAVNUM*R                                                       
      IF(ICAL.EQ.1) GO TO 1                                            
C                                                                      
C  ************  COULOMB FUNCTIONS.                                    
C                                                                      
      DELTA0=1.0D30                                                    
      CALL FCOUL(ETA,RLAMB,X,F,FP,G,GP,ERR)                            
      FP=FP*WAVNUM                                                     
      GP=GP*WAVNUM                                                     
      IF(DELTA0.LT.1.0D29) THEN                                        
        DELTA=DELTA0                                                   
      ELSE                                                             
        DELTA=DELTAC(ETA,RLAMB)                                        
      ENDIF                                                            
      RETURN                                                           
C                                                                      
C  ************  Z=0. SPHERICAL BESSEL FUNCTIONS.                      
C                                                                      
    1 CONTINUE                                                         
      F=X*BESLJN(1,L,X)                                                 
      G=-X*BESLJN(2,L,X)                                                
      FP=((L+1)*BESLJN(1,L,X)-X*BESLJN(1,L+1,X))*WAVNUM                  
      GP=-((L+1)*BESLJN(2,L,X)-X*BESLJN(2,L+1,X))*WAVNUM                 
      DELTA=0.0D0                                                      
      ERR=0.0D0                                                        
      RETURN                                                           
      END                
