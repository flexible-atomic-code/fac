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
