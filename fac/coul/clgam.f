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
