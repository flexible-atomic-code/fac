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
