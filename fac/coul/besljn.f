C  **************************************************************      
C                         FUNCION BESLJN                                
C  **************************************************************      
      FUNCTION BESLJN(JY,N,X)                                           
C                                                                      
C      THIS FUNCTION COMPUTES THE SPHERICAL BESSEL FUNCTIONS OF        
C   THE FIRST KIND AND SPHERICAL BESSEL FUNCTIONS OF THE SECOND        
C   KIND (ALSO KNOWN AS SPHERICAL NEUMANN FUNCTIONS) FOR REAL          
C   POSITIVE ARGUMENTS.                                                
C                                                                      
C      INPUT:                                                          
C         JY ...... KIND: 1(BESSEL) OR 2(NEUMANN).                     
C         N ....... ORDER (INTEGER).                                   
C         X ....... ARGUMENT (REAL AND POSITIVE).                      
C                                                                      
C   REF.: M. ABRAMOWITZ AND I.A. STEGUN, 'HANDBOOK OF MATHEMATI-       
C         CAL FUNCTIONS'. DOVER, NEW YORK (1974). PP 435-478.          
C                                                                      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                              
      PARAMETER (PI=3.1415926535897932D0,TPI=PI+PI)                    
      IF(X.LT.0) THEN                                                  
        WRITE(6,1000)                                                  
 1000   FORMAT(1X,'*** NEGATIVE ARGUMENT IN FUNCTION BESLJN.')          
        STOP                                                           
      ENDIF                                                            
C  ****  ORDER AND PHASE CORRECTION FOR NEUMANN FUNCTIONS.             
C        ABRAMOWITZ AND STEGUN, EQ. 10.1.15.                           
      IF(JY.EQ.2) THEN                                                 
        NL=-N-1                                                        
        IPH=2*MOD(IABS(N),2)-1                                         
      ELSE                                                             
        NL=N                                                           
        IPH=1                                                          
      ENDIF                                                            
C  ****  SELECTION OF CALCULATION MODE.                                
      IF(NL.LT.0) GO TO 10                                             
      IF(X.GT.1.0D0*NL) GO TO 7                                        
      XI=X*X                                                          
      IF(XI.GT.NL+NL+3.0D0) GO TO 4                                    
C  ****  POWER SERIES FOR SMALL ARGUMENTS AND POSITIVE ORDERS.         
C        ABRAMOWITZ AND STEGUN, EQ. 10.1.2.                            
      F1=1.0D0                                                         
      IP=1                                                             
      IF(NL.NE.0) THEN                                                 
        DO 1 I=1,NL                                                    
        IP=IP+2                                                        
    1   F1=F1*X/IP                                                     
      ENDIF                                                            
      XI=0.5D0*XI                                                      
      BESLJN=1.0D0                                                      
      PS=1.0D0                                                         
      DO 2 I=1,500                                                     
      IP=IP+2                                                          
      PS=-PS*XI/(I*IP)                                                 
      BESLJN=BESLJN+PS                                                   
      IF(DABS(PS).LT.1.0D-18*DABS(BESLJN)) GO TO 3                      
    2 CONTINUE                                                         
    3 BESLJN=IPH*F1*BESLJN                                               
      RETURN                                                           
C  ****  MILLER'S METHOD FOR POSITIVE ORDERS AND INTERMEDIATE          
C        ARGUMENTS. ABRAMOWITZ AND STEGUN, EQ. 10.1.19.                
    4 XI=1.0D0/X                                                       
      F2=0.0D0                                                         
      F3=1.0D-35                                                       
      IP=2*(NL+31)+3                                                   
      DO 5 I=1,31                                                      
      F1=F2                                                            
      F2=F3                                                            
      IP=IP-2                                                          
      F3=IP*XI*F2-F1                                                   
      IF(DABS(F3).GT.1.0D30) THEN                                      
        F2=F2/F3                                                       
        F3=1.0D0                                                       
      ENDIF                                                            
    5 CONTINUE                                                         
      BESLJN=1.0D0                                                      
      F2=F2/F3                                                         
      F3=1.0D0                                                         
      DO 6 I=1,NL                                                      
      F1=F2                                                            
      F2=F3                                                            
      IP=IP-2                                                          
      F3=IP*XI*F2-F1                                                   
      IF(DABS(F3).GT.1.0D30) THEN                                      
        BESLJN=BESLJN/F3                                                 
        F2=F2/F3                                                       
        F3=1.0D0                                                       
      ENDIF                                                            
    6 CONTINUE                                                         
      BESLJN=IPH*XI*DSIN(X)*BESLJN/F3                                    
      RETURN                                                           
C  ****  RECURRENCE RELATION FOR ARGUMENTS GREATER THAN ORDER.         
C        ABRAMOWITZ AND STEGUN, EQ. 10.1.19.                           
    7 XI=1.0D0/X                                                       
      F3=XI*DSIN(X)                                                    
      IF(NL.EQ.0) GO TO 9                                              
      F2=F3                                                            
      F3=XI*(F2-DCOS(X))                                               
      IF(NL.EQ.1) GO TO 9                                              
      IP=1                                                             
      DO 8 I=2,NL                                                      
      F1=F2                                                            
      F2=F3                                                            
      IP=IP+2                                                          
    8 F3=IP*XI*F2-F1                                                   
    9 BESLJN=IPH*F3                                                     
      RETURN                                                           
C  ****  RECURRENCE RELATION FOR NEGATIVE ORDERS.                      
C        ABRAMOWITZ AND STEGUN, EQ. 10.1.19.                           
   10 NL=IABS(NL)                                                      
      IF(X.LT.7.36D-1*(NL+1)*1.0D-35**(1.0D0/(NL+1))) THEN             
        BESLJN=-1.0D35                                                  
        RETURN                                                         
      ENDIF                                                            
      XI=1.0D0/X                                                       
      F3=XI*DSIN(X)                                                    
      F2=XI*(F3-DCOS(X))                                               
      IP=3                                                             
      DO 11 I=1,NL                                                     
      F1=F2                                                            
      F2=F3                                                            
      IP=IP-2                                                          
      F3=IP*XI*F2-F1                                                   
      IF(DABS(F3).GT.1.0D35) THEN                                      
        BESLJN=-1.0D35                                                  
        RETURN                                                         
      ENDIF                                                            
   11 CONTINUE                                                         
      BESLJN=IPH*F3                                                     
      RETURN                                                           
      END                                 
