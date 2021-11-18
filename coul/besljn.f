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
C     ABRAMOWITZ AND STEGUN, EQ. 10.1.15.
      IF(JY > 2) THEN
         JY=JY-2
         IM=1
      ELSE
         IM=0
      ENDIF
      IF(JY.EQ.2) THEN                                                 
        NL=-N-1                                                        
        IPH=2*MOD(IABS(N),2)-1                                         
      ELSE                                                             
        NL=N                                                           
        IPH=1                                                          
      ENDIF  
      X2 = X*X
C  ****  SELECTION OF CALCULATION MODE.                                
      IF(NL.LT.0) GO TO 10                                             
      IF(X.GT.1.0D0*NL) GO TO 7                                        
      XI=X2                                                          
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
      IF (IM .EQ. 0) THEN
         BESLJN=1.0D0
      ELSE
         IF (XI .LT. 1E-20) THEN
            BESLJN = 1D0/(IP+2)
            RETURN
         ENDIF
         BESLJN=0D0
      ENDIF
      PS=1.0D0                                                         
      DO 2 I=1,500                                                     
      IP=IP+2                                                          
      PS=-PS*XI/(I*IP)                                                 
      BESLJN=BESLJN+PS                                                   
      IF(DABS(PS).LT.1.0D-18*DABS(BESLJN)) GO TO 3                      
    2 CONTINUE                                                         
    3 CONTINUE
      IF (IM .EQ. 0) THEN
         BESLJN=IPH*F1*BESLJN
      ELSE
         BESLJN=-BESLJN/X2
      ENDIF
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
      IF (IM .EQ. 1) THEN                            
         F1=1.0D0                                                         
         IP=1      
         IF(NL.NE.0) THEN
            DO 71 I=1,NL                                                    
               IP=IP+2                                                        
 71            F1=F1*X/IP                                                     
         ENDIF           
         BESLJN = (1D0 - BESLJN/F1)/X2
      ENDIF
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
      IF (IM .EQ. 1) THEN                            
         F1=1.0D0                                                         
         IP=1      
         IF(NL.NE.0) THEN
            DO 81 I=1,NL                                                    
               IP=IP+2                                                        
 81            F1=F1*X/IP                                                     
         ENDIF           
         BESLJN = (1D0 - BESLJN/F1)/X2
      ENDIF
      RETURN                                                           
C  ****  RECURRENCE RELATION FOR NEGATIVE ORDERS.                      
C        ABRAMOWITZ AND STEGUN, EQ. 10.1.19.                           
 10   NL=IABS(NL)
C  ****  FOR SMALL ARG, NEED POWER SERIES EXPANSION
C  ****  POWER SERIES FOR SMALL ARGUMENTS.         
C     ABRAMOWITZ AND STEGUN, EQ. 10.1.2.
      XI=X2
      IF (N .GT. 0) THEN
         F1 = 0.1*(2*N-1)
      ELSE
         F1 = 0.1
      ENDIF
      IF ((N .EQ. 0) .OR. (XI .GT. F1)) THEN
         GOTO 30
      ENDIF
      IF (IM .EQ. 0) THEN
         F1=1/XI                                                         
         IP=1
         DO 21 I=1,N-1                                                    
            IP=IP+2                                                        
 21         F1=F1*IP/X                                                     
            F1 = -F1
      ENDIF
      XI=0.5D0*XI
      IF (IM .EQ. 0) THEN
         BESLJN=1.0D0
      ELSE
         IF (XI .LT. 1E-20) THEN
            BESLJN = 1D0/(1-2*N)
            RETURN
         ENDIF
         BESLJN=0D0
      ENDIF
      PS=1.0D0
      IP = -1 - 2*N
      DO 22 I=1,100                                                     
      IP=IP+2                                                          
      PS=-PS*XI/(I*IP)                                                 
      BESLJN=BESLJN+PS
      IF(DABS(PS).LT.1.0D-18*DABS(BESLJN)) GO TO 23                      
 22   CONTINUE                                                         
 23   CONTINUE
      IF (IM .EQ. 0) THEN
         BESLJN=F1*BESLJN
      ELSE
         BESLJN=-BESLJN/X2
      ENDIF
      RETURN
C     this block appears to cause problem
C      IF(X.LT.7.36D-1*(NL+1)*1.0D-35**(1.0D0/(NL+1))) THEN             
C        BESLJN=-1.0D35                                                  
C        RETURN                                                         
C      ENDIF
 30   XI=1.0D0/X           
      F3=XI*DSIN(X)                                                    
      F2=XI*(F3-DCOS(X))                                               
      IP=3                                                             
      DO 31 I=1,NL                                                     
      F1=F2                                                            
      F2=F3                                                            
      IP=IP-2                                                          
      F3=IP*XI*F2-F1                                                   
      IF(DABS(F3).GT.1.0D35) THEN                                      
        BESLJN=-1.0D35                                                  
        RETURN                                                         
      ENDIF                                                            
   31 CONTINUE                                                         
      BESLJN=IPH*F3
      IF (IM .EQ. 1) THEN
         XI=X2
         F1=1/XI                                                         
         IP=1
         DO 91 I=1,N-1                                                    
            IP=IP+2                                                        
 91         F1=F1*IP/X                                                     
            F1 = -F1
         BESLJN=(1-BESLJN/F1)/X2
      ENDIF                                                     
      RETURN                                                           
      END                                 
