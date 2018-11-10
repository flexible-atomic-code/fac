      DOUBLE PRECISION FUNCTION EXPINT(X,N)                        
C   FOR A DOUBLE PRECISION ARGUMENT X.GT.0.0 AND INTEGER N.GT.0         
C   THIS FUNCTION CALCULATES THE EXPONENTIAL-INTEGRAL FUNCTION          
C   E(N,X) = INTEGRAL (FROM 1 TO INFINITY) OF EXP(-X*T) * T**(-N) DT.   
C   CALLS DOUBLE PRECISION FUNCTION EIONE(IOPT,X,IOUT) TO CALCULATE     
C   THE CANONICAL EXPONENTIAL INTEGRAL FUNCTION E1(X) = E(1,X).         
C   IOUT IS THE LOGICAL UNIT NUMBER FOR WRITING (INPUT).                
C   ERROR FLAGS:                                                        
C     STOP 1  IF N .LE. 0,                                              
C     STOP 2  IF X .LE. 0.0                                             
C   S. BIENSTOCK, VERSION OF DECEMBER 5, 1982.                          
      IMPLICIT REAL*8(A-H,O-Z)                                          
      IF(N.LE.0) STOP 1                                                 
      IF(X.LE.0.D0) STOP 2                                              
      IER=0      
      if (X .GT. 10.0) then
         X2 = X*X
         X3 = X2*X
         X4 = X2*X2
         A1 = 8.5733287401
         A2 = 18.059016973
         A3 = 8.6347608925
         A4 = 0.2677737343
         B1 = 9.5733223454
         B2 = 25.6329561486
         B3 = 21.0996530827
         B4 = 3.9584969228
         EXPINT = X4 + A1*X3 + A2*X2 + A3*X1 + A4
         EXPINT = EXPINT/(X4 + B1*X3 + B2*X2 + B3*X1 + B4)
         EXPINT = EXPINT/X
      ELSE IF (X .LT. 1D-3) THEN
         A0 = -0.57721566
         A1 = 0.99999193
         A2 = -0.24991055
         A3 = 0.05519968
         A4 = -0.00976004
         A5 = 0.00107857
         X2 = X*X
         X3 = X2*X
         X4 = X3*X
         X5 = X4*X
         EXPINT = A0+A1*X+A2*X2+A3*X3+A4*X4+A5*X5-DLOG(X)
         EXPINT = EXPINT*DEXP(X)
      ELSE
         EXPINT=DEXP(DLOG(EIONE(2,X))+X)
      ENDIF
      IF(N.EQ.1) RETURN                                                 
      N1=N-1                                                            
      DO 10 I=1,N1                                                      
   10 EXPINT=(1D0-X*EXPINT)/DFLOAT(I)                              
      RETURN                                                            
      END
      
