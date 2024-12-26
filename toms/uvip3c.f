
      SUBROUTINE  UVIP3C(NP,ND,XD,YD,C1,C2,C3)

C Specification statement
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      DIMENSION   XD(*),YD(*),C1(*),C2(*),C3(*)
C Error check
C 10   IF (ND.LE.1)   GO TO 90
C      DO 11  ID=2,ND
C        IF (XD(ID).LE.XD(ID-1))     GO TO 92
C 11   CONTINUE

C Branches off special cases.
C Take calre of the linear case
      IF (ND .GE. 3 .AND. NP .EQ. 1) GOTO 30
C Take Care of the quadratic interpolation case.
      IF (ND .GE. 4 .AND. NP .EQ. 2) GOTO 30
C Now for higher polynomials.
      IF (ND.LE.4)   GO TO 50
C General case  --  Five data points of more
C Calculates some local variables.
      NP = MIN(3,NP)
 20   NP0=MAX(3,NP)
      NPM1=NP0-1
      RENPM1=NPM1
      RENNM2=NP0*(NP0-2)

C End of locating the interval of interest
C Interpolation or extrapolation in one of the three subcases
C Subcase 1  --  Linear extrapolation when the abscissa of the
C                desired point is equal to that of the first data
C                point or less.
C Estimates the first derivative when the interval is not the
C same as the one for the previous desired point.  --
C cf. Equation (8)
      X0=XD(1)
      X1=XD(2)-X0
      Y0=YD(1)
      Y1=YD(2)-Y0
      IF (NP .GE. 2) THEN
         X2=XD(3)-X0
         Y2=YD(3)-Y0
         IF (NP .GE. 3) THEN
            Y3=YD(4)-Y0    
            X3=XD(4)-X0
         ENDIF
      ENDIF
      IF (NP .EQ. 1) THEN
         A1 = Y1/X1
      ELSE IF (NP .EQ. 2) THEN
         DLT = X1*X2*(X2-X1)
         A1 = (X2*X2*Y1 - X1*X1*Y2)/DLT
      ELSE
         DLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
         A1=(((X2*X3)**2)*(X3-X2)*Y1
     1        +((X3*X1)**2)*(X1-X3)*Y2
     2        +((X1*X2)**2)*(X2-X1)*Y3)/DLT
      ENDIF
C Evaluates the YI value.
      C1(ND) = A1
C End of Subcase 1

C Subcase 2  --  Linear extrapolation when the abscissa of the
C                desired point is equal to that of the last data
C                point or greater.
C Estimates the first derivative when the interval is not the
C same as the one for the previous desired point.  --
C cf. Equation (8)
      X0=XD(ND)
      X1=XD(ND-1)-X0
      Y0=YD(ND)
      Y1=YD(ND-1)-Y0
      IF (NP .GE. 1) THEN
         X2=XD(ND-2)-X0
         Y2=YD(ND-2)-Y0
         IF (NP .GE. 2) THEN
            X3=XD(ND-3)-X0
            Y3=YD(ND-3)-Y0
         ENDIF
      ENDIF
      IF (NP .EQ. 1) THEN
         A1 = Y1/X1             
      ELSE IF (NP .EQ. 2) THEN
         DLT = X1*X2*(X2-X1)
         A1 = (X2*X2*Y1 - X1*X1*Y2)/DLT  
      ELSE
         DLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
         A1=(((X2*X3)**2)*(X3-X2)*Y1
     1        +((X3*X1)**2)*(X1-X3)*Y2
     2        +((X1*X2)**2)*(X2-X1)*Y3)/DLT
      ENDIF
C Evaluates the YI value.
      C2(ND) = A1

C Main calculation for the general case
C First (outermost) DO-loop with respect to the desired points
 30   DO 39  IINT=1,ND-1
C ADD QUADRATIC INTERPOLATION WHEN NP=2, WHICH ONLY 
C WORKS IF ND >= 4
         IF (NP .EQ. 2) THEN
            IF (IINT .EQ. 1) THEN
               X0 = XD(IINT)
               Y0 = YD(IINT)
               X1 = XD(IINT+1)-X0
               Y1 = YD(IINT+1)-Y0
               X2 = XD(IINT+2)-X0
               Y2 = YD(IINT+2)-Y0
               DLT = X1*X2*(X2-X1)
               A1=(X2*X2*Y1-X1*X1*Y2)/DLT
               A2=(X1*Y2-X2*Y1)/DLT
            ELSE IF (IINT .EQ. ND-1) THEN
               X0 = XD(IINT-1)
               Y0 = YD(IINT-1)
               X1 = XD(IINT)-X0
               Y1 = YD(IINT)-Y0
               X2 = XD(IINT+1)-X0
               Y2 = YD(IINT+1)-Y0
               DLT = X1*X2*(X2-X1)
               A1=(X2*X2*Y1-X1*X1*Y2)/DLT
               A2=(X1*Y2-X2*Y1)/DLT
            ELSE
               X0 = XD(IINT)
               Y0 = YD(IINT)
               X1 = XD(IINT+1)-X0
               Y1 = YD(IINT+1)-Y0
               X2 = XD(IINT+2)-X0
               Y2 = YD(IINT+2)-Y0
               DLT = X1*X2*(X2-X1)
               A1=(X2*X2*Y1-X1*X1*Y2)/DLT
               A2=(X1*Y2-X2*Y1)/DLT
               X00 = XD(IINT-1)
               Y00 = YD(IINT-1)
               X11 = XD(IINT)-X00
               Y11 = YD(IINT)-Y00
               X22 = XD(IINT+1)-X00
               Y22 = YD(IINT+1)-Y00
               DLT = X11*X22*(X22-X11)
               A11=(X22*X22*Y11-X11*X11*Y22)/DLT
               A22=(X11*Y22-X22*Y11)/DLT
            ENDIF

            IF (IINT .NE. 1 .AND. IINT .NE. ND-1) THEN
               XX = X0-X00
               C1(IINT) = 0.5*A1 + A11 + 2*A22
               C2(IINT) = 0.5*A2 + A22
            ELSE
               C1(IINT) = A1
               C2(IINT) = A2            
            ENDIF
         
            GOTO 39
         ENDIF
C Here is linear interpolation when NP=1
         IF (NP .EQ. 1) THEN
            X0 = XD(IINT)
            Y0 = YD(IINT)
            X1 = XD(IINT+1)-X0
            Y1 = YD(IINT+1)-Y0
            A1 = Y1/X1
            C1(IINT) = A1
            GOTO 39
         ENDIF
C Subcase 3  --  Interpolation when the abscissa of the desired
C                point is  between those of the first and last
C                data points.
C Calculates the coefficients of the third-degree polynomial (for
C NP.LE.3) or the factors for the higher-degree polynomials (for
C NP.GT.3), when the interval is not the same as the one for the
C previous desired point.
C The second DO-loop with respect to the two endpoints of the
C interval
         DO 37  IEPT=1,2
C     Calculates the estimate of the first derivative at an endpoint.
C Initial setting for calculation
            ID0=IINT+IEPT-1
            X0=XD(ID0)
            Y0=YD(ID0)
            SMPEF=0.0
            SMWTF=0.0
            SMPEI=0.0
            SMWTI=0.0
C The third (innermost) DO-loop with respect to the four primary
C estimate of the first derivative
            DO 36  IPE=1,4
C Selects point numbers of four consecutive data points for
C calculating the primary estimate of the first derivative.
               IF (IPE.EQ.1)  THEN
                  ID1=ID0-3
                  ID2=ID0-2
                  ID3=ID0-1
               ELSE IF (IPE.EQ.2)  THEN
                  ID1=ID0+1
               ELSE IF (IPE.EQ.3)  THEN
                  ID2=ID0+2
               ELSE
                  ID3=ID0+3
               END IF
C     Checks if any point number falls outside the legitimate range
C (between 1 and ND).  Skips calculation of the primary estimate
C if any does.
               IF (ID1.LT.1.OR.ID2.LT.1.OR.ID3.LT.1.OR.
     1              ID1.GT.ND.OR.ID2.GT.ND.OR.ID3.GT.ND)
     2              GO TO 36
C Calculates the primary estimate of the first derivative  --
C cf. Equation (8)
               X1=XD(ID1)-X0
               X2=XD(ID2)-X0
               X3=XD(ID3)-X0
               Y1=YD(ID1)-Y0
               Y2=YD(ID2)-Y0
               Y3=YD(ID3)-Y0
               DLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
               PE=(((X2*X3)**2)*(X3-X2)*Y1
     1              +((X3*X1)**2)*(X1-X3)*Y2
     2              +((X1*X2)**2)*(X2-X1)*Y3)/DLT
C Calculates the volatility factor, VOL, and distance factor,
C SXX, for the primary estimate.  --  cf. Equations (9) and (11)
               SX=X1+X2+X3
               SY=Y1+Y2+Y3
               SXX=X1*X1+X2*X2+X3*X3
               SXY=X1*Y1+X2*Y2+X3*Y3
               DNM=4.0*SXX-SX*SX
               B0=(SXX*SY-SX*SXY)/DNM
               B1=(4.0*SXY-SX*SY)/DNM
               DY0=-B0
               DY1=Y1-(B0+B1*X1)
               DY2=Y2-(B0+B1*X2)
               DY3=Y3-(B0+B1*X3)
               VOL=DY0*DY0+DY1*DY1+DY2*DY2+DY3*DY3
C Calculates the EPSLN value, which is used to decide whether or
C not the volatility factor, VOL, is essentially zero.
               EPSLN=(YD(ID0)**2+YD(ID1)**2
     1              +YD(ID2)**2+YD(ID3)**2)*1.0E-12
C Accumulates the weighted primary estimates.  --
C cf. Equations (13) and (14)
               IF (VOL.GT.EPSLN)  THEN
C - For finite weight.
                  WT=1.0/(VOL*SXX)
                  SMPEF=SMPEF+PE*WT
                  SMWTF=SMWTF+WT
               ELSE
C - For infinite weight.
                  SMPEI=SMPEI+PE
                  SMWTI=SMWTI+1.0
               END IF
 36         CONTINUE
C End of the third DO-loop
C Calculates the final estimate of the first derivative.  --
C cf. Equation (14)
            IF (SMWTI.LT.0.5)  THEN
C - When no infinite weights exist.
               YP=SMPEF/SMWTF
            ELSE
C - When infinite weights exist.
               YP=SMPEI/SMWTI
            END IF
            IF (IEPT.EQ.1)  THEN
               YP0=YP
            ELSE
               YP1=YP
            END IF
C End of the calculation of the estimate of the first derivative
C at an endpoint
 37      CONTINUE
C End of the second DO-loop

C Calculates the coefficients of the third-degree polynomial
C (when NP.LE.3).  --  cf. Equation (4)
         DX=XD(IINT+1)-XD(IINT)
         DY=YD(IINT+1)-YD(IINT)
         A0=YD(IINT)
         A1=YP0
         YP1=YP1-YP0
         YP0=YP0-DY/DX
         A2=-(3.0*YP0+YP1)/DX
         A3= (2.0*YP0+YP1)/(DX*DX)

         C1(IINT) = A1
         C2(IINT) = A2
         C3(IINT) = A3
C End of the calculation of the coefficients of the third-degree
C polynomial (when NP.LE.3) or the factors for the higher-degree
C polynomials (when NP.GT.3), when the interval is not the same
C as the one for the previous desired point.
C End of Subcase 3
 39   CONTINUE
C End of the first DO-loop
C End of general case
      RETURN
C Special cases  --  Four data points or less
C Preliminary processing for the special cases
 50   X0=XD(1)
      Y0=YD(1)
      X1=XD(2)-X0
      Y1=YD(2)-Y0
      IF (ND.EQ.2)   GO TO 60
      X2=XD(3)-X0
      Y2=YD(3)-Y0
      IF (ND.EQ.3)   GO TO 70
      X3=XD(4)-X0
      Y3=YD(4)-Y0
      GO TO 80
C     Special Case 1  --  Two data points
C     (Linear interpolation and extrapolation)
 60   A1=Y1/X1
      C1(1) = A1
      C1(2) = A1 
      C2(2) = A1
      IF (NP .GT. 1) THEN
         C2(1) = 0.0
         IF (NP .GT. 2) THEN
            C3(1) = 0.0
            C3(2) = 0.0
         ENDIF
      ENDIF
C End of Special Case 1
      RETURN
C Special Case 2  --  Three data points
C (Quadratic interpolation and linear extrapolation)
 70   DLT=X1*X2*(X2-X1)
      A1=(X2*X2*Y1-X1*X1*Y2)/DLT
      A2=(X1*Y2-X2*Y1)/DLT
      A12=2.0*A2*X2+A1
      C1(1) = A1
      C1(2) = A1      
      C1(3) = A1
      C2(3) = A12
      IF (NP .GT. 1) THEN
         C2(1) = A2
         C2(2) = A2
         IF (NP .GT. 2) THEN
            DO 71 IINT = 1, ND
               C3(IINT) = 0.0
 71         CONTINUE
         ENDIF
      ENDIF
C End of Special Case 2
      RETURN
C Special Case 3  --  Four data points
C (Cubic interpolation and linear extrapolation)
   80 DLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
      A1=(((X2*X3)**2)*(X3-X2)*Y1
     1   +((X3*X1)**2)*(X1-X3)*Y2
     2   +((X1*X2)**2)*(X2-X1)*Y3)/DLT
      A2=(X2*X3*(X2*X2-X3*X3)*Y1
     1   +X3*X1*(X3*X3-X1*X1)*Y2
     2   +X1*X2*(X1*X1-X2*X2)*Y3)/DLT
      A3=(X2*X3*(X3-X2)*Y1
     1   +X3*X1*(X1-X3)*Y2
     2   +X1*X2*(X2-X1)*Y3)/DLT
      A13=(3.0*A3*X3+2.0*A2)*X3+A1
      DO 81 IINT = 1, ND
         C1(IINT) = A1
 81   CONTINUE
      C2(ND) = A13
      IF (NP .GT. 1) THEN
         DO 82 IINT = 1, ND-1
            C2(IINT) = A2
 82      CONTINUE
         IF (NP .GT. 2) THEN
            DO 83 IINT = 1, ND
               C3(IINT) = A3
 83         CONTINUE
         ENDIF
      ENDIF
C End of Special Case 3
      RETURN
C Error exit
   90 WRITE (*,99090) ND
      GO TO 99
   92 WRITE (*,99092) ID,XD(ID-1),XD(ID)
   99 WRITE (*,99099)
      STOP
      RETURN
C Format statements for error messages
99090 FORMAT (1X/ ' ***   Insufficient data points.'
     1  7X,'ND =',I3)
99092 FORMAT (1X/ ' ***   Two data points identical or out of ',
     1  'sequence.'/
     2  7X,'ID, XD(ID-1), XD(ID) =',I5,2F10.3)
99099 FORMAT (' Error detected in the UVIP3C subroutine'/)
      END

