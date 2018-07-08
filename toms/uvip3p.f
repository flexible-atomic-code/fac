
      SUBROUTINE  UVIP3P(NP,ND,XD,YD,NI,XI, YI)
C
C Univariate Interpolation (Improved Akima Method)
C
C Hiroshi Akima
C U.S. Department of Commerce, NTIA/ITS
C Version of 89/07/04
C
C This subroutine performs univariate interpolation.  It is based
C on the improved A method developed by Hiroshi Akima, 'A method
C of univariate interpolation that has the accuracy of a third-
C degree polynomial,' ACM TOMS, vol. xx, pp. xxx-xxx, 19xx.  (The
C equation numbers referred to in the comments below are those in
C the paper.)
C
C In this method, the interpolating function is a piecewise
C function composed of a set of polynomials applicable to
C successive intervals of the given data points.  This method
C uses third-degree polynomials as the default, but the user has
C an option to use higher-degree polynomial to reduce undulations
C in resulting curves.
C
C This method has the accuracy of a third-degree polynomial if
C the degree of the polynomials for the interpolating function is
C set to three.
C
C The input arguments are
C   NP = degree of the polynomials for the interpolating
C        function,
C   ND = number of input data points
C        (must be equal to 2 or greater),
C   XD = array of dimension ND, containing the abscissas of
C        the input data points
C        (must be in a monotonic increasing order),
C   YD = array of dimension ND, containing the ordinates of
C        the input data points,
C   NI = number of points for which interpolation is desired
C        (must be equal to 1 or greater),
C   XI = array of dimension NI, containing the abscissas of
C        the desired points.
C
C The output argument is
C   YI = array of dimension NI, where the ordinates of the
C        desired points are to be stored.
C
C If an integer value smaller than 3 is given to the NP argument,
C this subroutine assumes NP = 3.
C
C The XI array elements need not be monotonic, but this
C subroutine interpolates faster if the XI array elements are
C given in a monotonic order.
C
C If the XI array element is less than XD(1) or greater than
C XD(ND), this subroutine linearly interpolates the YI value.
C
C
C Specification statement
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      DIMENSION   XD(*),YD(*),XI(*), YI(*)
      DOUBLE PRECISION XT(8), YT(8), AN(2)
C Error check
   10 IF (ND.LE.1)   GO TO 90
      IF (NI.LE.0)   GO TO 91
      DO 11  ID=2,ND
        IF (XD(ID).LE.XD(ID-1))     GO TO 92
 11   CONTINUE

c      call uvip3o(NP,ND,XD,YD,NI,XI, YI)
c      if (ni .gt. 1) then
c         return
c      endif
c      yio = yi(1)
      IF (ND .LE. 8) THEN
         CALL SCALEXY(AN, XT, YT, 1, ND, XD, YD)
         NDT = ND
      ENDIF
C Branches off special cases.
C Take calre of the linear case
      IF (ND .GE. 3 .AND. NP .EQ. 1) GOTO 30
C Take Care of the quadratic interpolation case.
      IF (ND .GE. 4 .AND. NP .EQ. 2) GOTO 30
C Now for higher polynomials.
      IF (ND.LE.4)   GO TO 50
C General case  --  Five data points of more
C Calculates some local variables.
   20 NP0=MAX(3,NP)
      NPM1=NP0-1
      RENPM1=NPM1
      RENNM2=NP0*(NP0-2)
C Main calculation for the general case
C First (outermost) DO-loop with respect to the desired points
   30 DO 39  II=1,NI
        IF (II.EQ.1)      IINTPV=-1
        XII=XI(II)
C Locates the interval that includes the desired point by binary
C search.
        IF (XII.LE.XD(1))  THEN
          IINT=0
        ELSE IF (XII.LT.XD(ND))  THEN
          IDMN=1
          IDMX=ND
          IDMD=(IDMN+IDMX)/2
   31     IF (XII.GE.XD(IDMD))  THEN
            IDMN=IDMD
          ELSE
            IDMX=IDMD
          END IF
          IDMD=(IDMN+IDMX)/2
          IF (IDMD.GT.IDMN)    GO TO 31
          IINT=IDMD
        ELSE
          IINT=ND
       END IF
       JINT = IINT
       IF (ND .GT. 8) THEN
          II0 = IINT - 3
          IF (IINT .GT. 3) THEN
             II0 = IINT - 3
             JINT = 4
          ELSE
             II0 = 1
             JINT = IINT
          ENDIF
          II1 = II0 + 7
          IF (II1 .GT. ND) THEN
             JINT = JINT + (II1-ND)
             II1 = ND
             II0 = ND-7
          ENDIF
          CALL SCALEXY(AN, XT, YT, II0, II1, XD, YD)
          NDT = 8
       ENDIF
       XII = XT(1) + (XII-XT(1))*AN(1)
C End of locating the interval of interest
C Interpolation or extrapolation in one of the three subcases
       IF (IINT.LE.0)  THEN
C Subcase 1  --  Linear extrapolation when the abscissa of the
C                desired point is equal to that of the first data
C                point or less.
C Estimates the first derivative when the interval is not the
C same as the one for the previous desired point.  --
C cf. Equation (8)
          IF (IINT.NE.IINTPV)  THEN
            IINTPV=IINT
            X0=XT(1)
            X1=XT(2)-X0
            Y0=YT(1)
            Y1=YT(2)-Y0
            IF (NP .GE. 2) THEN
               X2=XT(3)-X0
               Y2=YT(3)-Y0
               IF (NP .GE. 3) THEN
                  Y3=YT(4)-Y0    
                  X3=XT(4)-X0
               ENDIF
            ENDIF
            IF (NP .EQ. 1) THEN
               A1 = Y1/X1
            ELSE IF (NP .EQ. 2) THEN
               XDLT = X1*X2*(X2-X1)
               A1 = (X2*X2*Y1 - X1*X1*Y2)/XDLT
            ELSE
               XDLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
               A1=(((X2*X3)**2)*(X3-X2)*Y1
     1              +((X3*X1)**2)*(X1-X3)*Y2
     2              +((X1*X2)**2)*(X2-X1)*Y3)/XDLT
            ENDIF
         ENDIF
C Evaluates the YI value.
         YI(II)=(Y0+A1*(XII-X0))*AN(2)
C End of Subcase 1
      ELSE IF (IINT.GE.ND)  THEN
C Subcase 2  --  Linear extrapolation when the abscissa of the
C                desired point is equal to that of the last data
C                point or greater.
C Estimates the first derivative when the interval is not the
C same as the one for the previous desired point.  --
C cf. Equation (8)
         IF (IINT.NE.IINTPV)  THEN
            IINTPV=IINT
            X0=XT(NDT)
            X1=XT(NDT-1)-X0
            Y0=YT(NDT)
            Y1=YT(NDT-1)-Y0
            IF (NP .GE. 1) THEN
               X2=XT(NDT-2)-X0
               Y2=YT(NDT-2)-Y0
               IF (NP .GE. 2) THEN
                  X3=XT(NDT-3)-X0
                  Y3=YT(NDT-3)-Y0
               ENDIF
            ENDIF
            IF (NP .EQ. 1) THEN
               A1 = Y1/X1             
            ELSE IF (NP .EQ. 2) THEN
               XDLT = X1*X2*(X2-X1)
               A1 = (X2*X2*Y1 - X1*X1*Y2)/XDLT  
            ELSE
               XDLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
               A1=(((X2*X3)**2)*(X3-X2)*Y1
     1              +((X3*X1)**2)*(X1-X3)*Y2
     2              +((X1*X2)**2)*(X2-X1)*Y3)/XDLT
            ENDIF
          END IF
C Evaluates the YI value.
          YI(II)=(Y0+A1*(XII-X0))*AN(2)
C End of Subcase 2
        ELSE
C ADD QUADRATIC INTERPOLATION WHEN NP=2, WHICH ONLY 
C WORKS IF ND >= 4
           IF (NP .EQ. 2) THEN
              IF (IINT .NE. IINTPV) THEN
                 IINTPV = IINT
                 IF (IINT .EQ. 1) THEN
                    X0 = XT(JINT)
                    Y0 = YT(JINT)
                    X1 = XT(JINT+1)-X0
                    Y1 = YT(JINT+1)-Y0
                    X2 = XT(JINT+2)-X0
                    Y2 = YT(JINT+2)-Y0
                    XDLT = X1*X2*(X2-X1)
                    A1=(X2*X2*Y1-X1*X1*Y2)/XDLT
                    A2=(X1*Y2-X2*Y1)/XDLT
                 ELSE IF (IINT .EQ. ND-1) THEN
                    X0 = XT(JINT-1)
                    Y0 = YT(JINT-1)
                    X1 = XT(JINT)-X0
                    Y1 = YT(JINT)-Y0
                    X2 = XT(JINT+1)-X0
                    Y2 = YT(JINT+1)-Y0
                    XDLT = X1*X2*(X2-X1)
                    A1=(X2*X2*Y1-X1*X1*Y2)/XDLT
                    A2=(X1*Y2-X2*Y1)/XDLT
                 ELSE
                    X0 = XT(JINT)
                    Y0 = YT(JINT)
                    X1 = XT(JINT+1)-X0
                    Y1 = YT(JINT+1)-Y0
                    X2 = XT(JINT+2)-X0
                    Y2 = YT(JINT+2)-Y0
                    XDLT = X1*X2*(X2-X1)
                    A1=(X2*X2*Y1-X1*X1*Y2)/XDLT
                    A2=(X1*Y2-X2*Y1)/XDLT
                    X00 = XT(JINT-1)
                    Y00 = YT(JINT-1)
                    X11 = XT(JINT)-X00
                    Y11 = YT(JINT)-Y00
                    X22 = XT(JINT+1)-X00
                    Y22 = YT(JINT+1)-Y00
                    XDLT = X11*X22*(X22-X11)
                    A11=(X22*X22*Y11-X11*X11*Y22)/XDLT
                    A22=(X11*Y22-X22*Y11)/XDLT
                 ENDIF
              ENDIF
              XX = XII-X0              
              YI(II) = Y0 + XX*(A1+XX*A2)
              IF (IINT .NE. 1 .AND. IINT .NE. ND-1) THEN
                 XX = XII-X00
                 YI(II) = 0.5*(YI(II) + Y00+XX*(A11+XX*A22))
              ENDIF
              YI(II) = YI(II)*AN(2)
              GOTO 39
           ENDIF

C Here is linear interpolation when NP=1
           IF (NP .EQ. 1) THEN
              IF (IINT .NE. IINTPV) THEN
                 IINTPV = IINT
                 X0 = XT(JINT)
                 Y0 = YT(JINT)
                 X1 = XT(JINT+1)-X0
                 Y1 = YT(JINT+1)-Y0
                 A1 = Y1/X1
              ENDIF
              XX = XII-X0
              YI(II) = (Y0 + A1*XX)*AN(2)
              GOTO 39
           ENDIF

C Subcase 3  --  Interpolation when the abscissa of the desired
C                point is  between those of the first and last
C                data points.
C Calculates the coefficients of the third-degree polynomial (for
C NP.LE.3) or the factors for the higher-degree polynomials (for
C NP.GT.3), when the interval is not the same as the one for the
C previous desired point.
          IF (IINT.NE.IINTPV)  THEN
            IINTPV=IINT
C The second DO-loop with respect to the two endpoints of the
C interval
            DO 37  IEPT=1,2
C Calculates the estimate of the first derivative at an endpoint.
C Initial setting for calculation
              ID0=JINT+IEPT-1
              X0=XT(ID0)
              Y0=YT(ID0)
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
C Checks if any point number falls outside the legitimate range
C (between 1 and ND).  Skips calculation of the primary estimate
C if any does.
                IF (ID1.LT.1.OR.ID2.LT.1.OR.ID3.LT.1.OR.
     1              ID1.GT.NDT.OR.ID2.GT.NDT.OR.ID3.GT.NDT)
     2               GO TO 36
C Calculates the primary estimate of the first derivative  --
C cf. Equation (8)
                X1=XT(ID1)-X0
                X2=XT(ID2)-X0
                X3=XT(ID3)-X0
                Y1=YT(ID1)-Y0
                Y2=YT(ID2)-Y0
                Y3=YT(ID3)-Y0
                XDLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
                PE=(((X2*X3)**2)*(X3-X2)*Y1
     1             +((X3*X1)**2)*(X1-X3)*Y2
     2             +((X1*X2)**2)*(X2-X1)*Y3)/XDLT
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
                EPSLN=(YT(ID0)**2+YT(ID1)**2
     1                +YT(ID2)**2+YT(ID3)**2)*1.0E-12
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
   37       CONTINUE

C End of the second DO-loop
            IF (NP0.LE.3)  THEN
C Calculates the coefficients of the third-degree polynomial
C (when NP.LE.3).  --  cf. Equation (4)
              DX=XT(JINT+1)-XT(JINT)
              DY=YT(JINT+1)-YT(JINT)
              A0=YT(JINT)
              A1=YP0
              YP1=YP1-YP0
              YP0=YP0-DY/DX
              A2=-(3.0*YP0+YP1)/DX
              A3= (2.0*YP0+YP1)/(DX*DX)
            ELSE
C Calculates the factors for the higher-degree polynomials
C (when NP.GT.3).  --  cf. Equation (20)
              DX=XT(JINT+1)-XT(JINT)
              DY=YT(JINT+1)-YT(JINT)
              T0=YP0*DX-DY
              T1=YP1*DX-DY
              AA0= (T0+RENPM1*T1)/RENNM2
              AA1=-(RENPM1*T0+T1)/RENNM2
            END IF
          END IF
C End of the calculation of the coefficients of the third-degree
C polynomial (when NP.LE.3) or the factors for the higher-degree
C polynomials (when NP.GT.3), when the interval is not the same
C as the one for the previous desired point.
C Evaluates the YI value.
          IF (NP0.LE.3)  THEN
C     - With a third-degree polynomial.  --  cf. Equation (3)
             XX=XII-XT(JINT)
             YI(II)=A0+XX*(A1+XX*(A2+XX*A3))
          ELSE
C     - With a higher-degree polynomial.  --  cf. Equation (19)
             U=(XII-XT(JINT))/DX
             UC=1.0-U
             V=AA0*((U**NP0)-U)+AA1*((UC**NP0)-UC)
             YI(II)=YT(JINT)+DY*U+V
          END IF
          YI(II) = YI(II)*AN(2)
C End of Subcase 3
        END IF
   39 CONTINUE
C End of the first DO-loop
C     End of general case
      RETURN
C Special cases  --  Four data points or less
C Preliminary processing for the special cases
   50 X0=XT(1)
      Y0=YT(1)
      X1=XT(2)-X0
      Y1=YT(2)-Y0
      IF (ND.EQ.2)   GO TO 60
      X2=XT(3)-X0
      Y2=YT(3)-Y0
      IF (ND.EQ.3)   GO TO 70
      X3=XT(4)-X0
      Y3=YT(4)-Y0
      GO TO 80
C Special Case 1  --  Two data points
C (Linear interpolation and extrapolation)
   60 A1=Y1/X1
      DO 61  II=1,NI
         XII = XT(1) + (XI(II) - XT(1))*AN(1)
         YI(II)=(Y0+A1*(XII-X0))*AN(2)
   61 CONTINUE
C End of Special Case 1
      RETURN
C Special Case 2  --  Three data points
C (Quadratic interpolation and linear extrapolation)
 70   XDLT=X1*X2*(X2-X1)
      A1=(X2*X2*Y1-X1*X1*Y2)/XDLT
      A2=(X1*Y2-X2*Y1)/XDLT
      A12=2.0*A2*X2+A1
      DO 71  II=1,NI
         XII = XT(1) + (XI(II)-XT(1))*AN(1)
         XX=XII-X0
         IF (XX.LE.0.0)  THEN
            YI(II)=Y0+A1*XX
         ELSE IF (XX.LT.X2) THEN
            YI(II)=Y0+XX*(A1+XX*A2)
         ELSE
            YI(II)=Y0+Y2+A12*(XX-X2)
         END IF
         YI(II) = YI(II)*AN(2)
 71   CONTINUE
C End of Special Case 2
      RETURN
C Special Case 3  --  Four data points
C (Cubic interpolation and linear extrapolation)
   80 XDLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
      A1=(((X2*X3)**2)*(X3-X2)*Y1
     1   +((X3*X1)**2)*(X1-X3)*Y2
     2   +((X1*X2)**2)*(X2-X1)*Y3)/XDLT
      A2=(X2*X3*(X2*X2-X3*X3)*Y1
     1   +X3*X1*(X3*X3-X1*X1)*Y2
     2   +X1*X2*(X1*X1-X2*X2)*Y3)/XDLT
      A3=(X2*X3*(X3-X2)*Y1
     1   +X3*X1*(X1-X3)*Y2
     2   +X1*X2*(X2-X1)*Y3)/XDLT
      A13=(3.0*A3*X3+2.0*A2)*X3+A1
      DO 81  II=1,NI
         XII = XT(1) + (XI(II)-XT(1))*AN(1)
         XX=XII-X0
         IF (XX.LE.0.0)  THEN
            YI(II)=Y0+A1*XX
         ELSE IF (XX.LT.X3) THEN
            YI(II)=Y0+XX*(A1+XX*(A2+XX*A3))
         ELSE
            YI(II)=Y0+Y3+A13*(XX-X3)
         END IF
         YI(II) = YI(II)*AN(2)
   81 CONTINUE
C End of Special Case 3
      RETURN
C Error exit
   90 WRITE (*,99090) ND
      GO TO 99
   91 WRITE (*,99091) NI
      GO TO 99
   92 WRITE (*,99092) ID,XD(ID-1),XD(ID)
   99 WRITE (*,99099)
      STOP
      RETURN
C Format statements for error messages
99090 FORMAT (1X/ ' ***   Insufficient data points.'
     1  7X,'ND =',I3)
99091 FORMAT (1X/ ' ***   No desired points.'
     1  7X,'NI =',I3)
99092 FORMAT (1X/ ' ***   Two data points identical or out of ',
     1  'sequence.'/
     2  7X,'ID, XD(ID-1), XD(ID) =',I5,2F10.3)
99099 FORMAT (' Error detected in the UVIP3P subroutine'/)
      END

      SUBROUTINE SCALEXY(AN, XT, YT, II0, II1, XD, YD)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      DIMENSION   XD(*),YD(*),XT(*), YT(*), AN(*)

      YN=0D0
      XN=0D0
      DO II=II0,II1
         YN = MAX(DABS(YD(II)),YN)
         IF (II .GT. II0) THEN
            XN = XN + DABS(XD(II)-XD(II-1))
         ENDIF
      END DO
      IF (YN .LE. 0) YN = 1D0
      XN = XN / (II1-II0)
      DO II=II0,II1
         K = II-II0+1
         IF (II .GT. II0) THEN
            XT(K) = XD(II0) + (XD(II)-XD(II0))/XN
         ELSE
            XT(K) = XD(II0)
         ENDIF
         YT(K) = YD(II)/YN
      ENDDO
      AN(1) = 1D0/XN
      AN(2) = YN
      RETURN
      END
      

      SUBROUTINE  UVIP3O(NP,ND,XD,YD,NI,XI, YI)
C
C Univariate Interpolation (Improved Akima Method)
C
C Hiroshi Akima
C U.S. Department of Commerce, NTIA/ITS
C Version of 89/07/04
C
C This subroutine performs univariate interpolation.  It is based
C on the improved A method developed by Hiroshi Akima, 'A method
C of univariate interpolation that has the accuracy of a third-
C degree polynomial,' ACM TOMS, vol. xx, pp. xxx-xxx, 19xx.  (The
C equation numbers referred to in the comments below are those in
C the paper.)
C
C In this method, the interpolating function is a piecewise
C function composed of a set of polynomials applicable to
C successive intervals of the given data points.  This method
C uses third-degree polynomials as the default, but the user has
C an option to use higher-degree polynomial to reduce undulations
C in resulting curves.
C
C This method has the accuracy of a third-degree polynomial if
C the degree of the polynomials for the interpolating function is
C set to three.
C
C The input arguments are
C   NP = degree of the polynomials for the interpolating
C        function,
C   ND = number of input data points
C        (must be equal to 2 or greater),
C   XD = array of dimension ND, containing the abscissas of
C        the input data points
C        (must be in a monotonic increasing order),
C   YD = array of dimension ND, containing the ordinates of
C        the input data points,
C   NI = number of points for which interpolation is desired
C        (must be equal to 1 or greater),
C   XI = array of dimension NI, containing the abscissas of
C        the desired points.
C
C The output argument is
C   YI = array of dimension NI, where the ordinates of the
C        desired points are to be stored.
C
C If an integer value smaller than 3 is given to the NP argument,
C this subroutine assumes NP = 3.
C
C The XI array elements need not be monotonic, but this
C subroutine interpolates faster if the XI array elements are
C given in a monotonic order.
C
C If the XI array element is less than XD(1) or greater than
C XD(ND), this subroutine linearly interpolates the YI value.
C
C
C Specification statement
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      DIMENSION   XD(*),YD(*),XI(*), YI(*)
C Error check
   10 IF (ND.LE.1)   GO TO 90
      IF (NI.LE.0)   GO TO 91
      DO 11  ID=2,ND
        IF (XD(ID).LE.XD(ID-1))     GO TO 92
   11 CONTINUE
C Branches off special cases.
C Take calre of the linear case
      IF (ND .GE. 3 .AND. NP .EQ. 1) GOTO 30
C Take Care of the quadratic interpolation case.
      IF (ND .GE. 4 .AND. NP .EQ. 2) GOTO 30
C Now for higher polynomials.
      IF (ND.LE.4)   GO TO 50
C General case  --  Five data points of more
C Calculates some local variables.
   20 NP0=MAX(3,NP)
      NPM1=NP0-1
      RENPM1=NPM1
      RENNM2=NP0*(NP0-2)
C Main calculation for the general case
C First (outermost) DO-loop with respect to the desired points
   30 DO 39  II=1,NI
        IF (II.EQ.1)      IINTPV=-1
        XII=XI(II)
C Locates the interval that includes the desired point by binary
C search.
        IF (XII.LE.XD(1))  THEN
          IINT=0
        ELSE IF (XII.LT.XD(ND))  THEN
          IDMN=1
          IDMX=ND
          IDMD=(IDMN+IDMX)/2
   31     IF (XII.GE.XD(IDMD))  THEN
            IDMN=IDMD
          ELSE
            IDMX=IDMD
          END IF
          IDMD=(IDMN+IDMX)/2
          IF (IDMD.GT.IDMN)    GO TO 31
          IINT=IDMD
        ELSE
          IINT=ND
        END IF
C End of locating the interval of interest
C Interpolation or extrapolation in one of the three subcases
        IF (IINT.LE.0)  THEN
C Subcase 1  --  Linear extrapolation when the abscissa of the
C                desired point is equal to that of the first data
C                point or less.
C Estimates the first derivative when the interval is not the
C same as the one for the previous desired point.  --
C cf. Equation (8)
          IF (IINT.NE.IINTPV)  THEN
            IINTPV=IINT
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
               XDLT = X1*X2*(X2-X1)
               A1 = (X2*X2*Y1 - X1*X1*Y2)/XDLT
            ELSE
               XDLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
               A1=(((X2*X3)**2)*(X3-X2)*Y1
     1              +((X3*X1)**2)*(X1-X3)*Y2
     2              +((X1*X2)**2)*(X2-X1)*Y3)/XDLT
            ENDIF
         ENDIF
C Evaluates the YI value.
         YI(II)=Y0+A1*(XII-X0)
C End of Subcase 1
      ELSE IF (IINT.GE.ND)  THEN
C Subcase 2  --  Linear extrapolation when the abscissa of the
C                desired point is equal to that of the last data
C                point or greater.
C Estimates the first derivative when the interval is not the
C same as the one for the previous desired point.  --
C cf. Equation (8)
         IF (IINT.NE.IINTPV)  THEN
            IINTPV=IINT
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
               XDLT = X1*X2*(X2-X1)
               A1 = (X2*X2*Y1 - X1*X1*Y2)/XDLT  
            ELSE
               XDLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
               A1=(((X2*X3)**2)*(X3-X2)*Y1
     1              +((X3*X1)**2)*(X1-X3)*Y2
     2              +((X1*X2)**2)*(X2-X1)*Y3)/XDLT
            ENDIF
          END IF
C Evaluates the YI value.
          YI(II)=Y0+A1*(XII-X0)
C End of Subcase 2
        ELSE
C ADD QUADRATIC INTERPOLATION WHEN NP=2, WHICH ONLY 
C WORKS IF ND >= 4
           IF (NP .EQ. 2) THEN
              IF (IINT .NE. IINTPV) THEN
                 IINTPV = IINT
                 IF (IINT .EQ. 1) THEN
                    X0 = XD(IINT)
                    Y0 = YD(IINT)
                    X1 = XD(IINT+1)-X0
                    Y1 = YD(IINT+1)-Y0
                    X2 = XD(IINT+2)-X0
                    Y2 = YD(IINT+2)-Y0
                    XDLT = X1*X2*(X2-X1)
                    A1=(X2*X2*Y1-X1*X1*Y2)/XDLT
                    A2=(X1*Y2-X2*Y1)/XDLT
                 ELSE IF (IINT .EQ. ND-1) THEN
                    X0 = XD(IINT-1)
                    Y0 = YD(IINT-1)
                    X1 = XD(IINT)-X0
                    Y1 = YD(IINT)-Y0
                    X2 = XD(IINT+1)-X0
                    Y2 = YD(IINT+1)-Y0
                    XDLT = X1*X2*(X2-X1)
                    A1=(X2*X2*Y1-X1*X1*Y2)/XDLT
                    A2=(X1*Y2-X2*Y1)/XDLT
                 ELSE
                    X0 = XD(IINT)
                    Y0 = YD(IINT)
                    X1 = XD(IINT+1)-X0
                    Y1 = YD(IINT+1)-Y0
                    X2 = XD(IINT+2)-X0
                    Y2 = YD(IINT+2)-Y0
                    XDLT = X1*X2*(X2-X1)
                    A1=(X2*X2*Y1-X1*X1*Y2)/XDLT
                    A2=(X1*Y2-X2*Y1)/XDLT
                    X00 = XD(IINT-1)
                    Y00 = YD(IINT-1)
                    X11 = XD(IINT)-X00
                    Y11 = YD(IINT)-Y00
                    X22 = XD(IINT+1)-X00
                    Y22 = YD(IINT+1)-Y00
                    XDLT = X11*X22*(X22-X11)
                    A11=(X22*X22*Y11-X11*X11*Y22)/XDLT
                    A22=(X11*Y22-X22*Y11)/XDLT
                 ENDIF
              ENDIF
              XX = XI(II)-X0
              YI(II) = Y0 + XX*(A1+XX*A2)
              IF (IINT .NE. 1 .AND. IINT .NE. ND-1) THEN
                 XX = XI(II)-X00
                 YI(II) = 0.5*(YI(II) + Y00+XX*(A11+XX*A22))
              ENDIF
              GOTO 39
           ENDIF

C Here is linear interpolation when NP=1
           IF (NP .EQ. 1) THEN
              IF (IINT .NE. IINTPV) THEN
                 IINTPV = IINT
                 X0 = XD(IINT)
                 Y0 = YD(IINT)
                 X1 = XD(IINT+1)-X0
                 Y1 = YD(IINT+1)-Y0
                 A1 = Y1/X1
              ENDIF
              XX = XI(II)-X0
              YI(II) = Y0 + A1*XX
              GOTO 39
           ENDIF

C Subcase 3  --  Interpolation when the abscissa of the desired
C                point is  between those of the first and last
C                data points.
C Calculates the coefficients of the third-degree polynomial (for
C NP.LE.3) or the factors for the higher-degree polynomials (for
C NP.GT.3), when the interval is not the same as the one for the
C previous desired point.
          IF (IINT.NE.IINTPV)  THEN
            IINTPV=IINT
C The second DO-loop with respect to the two endpoints of the
C interval
            DO 37  IEPT=1,2
C Calculates the estimate of the first derivative at an endpoint.
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
C Checks if any point number falls outside the legitimate range
C (between 1 and ND).  Skips calculation of the primary estimate
C if any does.
                IF (ID1.LT.1.OR.ID2.LT.1.OR.ID3.LT.1.OR.
     1              ID1.GT.ND.OR.ID2.GT.ND.OR.ID3.GT.ND)
     2               GO TO 36
C Calculates the primary estimate of the first derivative  --
C cf. Equation (8)
                X1=XD(ID1)-X0
                X2=XD(ID2)-X0
                X3=XD(ID3)-X0
                Y1=YD(ID1)-Y0
                Y2=YD(ID2)-Y0
                Y3=YD(ID3)-Y0
                XDLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
                PE=(((X2*X3)**2)*(X3-X2)*Y1
     1             +((X3*X1)**2)*(X1-X3)*Y2
     2             +((X1*X2)**2)*(X2-X1)*Y3)/XDLT
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
     1                +YD(ID2)**2+YD(ID3)**2)*1.0E-12
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
   37       CONTINUE

C End of the second DO-loop
            IF (NP0.LE.3)  THEN
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
            ELSE
C Calculates the factors for the higher-degree polynomials
C (when NP.GT.3).  --  cf. Equation (20)
              DX=XD(IINT+1)-XD(IINT)
              DY=YD(IINT+1)-YD(IINT)
              T0=YP0*DX-DY
              T1=YP1*DX-DY
              AA0= (T0+RENPM1*T1)/RENNM2
              AA1=-(RENPM1*T0+T1)/RENNM2
            END IF
          END IF
C End of the calculation of the coefficients of the third-degree
C polynomial (when NP.LE.3) or the factors for the higher-degree
C polynomials (when NP.GT.3), when the interval is not the same
C as the one for the previous desired point.
C Evaluates the YI value.
          IF (NP0.LE.3)  THEN
C - With a third-degree polynomial.  --  cf. Equation (3)
            XX=XII-XD(IINT)
            YI(II)=A0+XX*(A1+XX*(A2+XX*A3))
          ELSE
C - With a higher-degree polynomial.  --  cf. Equation (19)
            U=(XII-XD(IINT))/DX
            UC=1.0-U
            V=AA0*((U**NP0)-U)+AA1*((UC**NP0)-UC)
            YI(II)=YD(IINT)+DY*U+V
          END IF
C End of Subcase 3
        END IF
   39 CONTINUE
C End of the first DO-loop
C End of general case
      RETURN
C Special cases  --  Four data points or less
C Preliminary processing for the special cases
   50 X0=XD(1)
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
C Special Case 1  --  Two data points
C (Linear interpolation and extrapolation)
   60 A1=Y1/X1
      DO 61  II=1,NI
        YI(II)=Y0+A1*(XI(II)-X0)
   61 CONTINUE
C End of Special Case 1
      RETURN
C Special Case 2  --  Three data points
C (Quadratic interpolation and linear extrapolation)
 70   XDLT=X1*X2*(X2-X1)
      A1=(X2*X2*Y1-X1*X1*Y2)/XDLT
      A2=(X1*Y2-X2*Y1)/XDLT
      A12=2.0*A2*X2+A1
      DO 71  II=1,NI
         XX=XI(II)-X0
         IF (XX.LE.0.0)  THEN
            YI(II)=Y0+A1*XX
         ELSE IF (XX.LT.X2) THEN
            YI(II)=Y0+XX*(A1+XX*A2)
         ELSE
            YI(II)=Y0+Y2+A12*(XX-X2)
         END IF
 71   CONTINUE
C End of Special Case 2
      RETURN
C Special Case 3  --  Four data points
C (Cubic interpolation and linear extrapolation)
   80 XDLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
      A1=(((X2*X3)**2)*(X3-X2)*Y1
     1   +((X3*X1)**2)*(X1-X3)*Y2
     2   +((X1*X2)**2)*(X2-X1)*Y3)/XDLT
      A2=(X2*X3*(X2*X2-X3*X3)*Y1
     1   +X3*X1*(X3*X3-X1*X1)*Y2
     2   +X1*X2*(X1*X1-X2*X2)*Y3)/XDLT
      A3=(X2*X3*(X3-X2)*Y1
     1   +X3*X1*(X1-X3)*Y2
     2   +X1*X2*(X2-X1)*Y3)/XDLT
      A13=(3.0*A3*X3+2.0*A2)*X3+A1
      DO 81  II=1,NI
        XX=XI(II)-X0
        IF (XX.LE.0.0)  THEN
          YI(II)=Y0+A1*XX
        ELSE IF (XX.LT.X3) THEN
          YI(II)=Y0+XX*(A1+XX*(A2+XX*A3))
        ELSE
          YI(II)=Y0+Y3+A13*(XX-X3)
        END IF
   81 CONTINUE
C End of Special Case 3
      RETURN
C Error exit
   90 WRITE (*,99090) ND
      GO TO 99
   91 WRITE (*,99091) NI
      GO TO 99
   92 WRITE (*,99092) ID,XD(ID-1),XD(ID)
   99 WRITE (*,99099)
      STOP
      RETURN
C Format statements for error messages
99090 FORMAT (1X/ ' ***   Insufficient data points.'
     1  7X,'ND =',I3)
99091 FORMAT (1X/ ' ***   No desired points.'
     1  7X,'NI =',I3)
99092 FORMAT (1X/ ' ***   Two data points identical or out of ',
     1  'sequence.'/
     2  7X,'ID, XD(ID-1), XD(ID) =',I5,2F10.3)
99099 FORMAT (' Error detected in the UVIP3P subroutine'/)
      END

