      FUNCTION F11(X,ETA,ZL,P,EPS,LIMIT,KIND,ERR,NITS,FPMAX,ACC8,ACC16)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 X,ETA,ZL,P,AA,BB,Z,F11,CDIGAM,CI
      COMPLEX*16 DD,G,F,AI2,BI2,T2
      LOGICAL ZLLIN

      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      INTEGER ND
      INTEGER MPDIM
      PARAMETER (ND = 5, MPDIM = ND+4)

      REAL AR(MPDIM), BR(MPDIM), GR(MPDIM), GI(MPDIM), DR(MPDIM)
      REAL DI(MPDIM), TR(MPDIM), TI(MPDIM), UR(MPDIM), UI(MPDIM)
      REAL FI(MPDIM), FI1(MPDIM), MPONE(MPDIM)
      REAL AI(MPDIM), BI(MPDIM), DEN(MPDIM), TR2(MPDIM)
      REAL TI2(MPDIM), MPTMP(MPDIM), ZI(MPDIM), ZR(MPDIM)

      DATA ZERO,ONE,TWO / 0D+0, 1D+0, 2D+0 /, CI / (0D+0, 1D+0) /

      ABSC(AA) = ABS(DBLE(AA)) + ABS(IMAG(AA))
      NINTC(AA) = NINT(REAL(REAL(AA)))
C
C *** evaluate the HYPERGEOMETRIC FUNCTION 1F1
C                                        i
C            F (AA;BB; Z) = SUM  (AA)   Z / ( (BB)  i] )
C           1 1              i       i            i
C
C     to accuracy EPS with at most LIMIT terms.
C  If KIND = 0 : using extended precision but real arithmetic only,
C            1 : using normal precision in complex arithmetic,
C   or       2 : using normal complex arithmetic, but with CDIGAM factor
C
C  where
         AA = ZL+ONE - ETA*P
         BB = TWO*(ZL+ONE)
C  and
         Z  = TWO*P*X
C
         ZLLIN = DBLE(BB).LE.ZERO .AND. ABS(BB-NINTC(BB)).LT.ACC8**0.25
             IF(.NOT.ZLLIN.OR.DBLE(BB)+LIMIT.LT.1.5) GO TO 10
                NITS = -1
                RETURN
   10 IF(LIMIT.LE.0) THEN
         F11 = ZERO
         ERR = ZERO
         NITS= 1
         RETURN
         ENDIF
      TA = ONE
      RK = ONE
      IF(KIND.LE.0.AND.ABSC(Z)*ABSC(AA).GT.ABSC(BB) * 1.0) THEN
         NW = ND
         CALL MPDMC(ONE, 0, DR)
         CALL MPDMC(ZERO, 0, DI)
         CALL MPDMC(ONE, 0, GR)
         CALL MPDMC(ZERO, 0, GI)
         CALL MPDMC(IMAG(AA), 0, AI)
         CALL MPDMC(DBLE(AA), 0, AR)
         CALL MPDMC(IMAG(BB), 0, BI)
         CALL MPDMC(DBLE(BB), 0, BR)
         CALL MPDMC(ZERO, 0, FI)
         CALL MPDMC(ONE, 0, MPONE)
         CALL MPDMC(DBLE(Z), 0, ZR)
         CALL MPDMC(IMAG(Z), 0, ZI)
         
      DO 20 I=2,LIMIT
         CALL MPADD(FI, MPONE, FI1)

         CALL MPMUL(BR, FI1, TR)
         CALL MPMUL(BI, FI1, TI)

         CALL MPNPWR(TR, 2, TR2)
         CALL MPNPWR(TI, 2, TI2)
         CALL MPADD(TR2, TI2, MPTMP)
         CALL MPDIV(MPONE, MPTMP, DEN)

         CALL MPMUL(AR, TR, TR2)
         CALL MPMUL(AI, TI, TI2)
         CALL MPADD(TR2, TI2, MPTMP)
         CALL MPMUL(MPTMP, DEN, UR)

         CALL MPMUL(AI, TR, TR2)
         CALL MPMUL(AR, TI, TI2)
         CALL MPSUB(TR2, TI2, MPTMP)
         CALL MPMUL(MPTMP, DEN, UI)

         CALL MPMUL(UR, GR, TR2)
         CALL MPMUL(UI, GI, TI2)
         CALL MPSUB(TR2, TI2, TR)

         CALL MPMUL(UR, GI, TR2)
         CALL MPMUL(UI, GR, TI2)
         CALL MPADD(TR2, TI2, TI)

         CALL MPMUL(ZR, TR, TR2)
         CALL MPMUL(ZI, TI, TI2)
         CALL MPSUB(TR2, TI2, GR)

         CALL MPMUL(ZR, TI, TR2)
         CALL MPMUL(ZI, TR, TI2)
         CALL MPADD(TR2, TI2, GI)

         CALL MPADD(DR, GR, MPTMP)
         CALL MPEQ(MPTMP, DR)
         CALL MPADD(DI, GI, MPTMP)
         CALL MPEQ(MPTMP, DI)
         
         CALL MPEQ(GR, TR2)
         TR2(1) = ABS(TR2(1))
         CALL MPEQ(GI, TI2)
         TI2(1) = ABS(TI2(1))
         CALL MPADD(TR2, TI2, MPTMP)
         CALL MPMDC(MPTMP, ERR, N)
         ERR = ERR * TWO**N
         IF(ERR.GT.FPMAX) GO TO 60

         CALL MPEQ(DR, TR2)
         TR2(1) = ABS(TR2(1))
         CALL MPEQ(DI, TI2)
         TI2(1) = ABS(TI2(1))
         CALL MPADD(TR2, TI2, MPTMP)
         CALL MPMDC(MPTMP, RK, N)
         RK = RK * TWO**N
         TA = MAX(TA,RK)

         IF(ERR.LT.RK*EPS .OR. I.GE.4.AND.ERR.LT.ACC16) GO TO 30

         CALL MPMDC(FI1, DFI1, N)
         CALL MPEQ(FI1, FI)
         CALL MPADD(AR, MPONE, MPTMP)
         CALL MPEQ(MPTMP, AR)
         CALL MPADD(BR, MPONE, MPTMP)
         CALL MPEQ(MPTMP, BR)
 20      CONTINUE
C
 30      CALL MPMDC(DR, F11R, N)
         F11R = F11R * TWO**N
         CALL MPMDC(DI, F11I, N)
         F11I = F11I * TWO**N
         F11 = DCMPLX(F11R, F11I)
         ERR = ACC16 * TA / RK
C
      ELSE
C* ---------------------------------- alternative code
C*    If REAL*16 arithmetic is not available, (or already using it]),
C*    then use KIND > 0
         G = ONE
          F = ONE
          IF(KIND.GE.2) F = CDIGAM(AA) - CDIGAM(BB) - CDIGAM(G)
         DD = F
         DO 40 I=2,LIMIT
            AI2 = AA + (I-2)
            BI2 = BB + (I-2)
            R  = I-ONE
         G = G * Z * AI2 / (BI2 * R)
         IF(KIND.GE.2)
C                              multiply by (psi(a+r)-psi(b+r)-psi(1+r))
     X        F = F + ONE/AI2 - ONE/BI2 - ONE/R
         T2  = G * F
         DD = DD + T2
            ERR = ABSC(T2)
               IF(ERR.GT.FPMAX) GO TO 60
            RK = ABSC(DD)
         TA = MAX(TA,RK)
         IF(ERR.LT.RK*EPS.OR.ERR.LT.ACC8.AND.I.GE.4) GO TO 50
   40    CONTINUE

   50    ERR = ACC8 * TA / RK
         F11 = DD
C* ------------------------------------------- end of alternative code
      ENDIF
   60    NITS = I
      RETURN
      END
