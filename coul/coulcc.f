      SUBROUTINE COULCC(XX,ETA1,ZLMIN,NL, FC,GC,FCP,GCP, SIG,
     X                  MODE1,KFN,IFAIL)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  COMPLEX COULOMB WAVEFUNCTION PROGRAM USING STEED'S METHOD           C
C                                                                      C
C  A. R. Barnett           Manchester  March   1981                    C
C  modified I.J. Thompson  Daresbury, Sept. 1983 for Complex Functions C
C                                                                      C
C  original program  RCWFN       in    CPC  8 (1974) 377-395           C
C                 +  RCWFF       in    CPC 11 (1976) 141-142           C
C                 +  COULFG      in    CPC 27 (1982) 147-166           C
C  description of real algorithm in    CPC 21 (1981) 297-314           C
C  description of complex algorithm    JCP XX (1985) YYY-ZZZ           C
C  this version written up       in    CPC XX (1985) YYY-ZZZ           C
C                                                                      C
C  COULCC returns F,G,G',G',SIG for complex XX, ETA1, and ZLMIN,       C
C   for NL integer-spaced lambda values ZLMIN to ZLMIN+NL-1 inclusive, C
C   thus giving  complex-energy solutions to the Coulomb Schrodinger   C
C   equation,to the Klein-Gordon equation and to suitable forms of     C
C   the Dirac equation ,also spherical & cylindrical Bessel equations  C
C                                                                      C
C  if /MODE1/= 1  get F,G,F',G'   for integer-spaced lambda values     C
C            = 2      F,G      unused arrays must be dimensioned in    C
C            = 3      F,  F'          call to at least length (1)      C
C            = 4      F                                                C
C            = 11 get F,H+,F',H+' ) if KFN=0, H+ = G + i.F        )    C
C            = 12     F,H+        )       >0, H+ = J + i.Y = H(1) ) in C
C            = 21 get F,H-,F',H-' ) if KFN=0, H- = G - i.F        ) GC C
C            = 22     F,H-        )       >0, H- = J - i.Y = H(2) )    C
C                                                                      C
C     if MODE1<0 then the values returned are scaled by an exponential C
C                factor (dependent only on XX) to bring nearer unity   C
C                the functions for large /XX/, small ETA & /ZL/ < /XX/ C
C        Define SCALE = (  0        if MODE1 > 0                       C
C                       (  IMAG(XX) if MODE1 < 0  &  KFN < 3           C
C                       (  DBLE(XX) if MODE1 < 0  &  KFN = 3           C
C        then FC = EXP(-ABS(SCALE)) * ( F, j, J, or I)                 C
C         and GC = EXP(-ABS(SCALE)) * ( G, y, or Y )                   C
C               or EXP(SCALE)       * ( H+, H(1), or K)                C
C               or EXP(-SCALE)      * ( H- or H(2) )                   C
C                                                                      C
C  if  KFN  =  0,-1  complex Coulomb functions are returned   F & G    C
C           =  1   spherical Bessel      "      "     "       j & y    C
C           =  2 cylindrical Bessel      "      "     "       J & Y    C
C           =  3 modified cyl. Bessel    "      "     "       I & K    C
C                                                                      C
C          and where Coulomb phase shifts put in SIG if KFN=0 (not -1) C
C                                                                      C
C  The use of MODE and KFN is independent                              C
C    (except that for KFN=3,  H(1) & H(2) are not given)               C
C                                                                      C
C  With negative orders lambda, COULCC can still be used but with      C
C    reduced accuracy as CF1 becomes unstable. The user is thus        C
C    strongly advised to use reflection formulae based on              C
C    H+-(ZL,,) = H+-(-ZL-1,,) * exp +-i(sig(ZL)-sig(-ZL-1)-(ZL+1/2)pi) C
C                                                                      C
C  Precision:  results to within 2-3 decimals of 'machine accuracy',   C
C               but if CF1A fails because X too small or ETA too large C
C               the F solution  is less accurate if it decreases with  C
C               decreasing lambda (e.g. for lambda.LE.-1 & ETA.NE.0)   C
C              RERR in COMMON/STEED/ traces the main roundoff errors.  C
C                                                                      C
C   COULCC is coded for real*8 on IBM or equivalent  ACCUR >= 10**-14  C
C          with a section of doubled REAL*16 for less roundoff errors. C
C          (If no doubled precision available, increase JMAX to eg 100)C
C   Use IMPLICIT COMPLEX*32 & REAL*16 on VS compiler ACCUR >= 10**-32  C
C   For single precision CDC (48 bits) reassign REAL*8=REAL etc.       C
C                                                                      C
C   IFAIL  on input   = 0 : no printing of error messages              C
C                    ne 0 : print error messages on file 6             C
C   IFAIL  in output = -2 : argument out of range                      C
C                    = -1 : one of the continued fractions failed,     C
C                           or arithmetic check before final recursion C
C                    =  0 : All Calculations satisfactory              C
C                    ge 0 : results available for orders up to & at    C
C                             position NL-IFAIL in the output arrays.  C
C                    = -3 : values at ZLMIN not found as over/underflowC
C                    = -4 : roundoff errors make results meaningless   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     Machine dependent constants :                                    C
C                                                                      C
C     ACCUR    target bound on relative error (except near 0 crossings)C
C               (ACCUR should be at least 100 * ACC8)                  C
C     ACC8     smallest number with 1+ACC8 .ne.1 in REAL*8  arithmetic C
C     ACC16    smallest number with 1+ACC16.ne.1 in REAL*16 arithmetic C
C     FPMAX    magnitude of largest floating point number * ACC8       C
C     FPMIN    magnitude of smallest floating point number / ACC8      C
C     FPLMAX   LOG(FPMAX)                                              C
C     FPLMIN   LOG(FPMIN)                                              C
C                                                                      C
C     ROUTINES CALLED :       LOGAM/CLOGAM/CDIGAM,                     C
C                             F20, CF1A, RCF, CF1C, CF2, F11, CF1R     C
C     Intrinsic functions :   MIN, MAX, SQRT, REAL, IMAG, ABS, LOG, EXP,
C      (Generic names)        NINT, MOD, ATAN, ATAN2, COS, SIN, DCMPLX,
C                             SIGN, CONJG, INT, TANH                   C
C     Note: Statement fntn.   NINTC = integer nearest to a complex no. C
C                                                                      C
C     Parameters determining region of calculations :                  C
C                                                                      C
C        R20      estimate of (2F0 iterations)/(CF2 iterations)        C
C        ASYM     minimum X/(ETA**2+L) for CF1A to converge easily     C
C        XNEAR    minimum ABS(X) for CF2 to converge accurately        C
C        LIMIT    maximum no. iterations for CF1, CF2, and 1F1 series  C
C        JMAX     size of work arrays for Pade accelerations           C
C        NDROP    number of successive decrements to define instabilityC
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT COMPLEX*16 (A-H,O-Z)
      PARAMETER(JMAX=50)
      DIMENSION FC(NL),GC(NL),FCP(NL),GCP(NL),SIG(NL),XRCF(JMAX,4)
      LOGICAL PR,ETANE0,IFCP,RLEL,DONEM,UNSTAB,ZLNEG,AXIAL,NOCF2,NPINT
      REAL*8 ERR,RERR,ABSC,ACCUR,ACCT,ACC8,ACCH,ACC16,ACCB, XNEAR,CF1R,
     X       ZERO,ONE,TWO,HALF,HPI,TLOG,FPMAX,FPMIN,FPLMIN,FPLMAX,
     X       PACCQ,EPS,OFF,SCALE,SF,SFSH,TA,RK,OMEGA,R20,ASYM,ABSX
      COMPLEX*16 LOGAM, ZTMP

C
      COMMON       /CSTEED/ RERR,NFP,N11,NPQ(2),N20,KAS(2)
C***  common blocks are for information & storage only.
C     (they are not essential to working of the code)
C      COMMON /RCFCM1/ PK,EK,CLGAA,CLGAB,CLGBB,DSIG,TPK1,W,RL,FCL1,Q,GAM,
C     X                HCL,HPL,FCM,HCL1,ALPHA,BETA,PL
C      EQUIVALENCE            (PK,XRCF(1,1))
C
      DATA ZERO,ONE,TWO,LIMIT /0.0D+0, 1.0D+0, 2.0D+0, 20000 /,
     X     HALF, CI / 0.5D+0, (0D+0, 1D+0) /,
     X     FPMAX,FPMIN,FPLMAX,FPLMIN / 1D+60,1D-60 ,140D+0, -140D+0 /,
     X     R20,ASYM,XNEAR,NDROP / 3., 3., .5, 5 /,
     X     ACCUR, ACC8, ACC16 / 1D-14, 2D-16, 3D-33 /      
!$OMP THREADPRIVATE(/CSTEED/)
      NINTC(W) = NINT(REAL(REAL(W)))
      ABSC(W) = ABS(DBLE(W)) + ABS(IMAG(W))
      NPINT(W,ACCB) = ABSC(NINTC(W)-W).LT.ACCB .AND. DBLE(W).LT.HALF
C
      MODE = MOD(ABS(MODE1),10)
      IFCP = MOD(MODE,2).EQ.1
      PR = IFAIL.NE.0
      IFAIL = -2
      N11   = 0
      NFP   = 0
      KAS(1)   = 0
      KAS(2)   = 0
      NPQ(1)   = 0
      NPQ(2)   = 0
      N20 = 0
      HPI = TWO*ATAN(ONE)
      TLOG = LOG(TWO)
      ACCUR = MAX(ACCUR, 50*ACC8)
      ACCT = ACCUR * .5
C                       initialise the log-gamma function :
      ZTMP = LOGAM(ACC8)
      ACCH  = SQRT(ACCUR)
      ACCB  = SQRT(ACCH)
      RERR = ACCT
C
      CIK = ONE
         IF(KFN.GE.3) CIK = CI * SIGN(ONE,ACC8-IMAG(XX))
      X     = XX * CIK
      ETA   = ETA1
      IF(KFN .GT. 0) ETA = ZERO
         ETANE0  = ABSC(ETA).GT.ACC8
         ETAI = ETA*CI
      DELL  = ZERO
      IF(KFN .GE. 2)  DELL = HALF
      ZM1   = ZLMIN - DELL
      SCALE = ZERO
      IF(MODE1.LT.0) SCALE = IMAG(X)
C
      M1 = 1
      L1  = M1 + NL - 1
      RLEL = ABS(IMAG(ETA)) + ABS(IMAG(ZM1)) .LT. ACC8
      ABSX = ABS(X)
      AXIAL = RLEL .AND. ABS(IMAG(X)) .LT. ACC8 * ABSX
      IF(MODE.LE.2 .AND. ABSX.LT.FPMIN) GO TO 310
      XI  = ONE/X
      XLOG = LOG(X)
C            log with cut along the negative real axis] see also OMEGA
      ID = 1
      DONEM = .FALSE.
         UNSTAB = .FALSE.
      LF = M1
      IFAIL = -1
   10    ZLM = ZM1 + LF - M1
         ZLL = ZM1 + L1 - M1
C
C ***       ZLL  is final lambda value, or 0.5 smaller for J,Y Bessels
C
              Z11 = ZLL
              IF(ID.LT.0) Z11 = ZLM
              P11 = CI*SIGN(ONE,ACC8-IMAG(ETA))
      LAST = L1
C
C ***       Find phase shifts and Gamow factor at lambda = ZLL
C
      PK = ZLL + ONE
      AA = PK - ETAI
      AB = PK + ETAI
      BB = TWO*PK
         ZLNEG = NPINT(BB,ACCB)
                     CLGAA = CLOGAM(AA)
                     CLGAB = CLGAA
         IF(ETANE0.AND..NOT.RLEL)  CLGAB = CLOGAM(AB)
         IF(ETANE0.AND.     RLEL)  CLGAB = CONJG(CLGAA)
         SIGMA = (CLGAA - CLGAB) * CI*HALF
         IF(KFN.EQ.0) SIG(L1) = SIGMA
         IF(.NOT.ZLNEG) CLL = ZLL*TLOG- HPI*ETA - CLOGAM(BB)
     X                                          + (CLGAA+CLGAB)*HALF
              THETA  = X - ETA*(XLOG+TLOG) - ZLL*HPI + SIGMA
C
        TA = (IMAG(AA)**2+IMAG(AB)**2+ABS(DBLE(AA))+ABS(DBLE(AB)))*HALF
      IF(ID.GT.0 .AND. ABSX .LT. TA*ASYM .AND. .NOT.ZLNEG) GO TO 20
C
C ***         use CF1 instead of CF1A, if predicted to converge faster,
C                 (otherwise using CF1A as it treats negative lambda &
C                  recurrence-unstable cases properly)
C
           RK = SIGN(ONE, DBLE(X) + ACC8)
           P =  THETA
           IF(RK.LT.0) P = -X + ETA*(LOG(-X)+TLOG)-ZLL*HPI-SIGMA
      F = RK * CF1A(X*RK,ETA*RK,ZLL,P,ACCT,JMAX,NFP,FEST,ERR,FPMAX,XRCF,
     X                                      XRCF(1,3), XRCF(1,4))
      FESL = LOG(FEST) + ABS(IMAG(X))
         NFP = - NFP
      IF(NFP.LT.0   .OR.(UNSTAB.AND.ERR.LT.ACCB)) GO TO 40
      IF(.NOT.ZLNEG .OR. UNSTAB.AND.ERR.GT.ACCB)  GO TO 20
c         IF(PR) WRITE(6,1060) '-L',ERR
         IF(ERR.GT.ACCB) GO TO 280
         GO TO 40
C
C ***    evaluate CF1  =  f   =  F'(ZLL,ETA,X)/F(ZLL,ETA,X)
C
   20 IF(AXIAL) THEN
C                                                        REAL VERSION
      F = CF1R(X,ETA,ZLL,ACC8,SF ,RK,  ETANE0,LIMIT,ERR,NFP,
     X         ACCH,FPMIN,FPMAX,PR,'COULCC')
          FCL = SF
          TPK1= RK
         ELSE
C                                                        COMPLEX VERSION
      F = CF1C(X,ETA,ZLL,ACC8,FCL,TPK1,ETANE0,LIMIT,ERR,NFP,
     X         ACCH,FPMIN,FPMAX,PR,'COULCC')
         ENDIF
      IF(ERR.GT.ONE) GO TO 390
C
C ***  Make a simple check for CF1 being badly unstable:
C
      IF(ID.LT.0) GO TO 30
      UNSTAB = DBLE((ONE-ETA*XI)*CI*IMAG(THETA)/F).GT.ZERO
     X .AND..NOT.AXIAL .AND. ABS(IMAG(THETA)).GT.-LOG(ACC8)*.5
     X .AND. ABSC(ETA)+ABSC(ZLL).LT.ABSC(X)
      IF(UNSTAB) GO TO 60
C
C *** compare accumulated phase FCL with asymptotic phase for G(k+1) :
C     to determine estimate of F(ZLL) (with correct sign) to start recur
C
   30 W   =  X*X  *(HALF/TPK1 + ONE/TPK1**2) + ETA*(ETA-TWO*X)/TPK1
      FESL   = (ZLL+ONE) * XLOG + CLL - W - LOG(FCL)
   40 FESL = FESL - ABS(SCALE)
          RK   =        MAX(DBLE(FESL), FPLMIN*HALF)
          FESL = DCMPLX(MIN(RK,   FPLMAX*HALF ) , IMAG(FESL))
      FEST= EXP(FESL)
C
           RERR = MAX(RERR, ERR, ACC8 * ABS(DBLE(THETA)) )
C
      FCL = FEST
      FPL = FCL*F
      IF(IFCP) FCP(L1) = FPL
               FC (L1) = FCL
C
C *** downward recurrence to lambda = ZLM. array GC,if present,stores RL
C
      I  = MAX(-ID, 0)
      ZL  = ZLL + I
         MONO = 0
        OFF = ABS(FCL)
         TA = ABSC(SIGMA)
      DO 70  L  = L1-ID,LF,-ID
         IF(ETANE0) THEN
               IF(RLEL) THEN
                    DSIG = ATAN2(DBLE(ETA),DBLE(ZL))
                    RL = SQRT(DBLE(ZL)**2 + DBLE(ETA)**2)
                  ELSE
                    AA = ZL - ETAI
                    BB = ZL + ETAI
                    IF(ABSC(AA).LT.ACCH.OR.ABSC(BB).LT.ACCH) GOTO 50
                    DSIG = (LOG(AA) - LOG(BB)) * CI*HALF
                    RL = AA * EXP(CI*DSIG)
                 ENDIF
             IF(ABSC(SIGMA).LT.TA*HALF) THEN
C               re-calculate SIGMA because of accumulating roundoffs:
                SL =(CLOGAM(ZL+I-ETAI)-CLOGAM(ZL+I+ETAI))*CI*HALF
                RL = (ZL - ETAI) * EXP(CI*ID*(SIGMA - SL))
                SIGMA = SL
                TA = ZERO
              ELSE
                SIGMA = SIGMA - DSIG * ID
              ENDIF
                TA = MAX(TA, ABSC(SIGMA))
             SL    =  ETA  + ZL*ZL*XI
                PL = ZERO
                IF(ABSC(ZL).GT.ACCH) PL = (SL*SL - RL*RL)/ZL
             FCL1  = (FCL *SL + ID*ZL*FPL)/RL
              SF = ABS(FCL1)
                       IF(SF.GT.FPMAX) GO TO 350
             FPL   = (FPL *SL + ID*PL*FCL)/RL
             IF(MODE .LE. 1) GCP(L+ID)= PL * ID
        ELSE
C                               ETA = 0, including Bessels.  NB RL==SL
           RL = ZL* XI
           FCL1 = FCL * RL + FPL*ID
              SF = ABS(FCL1)
                      IF(SF.GT.FPMAX) GO TO 350
           FPL  =(FCL1* RL - FCL) * ID
        ENDIF
C             IF(ABSC(FCL1).LT.ABSC(FCL)) THEN
              IF(SF.LT.OFF) THEN
                 MONO = MONO + 1
                ELSE
                 MONO = 0
                ENDIF
         FCL   =  FCL1
           OFF = SF
         FC(L) =  FCL
         IF(IFCP) FCP(L)  = FPL
           IF(KFN.EQ.0) SIG(L) = SIGMA
           IF(MODE .LE. 2) GC(L+ID) = RL
      ZL = ZL - ID
      IF(MONO.LT.NDROP) GO TO 70
      IF(AXIAL .OR. DBLE(ZLM)*ID.GT.-NDROP.AND..NOT.ETANE0) GO TO 70
         UNSTAB = .TRUE.
C
C ***    take action if cannot or should not recur below this ZL:
   50    ZLM = ZL
         LF = L
            IF(ID.LT.0) GO TO 380
         IF(.NOT.UNSTAB) LF = L + 1
         IF(L+MONO.LT.L1-2 .OR. ID.LT.0 .OR. .NOT.UNSTAB) GO TO 80
C             otherwise, all L values (for stability) should be done
C                        in the reverse direction:
         GO TO 60
   70 CONTINUE
      GO TO 80
   60       ID = -1
            LF = L1
            L1 = M1
            RERR = ACCT
            GO TO 10
   80 IF(FCL .EQ. ZERO) FCL = + ACC8
      F  = FPL/FCL
C
C *** Check, if second time around, that the 'f' values agree]
C
      IF(ID.GT.0) FIRST = F
      IF(DONEM) RERR = MAX(RERR, ABSC(F-FIRST)/ABSC(F))
      IF(DONEM) GO TO 90
C
       NOCF2 = .FALSE.
      THETAM  = X - ETA*(XLOG+TLOG) - ZLM*HPI + SIGMA
C
C *** on left x-plane, determine OMEGA by requiring cut on -x axis
C     on right x-plane, choose OMEGA (using estimate based on THETAM)
C       so H(omega) is smaller and recurs upwards accurately.
C     (x-plane boundary is shifted to give CF2(LH) a chance to converge)
C
                           OMEGA = SIGN(ONE,IMAG(X)+ACC8)
      IF(DBLE(X).GE.XNEAR) OMEGA = SIGN(ONE,IMAG(THETAM)+ACC8)
C     correction from erratum.
      IF (AXIAL) OMEGA = ONE
C
         SFSH = EXP(OMEGA*SCALE - ABS(SCALE))
         OFF=EXP(MIN(TWO * MAX(ABS(IMAG(X)),ABS(IMAG(THETAM)),
     X                         ABS(IMAG(ZLM))*3 ) , FPLMAX) )
          EPS = MAX(ACC8 , ACCT * HALF / OFF)
C
C ***    Try first estimated omega, then its opposite,
C        to find the H(omega) linearly independent of F
C        i.e. maximise  CF1-CF2 = 1/(F H(omega)) , to minimise H(omega)
C
   90 DO 100 L=1,2
         LH = 1
         IF(OMEGA.LT.ZERO) LH = 2
      PM = CI*OMEGA
      ETAP = ETA * PM
         IF(DONEM) GO TO 130
         PQ1 = ZERO
         PACCQ = ONE
         KASE = 0
C
C ***            Check for small X, i.e. whether to avoid CF2 :
C
      IF(MODE.GE.3 .AND. ABSX.LT.ONE ) GO TO 190
      IF(MODE.LT.3 .AND. (NOCF2 .OR. ABSX.LT.XNEAR .AND.
     X   ABSC(ETA)*ABSX .LT. 5 .AND. ABSC(ZLM).LT.4)) THEN
        KASE = 5
        GO TO 120
        ENDIF
C
C ***  Evaluate   CF2 : PQ1 = p + i.omega.q  at lambda = ZLM
C
         PQ1 = CF2(X,ETA,ZLM,PM,EPS,LIMIT,ERR,NPQ(LH),ACC8,ACCH,
     X             PR,ACCUR,DELL,'COULCC')
C
       ERR = ERR * MAX(ONE,ABSC(PQ1)/MAX(ABSC(F-PQ1),ACC8) )
       IF(ERR.LT.ACCH)       GO TO 110
C
C *** check if impossible to get F-PQ accurately because of cancellation
               NOCF2 = DBLE(X).LT.XNEAR .AND. ABS(IMAG(X)).LT.-LOG(ACC8)
C                original guess for OMEGA (based on THETAM) was wrong
C                Use KASE 5 or 6 if necessary if Re(X) < XNEAR
  100            OMEGA = - OMEGA
                IF(UNSTAB) GO TO 360
c                IF(DBLE(X).LT.-XNEAR .AND. PR) WRITE(6,1060) '-X',ERR
  110     RERR = MAX(RERR,ERR)
C
C ***  establish case of calculation required for irregular solution
C
  120 IF(KASE.GE.5) GO TO 130
      IF(DBLE(X) .GT. XNEAR) THEN
C          estimate errors if KASE 2 or 3 were to be used:
         PACCQ = EPS * OFF * ABSC(PQ1) / MAX(ABS(IMAG(PQ1)),ACC8)
        ENDIF
      IF(PACCQ .LT. ACCUR) THEN
          KASE = 2
          IF(AXIAL) KASE = 3
      ELSE
          KASE = 1
          IF(NPQ(1) * R20 .LT. JMAX)     KASE = 4
C             i.e. change to kase=4 if the 2F0 predicted to converge
      ENDIF
  130 GO TO (190,140,150,170,190,190),  ABS(KASE)
  140    IF(.NOT.DONEM)
C
C ***  Evaluate   CF2 : PQ2 = p - i.omega.q  at lambda = ZLM   (Kase 2)
C
     X  PQ2 = CF2(X,ETA,ZLM,-PM,EPS,LIMIT,ERR,NPQ(3-LH),ACC8,ACCH,
     X             PR,ACCUR,DELL,'COULCC')
C
        P     = (PQ2 + PQ1) * HALF
        Q     = (PQ2 - PQ1) * HALF*PM
      GO TO 160
  150   P     = DBLE(PQ1)
        Q     = IMAG(PQ1)
C
C ***   With Kase = 3 on the real axes, P and Q are real & PQ2 = PQ1*
C
        PQ2 = CONJG(PQ1)
C
C *** solve for FCM = F at lambda = ZLM,then find norm factor W=FCM/FCL
C
  160 W   = (PQ1 - F) * (PQ2 - F)
         SF = EXP(-ABS(SCALE))
      FCM = SQRT(Q / W) * SF
C                  any SQRT given here is corrected by
C                  using sign for FCM nearest to phase of FCL
      IF(DBLE(FCM/FCL).LT.ZERO) FCM  = - FCM
      GAM = (F - P)/Q
         TA = ABSC(GAM + PM)
         PACCQ= EPS * MAX(TA,ONE/TA)
      HCL = FCM * (GAM + PM) * (SFSH/(SF*SF))
C
      IF(PACCQ.GT.ACCUR .AND. KASE.GT.0) THEN
C                                    Consider a KASE = 1 Calculation
          F11V= F11(X,ETA,Z11,P11,ACCT,LIMIT,0,ERR,N11,FPMAX,ACC8,ACC16)
          IF(ERR.LT.PACCQ) GO TO 200
          ENDIF
      RERR=MAX(RERR,PACCQ)
      GO TO 230
C
C *** Arrive here if KASE = 4
C     to evaluate the exponentially decreasing H(LH) directly.
C
  170  IF(DONEM) GO TO 180
      AA = ETAP - ZLM
      BB = ETAP + ZLM + ONE
      F20V = F20(AA,BB,-HALF*PM*XI, ACCT,JMAX,ERR,FPMAX,N20,XRCF)
        IF(N20.LE.0) GO TO 190
        RERR = MAX(RERR,ERR)
         HCL = FPMIN
         IF(ABS(DBLE(PM*THETAM)+OMEGA*SCALE).GT.FPLMAX) GO TO 330
  180 HCL = F20V * EXP(PM * THETAM + OMEGA*SCALE)
      FCM = SFSH / ((F - PQ1) * HCL )
      GO TO 230
C
C *** Arrive here if KASE=1   (or if 2F0 tried mistakenly & failed)
C
C           for small values of X, calculate F(X,SL) directly from 1F1
C               using REAL*16 arithmetic if possible.
C           where Z11 = ZLL if ID>0, or = ZLM if ID<0
C
  190 F11V = F11(X,ETA,Z11,P11,ACCT,LIMIT,0,ERR,N11,FPMAX,ACC8,ACC16)
C
  200       IF(N11.LT.0) THEN
C                               F11 failed from BB = negative integer
c               WRITE(6,1060) '-L',ONE
               GO TO 390
               ENDIF
            IF(ERR.GT.PACCQ .AND. PACCQ.LT.ACCB) THEN
C                               Consider a KASE 2 or 3 calculation :
                KASE = -2
                IF(AXIAL) KASE = -3
                GO TO 130
                ENDIF
         RERR = MAX(RERR, ERR)
         IF(ERR.GT.FPMAX) GO TO 370
         IF(ID.LT.0) CLL = Z11*TLOG- HPI*ETA - CLOGAM(BB)
     X                       + CLOGAM(Z11 + ONE + P11*ETA) - P11*SIGMA
      EK   = (Z11+ONE)*XLOG - P11*X + CLL  - ABS(SCALE)
      IF(ID.GT.0) EK = EK - FESL + LOG(FCL)
         IF(DBLE(EK).GT.FPLMAX) GO TO 350
         IF(DBLE(EK).LT.FPLMIN) GO TO 340
      FCM = F11V * EXP(EK)
C
      IF(KASE.GE.5) THEN
        IF(ABSC(ZLM+ZLM-NINTC(ZLM+ZLM)).LT.ACCH) KASE = 6
C
C ***  For abs(X) < XNEAR, then CF2 may not converge accurately, so
C ***      use an expansion for irregular soln from origin :
C
         SL = ZLM
            ZLNEG = DBLE(ZLM) .LT. -ONE + ACCB
         IF(KASE.EQ.5 .OR. ZLNEG) SL = - ZLM - ONE
         PK = SL + ONE
            AA = PK - ETAP
            AB = PK + ETAP
            BB = TWO*PK
                     CLGAA = CLOGAM(AA)
                     CLGAB = CLGAA
         IF(ETANE0)  CLGAB = CLOGAM(AB)
                     CLGBB = CLOGAM(BB)
           IF(KASE.EQ.6 .AND. .NOT.ZLNEG) THEN
              IF(NPINT(AA,ACCUR)) CLGAA = CLGAB - TWO*PM*SIGMA
              IF(NPINT(AB,ACCUR)) CLGAB = CLGAA + TWO*PM*SIGMA
             ENDIF
          CLL = SL*TLOG- HPI*ETA - CLGBB + (CLGAA + CLGAB) * HALF
          DSIG = (CLGAA - CLGAB) * PM*HALF
             IF(KASE.EQ.6) P11 = - PM
          EK  = PK * XLOG - P11*X + CLL  - ABS(SCALE)
                     SF = EXP(-ABS(SCALE))
                     CHI = ZERO
       IF(.NOT.( KASE.EQ.5 .OR. ZLNEG ) ) GO TO 210
C
C *** Use  G(l)  =  (cos(CHI) * F(l) - F(-l-1)) /  sin(CHI)
C
C      where CHI = sig(l) - sig(-l-1) - (2l+1)*pi/2
C
         CHI = SIGMA - DSIG - (ZLM-SL) * HPI
         F11V=F11(X,ETA,SL,P11,ACCT,LIMIT,0,ERR,NPQ(1),FPMAX,ACC8,ACC16)
                    RERR = MAX(RERR,ERR)
            IF(KASE.EQ.6) GO TO 210
         FESL = F11V * EXP( EK )
         FCL1 = EXP(PM*CHI) * FCM
         HCL = FCL1 - FESL
               RERR=MAX(RERR,ACCT*MAX(ABSC(FCL1),ABSC(FESL))/ABSC(HCL))
         HCL = HCL / SIN(CHI) * (SFSH/(SF*SF))
       GO TO 220
C
C *** Use the logarithmic expansion for the irregular solution (KASE 6)
C        for the case that BB is integral so sin(CHI) would be zero.
C
  210    RL = BB - ONE
         N  = NINTC(RL)
         ZLOG = XLOG + TLOG - PM*HPI
         CHI = CHI + PM * THETAM + OMEGA * SCALE + AB * ZLOG
            AA  = ONE - AA
         IF(NPINT(AA,ACCUR)) THEN
            HCL = ZERO
         ELSE
               IF(ID.GT.0 .AND. .NOT.ZLNEG) F11V = FCM * EXP(-EK)
            HCL = EXP(CHI - CLGBB - CLOGAM(AA)) * (-1)**(N+1)
     X              * ( F11V * ZLOG +
     X      F11(X,ETA,SL,-PM,ACCT,LIMIT,2,ERR,NPQ(2),FPMAX,ACC8,ACC16))
                RERR = MAX(RERR,ERR)
            ENDIF
         IF(N.GT.0) THEN
             EK = CHI + CLOGAM(RL) - CLGAB - RL*ZLOG
             DF =F11(X,ETA,-SL-ONE,-PM,ZERO,N,0,ERR,L,FPMAX,ACC8,ACC16)
             HCL = HCL + EXP(EK) * DF
            ENDIF
C
  220    PQ1 = F - SFSH/(FCM * HCL)
      ELSE
           IF(MODE.LE.2) HCL = SFSH/((F - PQ1) * FCM)
           KASE = 1
      ENDIF
C
C ***  Now have absolute normalisations for Coulomb Functions
C          FCM & HCL  at lambda = ZLM
C      so determine linear transformations for Functions required :
C
  230 IH = ABS(MODE1) / 10
        IF(KFN.EQ.3) IH = (3-IMAG(CIK))/2  + HALF
      P11 = ONE
      IF(IH.EQ.1) P11 = CI
      IF(IH.EQ.2) P11 = -CI
                  DF = - PM
      IF(IH.GE.1) DF = - PM + P11
          IF(ABSC(DF).LT.ACCH) DF = ZERO
C
C *** Normalisations for spherical or cylindrical Bessel functions
C
                          ALPHA = ZERO
          IF(KFN  .EQ. 1) ALPHA = XI
          IF(KFN  .GE. 2) ALPHA = XI*HALF
                          BETA  = ONE
          IF(KFN  .EQ. 1) BETA  = XI
          IF(KFN  .GE. 2) BETA  = SQRT(XI/HPI)
          IF(KFN  .GE. 2 .AND. DBLE(BETA).LT.ZERO) BETA  = - BETA
C
      AA = ONE
      IF(KFN.GT.0) AA = -P11 * BETA
      IF(KFN.GE.3) THEN
C                        Calculate rescaling factors for I & K output
         P = EXP((ZLM+DELL) * HPI * CIK)
         AA= BETA * HPI * P
         BETA = BETA / P
         Q = CIK * ID
        ENDIF
C                        Calculate rescaling factors for GC output
      IF(IH.EQ.0) THEN
         TA = ABS(SCALE) + IMAG(PM)*SCALE
         RK = ZERO
         IF(TA.LT.FPLMAX) RK = EXP(-TA)
       ELSE
         TA = ABS(SCALE) + IMAG(P11)*SCALE
C
         IF(ABSC(DF).GT.ACCH .AND. TA.GT.FPLMAX) GO TO 320
         IF(ABSC(DF).GT.ACCH) DF = DF * EXP(TA)
         SF = TWO * (LH-IH) * SCALE
         RK = ZERO
         IF(SF.GT.FPLMAX) GO TO 320
         IF(SF.GT.FPLMIN) RK = EXP(SF)
      ENDIF
C
         KAS((3-ID)/2) = KASE
      W = FCM / FCL
         IF(LOG(ABSC(W))+LOG(ABSC(FC(LF))) .LT. FPLMIN) GO TO 340
         IF(MODE.GE.3) GO TO 240
c         IF(ABSC(F-PQ1) .LT. ACCH*ABSC(F) .AND. PR)
c     X        WRITE(6,1020) LH,ZLM+DELL
      HPL = HCL * PQ1
         IF(ABSC(HPL).LT.FPMIN.OR.ABSC(HCL).LT.FPMIN) GO TO 330
C
C *** IDward recurrence from HCL,HPL(LF) (stored GC(L) is RL if reqd)
C *** renormalise FC,FCP at each lambda
C ***    ZL   = ZLM - MIN(ID,0) here
C
  240 DO 270 L = LF,L1,ID
                     FCL = W* FC(L)
                      IF(ABSC(FCL).LT.FPMIN) GO TO 340
            IF(IFCP) FPL = W*FCP(L)
                     FC(L)  = BETA * FCL
            IF(IFCP) FCP(L) = BETA * (FPL - ALPHA * FCL) * CIK
                     FC(L)  = TIDY(FC(L),ACCUR)
            IF(IFCP) FCP(L) = TIDY(FCP(L),ACCUR)
       IF(MODE .GE. 3) GO TO 260
       IF(L.EQ.LF)  GO TO 250
                      ZL = ZL + ID
                      ZID= ZL * ID
                      RL = GC(L)
         IF(ETANE0)   THEN
                      SL = ETA + ZL*ZL*XI
            IF(MODE.EQ.1) THEN
              PL = GCP(L)
            ELSE
              PL = ZERO
              IF(ABSC(ZL).GT.ACCH) PL = (SL*SL - RL*RL)/ZID
            ENDIF
           HCL1     = (SL*HCL - ZID*HPL) / RL
           HPL      = (SL*HPL - PL *HCL) / RL
         ELSE
           HCL1 = RL * HCL - HPL * ID
           HPL  = (HCL - RL * HCL1) * ID
         ENDIF
         HCL      = HCL1
         IF(ABSC(HCL).GT.FPMAX) GO TO 320
  250    GC(L) = AA * (RK * HCL + DF * FCL)
      IF(MODE.EQ.1) GCP(L) = (AA *(RK*HPL +DF*FPL) - ALPHA * GC(L)) *CIK
         GC(L) = TIDY(GC(L),ACCUR)
      IF(MODE.EQ.1) GCP(L) = TIDY(GCP(L),ACCUR)
         IF(KFN.GE.3) AA = AA * Q
  260    IF(KFN.GE.3) BETA = - BETA * Q
  270  LAST = MIN(LAST,(L1 - L)*ID)
C
C *** Come here after all soft errors to determine how many L values ok
C
  280  IF(ID.GT.0 .OR.  LAST.EQ.0) IFAIL = LAST
       IF(ID.LT.0 .AND. LAST.NE.0) IFAIL = -3
C
C *** Come here after ALL errors for this L range (ZLM,ZLL)
C
  290 IF(ID.GT.0 .AND. LF.NE.M1) GO TO 300
         IF(IFAIL.LT.0) RETURN
c         IF(RERR.GT.ACCB) WRITE(6,1070) RERR
         IF(RERR.GT.0.1) IFAIL = -4
         RETURN
C
C *** so on first block, 'F' started decreasing monotonically,
C                        or hit bound states for low ZL.
C     thus redo M1 to LF-1 in reverse direction
C      i.e. do CF1A at ZLMIN & CF2 at ZLM (midway between ZLMIN & ZLMAX)
C
  300 ID = -1
      IF(.NOT.UNSTAB) LF = LF - 1
      DONEM = UNSTAB
      LF = MIN(LF,L1)
      L1 = M1
      GO TO 10
C
C ***    error messages
C
 310  return
 320  goto 280
 330  goto 280
 340  goto 280
 350  goto 280
 360  goto 280
 370  goto 390
 380  goto 390
 390  IFAIL = -1
      goto 290

c$$$
c$$$  310 IF(PR) WRITE (6,1000) XX
c$$$ 1000 FORMAT(/' COULCC: CANNOT CALCULATE IRREGULAR SOLUTIONS FOR X =',
c$$$     X 1P,2D10.2,', AS ABS(X) IS TOO SMALL'/)
c$$$      RETURN
c$$$  320 IF(PR) WRITE(6,1010) ZL+DELL,'IR',HCL,'MORE',FPMAX
c$$$ 1010 FORMAT(' COULCC: AT ZL =',2F8.3,' ',A2,'REGULAR SOLUTION (',1P,
c$$$     X 2E10.1,') WILL BE ',A4,' THAN',E10.1)
c$$$      GO TO 280
c$$$  330 IF(PR) WRITE(6,1010) ZL+DELL,'IR',HCL,'LESS',FPMIN
c$$$      GO TO 280
c$$$  340 IF(PR) WRITE(6,1010) ZL+DELL,'  ',FCL,'LESS',FPMIN
c$$$      GO TO 280
c$$$  350 IF(PR) WRITE(6,1010) ZL+DELL,'  ',FCL,'MORE',FPMAX
c$$$      GO TO 280
c$$$ 1020 FORMAT('0COULCC WARNING: LINEAR INDEPENDENCE BETWEEN ''F'' AND ''H
c$$$     X(',I1,')'' IS LOST AT ZL =',2F7.2,' (EG. COULOMB EIGENSTATE, OR CF
c$$$     X1 UNSTABLE)'/)
c$$$  360 IF(PR) WRITE(6,1030) ZLL+DELL
c$$$ 1030 FORMAT(' COULCC: (ETA&L)/X TOO LARGE FOR CF1A, AND CF1 UNSTABLE AT
c$$$     X L =',2F8.2)
c$$$      GO TO 280
c$$$  370 IF(PR) WRITE(6,1040) Z11,I
c$$$ 1040 FORMAT(' COULCC: OVERFLOW IN 1F1 SERIES AT ZL =',2F8.3,' AT TERM',
c$$$     X I5)
c$$$      GO TO 390
c$$$  380 IF(PR) WRITE(6,1050) ZLMIN,ZLM,ZLM+ONE,ZLMIN+NL-ONE
c$$$ 1050 FORMAT(' COULCC: BOTH BOUND-STATE POLES AND F-INSTABILITIES OCCUR'
c$$$     X ,', OR MULTIPLE INSTABILITIES PRESENT.'
c$$$     X,/,' TRY CALLING TWICE,  FIRST FOR ZL FROM',2F8.3,' TO',2F8.3,
c$$$     X ' (INCL.)',/,20X,     'SECOND FOR ZL FROM',2F8.3,' TO',2F8.3)
c$$$      GO TO 390
c$$$  390 IFAIL = -1
c$$$      GO TO 290
c$$$ 1060 FORMAT('0COULCC WARNING: AS ''',A2,''' REFLECTION RULES NOT USED,
c$$$     #ERRORS CAN BE UP TO',1P,D12.2/)
c$$$ 1070 FORMAT('0COULCC WARNING: OVERALL ROUNDOFF ERROR APPROX.',1P,E11.1)
      END
      FUNCTION CF1C(X,ETA,ZL,EPS,FCL,TPK1,ETANE0,LIMIT,ERR,NFP,
     X              ACCH,FPMIN,FPMAX,PR,CALLER)
      IMPLICIT COMPLEX*16(A-H,O-Z)
      LOGICAL PR,ETANE0
      REAL*8 ONE,TWO,EPS,ERR,ACCH,FPMIN,FPMAX,ABSC,SMALL,RK,PX
      CHARACTER*6 CALLER
      DATA ONE,TWO / 1D+0, 2D+0 /
      ABSC(W) = ABS(DBLE(W)) + ABS(IMAG(W))
C
C
C ***    Evaluate CF1  =  F   =  F'(ZL,ETA,X)/F(ZL,ETA,X)
C
C        using complex arithmetic
C
      FCL = ONE
      XI = ONE/X
      PK  = ZL + ONE
      PX  = PK  + LIMIT
   10 EK  = ETA / PK
        RK2 =          ONE + EK*EK
      F   = (EK + PK*XI)*FCL + (FCL - ONE)*XI
      PK1 =  PK + ONE
         TPK1 = PK + PK1
      TK  = TPK1*(XI + EK/PK1)
      IF(ETANE0) THEN
C ***   test ensures b1 .ne. zero for negative ETA etc.; fixup is exact.
             IF(ABSC(TK) .GT. ACCH)  GO TO 20
             FCL  = RK2/(ONE + (ETA/PK1)**2)
             SL   = TPK1*XI * (TPK1+TWO)*XI
             PK   =  TWO + PK
             GO TO 10
         ENDIF
   20 D   =  ONE/TK
      DF  = -FCL*RK2*D
            IF(DBLE(PK).GT.DBLE(ZL)+TWO) FCL = - RK2 * SL
            FCL = FCL * D * TPK1 * XI
      F   =  F  + DF
C
C ***   begin CF1 loop on PK = k = lambda + 1
C
      RK    = ONE
      SMALL    = SQRT(FPMIN)
   30 PK    = PK1
        PK1 = PK1 + ONE
         TPK1 = PK + PK1
         IF(ETANE0) THEN
           EK  = ETA / PK
           RK2 =          ONE + EK*EK
          ENDIF
        TK  = TPK1*(XI + EK/PK1)
        D   =  TK - D*RK2
              IF(ABSC(D) .GT. ACCH)             GO TO 40
C              IF(PR) WRITE (6,1000) CALLER,D,DF,ACCH,PK,EK,ETA,X
              RK= RK +   ONE
              IF( RK .GT. TWO )                  GO TO 50
   40 D     = ONE/D
            FCL = FCL * D * TPK1*XI
            IF(ABSC(FCL).LT.SMALL) FCL = FCL / SMALL
            IF(ABSC(FCL).GT.FPMAX) FCL = FCL / FPMAX
        DF  = DF*(D*TK - ONE)
        F   = F  + DF
              IF( DBLE(PK) .GT. PX ) GO TO 50
      IF(ABSC(DF) .GE. ABSC(F)*EPS)             GO TO 30
                NFP = PK - ZL - 1
                  ERR = EPS * SQRT(DBLE(NFP))
      CF1C = F
      RETURN
 1000 FORMAT(/' ',A6,': CF1 ACCURACY LOSS: D,DF,ACCH,K,ETA/K,ETA,X = ',
     X    /1X,1P,13D9.2/)
   50 IF(PR) WRITE (6,1010) CALLER,LIMIT,ABS(X)
 1010 FORMAT(' ',A6,': CF1 HAS FAILED TO CONVERGE AFTER ',I10  ,' ITERAT
     XIONS AS ABS(X) =',F15.0)
      ERR = TWO
      RETURN
      END
      FUNCTION CF2(X,ETA,ZL,PM,EPS,LIMIT,ERR,NPQ,ACC8,ACCH,
     X             PR,ACCUR,DELL,CALLER)
      IMPLICIT COMPLEX*16(A-H,O-Z)
      LOGICAL PR
      REAL*8 EPS,ERR,ACC8,ACCH,ACCUR,TA,RK,
     X       ABSC,ZERO,HALF,ONE,TWO
      CHARACTER*6 CALLER
      DATA ZERO,HALF,ONE,TWO / 0D+0, .5D+0, 1D+0, 2D+0 /
      ABSC(W) = ABS(DBLE(W)) + ABS(IMAG(W))
C
C                                    (omega)        (omega)
C *** Evaluate  CF2  = p + PM.q  =  H   (ETA,X)' / H   (ETA,X)
C                                    ZL             ZL
C     where PM = omega.i
C
      TA = TWO*LIMIT
      E2MM1 = ETA*ETA + ZL*ZL + ZL
      ETAP = ETA * PM
      XI = ONE/X
      WI = TWO*ETAP
      RK = ZERO
      PQ = (ONE - ETA*XI) * PM
      AA = -E2MM1 + ETAP
      BB = TWO*(X - ETA + PM)
         RL = XI * PM
      IF(ABSC(BB).LT.ACCH) THEN
         RL = RL * AA / (AA + RK + WI)
         PQ = PQ + RL * (BB + TWO*PM)
            AA = AA + TWO*(RK+ONE+WI)
            BB = BB + (TWO+TWO)*PM
            RK = RK + (TWO+TWO)
         ENDIF
      DD = ONE/BB
      DL = AA*DD* RL
   10 PQ    = PQ + DL
         RK = RK + TWO
         AA = AA + RK + WI
         BB = BB + TWO*PM
         DD = ONE/(AA*DD + BB)
         DL = DL*(BB*DD - ONE)
            ERR = ABSC(DL)/ABSC(PQ)
         IF(ERR.GE.MAX(EPS,ACC8*RK*HALF) .AND. RK.LE.TA) GO TO 10
C
         NPQ   = RK/TWO
         PQ    = PQ + DL
C           IF(PR.AND.NPQ.GE.LIMIT-1 .AND. ERR.GT.ACCUR)
C     X             WRITE(6,1000) CALLER,INT(IMAG(PM)),NPQ,ERR,ZL+DELL
C 1000 FORMAT(' ',A6,': CF2(',I2,') NOT CONVERGED FULLY IN ',I7,
C     X' ITERATIONS, SO ERROR IN IRREGULAR SOLUTION =',1P,D11.2,' AT ZL
C     X=', 0P,2F8.3)
      CF2 = PQ
      RETURN
      END

      FUNCTION CF1R(X,ETA,ZL,EPS,FCL,TPK1,ETANE0,LIMIT,ERR,NFP,
     X              ACCH,FPMIN,FPMAX,PR,CALLER)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL PR,ETANE0
      CHARACTER*6 CALLER
      DATA ONE,TWO / 1D+0, 2D+0 /
C
C
C ***    Evaluate CF1  =  F   =  F'(ZL,ETA,X)/F(ZL,ETA,X)
C
C        using real arithmetic
C
      FCL = ONE
      XI = ONE/X
      PK  = ZL + ONE
      PX  = PK  + LIMIT
   10 EK  = ETA / PK
        RK2 =          ONE + EK*EK
      F   = (EK + PK*XI)*FCL + (FCL - ONE)*XI
      PK1 =  PK + ONE
         TPK1 = PK + PK1
      TK  = TPK1*(XI + EK/PK1)
      IF(ETANE0) THEN
C ***   test ensures b1 .ne. zero for negative ETA etc.; fixup is exact.
             IF(ABS(TK) .GT. ACCH)  GO TO 20
             FCL  = RK2/(ONE + (ETA/PK1)**2)
             SL   = TPK1*XI * (TPK1+TWO)*XI
             PK   =  TWO + PK
             GO TO 10
         ENDIF
   20 D   =  ONE/TK
      DF  = -FCL*RK2*D
            IF(PK.GT.ZL+TWO) FCL = - RK2 * SL
            FCL = FCL * D * TPK1 * XI
      F   =  F  + DF
C
C ***   begin CF1 loop on PK = k = lambda + 1
C
      RK    = ONE
      SMALL    = SQRT(FPMIN)
   30 PK    = PK1
        PK1 = PK1 + ONE
         TPK1 = PK + PK1
         IF(ETANE0) THEN
           EK  = ETA / PK
           RK2 =          ONE + EK*EK
          ENDIF
        TK  = TPK1*(XI + EK/PK1)
        D   =  TK - D*RK2
              IF(ABS(D) .GT. ACCH)             GO TO 40
C              IF(PR) WRITE (6,1000) CALLER,D,DF,ACCH,PK,EK,ETA,X
              RK= RK +   ONE
              IF( RK .GT. TWO )                  GO TO 50
   40 D     = ONE/D
            FCL = FCL * D * TPK1*XI
            IF(ABS(FCL).LT.SMALL) FCL = FCL / SMALL
            IF(ABS(FCL).GT.FPMAX) FCL = FCL / FPMAX
        DF  = DF*(D*TK - ONE)
        F   = F  + DF
              IF( PK .GT. PX ) GO TO 50
      IF(ABS(DF) .GE. ABS(F)*EPS)             GO TO 30
                NFP = PK - ZL - 1
                  ERR = EPS * SQRT(DBLE(NFP))
      CF1R = F
      RETURN
 1000 FORMAT(/' ',A6,': CF1 ACCURACY LOSS: D,DF,ACCH,K,ETA/K,ETA,X = ',
     X    /1X,1P,7D9.2/)
   50 IF(PR) WRITE (6,1010) CALLER,LIMIT,ABS(X)
 1010 FORMAT(' ',A6,': CF1 HAS FAILED TO CONVERGE AFTER ',I10  ,' ITERAT
     XIONS AS ABS(X) =',F15.0)
      ERR = TWO
      RETURN
      END
      FUNCTION F20(AA,BB,Z,EPS,JMAX,RE,FPMAX,N,X)
C
C     evaluate the HYPERGEOMETRIC FUNCTION 2F0
C                                             i
C            F (AA,BB;;Z) = SUM  (AA)  (BB)  Z / i]
C           2 0              i       i     i
C
C     to accuracy EPS with at most JMAX terms.
C
C     if the terms start diverging,
C     the corresponding continued fraction is found by RCF
C     & evaluated progressively by Steed's method to obtain convergence.
C
C      useful number also input:  FPMAX = near-largest f.p. number
C
      IMPLICIT COMPLEX*16(A-H,O-Z)
      DIMENSION X(JMAX,4)
      LOGICAL FINITE
      REAL*8 EP,EPS,AT,ATL,ABSC,RE,FPMAX
      DATA ONE,ZERO / (1D+0,0D+0), (0D+0,0D+0) /
      ABSC(W) = ABS(DBLE(W)) + ABS(IMAG(W))
      NINTC(W) = NINT(REAL(REAL(W)))
C
      RE = 0.0
      X(1,1) = ONE
      SUM = X(1,1)
      ATL = ABSC(X(1,1))
         F    = SUM
         D = ONE
         DF   = SUM
      J = 0
      EP = EPS * JMAX *10.
      MA = - NINTC(AA)
      MB = - NINTC(BB)
      FINITE = ABS(ABS(DBLE(AA))-MA).LT.EP .AND. ABS(IMAG(AA)).LT.EP
     X    .OR. ABS(ABS(DBLE(BB))-MB).LT.EP .AND. ABS(IMAG(BB)).LT.EP
      IMAX = JMAX
      IF(FINITE.AND.MA.GE.0) IMAX = MIN(MA+1,IMAX)
      IF(FINITE.AND.MB.GE.0) IMAX = MIN(MB+1,IMAX)
      DO 10 I=2,IMAX
      X(I,1) = X(I-1,1) * Z * (AA+I-2) * (BB+I-2) / (I-1)
         IF(ABSC(X(I,1)).GT.FPMAX) GO TO 40
      AT = ABSC(X(I,1))
         IF(J.EQ.0) THEN
                 SUM = SUM + X(I,1)
                 IF(AT .LT. ABSC(SUM)*EPS) GO TO 20
               ENDIF
      IF(FINITE) GO TO 10
      IF(J.GT.0 .OR. AT.GT.ATL .OR. I.GE.JMAX-2) J = J + 1
         IF(J.EQ.0) GO TO 10
         CALL RCF(X(1,1),X(1,2),J,I,X(1,3),EPS)
              IF(I.LT.0) GO TO 40
            DO 50 K=MAX(J,2),I
            D = ONE/(D*X(K,2) + ONE)
            DF = DF*(D - ONE)
            F = F + DF
            IF(ABSC(DF) .LT. ABSC(F)*EPS) GO TO 30
            IF(DF.EQ.ZERO.AND.F.EQ.ZERO.AND.I.GE.4) GO TO 30
   50       CONTINUE
         J = I
   10 ATL = AT
      IF(.NOT.FINITE) I = -JMAX
   20 N = I
       F20 = SUM
       IF(.NOT.FINITE) RE  = AT / ABSC(SUM)
       RETURN
   30 F20 = F
      RE = ABSC(DF) / ABSC(F)
      N = K
      RETURN
   40 I = 0
      GO TO 20
      END
      FUNCTION CF1A(RHO,ETA,XL,PSI,EPS,NMAX,NUSED,FCL,RE,FPMAX,XX,G,C)
C
C     evaluate the ASYMPTOTIC EXPANSION for the
C            LOGARITHMIC DERIVATIVE OF THE REGULAR SOLUTION
C
C ***        CF1A  =  f   =  F'(XL,ETA,RHO)/F(XL,ETA,RHO)
C
C      that is valid for DBLE(RHO)>0, and best for RHO >> ETA**2, XL,
C      and is derived from the 2F0 expansions for H+ and H-
C      e.g. by Froeberg (Rev. Mod. Physics Vol 27, p399 , 1955)
C      Some lines of this subprogram are for convenience copied from
C           Takemasa, Tamura & Wolter CPC 17 (1979) 351.
C
C     Evaluate to accuracy EPS with at most NMAX terms.
C
C     If the terms start diverging,
C     the corresponding continued fraction is found by RCF
C     & evaluated progressively by Steed's method to obtain convergence.
C
C      useful number also input:  FPMAX = near-largest f.p. number
C
      IMPLICIT COMPLEX*16(A-H,O-Z)
      DIMENSION XX(2,NMAX),G(NMAX),C(NMAX)
      REAL*8 RE,EPS,T1,T2,T3,ZERO,ONE,TWO,AT,ATL,ABSC,FPMAX
      DATA ZERO,ONE,TWO,CI / 0D+0, 1D+0, 2D+0, (0D+0,1D+0) /
      ABSC(W) = ABS(DBLE(W)) + ABS(IMAG(W))
C
      HPI = TWO*ATAN(ONE)
      T1 = SIN(DBLE(PSI))
      T2 = COS(DBLE(PSI))
      ATL= TANH(IMAG(PSI))
C             GIVE COS(PSI)/COSH(IM(PSI)), WHICH ALWAYS HAS CORRECT SIGN
          COSL = DCMPLX( T2 , -T1 * ATL )
      TANL = DCMPLX(T1,T2*ATL) / COSL
      RE = ZERO
      XLL1= XL*(XL+ONE)
      ETASQ = ETA*ETA
      SL1=ONE
      SL=SL1
      SC1=ZERO
      SC=SC1
      TL1=SC
      TL=TL1
      TC1=ONE-ETA/RHO
      TC=TC1
      FCL  = TL + SL*TANL
      G(1) = (TC + SC*TANL) / FCL
      GLAST = G(1)
      ATL = ABSC(GLAST)
         F    = GLAST
         D = ONE
         DF   = GLAST
      J = 0
      DO 10 N=2,NMAX
      T1=N-1
      T2=TWO*T1-ONE
      T3=T1*(T1-ONE)
      DENOM=TWO*RHO*T1
      C1=(ETA*T2)/DENOM
      C2=(ETASQ+XLL1-T3)/DENOM
      SL2=C1*SL1-C2*TL1
      TL2=C1*TL1+C2*SL1
      SC2=C1*SC1-C2*TC1-SL2/RHO
      TC2=C1*TC1+C2*SC1-TL2/RHO
      SL=SL+SL2
      TL=TL+TL2
      SC=SC+SC2
      TC=TC+TC2
      SL1=SL2
      TL1=TL2
      SC1=SC2
      TC1=TC2
      FCL  =  TL + SL*TANL
         IF(ABSC(FCL).GT.FPMAX .OR. ABSC(FCL).LT.1./FPMAX) GO TO 40
      GSUM = (TC + SC*TANL) / FCL
      G(N) = GSUM - GLAST
      GLAST = GSUM
         AT = ABSC(G(N))
         IF(AT.LT.ABSC(GSUM)*EPS) GO TO 20
      IF(J.GT.0 .OR. AT.GT.ATL .OR. N.GE.NMAX-2) J = J + 1
         IF(J.EQ.0) GO TO 10
            CALL RCF(G,C,J,N,XX,EPS)
              IF(N.LT.0) GO TO 40
            DO 60 K=MAX(J,2),N
               D = ONE/(D*C(K) + ONE)
               DF = DF*(D - ONE)
               F = F + DF
         IF(ABSC(DF) .LT. ABSC(F)*EPS) GO TO 30
         IF(DF.EQ.ZERO.AND.F.EQ.ZERO.AND.N.GE.4) GO TO 30
   60         CONTINUE
         J = N
   10    ATL = AT
      K = -NMAX
      GO TO 30
   20 FCL = FCL * COSL
         CF1A = GSUM
         RE = AT / ABSC(GSUM)
         NUSED = N
         RETURN
   30 CF1A = F
      FCL = FCL * COSL
         RE = ABSC(DF) / ABSC(F)
         NUSED = K
      RETURN
   40 CF1A = G(1)
      FCL = 1.0
      RE = 1.0
      NUSED = 0
      RETURN
      END
      SUBROUTINE RCF(A,B,IBEG,INUM,XX,EPS)
C
C*******************************************************************
C
C  RCF converts polynomial A to the corresponding continued
C         fraction, in 'normal'  form with coefficients B
C         by the 'P algorithmn' of Patry & Gupta
C
C   A(z) = A1/z + A2/z**3 + A3/z**5 + ... + An/z**(2n-1)
C
C   B(z) = B1/z+ B2/z+ B3/z+ .../(z+ Bn/z)
C
C  data:
C   A     vector A(k), k=1,INUM         input
C   B     vector B(k), k=IBEG,INUM      output
C   IBEG  order of first coef. calc.    input
C   INUM  order of A, even or odd       input
C   XX    auxiliary vector of length .ge. length of vector B
C         caller provides space for A,B,XX
C     Note that neither of the first two terms A(1) A(2) should be zero
C             & the user can start the calculation with any value of
C                IBEG provided the c.f. coefs have been already
C                calculated up to INUM = IBEG-1
C             & the method breaks down as soon as the absolute value
C                of a c.f. coef. is less than EPS.    At the time of the
C                break up XX(1) has been replaced by 1E-50, and INUM has
C                been replaced by minus times the number of this coef.
C   algorithm: J.Patry & S.Gupta,
C              EIR-bericht nr. 247,
C              Eidg. Institut fur Reaktorforschung Wuerenlingen
C              Wueringlingen, Schweiz.
C              November 1973
C   see also:  Haenggi,Roesel & Trautmann,
C              Jnl. Computational Physics, vol 137, pp242-258 (1980)
C   note:      restart procedure modified by I.J.Thompson
C
C*******************************************************************
C
      IMPLICIT COMPLEX*16(A-H,O-Z)
      DIMENSION A(100),B(100),XX(2,100)
      LOGICAL EVEN
      REAL*8 EPS
      COMMON /RCFCM2/ X1,M2M1,MP12,EVEN,M      
!$OMP THREADPRIVATE(/RCFCM2/)
C     ibn = ibeg + inum - 1
      IBN = INUM
C                             B(IBN) is last value set on this call
      IF(IBEG.GT.4 .AND. M .NE. IBEG-1) GO TO 90
C                             B(M) is last value set in previous call
      IF(IBEG.GT.4) GO TO 50
      IF(IBEG.EQ.4) GO TO 20
      B(1) = A(1)
      IF(IBN.GE.2) B(2) = - A(2)/A(1)
      IF(IBN.LT.3) GO TO 10
      X0 = A(3) / A(2)
      XX(2,1) = B(2)
      XX(1,1) = - X0
      XX(1,2) = 0.
      B(3) = -X0 - B(2)
      X0 = -B(3) * A(2)
      M = 3
      MP12 = 2
      EVEN = .TRUE.
      IF(IBN.GT.3) GO TO 20
   10 RETURN
   20 IF(ABS(B(3)) .LT. EPS*ABS(X0)) GOTO 80
      M = 4
   30 X1 = A(M)
      M2M1 = MP12
      MP12 = M2M1 + 1
      IF(EVEN) MP12 = M2M1
      DO 40 K=2,MP12
   40 X1 = X1 + A(M-K+1) * XX(1,K-1)
      B(M) = - X1/X0
      IF(M.GE.IBN) RETURN
   50 IF(ABS(B(M)).LT.EPS*ABS(X0)) GO TO 80
      K = M2M1
   60 XX(2,K) = XX(1,K) + B(M) * XX(2,K-1)
      K = K-1
      IF(K.GT.1) GO TO 60
      XX(2,1) = XX(1,1) + B(M)
      DO 70 K=1,M2M1
      X0 = XX(2,K)
      XX(2,K) = XX(1,K)
   70 XX(1,K) = X0
      X0 = X1
      XX(1,M2M1+1) = 0.
      M = M+1
      EVEN = .NOT.EVEN
      GO TO 30
   80 INUM = -M
C     XX(1,1) = 1.E-50
C     PRINT 1000,M
C1000 FORMAT('0RCF: ZERO CF COEFFICIENT AT POSITION ',I4/)
      RETURN
   90 PRINT 1000,M,IBEG-1
 1000 FORMAT('0RCF: LAST CALL SET M =',I4,', BUT RESTART REQUIRES',I4)
      STOP
      END

      FUNCTION TIDY(Z,ACC)
C                     TIDY A COMPLEX NUMBER
      REAL*8 X,Y,ACC,AZ
      COMPLEX*16 Z,TIDY
C
      X = DBLE(Z)
      Y = IMAG(Z)
      AZ= (ABS(X) + ABS(Y)) * ACC * 5
      IF(ABS(X) .LT. AZ) X = 0D+0
      IF(ABS(Y) .LT. AZ) Y = 0D+0
      TIDY = DCMPLX(X,Y)
      RETURN
      END


 
