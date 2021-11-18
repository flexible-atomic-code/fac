C      ALGORITHM 745, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 21, NO. 3, September, 1995, P.  221-232.
C
C      Incorporates remark by Goano
C
C This file contains 3 files separated by lines of the form
C         C*** filename
C
* **********************************************************************
*
        SUBROUTINE FERMID(ORD, X, RELERR, FD, IERR)
*
* **********************************************************************
* FERMID returns in FD the value of the Fermi-Dirac integral of real
*        order ORD and real argument X, approximated with a relative
*        error RELERR.  FERMID is a driver routine that selects FDNINT
*        for integer ORD .LE. 0, FDNEG for X .LE. 0, and FDETA, FDPOS or
*        FDASYM for X .GT. 0.  A nonzero value is assigned to the error
*        flag IERR when an error condition occurs:
*           IERR = 1:  on input, the requested relative error RELERR is
*                      smaller than the machine precision;
*           IERR = 3:  an integral of large negative integer order could
*                      not be evaluated:  increase the parameter NMAX
*                      in subroutine FDNINT.
*           IERR = 4:  an integral (probably of small argument and large
*                      negative order) could not be evaluated with the
*                      requested accuracy after the inclusion of ITMAX
*                      terms of the series expansion:  increase the
*                      parameter ITMAX in the routine which produced the
*                      error message and in its subroutines.
*        When an error occurs, a message is also printed on the standard
*        output unit by the subroutine FERERR, and the execution of the
*        program is not interrupted; to change/suppress the output unit
*        or to stop the program when an error occurs, only FERERR should
*        be modified.
*
* References:
*
*   [1] M. Goano, "Series expansion of the Fermi-Dirac integral F_j(x)
*       over the entire domain of real j and x", Solid-State
*       Electronics, vol. 36, no. 2, pp. 217-221, 1993.
*
*   [2] J. S. Blakemore, "Approximation for Fermi-Dirac integrals,
*       especially the function F_1/2(eta) used to describe electron
*       density in a semiconductor", Solid-State Electronics, vol. 25,
*       no. 11, pp. 1067-1076, 1982.
*
* If a single precision version is desired, change all occurrences of
* *SP in columns 1 to 3 to blanks and comment the corresponding double
* precision statements.
*
* Michele Goano, Politecnico di Torino  (goano@polito.it).
* Latest revision:  March 23, 1994.
* **********************************************************************
*   Parameters
*SP     REAL             ONE, TEN, THREE, TWO, ZERO
        DOUBLE PRECISION ONE, TEN, THREE, TWO, ZERO
*SP     PARAMETER (ONE = 1.0E+0, TEN = 10.0E+0, THREE = 3.0E+0,
*SP  &             TWO = 2.0E+0, ZERO = 0.0E+0)
        PARAMETER (ONE = 1.0D+0, TEN = 10.0D+0, THREE = 3.0D+0,
     &             TWO = 2.0D+0, ZERO = 0.0D+0)
*   Scalar arguments
        INTEGER IERR
*SP     REAL             ORD, X, RELERR, FD
        DOUBLE PRECISION ORD, X, RELERR, FD
*   Local scalars
        LOGICAL INTORD, TRYASY
        INTEGER NORD
*SP     REAL             RKDIV, XASYMP
        DOUBLE PRECISION RKDIV, XASYMP
*   External subroutines
        EXTERNAL FDASYM, FDETA, FDNEG, FDNINT, FDPOS, FERERR
*   Intrinsic functions
*SP     INTRINSIC ABS, ANINT, EXP, LOG, LOG10, MAX, NINT, SQRT
        INTRINSIC ABS, ANINT, EXP, LOG, LOG10, MAX, NINT, SQRT
* ----------------------------------------------------------------------
*   Parameters of the floating-point arithmetic system.  Only the values
*   for very common machines are provided:  the subroutine MACHAR [1]
*   can be used to determine the machine constants of any other system.
*
*   [1] W. J. Cody,"Algorithm 665. MACHAR: A subroutine to dynamically
*       determine machine parameters", ACM Transactions on Mathematical
*       Software, vol. 14, no. 4, pp. 303-311, 1988.
*
        INTEGER MACHEP, MINEXP, MAXEXP, NEGEXP
*   ANSI/IEEE standard 745-1985:  IBM RISC 6000, DEC Alpha (S_floating
*   and T_floating), Apple Macintosh, SunSparc, most IBM PC compilers...
*SP     PARAMETER (MACHEP = -23, MINEXP = -126, MAXEXP = 128,
*SP  &             NEGEXP = -24)
        PARAMETER (MACHEP = -52, MINEXP = -1022, MAXEXP = 1024,
     &             NEGEXP = -53)
*   DEC VAX (F_floating and D_floating)
*SP     PARAMETER (MACHEP = -24, MINEXP = -128, MAXEXP = 127,
*SP  &             NEGEXP = -24)
*DP     PARAMETER (MACHEP = -56, MINEXP = -128, MAXEXP = 127,
*DP  &             NEGEXP = -56)
*   CRAY
*SP     PARAMETER (MACHEP = -47, MINEXP = -8193, MAXEXP = 8191,
*SP  &             NEGEXP = -47)
*DP     PARAMETER (MACHEP = -95, MINEXP = -8193, MAXEXP = 8191,
*DP  &             NEGEXP = -95)
*
*SP     REAL             BETA, EPS, XBIG, XMIN, XMAX
        DOUBLE PRECISION BETA, EPS, XBIG, XMIN, XMAX
        PARAMETER (BETA = TWO, EPS = BETA**MACHEP, XMIN = BETA**MINEXP,
     &             XMAX = (BETA**(MAXEXP-1) -
     &                                    BETA**(MAXEXP+NEGEXP-1))*BETA)
        XBIG = LOG(XMAX)
* ----------------------------------------------------------------------
        IERR = 0
        FD = ZERO
        INTORD = ABS(ORD - ANINT(ORD)).LE.ABS(ORD)*EPS
        IF (RELERR.LT.EPS) THEN
*   Test on the accuracy requested by the user
          IERR = 1
          CALL FERERR(' FERMID:  Input error: ' //
     &                ' RELERR is smaller than the machine precision')
        ELSE IF (INTORD .AND. ORD.LE.ZERO) THEN
*   Analytic expression for integer ORD .le. 0
          NORD = NINT(ORD)
          CALL FDNINT(NORD, X, FD, IERR)
        ELSE IF (X.LE.ZERO) THEN
*   Series expansion for negative argument
          CALL FDNEG(ORD, X, XMIN, RELERR, FD, IERR)
        ELSE
*   Positive argument:  approximations for k_div - 1 and x_min (RKDIV
*   and XASYMP)
          RKDIV = -LOG10(RELERR)
          XASYMP = MAX(ORD - ONE,
     &                 TWO*RKDIV - ORD*(TWO+RKDIV/TEN),
     &                 SQRT(ABS((2*RKDIV-ONE-ORD)*(2*RKDIV-ORD))))
          IF (X.GT.XASYMP .OR. INTORD) THEN
*   Asymptotic expansion, used also for positive integer order
            TRYASY = .TRUE.
            CALL FDASYM(ORD, X, EPS, XMAX, XMIN, RELERR, FD, IERR)
          ELSE
            TRYASY = .FALSE.
          END IF
          IF (.NOT.TRYASY .OR. IERR.NE.0) THEN
            IF (ORD.GT.-TWO .AND. X.LT.TWO/THREE) THEN
*   Taylor series expansion, involving eta function
              CALL FDETA(ORD, X, EPS, XMAX, RELERR, FD, IERR)
            ELSE
*   Series expansion for positive argument, involving confluent
*   hypergeometric functions
              CALL FDPOS(ORD, X, EPS, XMAX, RELERR, FD, IERR)
            END IF
          END IF
        END IF
        RETURN
        END

* **********************************************************************
*
        SUBROUTINE FERINC(ORD, X, B, RELERR, FDI, IERR)
*
* **********************************************************************
* FERINC returns in FDI the value of the incomplete Fermi-Dirac integral
*        of real order ORD and real arguments X and B, approximated with
*        a relative error RELERR.  Levin's u transform [2] is used to
*        sum the alternating series (21) of [1].  A nonzero value is
*        assigned to the error flag IERR when an error condition occurs:
*           IERR = 1:  on input, the requested relative error RELERR is
*                      smaller than the machine precision;
*           IERR = 2:  on input, the lower bound B of the incomplete
*                      integral is lower than zero;
*           IERR = 3:  a complete integral of very large negative
*                      integer order could not be evaluated:  increase
*                      the parameter NMAX in subroutine FDNINT.
*           IERR = 4:  an integral (probably of small argument and large
*                      negative order) could not be evaluated with the
*                      requested accuracy after the inclusion of ITMAX
*                      terms of the series expansion:  increase the
*                      parameter ITMAX in the routine which produced the
*                      error message and in its subroutines.
*        When an error occurs, a message is also printed on the standard
*        output unit by the subroutine FERERR, and the execution of the
*        program is not interrupted; to change/suppress the output unit
*        and/or to stop the program when an error occurs, only FERERR
*        should be modified.
*
* References:
*
*   [1] M. Goano, "Series expansion of the Fermi-Dirac integral F_j(x)
*       over the entire domain of real j and x", Solid-State
*       Electronics, vol. 36, no. 2, pp. 217-221, 1993.
*
*   [2] T. Fessler, W. F. Ford, D. A. Smith, "ALGORITHM 602. HURRY: An
*       acceleration algorithm for scalar sequences and series", ACM
*       Transactions on Mathematical Software, vol. 9, no. 3,
*       pp. 355-357, September  1983.
*
* If a single precision version is desired, change all occurrences of
* *SP in columns 1 to 3 to blanks and comment the corresponding double
* precision statements.
*
* Michele Goano, Politecnico di Torino  (goano@polito.it).
* Latest revision:  March 22, 1994.
* **********************************************************************
*   Parameters
        INTEGER ITMAX
        PARAMETER (ITMAX = 100)
*SP     REAL             ONE, TWO, ZERO
        DOUBLE PRECISION ONE, TWO, ZERO
*SP     PARAMETER (ONE = 1.0E+0, TWO = 2.0E+0, ZERO = 0.0E+0)
        PARAMETER (ONE = 1.0D+0, TWO = 2.0D+0, ZERO = 0.0D+0)
*   Scalar arguments
        INTEGER IERR
*SP     REAL             ORD, X, B, RELERR, FDI
        DOUBLE PRECISION ORD, X, B, RELERR, FDI
*   Local scalars
        LOGICAL LOGGAM
        INTEGER JTERM
*SP     REAL             BMX, BMXN, BN, EBMX, ENBMX, ENXMB,
*SP  &                   EXMB, FD, FDOLD, GAMMA, M, S, TERM, U,
*SP  &                   XMB, XMBN
        DOUBLE PRECISION BMX, BMXN, BN, EBMX, ENBMX, ENXMB,
     &                   EXMB, FD, FDOLD, GAMMA, M, S, TERM, U,
     &                   XMB, XMBN
*   Local arrays
*SP     REAL             QNUM(ITMAX), QDEN(ITMAX)
        DOUBLE PRECISION QNUM(ITMAX), QDEN(ITMAX)
*   External subroutines
        EXTERNAL FERERR, FERMID, GAMMAC, M1KUMM, U1KUMM, WHIZ
*   Intrinsic functions
        INTRINSIC ABS, ANINT, EXP, LOG
* ----------------------------------------------------------------------
*   Parameters of the floating-point arithmetic system.  Only the values
*   for very common machines are provided:  the subroutine MACHAR [1]
*   can be used to determine the machine constants of any other system.
*
*   [1] W. J. Cody,"Algorithm 665. MACHAR: A subroutine to dynamically
*       determine machine parameters", ACM Transactions on Mathematical
*       Software, vol. 14, no. 4, pp. 303-311, 1988.
*
        INTEGER MACHEP, MINEXP, MAXEXP, NEGEXP
*   ANSI/IEEE standard 745-1985:  IBM RISC 6000, DEC Alpha (S_floating
*   and T_floating), Apple Macintosh, SunSparc, most IBM PC compilers...
*SP     PARAMETER (MACHEP = -23, MINEXP = -126, MAXEXP = 128,
*SP  &             NEGEXP = -24)
        PARAMETER (MACHEP = -52, MINEXP = -1022, MAXEXP = 1024,
     &             NEGEXP = -53)
*   DEC VAX (F_floating and D_floating)
*SP     PARAMETER (MACHEP = -24, MINEXP = -128, MAXEXP = 127,
*SP  &             NEGEXP = -24)
*        PARAMETER (MACHEP = -56, MINEXP = -128, MAXEXP = 127,
*     &       NEGEXP = -56)
*   CRAY
*SP     PARAMETER (MACHEP = -47, MINEXP = -8193, MAXEXP = 8191,
*SP  &             NEGEXP = -47)
*DP     PARAMETER (MACHEP = -95, MINEXP = -8193, MAXEXP = 8191,
*DP  &             NEGEXP = -95)
*
*SP     REAL             EPS, XMIN, XMAX, XTINY
        DOUBLE PRECISION EPS, XMIN, XMAX, XTINY
        PARAMETER (EPS = TWO**MACHEP, XMIN = TWO**MINEXP, XMAX =
     &       (TWO**(MAXEXP-1) - TWO**(MAXEXP+NEGEXP-1))*TWO)
        
        XTINY = LOG(XMIN)
* ----------------------------------------------------------------------
        IERR = 0
        FDI = ZERO
        IF (RELERR.LT.EPS) THEN
*   Test on the accuracy requested by the user
          IERR = 1
          CALL FERERR(
     &  ' FERINC:  Input error:  RELERR smaller than machine precision')
        ELSE IF (B.LT.ZERO) THEN
*   Error in the argument B
          IERR = 2
          CALL FERERR(' FERINC:  Input error:  B is lower than zero')
        ELSE IF (B.EQ.ZERO) THEN
*   Complete integral
          CALL FERMID(ORD, X, RELERR, FDI, IERR)
        ELSE IF (ORD.LE.ZERO .AND.
     &                         ABS(ORD-ANINT(ORD)).LE.ABS(ORD)*EPS) THEN
*   Analytic expression for integer ORD .le. 0
          IF (NINT(ORD).EQ.0) THEN
            XMB = X - B
            IF (XMB.GE.ZERO) THEN
              FDI = XMB + LOG(ONE + EXP(-XMB))
            ELSE
              FDI = LOG(ONE + EXP(XMB))
            END IF
          ELSE
            FDI = ZERO
          END IF
        ELSE IF (B.LT.X) THEN
*   Series involving Kummer's function M
          CALL FERMID(ORD, X, RELERR, FD, IERR)
          CALL GAMMAC(ORD + TWO, EPS, XMAX, GAMMA, IERR)
          IF (IERR.EQ.-1) THEN
            LOGGAM = .TRUE.
            IERR = 0
          END IF
          BMX = B - X
          BMXN = BMX
          EBMX = -EXP(BMX)
          ENBMX = -EBMX
          BN = B
          FDI = XMAX
          DO 10 JTERM = 1, ITMAX
            FDOLD = FDI
            CALL M1KUMM(ORD, BN, EPS, XMAX, RELERR, M)
            TERM = ENBMX*M
            CALL WHIZ(TERM, JTERM, QNUM, QDEN, FDI, S)
*   Check truncation error and convergence
            BMXN = BMXN + BMX
            IF (ABS(FDI-FDOLD).LE.ABS(ONE-FDI)*RELERR .OR.
     &                                           BMXN.LT.XTINY) GO TO 20
            ENBMX = ENBMX*EBMX
            BN = BN + B
   10     CONTINUE
          IERR = 4
          CALL FERERR(
     &       ' FERINC:  RELERR not achieved:  increase parameter ITMAX')
   20     CONTINUE
          IF (LOGGAM) THEN
            FDI = FD - EXP((ORD+ONE)*LOG(B) - GAMMA)*(ONE - FDI)
          ELSE
            FDI = FD - B**(ORD + ONE)/GAMMA*(ONE - FDI)
          END IF
        ELSE
*   Series involving Kummer's function U
          CALL GAMMAC(ORD + ONE, EPS, XMAX, GAMMA, IERR)
          IF (IERR.EQ.-1) THEN
            LOGGAM = .TRUE.
            IERR = 0
          END IF
          XMB = X - B
          XMBN = XMB
          EXMB = -EXP(XMB)
          ENXMB = -EXMB
          BN = B
          FDI = XMAX
          DO 30 JTERM = 1, ITMAX
            FDOLD = FDI
            CALL U1KUMM(ORD + ONE, BN, EPS, XMAX, RELERR, U)
            TERM = ENXMB*U
            CALL WHIZ(TERM, JTERM, QNUM, QDEN, FDI, S)
c            write(*,*) term,jterm,enxmb,u,fdi
*   Check truncation error and convergence
            XMBN = XMBN + XMB
            IF (ABS(FDI-FDOLD).LE.ABS(FDI)*RELERR .OR.
     &                                           XMBN.LT.XTINY) GO TO 40
            ENXMB = ENXMB*EXMB
            BN = BN + B
   30     CONTINUE
          IERR = 4
          CALL FERERR(
     &       ' FERINC:  RELERR not achieved:  increase parameter ITMAX')
 40       CONTINUE
c          write(*,*) loggam,ord,x,b,gamma,xmax,fdi
          IF (LOGGAM) THEN
            FDI = EXP((ORD+ONE)*LOG(B) - GAMMA)*FDI
          ELSE
            FDI = B**(ORD + ONE)/GAMMA*FDI
          END IF
        END IF
        RETURN
        END

* **********************************************************************
*
        SUBROUTINE FDNINT(NORD, X, FD, IERR)
*
* **********************************************************************
* FDNINT returns in FD the value of the Fermi-Dirac integral of integer
*        order NORD (-NMAX-1 .LE. NORD .LE. 0) and argument X, for
*        which an analytical expression is available.  A nonzero value
*        is assigned to the error flag IERR when ABS(NORD).GT.NMAX+1:
*        to remedy, increase the parameter NMAX.
*
* If a single precision version is desired, change all occurrences of
* *SP in columns 1 to 3 to blanks and comment the corresponding double
* precision statements.
*
* Michele Goano, Politecnico di Torino  (goano@polito.it).
* Latest revision:  March 22, 1994.
* **********************************************************************
*   Parameters
        INTEGER NMAX
        PARAMETER (NMAX = 100)
*SP     REAL             ONE
        DOUBLE PRECISION ONE, ZERO
*SP     PARAMETER (ONE = 1.0E+0, ZERO = 0.0E+0)
        PARAMETER (ONE = 1.0D+0, ZERO = 0.0D+0)
*   Scalar arguments
        INTEGER NORD, IERR
*SP     REAL             X, FD
        DOUBLE PRECISION X, FD
*   Local scalars
        INTEGER I, K, N
*SP     REAL             A
        DOUBLE PRECISION A
*   Local arrays
*SP     REAL             QCOEF(NMAX)
        DOUBLE PRECISION QCOEF(NMAX)
*   External subroutines
        EXTERNAL FERERR
*   Intrinsic functions
*SP     INTRINSIC EXP, LOG, REAL
        INTRINSIC DBLE, EXP, LOG
* ----------------------------------------------------------------------
        IERR = 0
        FD = ZERO
*   Test on the order, whose absolute value must be lower or equal than
*   NMAX+1
        IF (NORD.LT.-NMAX - 1) THEN
          IERR = 3
          CALL FERERR(
     &            ' FDNINT:  order too large:  increase parameter NMAX')
        ELSE IF (NORD.EQ.0) THEN
*   Analytic expression for NORD .eq. 0
          IF (X.GE.ZERO) THEN
            FD = X + LOG(ONE + EXP(-X))
          ELSE
            FD = LOG(ONE + EXP(X))
          END IF
        ELSE IF (NORD.EQ.-1) THEN
*   Analytic expression for NORD .eq. -1
          IF (X.GE.ZERO) THEN
            FD = ONE/(ONE + EXP(-X))
          ELSE
            A = EXP(X)
            FD = A/(ONE + A)
          END IF
        ELSE
*   Evaluation of the coefficients of the polynomial P(a), having degree
*   (-NORD - 2), appearing at the numerator of the analytic expression
*   for NORD .le. -2
          N = -NORD - 1
          QCOEF(1) = ONE
          DO 20 K = 2, N
            QCOEF(K) = -QCOEF(K - 1)
            DO 10 I = K - 1, 2, -1
*SP           QCOEF(I) = REAL(I)*QCOEF(I) - REAL(K - (I-1))*QCOEF(I - 1)
              QCOEF(I) = DBLE(I)*QCOEF(I) - DBLE(K - (I-1))*QCOEF(I - 1)
   10       CONTINUE
   20     CONTINUE
*   Computation of P(a)
          IF (X.GE.ZERO) THEN
            A = EXP(-X)
            FD = QCOEF(1)
            DO 30 I = 2, N
              FD = FD*A + QCOEF(I)
   30       CONTINUE
          ELSE
            A = EXP(X)
            FD = QCOEF(N)
            DO 40 I = N - 1, 1, -1
              FD = FD*A + QCOEF(I)
   40       CONTINUE
          END IF
*   Evaluation of the Fermi-Dirac integral
          FD = FD*A*(ONE + A)**NORD
        END IF
        RETURN
        END


* **********************************************************************
*
        SUBROUTINE FDNEG(ORD, X, XMIN, RELERR, FD, IERR)
*
* **********************************************************************
* FDNEG returns in FD the value of the Fermi-Dirac integral of real
*       order ORD and negative argument X, approximated with a relative
*       error RELERR.  XMIN represent the smallest non-vanishing
*       floating-point number.  Levin's u transform [2] is used to sum
*       the alternating series (13) of [1].
*
* References:
*
*   [1] J. S. Blakemore, "Approximation for Fermi-Dirac integrals,
*       especially the function F_1/2(eta) used to describe electron
*       density in a semiconductor", Solid-State Electronics, vol. 25,
*       no. 11, pp. 1067-1076, 1982.
*
*   [2] T. Fessler, W. F. Ford, D. A. Smith, "ALGORITHM 602. HURRY: An
*       acceleration algorithm for scalar sequences and series", ACM
*       Transactions on Mathematical Software, vol. 9, no. 3,
*       pp. 355-357, September  1983.
*
* If a single precision version is desired, change all occurrences of
* *SP in columns 1 to 3 to blanks and comment the corresponding double
* precision statements.
*
* Michele Goano, Politecnico di Torino  (goano@polito.it).
* Latest revision:  February 5, 1996.
* **********************************************************************
*   Parameters
        INTEGER ITMAX
        PARAMETER (ITMAX = 100)
*SP     REAL             ONE, ZERO
        DOUBLE PRECISION ONE, ZERO
*SP     PARAMETER (ONE = 1.0E+0, ZERO = 0.0E+0)
        PARAMETER (ONE = 1.0D+0, ZERO = 0.0D+0)
*   Scalar arguments
        INTEGER IERR
*SP     REAL             ORD, X, XMIN, RELERR, FD
        DOUBLE PRECISION ORD, X, XMIN, RELERR, FD
*   Local scalars
        INTEGER JTERM
*SP     REAL             EX, ENX, FDOLD, S, TERM, XN, XTINY
        DOUBLE PRECISION EX, ENX, FDOLD, S, TERM, XN, XTINY
*   Local arrays
*SP     REAL             QNUM(ITMAX), QDEN(ITMAX)
        DOUBLE PRECISION QNUM(ITMAX), QDEN(ITMAX)
*   External subroutines
        EXTERNAL FERERR, WHIZ
*   Intrinsic functions
*SP     INTRINSIC ABS, EXP, LOG, REAL
        INTRINSIC ABS, DBLE, EXP, LOG
* ----------------------------------------------------------------------
        IERR = 0
        FD = ZERO
        XTINY = LOG(XMIN)
*
        IF (X.GT.XTINY) THEN
          XN = X
          EX = -EXP(X)
          ENX = -EX
          DO 10 JTERM = 1, ITMAX
            FDOLD = FD
*SP         TERM = ENX/REAL(JTERM)**(ORD + ONE)
            TERM = ENX/DBLE(JTERM)**(ORD + ONE)
            CALL WHIZ(TERM, JTERM, QNUM, QDEN, FD, S)
*   Check truncation error and convergence
            XN = XN + X
            IF (ABS(FD-FDOLD).LE.ABS(FD)*RELERR .OR. XN.LT.XTINY) RETURN
            ENX = ENX*EX
   10     CONTINUE
          IERR = 4
          CALL FERERR(
     &        ' FDNEG:  RELERR not achieved:  increase parameter ITMAX')
        END IF
        RETURN
        END
* **********************************************************************
*
        SUBROUTINE FDPOS(ORD, X, EPS, XMAX, RELERR, FD, IERR)
*
* **********************************************************************
* FDPOS returns in FD the value of the Fermi-Dirac integral of real
*       order ORD and argument X .GT. 0, approximated with a relative
*       error RELERR.  EPS and XMAX represent the smallest positive
*       floating-point number such that 1.0+EPS .NE. 1.0, and the
*       largest finite floating-point number, respectively.  Levin's u
*       transform [2] is used to sum the alternating series (11) of [1].
*
* References:
*
*   [1] M. Goano, "Series expansion of the Fermi-Dirac integral F_j(x)
*       over the entire domain of real j and x", Solid-State
*       Electronics, vol. 36, no. 2, pp. 217-221, 1993.
*
*   [2] T. Fessler, W. F. Ford, D. A. Smith, "ALGORITHM 602. HURRY: An
*       acceleration algorithm for scalar sequences and series", ACM
*       Transactions on Mathematical Software, vol. 9, no. 3,
*       pp. 355-357, September  1983.
*
* If a single precision version is desired, change all occurrences of
* *SP in columns 1 to 3 to blanks and comment the corresponding double
* precision statements.
*
* Michele Goano, Politecnico di Torino  (goano@polito.it).
* Latest revision:  March 22, 1994.
* **********************************************************************
*   Parameters
        INTEGER ITMAX
        PARAMETER (ITMAX = 100)
*SP     REAL             ONE, TWO, ZERO
        DOUBLE PRECISION ONE, TWO, ZERO
*SP     PARAMETER (ONE = 1.0E+0, TWO = 2.0E+0, ZERO = 0.0E+0)
        PARAMETER (ONE = 1.0D+0, TWO = 2.0D+0, ZERO = 0.0D+0)
*   Scalar arguments
        INTEGER IERR
*SP     REAL             ORD, X, EPS, XMAX, RELERR, FD
        DOUBLE PRECISION ORD, X, EPS, XMAX, RELERR, FD
*   Local scalars
        LOGICAL LOGGAM
        INTEGER JTERM
*SP     REAL             FDOLD, GAMMA, M, S, SEGN, TERM, U, XN
        DOUBLE PRECISION FDOLD, GAMMA, M, S, SEGN, TERM, U, XN
*   Local arrays
*SP     REAL             QNUM(ITMAX), QDEN(ITMAX)
        DOUBLE PRECISION QNUM(ITMAX), QDEN(ITMAX)
*   External subroutines
        EXTERNAL FERERR, GAMMAC, M1KUMM, U1KUMM, WHIZ
*   Intrinsic functions
        INTRINSIC ABS, EXP, LOG
* ----------------------------------------------------------------------
        IERR = 0
        FD = ZERO
*
        CALL GAMMAC(ORD + TWO, EPS, XMAX, GAMMA, IERR)
        IF (IERR.EQ.-1) LOGGAM = .TRUE.
        SEGN = ONE
        XN = X
        FD = XMAX
        DO 10 JTERM = 1, ITMAX
          FDOLD = FD
          CALL U1KUMM(ORD + ONE, XN, EPS, XMAX, RELERR, U)
          CALL M1KUMM(ORD, XN, EPS, XMAX, RELERR, M)
          TERM = SEGN*((ORD+ONE)*U - M)
          CALL WHIZ(TERM, JTERM, QNUM, QDEN, FD, S)
*   Check truncation error and convergence
          IF (ABS(FD-FDOLD).LE.ABS(FD+ONE)*RELERR) GO TO 20
          SEGN = -SEGN
          XN = XN + X
   10   CONTINUE
        IERR = 4
        CALL FERERR(
     &        ' FDPOS:  RELERR not achieved:  increase parameter ITMAX')
   20   CONTINUE
        IF (LOGGAM) THEN
          FD = EXP((ORD+ONE)*LOG(X) - GAMMA)*(ONE + FD)
        ELSE
          FD = X**(ORD + ONE)/GAMMA*(ONE + FD)
        END IF
        RETURN
        END

* **********************************************************************
*
        SUBROUTINE FDETA(ORD, X, EPS, XMAX, RELERR, FD, IERR)
*
* **********************************************************************
* FDETA returns in FD the value of the Fermi-Dirac integral of real
*       order ORD and argument X such that ABS(X) .LE. PI, approximated
*       with a relative error RELERR.  EPS and XMAX represent the
*       smallest positive floating-point number such that
*       1.0+EPS .NE. 1.0, and the largest finite floating-point number,
*       respectively.  Taylor series expansion (4) of [1] is used,
*       involving eta function defined in (23.2.19) of [2].
*
*
* References:
*
*   [1] W. J. Cody and H. C. Thacher, Jr., "Rational Chebyshev
*       approximations for Fermi-Dirac integrals of orders -1/2, 1/2 and
*       3/2", Mathematics of Computation, vol. 21, no. 97, pp. 30-40,
*       1967.
*
*   [2] E. V. Haynsworth and K. Goldberg, "Bernoulli and Euler
*       Polynomials - Riemann Zeta Function", in "Handbook of
*       Mathematical Functions with Formulas, Graphs and Mathematical
*       Tables" (M. Abramowitz and I. A. Stegun, eds.), no. 55 in
*       National Bureau of Standards Applied Mathematics Series, ch. 23,
*       pp. 803-819, Washington, D.C.:  U.S. Government Printing Office,
*       1964.
*
* If a single precision version is desired, change all occurrences of
* *SP in columns 1 to 3 to blanks and comment the corresponding double
* precision statements.
*
* Michele Goano, Politecnico di Torino  (goano@polito.it).
* Latest revision:  March 22, 1994.
* **********************************************************************
*   Parameters
        INTEGER ITMAX
        PARAMETER (ITMAX = 100)
*SP     REAL             ONE, PI, TWO, ZERO
        DOUBLE PRECISION ONE, PI, TWO, ZERO
*SP     PARAMETER (ONE = 1.0E+0, PI = 3.141592653589793238462643E+0,
*SP  &             TWO = 2.0E+0, ZERO = 0.0E+0)
        PARAMETER (ONE = 1.0D+0, PI = 3.141592653589793238462643D+0,
     &             TWO = 2.0D+0, ZERO = 0.0D+0)
*   Scalar arguments
        INTEGER IERR
*SP     REAL             ORD, X, EPS, XMAX, RELERR, FD
        DOUBLE PRECISION ORD, X, EPS, XMAX, RELERR, FD
*   Local scalars
        LOGICAL OKJM1, OKJM2
        INTEGER JTERM
*SP     REAL             ETA, RJTERM, TERM, XNOFAC
        DOUBLE PRECISION ETA, RJTERM, TERM, XNOFAC
*   External subroutines
        EXTERNAL ETARIE, FERERR
*   Intrinsic functions
*SP     INTRINSIC ABS, REAL
        INTRINSIC ABS, DBLE
* ----------------------------------------------------------------------
        IERR = 0
        FD = ZERO
*
        OKJM1 = .FALSE.
        OKJM2 = .FALSE.
        XNOFAC = ONE
        DO 10 JTERM = 1, ITMAX
*SP       RJTERM = REAL(JTERM)
          RJTERM = DBLE(JTERM)
          CALL ETARIE(ORD + TWO - RJTERM, EPS, XMAX, RELERR, ETA)
          TERM = ETA*XNOFAC
          FD = FD + TERM
*   Check truncation error and convergence.  The summation is terminated
*   when three consecutive terms of the series satisfy the bound on the
*   relative error
          IF (ABS(TERM).GT.ABS(FD)*RELERR) THEN
            OKJM1 = .FALSE.
            OKJM2 = .FALSE.
          ELSE IF (.NOT.OKJM1) THEN
            OKJM1 = .TRUE.
          ELSE IF (OKJM2) THEN
            RETURN
          ELSE
            OKJM2 = .TRUE.
          END IF
          XNOFAC = XNOFAC*X/RJTERM
   10   CONTINUE
        IERR = 4
        CALL FERERR(
     &        ' FDETA:  RELERR not achieved:  increase parameter ITMAX')
        RETURN
        END

* **********************************************************************
*
        SUBROUTINE FDASYM(ORD, X, EPS, XMAX, XMIN, RELERR, FD, IERR)
*
* **********************************************************************
* FDASYM returns in FD the value of the Fermi-Dirac integral of real
*        order ORD and argument X .GT. 0, approximated with a relative
*        error RELERR by means of an asymptotic expansion.  EPS, XMAX
*        and XMIN represent the smallest positive floating-point number
*        such that 1.0+EPS .NE. 1.0, the largest finite floating-point
*        number, and the smallest non-vanishing floating-point number,
*        respectively.  A nonzero value is assigned to the error flag
*        IERR when the series does not converge.  The expansion always
*        terminates after a finite number of steps in case of integer
*        ORD.
*
* References:
*
*   [1] P. Rhodes, "Fermi-Dirac function of integral order", Proceedings
*       of the Royal Society of London. Series A - Mathematical and
*       Physical Sciences, vol. 204, pp. 396-405, 1950.
*
*   [2] R. B. Dingle, "Asymptotic Expansions: Their Derivation and
*       Interpretation", London and New York:  Academic Press, 1973.
*
* If a single precision version is desired, change all occurrences of
* *SP in columns 1 to 3 to blanks and comment the corresponding double
* precision statements.
*
* Michele Goano, Politecnico di Torino  (goano@polito.it).
* Latest revision:  March 22, 1994.
* **********************************************************************
*   Parameters
        INTEGER ITMAX
        PARAMETER (ITMAX = 100)
*SP     REAL             HALF, ONE, PI, TWO, ZERO
        DOUBLE PRECISION HALF, ONE, PI, TWO, ZERO
*SP     PARAMETER (HALF = 0.5E+0, ONE = 1.0E+0,
*SP  &             PI = 3.141592653589793238462643E+0, TWO = 2.0E+0,
*SP  &             ZERO = 0.0E+0)
        PARAMETER (HALF = 0.5D+0, ONE = 1.0D+0,
     &             PI = 3.141592653589793238462643D+0, TWO = 2.0D+0,
     &             ZERO = 0.0D+0)
*   Scalar arguments
        INTEGER IERR
*SP     REAL             ORD, X, EPS, XMAX, XMIN, RELERR, FD
        DOUBLE PRECISION ORD, X, EPS, XMAX, XMIN, RELERR, FD
*   Local scalars
        LOGICAL LOGGAM
        INTEGER N
*SP     REAL             ADD, ADDOLD, ETA, GAMMA, SEQN, XGAM, XM2
        DOUBLE PRECISION ADD, ADDOLD, ETA, GAMMA, SEQN, XGAM, XM2
*   External subroutines
        EXTERNAL ETAN, FDNEG, FERERR, GAMMAC
*   Intrinsic functions
*SP     INTRINSIC ABS, ANINT, COS, EXP, LOG, REAL
        INTRINSIC ABS, ANINT, COS, DBLE, EXP, LOG
* ----------------------------------------------------------------------
        IERR = 0
        FD = ZERO
*
        CALL GAMMAC(ORD + TWO, EPS, XMAX, GAMMA, IERR)
        IF (IERR.EQ.-1) THEN
          LOGGAM = .TRUE.
          IERR = 0
        END IF
        SEQN = HALF
        XM2 = X**(-2)
        XGAM = ONE
        ADD = XMAX
        DO 10 N = 1, ITMAX
          ADDOLD = ADD
*SP       XGAM = XGAM*XM2*(ORD + ONE - REAL(2*N-2))*
*SP  &                    (ORD + ONE - REAL(2*N-1))
          XGAM = XGAM*XM2*(ORD + ONE - DBLE(2*N-2))*
     &                    (ORD + ONE - DBLE(2*N-1))
          CALL ETAN(2*N, ETA)
          ADD = ETA*XGAM
          IF (ABS(ADD).GE.ABS(ADDOLD) .AND.
     &                       ABS(ORD - ANINT(ORD)).GT.ABS(ORD)*EPS) THEN
*   Asymptotic series is diverging
            IERR = 1
            RETURN
          END IF
          SEQN = SEQN + ADD
*   Check truncation error and convergence
          IF (ABS(ADD).LE.ABS(SEQN)*RELERR) GO TO 20
   10   CONTINUE
        IERR = 4
        CALL FERERR(
     &       ' FDASYM:  RELERR not achieved:  increase parameter ITMAX')
   20   CONTINUE
        CALL FDNEG(ORD, -X, XMIN, RELERR, FD, IERR)
        IF (LOGGAM) THEN
          FD = COS(ORD*PI)*FD + TWO*SEQN*EXP((ORD + ONE)*LOG(X) - GAMMA)
        ELSE
          FD = COS(ORD*PI)*FD + X**(ORD + ONE)*TWO*SEQN/GAMMA
        END IF
        RETURN
        END

* **********************************************************************
*
        SUBROUTINE M1KUMM(A, X, EPS, XMAX, RELERR, M)
*
* **********************************************************************
* M1KUMM returns in M the value of Kummer's confluent hypergeometric
*        function M(1,2+A,-X), defined in (13.1.2) of [1], for real
*        arguments A and X, approximated with a relative error RELERR.
*        EPS and XMAX represent the smallest positive floating-point
*        number such that 1.0+EPS .NE. 1.0, and the largest finite
*        floating-point number, respectively.  Asymptotic expansion [1]
*        or continued fraction representation [2] is used.
*        Renormalization is carried out as proposed in [3].
*
* References:
*
*   [1] L. J. Slater, "Confluent Hypergeometric Functions", in "Handbook
*       of Mathematical Functions with Formulas, Graphs and Mathematical
*       Tables" (M. Abramowitz and I. A. Stegun, eds.), no. 55 in
*       National Bureau of Standards Applied Mathematics Series, ch. 13,
*       pp. 503-535, Washington, D.C.:  U.S. Government Printing Office,
*       1964.
*
*   [2] P. Henrici, "Applied and Computational Complex Analysis.
*       Volume 2.  Special Functions-Integral Transforms-Asymptotics-
*       Continued Fractions", New York:  John Wiley & Sons, 1977.
*
*   [3] W. H. Press, B. P. Flannery, S. A. Teukolsky, W. T. Vetterling,
*       "Numerical Recipes. The Art of Scientific Computing", Cambridge:
*       Cambridge University Press, 1986.
*
* If a single precision version is desired, change all occurrences of
* *SP in columns 1 to 3 to blanks and comment the corresponding double
* precision statements.
*
* Michele Goano, Politecnico di Torino  (goano@polito.it).
* Latest revision:  March 22, 1994.
* **********************************************************************
*   Parameters
        INTEGER ITMAX
        PARAMETER (ITMAX = 100)
*SP     REAL             ONE, PI, TEN, THREE, TWO, ZERO
        DOUBLE PRECISION ONE, PI, TEN, THREE, TWO, ZERO
*SP     PARAMETER (ONE = 1.0E+0, PI = 3.141592653589793238462643E+0,
*SP  &             TEN = 10.0E+0, THREE = 3.0E+0, TWO = 2.0E+0,
*SP  &             ZERO = 0.0E+0)
        PARAMETER (ONE = 1.0D+0, PI = 3.141592653589793238462643D+0,
     &             TEN = 10.0D+0, THREE = 3.0D+0, TWO = 2.0D+0,
     &             ZERO = 0.0D+0)
*   Scalar arguments
*SP     REAL             A, X, EPS, XMAX, RELERR, M
        DOUBLE PRECISION A, X, EPS, XMAX, RELERR, M
*   Local scalars
        LOGICAL OKASYM
        INTEGER IERR, N
*SP     REAL             AA, ADD, ADDOLD, BB, FAC, GAMMA, GOLD, MLOG,
*SP  &                   P1, P2, Q1, Q2, RKDIV, RN, XASYMP, XBIG
        DOUBLE PRECISION AA, ADD, ADDOLD, BB, FAC, GAMMA, GOLD, MLOG,
     &                   P1, P2, Q1, Q2, RKDIV, RN, XASYMP, XBIG
*   External subroutines
        EXTERNAL GAMMAC
*   Intrinsic functions
*SP     INTRINSIC ABS, COS, EXP, LOG, LOG10, REAL
        INTRINSIC ABS, COS, EXP, LOG, LOG10, DBLE
* ----------------------------------------------------------------------
        XBIG = LOG(XMAX)
        OKASYM = .TRUE.
        M = ZERO
*   Special cases
        IF (X.EQ.ZERO) THEN
          M = ONE
        ELSE
*   Approximations for k_div - 1 and x_min (RKDIV and XASYMP)
          RKDIV = -TWO/THREE*LOG10(RELERR)
          XASYMP = MAX(A - ONE,
     &                 TWO + RKDIV - A*(TWO+RKDIV/TEN),
     &                 ABS(RKDIV-A))
          IF (X.GT.XASYMP) THEN
*   Asymptotic expansion
            CALL GAMMAC(A+ONE, EPS, XMAX, GAMMA, IERR)
            IF (IERR.EQ.-1) THEN
*   Handling of the logarithm of the gamma function to avoid overflow
              MLOG = GAMMA - X - A*LOG(X)
              IF (MLOG.LT.XBIG) THEN
                M = ONE - COS(PI*A)*EXP(MLOG)
              ELSE
                OKASYM = .FALSE.
                GO TO 20
              END IF
            ELSE
              M = ONE - COS(PI*A)*GAMMA*EXP(-X)/X**A
            END IF
            ADDOLD = XMAX
            ADD = -A/X
            DO 10 N = 1, ITMAX
*   Divergence
              IF (ABS(ADD).GE.ABS(ADDOLD)) THEN
                OKASYM = .FALSE.
                GO TO 20
              END IF
              M = M + ADD
*   Check truncation error and convergence
              IF (ABS(ADD).LE.ABS(M)*RELERR) THEN
                M = M*(A + ONE)/X
                RETURN
              END IF
              ADDOLD = ADD
*SP           ADD = -ADD*(A - REAL(N))/X
              ADD = -ADD*(A - DBLE(N))/X
   10       CONTINUE
          END IF
*   Continued fraction:  initial conditions
   20     CONTINUE
          GOLD = ZERO
          P1 = ONE
          Q1 = ONE
          P2 = A + TWO
          Q2 = X + A + TWO
          BB = A + TWO
*   Initial value of the normalization factor
          FAC = ONE
          DO 30 N = 1, ITMAX
*   Evaluation of a_(2N+1) and b_(2N+1)
*SP         RN = REAL(N)
            RN = DBLE(N)
            AA = -RN*X
            BB = BB + ONE
            P1 = (AA*P1 + BB*P2)*FAC
            Q1 = (AA*Q1 + BB*Q2)*FAC
*   Evaluation of a_(2N+2) and b_(2N+2)
            AA = (A + RN + ONE)*X
            BB = BB + ONE
            P2 = BB*P1 + AA*P2*FAC
            Q2 = BB*Q1 + AA*Q2*FAC
            IF (Q2.NE.ZERO) THEN
*   Renormalization and evaluation of w_(2N+2)
              FAC = ONE/Q2
              M = P2*FAC
*   Check truncation error and convergence
              IF (ABS(M-GOLD).LT.ABS(M)*RELERR) RETURN
              GOLD = M
            END IF
   30     CONTINUE
        END IF
        RETURN
        END

* **********************************************************************
*
        SUBROUTINE U1KUMM(A, X, EPS, XMAX, RELERR, U)
*
* **********************************************************************
* U1KUMM returns in U the value of Kummer's confluent hypergeometric
*        function U(1,1+A,X), defined in (13.1.3) of [1], for real
*        arguments A and X, approximated with a relative error RELERR.
*        EPS and XMAX represent the smallest positive floating-point
*        number such that 1.0+EPS .NE. 1.0, and the largest finite
*        floating-point number, respectively.  The relation with the
*        incomplete gamma function is exploited, by means of (13.6.28)
*        and (13.1.29) of [1].  For A .LE. 0 an expansion in terms of
*        Laguerre polynomials is used [3].  Otherwise the recipe of [4]
*        is followed: series expansion (6.5.29) of [2] if X .LT. A+1,
*        continued fraction (6.5.31) of [2] if X .GE. A+1.
*
* References:
*
*   [1] L. J. Slater, "Confluent Hypergeometric Functions", ch. 13 in
*       [5], pp. 503-535.
*
*   [2] P. J. Davis, "Gamma Function and Related Functions", ch. 6 in
*       [5], pp. 253-293.
*
*   [3] P. Henrici, "Computational Analysis with the HP-25 Pocket
*       Calculator", New York:  John Wiley & Sons, 1977.
*
*   [4] W. H. Press, B. P. Flannery, S. A. Teukolsky, W. T. Vetterling,
*       "Numerical Recipes. The Art of Scientific Computing", Cambridge:
*       Cambridge University Press, 1986.
*
*   [5] M. Abramowitz and I. A. Stegun (eds.), "Handbook of Mathematical
*       Functions with Formulas, Graphs and Mathematical Tables", no. 55
*       in National Bureau of Standards Applied Mathematics Series,
*       Washington, D.C.:  U.S. Government Printing Office, 1964.
*
* If a single precision version is desired, change all occurrences of
* *SP in columns 1 to 3 to blanks and comment the corresponding double
* precision statements.
*
* Michele Goano, Politecnico di Torino  (goano@polito.it).
* Latest revision:  March 22, 1994.
* **********************************************************************
*   Parameters
        INTEGER ITMAX
        PARAMETER (ITMAX = 100)
*SP     REAL             ONE, ZERO
        DOUBLE PRECISION ONE, ZERO
*SP     PARAMETER (ONE = 1.0E+0, ZERO = 0.0E+0)
        PARAMETER (ONE = 1.0D+0, ZERO = 0.0D+0)
*   Scalar arguments
*SP     REAL             A, X, EPS, XMAX, RELERR, U
        DOUBLE PRECISION A, X, EPS, XMAX, RELERR, U
*   Local scalars
        LOGICAL LOGGAM
        INTEGER IERR, N
*SP     REAL             A0, A1, ANA, ANF, AP, B0, B1, DEL, FAC, G,
*SP  &                   GAMMA, GOLD, PLAGN, PLAGN1, PLAGN2, RN, T,
*SP  &                   ULOG, XBIG
        DOUBLE PRECISION A0, A1, ANA, ANF, AP, B0, B1, DEL, FAC, G,
     &                   GAMMA, GOLD, PLAGN, PLAGN1, PLAGN2, RN, T,
     &                   ULOG, XBIG
*   External subroutines
        EXTERNAL GAMMAC
*   Intrinsic functions
*SP     INTRINSIC ABS, EXP, LOG, REAL
        INTRINSIC ABS, DBLE, EXP, LOG
* ----------------------------------------------------------------------
        XBIG = LOG(XMAX)
        U = ZERO
*   Special cases
        IF (X.EQ.ZERO) THEN
          U = -ONE/A
*   Laguerre polynomials
        ELSE IF (A.LE.ZERO) THEN
          U = ZERO
          PLAGN2 = ZERO
          PLAGN1 = ONE
          G = ONE
          DO 10 N = 1, ITMAX
*SP         RN = REAL(N)
            RN = DBLE(N)
            PLAGN = ((RN-A-ONE)*(PLAGN1-PLAGN2) + (RN+X)*PLAGN1)/RN
            T = G/(PLAGN1*PLAGN)
            U = U + T
            IF (ABS(T).LT.ABS(U)*RELERR) RETURN
            G = G*(RN - A)/(RN + ONE)
            PLAGN2 = PLAGN1
            PLAGN1 = PLAGN
   10     CONTINUE
*   Series expansion
        ELSE IF (X.LT.A+ONE) THEN
          CALL GAMMAC(A, EPS, XMAX, GAMMA, IERR)
          IF (IERR.EQ.-1) LOGGAM = .TRUE.
          AP = A
          U = ONE/A
          DEL = U
          DO 20 N = 1, ITMAX
            AP = AP + ONE
            DEL = DEL*X/AP
            U = U + DEL
            IF (ABS(DEL).LT.ABS(U)*RELERR) THEN
              IF (LOGGAM) THEN
                ULOG = GAMMA + X + A*LOG(X)
                U = EXP(ULOG) - U
              ELSE
                U = GAMMA*EXP(X)/X**A - U
              END IF
              RETURN
            END IF
   20     CONTINUE
*   Continued fraction
        ELSE
          GOLD = ZERO
          A0 = ONE
          A1 = X
          B0 = ZERO
          B1 = ONE
          FAC = ONE
          DO 30 N = 1, ITMAX
*SP         RN = REAL(N)
            RN = DBLE(N)
            ANA = RN - A
            A0 = (A1 + A0*ANA)*FAC
            B0 = (B1 + B0*ANA)*FAC
            ANF = RN*FAC
            A1 = X*A0 + ANF*A1
            B1 = X*B0 + ANF*B1
            IF (A1.NE.ZERO) THEN
              FAC = ONE/A1
              U = B1*FAC
              IF (ABS(U-GOLD).LT.ABS(U)*RELERR) RETURN
              GOLD = U
            END IF
   30     CONTINUE
        END IF
        RETURN
        END

* **********************************************************************
*
        SUBROUTINE ETARIE(S, EPS, XMAX, RELERR, ETA)
*
* **********************************************************************
* ETARIE returns in ETA the value of the eta function, for real argument
*        S, approximated with a relative error RELERR.  EPS and XMAX
*        represent the smallest positive floating-point number such that
*        1.0+EPS .NE. 1.0, and the largest finite floating-point number,
*        respectively.  For S .GT. -1 Levin's u transform [2] is used to
*        sum the alternating series (23.2.19) of [1], except when S is a
*        positive integer.  Otherwise the reflection formula (23.2.6) of
*        [1] is employed, involving gamma function evaluation, except in
*        the trivial zeros S = -2N.
*
* References:
*
*   [1] E. V. Haynsworth and K. Goldberg, "Bernoulli and Euler
*       Polynomials - Riemann Zeta Function", in "Handbook of
*       Mathematical Functions with Formulas, Graphs and Mathematical
*       Tables" (M. Abramowitz and I. A. Stegun, eds.), no. 55 in
*       National Bureau of Standards Applied Mathematics Series, ch. 23,
*       pp. 803-819, Washington, D.C.:  U.S. Government Printing Office,
*       1964.
*
*   [2] T. Fessler, W. F. Ford, D. A. Smith, "ALGORITHM 602. HURRY: An
*       acceleration algorithm for scalar sequences and series", ACM
*       Transactions on Mathematical Software, vol. 9, no. 3,
*       pp. 355-357, September  1983.
*
* If a single precision version is desired, change all occurrences of
* *SP in columns 1 to 3 to blanks and comment the corresponding double
* precision statements.
*
* Michele Goano, Politecnico di Torino  (goano@polito.it).
* Latest revision:  March 22, 1994.
* **********************************************************************
*   Parameters
*SP     REAL             ONE, PI, PILOG, TWO, ZERO
        DOUBLE PRECISION ONE, PI, PILOG, TWO, ZERO
*SP     PARAMETER (ONE = 1.0E+0, PI = 3.141592653589793238462643E+0,
*SP  &             PILOG = 1.144729885849400174143427E+0, TWO = 2.0E+0,
*SP  &             ZERO = 0.0E+0)
        PARAMETER (ONE = 1.0D+0, PI = 3.141592653589793238462643D+0,
     &             PILOG = 1.144729885849400174143427D+0, TWO = 2.0D+0,
     &             ZERO = 0.0D+0)
*   Scalar arguments
*SP     REAL             S, EPS, XMAX, RELERR, ETA
        DOUBLE PRECISION S, EPS, XMAX, RELERR, ETA
*   Local scalars
        LOGICAL LOGGAM
        INTEGER IERR
*SP     REAL             ETALOG, GAMMA, TWOTOS, XBIG
        DOUBLE PRECISION ETALOG, GAMMA, TWOTOS, XBIG
*   External subroutines
        EXTERNAL ETALEV, ETAN, GAMMAC
*   Intrinsic functions
        INTRINSIC ANINT, EXP, LOG, MOD, SIN
* ----------------------------------------------------------------------
        XBIG = LOG(XMAX)
        ETA = ZERO
*
        IF (S.EQ.ZERO) THEN
          ETA = ONE/TWO
        ELSE IF (S.LT.ZERO .AND. MOD(S, TWO).EQ.ZERO) THEN
          ETA = ZERO
        ELSE IF (S.GT.-ONE) THEN
          IF (ABS(S-ANINT(S)).LE.ABS(S)*EPS) THEN
            CALL ETAN(NINT(S), ETA)
          ELSE
            CALL ETALEV(S, RELERR, ETA)
          END IF
        ELSE
          TWOTOS = TWO**S
          CALL GAMMAC(ONE - S, EPS, XMAX, GAMMA, IERR)
          IF (IERR.EQ.-1) LOGGAM = .TRUE.
          CALL ETALEV(ONE - S, RELERR, ETA)
          IF (LOGGAM) THEN
            ETALOG = (S - ONE)*PILOG+GAMMA+LOG(ETA)
            ETA = (TWOTOS - TWO)/(ONE - TWOTOS)*SIN(S*PI/TWO)*
     &                                                       EXP(ETALOG)
          ELSE
            ETA = (TWOTOS - TWO)/(ONE - TWOTOS)*SIN(S*PI/TWO)*
     &                                           PI**(S - ONE)*GAMMA*ETA
          END IF
        END IF
        RETURN
        END

* **********************************************************************
*
        SUBROUTINE ETALEV(S, RELERR, ETA)
*
* **********************************************************************
* ETALEV returns in ETA the value of the eta function, for real argument
*        S, approximated with a relative error RELERR.  Levin's u
*        transform [2] is used to sum the alternating series (23.2.19)
*        of [1].
*
* References:
*
*   [1] E. V. Haynsworth and K. Goldberg, "Bernoulli and Euler
*       Polynomials - Riemann Zeta Function", in "Handbook of
*       Mathematical Functions with Formulas, Graphs and Mathematical
*       Tables" (M. Abramowitz and I. A. Stegun, eds.), no. 55 in
*       National Bureau of Standards Applied Mathematics Series, ch. 23,
*       pp. 803-819, Washington, D.C.:  U.S. Government Printing Office,
*       1964.
*
*   [2] T. Fessler, W. F. Ford, D. A. Smith, "ALGORITHM 602. HURRY: An
*       acceleration algorithm for scalar sequences and series", ACM
*       Transactions on Mathematical Software, vol. 9, no. 3,
*       pp. 355-357, September  1983.
*
* If a single precision version is desired, change all occurrences of
* *SP in columns 1 to 3 to blanks and comment the corresponding double
* precision statements.
*
* Michele Goano, Politecnico di Torino  (goano@polito.it).
* Latest revision:  March 22, 1994.
* **********************************************************************
*   Parameters
        INTEGER ITMAX
        PARAMETER (ITMAX = 100)
*SP     REAL             ONE, ZERO
        DOUBLE PRECISION ONE, ZERO
*SP     PARAMETER (ONE = 1.0E+0, ZERO = 0.0E+0)
        PARAMETER (ONE = 1.0D+0, ZERO = 0.0D+0)
*   Scalar arguments
*SP     REAL             S, RELERR, ETA
        DOUBLE PRECISION S, RELERR, ETA
*   Local scalars
        INTEGER JTERM
*SP     REAL             ETAOLD, SEGN, SUM, TERM
        DOUBLE PRECISION ETAOLD, SEGN, SUM, TERM
*   Local arrays
*SP     REAL             QNUM(ITMAX), QDEN(ITMAX)
        DOUBLE PRECISION QNUM(ITMAX), QDEN(ITMAX)
*   External subroutines
        EXTERNAL WHIZ
*   Intrinsic functions
*SP     INTRINSIC ABS, REAL
        INTRINSIC ABS, DBLE
* ----------------------------------------------------------------------
        ETA = ZERO
*
        SEGN = ONE
        DO 10 JTERM = 1, ITMAX
          ETAOLD = ETA
*SP       TERM = SEGN/REAL(JTERM)**S
          TERM = SEGN/DBLE(JTERM)**S
          CALL WHIZ(TERM, JTERM, QNUM, QDEN, ETA, SUM)
*   Check truncation error and convergence
          IF (ABS(ETA-ETAOLD).LE.ABS(ETA)*RELERR) RETURN
          SEGN = -SEGN
   10   CONTINUE
        END

* **********************************************************************
*
        SUBROUTINE ETAN(N, ETA)
*
* **********************************************************************
* ETAN returns in ETA the value of the eta function for integer
*      nonnegative argument N, approximated to 25 significant decimal
*      digits.
*
* Reference:
*
*   E. V. Haynsworth and K. Goldberg, "Bernoulli and Euler Polynomials -
*   Riemann Zeta Function", in "Handbook of Mathematical Functions with
*   Formulas, Graphs and Mathematical Tables" (M. Abramowitz and
*   I. A. Stegun, eds.), no. 55 in National Bureau of Standards Applied
*   Mathematics Series, ch. 23, pp. 803-819, Washington, D.C.:  U.S.
*   Government Printing Office, 1964.
*
* If a single precision version is desired, change all occurrences of
* *SP in columns 1 to 3 to blanks and comment the corresponding double
* precision statements.
*
* Michele Goano, Politecnico di Torino  (goano@polito.it).
* Latest revision:  March 22, 1994.
* **********************************************************************
*   Parameters
*SP     REAL             HALF, ONE, ZERO
        DOUBLE PRECISION HALF, ONE, ZERO
*SP     PARAMETER (HALF = 0.5E+0, ONE = 1.0E+0, ZERO = 0.0E+0)
        PARAMETER (HALF = 0.5D+0, ONE = 1.0D+0, ZERO = 0.0D+0)
*   Scalar arguments
        INTEGER N
*SP     REAL             ETA
        DOUBLE PRECISION ETA
*   Local arrays
*SP     REAL             ETABLE(84)
        DOUBLE PRECISION ETABLE(84)
* ----------------------------------------------------------------------
        SAVE ETABLE
*SP     DATA ETABLE(1),  ETABLE(2),  ETABLE(3),  ETABLE(4),
*SP  &       ETABLE(5),  ETABLE(6),  ETABLE(7),  ETABLE(8),
*SP  &       ETABLE(9),  ETABLE(10), ETABLE(11), ETABLE(12),
*SP  &       ETABLE(13), ETABLE(14), ETABLE(15), ETABLE(16) /
*SP  &  0.6931471805599453094172321E+0, 0.8224670334241132182362076E+0,
*SP  &  0.9015426773696957140498036E+0, 0.9470328294972459175765032E+0,
*SP  &  0.9721197704469093059356551E+0, 0.9855510912974351040984392E+0,
*SP  &  0.9925938199228302826704257E+0, 0.9962330018526478992272893E+0,
*SP  &  0.9980942975416053307677830E+0, 0.9990395075982715656392218E+0,
*SP  &  0.9995171434980607541440942E+0, 0.9997576851438581908531797E+0,
*SP  &  0.9998785427632651154921750E+0, 0.9999391703459797181709542E+0,
*SP  &  0.9999695512130992380826329E+0, 0.9999847642149061064416828E+0 /
*SP     DATA ETABLE(17), ETABLE(18), ETABLE(19), ETABLE(20),
*SP  &       ETABLE(21), ETABLE(22), ETABLE(23), ETABLE(24),
*SP  &       ETABLE(25), ETABLE(26), ETABLE(27), ETABLE(28),
*SP  &       ETABLE(29), ETABLE(30), ETABLE(31), ETABLE(32) /
*SP  &  0.9999923782920410119769379E+0, 0.9999961878696101134796892E+0,
*SP  &  0.9999980935081716751068565E+0, 0.9999990466115815221150508E+0,
*SP  &  0.9999995232582155428163167E+0, 0.9999997616132308225478972E+0,
*SP  &  0.9999998808013184395032238E+0, 0.9999999403988923946283614E+0,
*SP  &  0.9999999701988569628344151E+0, 0.9999999850992319965687877E+0,
*SP  &  0.9999999925495504849635159E+0, 0.9999999962747534001087275E+0,
*SP  &  0.9999999981373694181121867E+0, 0.9999999990686822814539786E+0,
*SP  &  0.9999999995343403314542175E+0, 0.9999999997671698959514908E+0 /
*SP     DATA ETABLE(33), ETABLE(34), ETABLE(35), ETABLE(36),
*SP  &       ETABLE(37), ETABLE(38), ETABLE(39), ETABLE(40),
*SP  &       ETABLE(41), ETABLE(42), ETABLE(43), ETABLE(44),
*SP  &       ETABLE(45), ETABLE(46), ETABLE(47), ETABLE(48) /
*SP  &  0.9999999998835848580460305E+0, 0.9999999999417923990453159E+0,
*SP  &  0.9999999999708961895298095E+0, 0.9999999999854480914338848E+0,
*SP  &  0.9999999999927240446065848E+0, 0.9999999999963620219331688E+0,
*SP  &  0.9999999999981810108432087E+0, 0.9999999999990905053804789E+0,
*SP  &  0.9999999999995452526765309E+0, 0.9999999999997726263336959E+0,
*SP  &  0.9999999999998863131653248E+0, 0.9999999999999431565821547E+0,
*SP  &  0.9999999999999715782909081E+0, 0.9999999999999857891453976E+0,
*SP  &  0.9999999999999928945726800E+0, 0.9999999999999964472863337E+0 /
*SP     DATA ETABLE(49), ETABLE(50), ETABLE(51), ETABLE(52),
*SP  &       ETABLE(53), ETABLE(54), ETABLE(55), ETABLE(56),
*SP  &       ETABLE(57), ETABLE(58), ETABLE(59), ETABLE(60),
*SP  &       ETABLE(61), ETABLE(62), ETABLE(63), ETABLE(64) /
*SP  &  0.9999999999999982236431648E+0, 0.9999999999999991118215817E+0,
*SP  &  0.9999999999999995559107906E+0, 0.9999999999999997779553952E+0,
*SP  &  0.9999999999999998889776976E+0, 0.9999999999999999444888488E+0,
*SP  &  0.9999999999999999722444244E+0, 0.9999999999999999861222122E+0,
*SP  &  0.9999999999999999930611061E+0, 0.9999999999999999965305530E+0,
*SP  &  0.9999999999999999982652765E+0, 0.9999999999999999991326383E+0,
*SP  &  0.9999999999999999995663191E+0, 0.9999999999999999997831596E+0,
*SP  &  0.9999999999999999998915798E+0, 0.9999999999999999999457899E+0 /
*SP     DATA ETABLE(65), ETABLE(66), ETABLE(67), ETABLE(68),
*SP  &       ETABLE(69), ETABLE(70), ETABLE(71), ETABLE(72),
*SP  &       ETABLE(73), ETABLE(74), ETABLE(75), ETABLE(76),
*SP  &       ETABLE(77), ETABLE(78), ETABLE(79), ETABLE(80) /
*SP  &  0.9999999999999999999728949E+0, 0.9999999999999999999864475E+0,
*SP  &  0.9999999999999999999932237E+0, 0.9999999999999999999966119E+0,
*SP  &  0.9999999999999999999983059E+0, 0.9999999999999999999991530E+0,
*SP  &  0.9999999999999999999995765E+0, 0.9999999999999999999997882E+0,
*SP  &  0.9999999999999999999998941E+0, 0.9999999999999999999999471E+0,
*SP  &  0.9999999999999999999999735E+0, 0.9999999999999999999999868E+0,
*SP  &  0.9999999999999999999999934E+0, 0.9999999999999999999999967E+0,
*SP  &  0.9999999999999999999999983E+0, 0.9999999999999999999999992E+0 /
*SP     DATA ETABLE(81), ETABLE(82), ETABLE(83), ETABLE(84) /
*SP  &  0.9999999999999999999999996E+0, 0.9999999999999999999999998E+0,
*SP  &  0.9999999999999999999999999E+0, 0.9999999999999999999999999E+0 /
        DATA ETABLE(1),  ETABLE(2),  ETABLE(3),  ETABLE(4),
     &       ETABLE(5),  ETABLE(6),  ETABLE(7),  ETABLE(8),
     &       ETABLE(9),  ETABLE(10), ETABLE(11), ETABLE(12),
     &       ETABLE(13), ETABLE(14), ETABLE(15), ETABLE(16) /
     &  0.6931471805599453094172321D+0, 0.8224670334241132182362076D+0,
     &  0.9015426773696957140498036D+0, 0.9470328294972459175765032D+0,
     &  0.9721197704469093059356551D+0, 0.9855510912974351040984392D+0,
     &  0.9925938199228302826704257D+0, 0.9962330018526478992272893D+0,
     &  0.9980942975416053307677830D+0, 0.9990395075982715656392218D+0,
     &  0.9995171434980607541440942D+0, 0.9997576851438581908531797D+0,
     &  0.9998785427632651154921750D+0, 0.9999391703459797181709542D+0,
     &  0.9999695512130992380826329D+0, 0.9999847642149061064416828D+0 /
        DATA ETABLE(17), ETABLE(18), ETABLE(19), ETABLE(20),
     &       ETABLE(21), ETABLE(22), ETABLE(23), ETABLE(24),
     &       ETABLE(25), ETABLE(26), ETABLE(27), ETABLE(28),
     &       ETABLE(29), ETABLE(30), ETABLE(31), ETABLE(32) /
     &  0.9999923782920410119769379D+0, 0.9999961878696101134796892D+0,
     &  0.9999980935081716751068565D+0, 0.9999990466115815221150508D+0,
     &  0.9999995232582155428163167D+0, 0.9999997616132308225478972D+0,
     &  0.9999998808013184395032238D+0, 0.9999999403988923946283614D+0,
     &  0.9999999701988569628344151D+0, 0.9999999850992319965687877D+0,
     &  0.9999999925495504849635159D+0, 0.9999999962747534001087275D+0,
     &  0.9999999981373694181121867D+0, 0.9999999990686822814539786D+0,
     &  0.9999999995343403314542175D+0, 0.9999999997671698959514908D+0 /
        DATA ETABLE(33), ETABLE(34), ETABLE(35), ETABLE(36),
     &       ETABLE(37), ETABLE(38), ETABLE(39), ETABLE(40),
     &       ETABLE(41), ETABLE(42), ETABLE(43), ETABLE(44),
     &       ETABLE(45), ETABLE(46), ETABLE(47), ETABLE(48) /
     &  0.9999999998835848580460305D+0, 0.9999999999417923990453159D+0,
     &  0.9999999999708961895298095D+0, 0.9999999999854480914338848D+0,
     &  0.9999999999927240446065848D+0, 0.9999999999963620219331688D+0,
     &  0.9999999999981810108432087D+0, 0.9999999999990905053804789D+0,
     &  0.9999999999995452526765309D+0, 0.9999999999997726263336959D+0,
     &  0.9999999999998863131653248D+0, 0.9999999999999431565821547D+0,
     &  0.9999999999999715782909081D+0, 0.9999999999999857891453976D+0,
     &  0.9999999999999928945726800D+0, 0.9999999999999964472863337D+0 /
        DATA ETABLE(49), ETABLE(50), ETABLE(51), ETABLE(52),
     &       ETABLE(53), ETABLE(54), ETABLE(55), ETABLE(56),
     &       ETABLE(57), ETABLE(58), ETABLE(59), ETABLE(60),
     &       ETABLE(61), ETABLE(62), ETABLE(63), ETABLE(64) /
     &  0.9999999999999982236431648D+0, 0.9999999999999991118215817D+0,
     &  0.9999999999999995559107906D+0, 0.9999999999999997779553952D+0,
     &  0.9999999999999998889776976D+0, 0.9999999999999999444888488D+0,
     &  0.9999999999999999722444244D+0, 0.9999999999999999861222122D+0,
     &  0.9999999999999999930611061D+0, 0.9999999999999999965305530D+0,
     &  0.9999999999999999982652765D+0, 0.9999999999999999991326383D+0,
     &  0.9999999999999999995663191D+0, 0.9999999999999999997831596D+0,
     &  0.9999999999999999998915798D+0, 0.9999999999999999999457899D+0 /
        DATA ETABLE(65), ETABLE(66), ETABLE(67), ETABLE(68),
     &       ETABLE(69), ETABLE(70), ETABLE(71), ETABLE(72),
     &       ETABLE(73), ETABLE(74), ETABLE(75), ETABLE(76),
     &       ETABLE(77), ETABLE(78), ETABLE(79), ETABLE(80) /
     &  0.9999999999999999999728949D+0, 0.9999999999999999999864475D+0,
     &  0.9999999999999999999932237D+0, 0.9999999999999999999966119D+0,
     &  0.9999999999999999999983059D+0, 0.9999999999999999999991530D+0,
     &  0.9999999999999999999995765D+0, 0.9999999999999999999997882D+0,
     &  0.9999999999999999999998941D+0, 0.9999999999999999999999471D+0,
     &  0.9999999999999999999999735D+0, 0.9999999999999999999999868D+0,
     &  0.9999999999999999999999934D+0, 0.9999999999999999999999967D+0,
     &  0.9999999999999999999999983D+0, 0.9999999999999999999999992D+0 /
        DATA ETABLE(81), ETABLE(82), ETABLE(83), ETABLE(84) /
     &  0.9999999999999999999999996D+0, 0.9999999999999999999999998D+0,
     &  0.9999999999999999999999999D+0, 0.9999999999999999999999999D+0 /
* ----------------------------------------------------------------------
        ETA = ZERO
        IF (N.EQ.0) THEN
          ETA = HALF
        ELSE IF (N.LE.84) THEN
          ETA = ETABLE(N)
        ELSE IF (N.GT.84) THEN
          ETA = ONE
        END IF
        RETURN
        END

* **********************************************************************
*
        SUBROUTINE FERERR(ERRMSG)
*
* **********************************************************************
* FERERR prints on the standard output unit an explanatory message of
*        the error condition occured in the package which approximates
*        the complete and incomplete Fermi-Dirac integral.
*
* Michele Goano, Politecnico di Torino  (goano@polito.it).
* Latest revision:  March 22, 1994.
* **********************************************************************
*   Scalar arguments
        CHARACTER*(*) ERRMSG
* ----------------------------------------------------------------------
        WRITE (*, FMT = 99999) ERRMSG
*   If you want to interrupt the execution after an error has occurred,
*   replace the RETURN statement with a STOP
        RETURN
99999   FORMAT (A)
        END
*
* **********************************************************************
*
        SUBROUTINE GAMMAC(X, EPS, XINF, GAMMA, IERR)
C-----------------------------------------------------------------------
C This routine calculates the gamma function for a real argument X.  The
C logarithm of the gamma function is computed, and the error flag IERR
C is set to -1, whenever the result would be too large to be represented
C on the floating-point arithmetic system.  Computation is based on an
C algorithm outlined in W. J. Cody, 'An overview of software development
C for special functions', Lecture Notes in Mathematics, 506, Numerical
C Analysis Dundee, 1975, G. A. Watson (ed.), Springer Verlag, Berlin,
C 1976.  The program uses rational functions that approximate the gamma
C function to at least 20 significant decimal digits.  Coefficients for
C the approximation over the interval (1,2) are unpublished.  Those for
C the approximation for X .GE. 12 are from Hart et al., Computer
C Approximations, Wiley and Sons, New York, 1968.
C
C If a single precision version is desired, change all occurrences of CS
C in columns 1 and 2 to blanks and comment the corresponding double
C precision statements.
C
C Explanation of machine-dependent variables
C
C EPS    - the smallest positive floating-point number such that
C          1.0 + EPS .GT. 1.0
C XINF   - the largest machine representable floating-point number.
C XBIG   - the largest floating-point number such that EXP(XBIG) is
C          machine representable.
C
C Error returns
C
C  The program returns LOG(GAMMA) and sets IERR = -1 when overflow would
C  occur.
C
C Author: W. J. Cody
C         Argonne National Laboratory
C
C Revised by M. Goano, Politecnico di Torino, to take advantage of
C Fortran 77 control structures.
C
C Latest modification of the original version: May 18, 1982
C                     of the revised version:  March 21, 1994
C-----------------------------------------------------------------------
        INTEGER I, IERR, J, N
CS      REAL             C, EPS, FACT, GAMMA, HALF, ONE, P, PI, Q, RES,
CS   &                   SQRTPI, SUM, TWELVE, X, XBIG, XDEN, XINF,
CS   &                   XNUM, Y, Y1, YSQ, Z, ZERO
        DOUBLE PRECISION C, EPS, FACT, GAMMA, HALF, ONE, P, PI, Q, RES,
     &                   SQRTPI, SUM, TWELVE, X, XBIG, XDEN, XINF,
     &                   XNUM, Y, Y1, YSQ, Z, ZERO
        LOGICAL PARITY
        DIMENSION C(7), P(8), Q(8)
CS      INTRINSIC ALOG, EXP, FLOAT, IFIX, SIN
        INTRINSIC DBLE, DEXP, DLOG, DSIN, FLOAT, IFIX, SNGL
C-----------------------------------------------------------------------
C Mathematical constants
C-----------------------------------------------------------------------
CS      PARAMETER (ONE = 1.0E+0, HALF = 0.5E+0, TWELVE = 12.0E+0,
CS   &             ZERO = 0.0E+0, PI = 3.1415926535897932384626434E+0,
CS   &             SQRTPI = 0.9189385332046727417803297E+0)
        PARAMETER (ONE = 1.0D+0, HALF = 0.5D+0, TWELVE = 12.0D+0,
     &             ZERO = 0.0D+0, PI = 3.1415926535897932384626434D+0,
     &             SQRTPI = 0.9189385332046727417803297D+0)
C-----------------------------------------------------------------------
C SAVE declaration for the arrays of the coefficients
C-----------------------------------------------------------------------
        SAVE C, P, Q
C-----------------------------------------------------------------------
C Numerator and denominator coefficients for rational minimax
C approximation over (1,2)
C-----------------------------------------------------------------------
CS      DATA P /-1.71618513886549492533811E+0,
CS   &           2.47656508055759199108314E+1,
CS   &          -3.79804256470945635097577E+2,
CS   &           6.29331155312818442661052E+2,
CS   &           8.66966202790413211295064E+2,
CS   &          -3.14512729688483675254357E+4,
CS   &          -3.61444134186911729807069E+4,
CS   &           6.64561438202405440627855E+4/
        DATA P /-1.71618513886549492533811D+0,
     &           2.47656508055759199108314D+1,
     &          -3.79804256470945635097577D+2,
     &           6.29331155312818442661052D+2,
     &           8.66966202790413211295064D+2,
     &          -3.14512729688483675254357D+4,
     &          -3.61444134186911729807069D+4,
     &           6.64561438202405440627855D+4/
CS      DATA Q /-3.08402300119738975254353E+1,
CS   &           3.15350626979604161529144E+2,
CS   &          -1.01515636749021914166146E+3,
CS   &          -3.10777167157231109440444E+3,
CS   &           2.25381184209801510330112E+4,
CS   &           4.75584627752788110767815E+3,
CS   &          -1.34659959864969306392456E+5,
CS   &          -1.15132259675553483497211E+5/
        DATA Q /-3.08402300119738975254353D+1,
     &           3.15350626979604161529144D+2,
     &          -1.01515636749021914166146D+3,
     &          -3.10777167157231109440444D+3,
     &           2.25381184209801510330112D+4,
     &           4.75584627752788110767815D+3,
     &          -1.34659959864969306392456D+5,
     &          -1.15132259675553483497211D+5/
C-----------------------------------------------------------------------
C Coefficients for minimax approximation over (12, INF)
C-----------------------------------------------------------------------
CS      DATA C /-1.910444077728E-03,
CS   &           8.4171387781295E-04,
CS   &          -5.952379913043012E-04,
CS   &           7.93650793500350248E-04,
CS   &          -2.777777777777681622553E-03,
CS   &           8.333333333333333331554247E-02,
CS   &           5.7083835261E-03/
        DATA C /-1.910444077728D-03,
     &           8.4171387781295D-04,
     &          -5.952379913043012D-04,
     &           7.93650793500350248D-04,
     &          -2.777777777777681622553D-03,
     &           8.333333333333333331554247D-02,
     &           5.7083835261D-03/
C-----------------------------------------------------------------------
C Machine dependent local variables
C-----------------------------------------------------------------------
CS      XBIG = ALOG(XINF)
        XBIG = DLOG(XINF)
C-----------------------------------------------------------------------
        IERR = 0
        PARITY = .FALSE.
        FACT = ONE
        N = 0
        Y = X
        IF (Y.LE.ZERO) THEN
C-----------------------------------------------------------------------
C Argument is negative
C-----------------------------------------------------------------------
          Y = -X
CS        J = IFIX(Y)
          J = IFIX(SNGL(Y))
CS        RES = Y - FLOAT(J)
          RES = Y - DBLE(FLOAT(J))
          IF (J.NE.(J/2)*2) PARITY = .TRUE.
CS        FACT = -PI/SIN(PI*RES)
          FACT = -PI/DSIN(PI*RES)
          Y = Y + ONE
        END IF
C-----------------------------------------------------------------------
C Argument is positive
C-----------------------------------------------------------------------
        IF (Y.LT.EPS) THEN
C-----------------------------------------------------------------------
C Argument .LT. EPS
C-----------------------------------------------------------------------
          RES = ONE/Y
        ELSE IF (Y.GE.TWELVE) THEN
C-----------------------------------------------------------------------
C Evaluate for argument .GE. 12.0
C-----------------------------------------------------------------------
          YSQ = Y*Y
          SUM = C(7)
          DO 10 I = 1, 6
            SUM = SUM/YSQ + C(I)
   10     CONTINUE
CS        SUM = SUM/Y + (Y - HALF)*ALOG(Y) - Y + SQRTPI
          SUM = SUM/Y + (Y - HALF)*DLOG(Y) - Y + SQRTPI
          IF (SUM.GT.XBIG) THEN
C-----------------------------------------------------------------------
C Return the logarithm to avoid overflow
C-----------------------------------------------------------------------
            RES = SUM
            IERR = -1
          ELSE
CS          RES = EXP(SUM)
            RES = DEXP(SUM)
          END IF
        ELSE
          Y1 = Y
          IF (Y.GE.ONE) THEN
C-----------------------------------------------------------------------
C 1.0 .LT. argument .LT. 12.0, reduce argument if necessary
C-----------------------------------------------------------------------
CS          N = IFIX(Y) - 1
            N = IFIX(SNGL(Y)) - 1
CS          Y = Y - FLOAT(N)
            Y = Y - DBLE(FLOAT(N))
            Z = Y - ONE
          ELSE
C-----------------------------------------------------------------------
C 0.0 .LT. argument .LT. 1.0
C-----------------------------------------------------------------------
            Z = Y
            Y = Y + ONE
          END IF
C-----------------------------------------------------------------------
C Evaluate approximation for 1.0 .LT. argument .LT. 2.0
C-----------------------------------------------------------------------
          XNUM = ZERO
          XDEN = ONE
          DO 20 I = 1, 8
            XNUM = (XNUM + P(I))*Z
            XDEN = XDEN*Z + Q(I)
   20     CONTINUE
          RES = XNUM/XDEN + ONE
          IF (Y.NE.Y1) THEN
            IF (Y1.GT.Y) THEN
C-----------------------------------------------------------------------
C Adjust result for case  2.0 .LT. argument .LT. 12.0
C-----------------------------------------------------------------------
              DO 30 I = 1, N
                RES = RES*Y
                Y = Y + ONE
   30         CONTINUE
            ELSE
C-----------------------------------------------------------------------
C Adjust result for case  0.0 .LT. argument .LT. 1.0
C-----------------------------------------------------------------------
              RES = RES/Y1
            END IF
          END IF
        END IF
C-----------------------------------------------------------------------
C Final adjustments and return
C-----------------------------------------------------------------------
        IF (PARITY) RES = -RES
        IF (FACT.NE.ONE) RES = FACT/RES
        GAMMA = RES
   40   CONTINUE
        RETURN
        END
*
* **********************************************************************
*
        SUBROUTINE WHIZ(TERM, ITERM, QNUM, QDEN, RESULT, S)
************************************************************************
*     ALGORITHM 602, COLLECTED ALGORITHMS FROM ACM.
*     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.9, NO. 3,
*     SEP., 1983, P. 355-357.
*
* The u algorithm for accelerating a series.
*
* Arguments:
*    TERM   = last element of series
*    ITERM   = order of TERM in the series = number of calls to WHIZ
*    QNUM   = backward diagonal of numerator array, at least N long
*    QDEN   = backward diagonal of denominator array, at least N long
*    RESULT = accelerated value of the sum
*    S      = simple sum of the series
*
* Inputs:  TERM, ITERM
*
* Outputs:  RESULT, S
*
* If a single precision version is desired, change all occurrences of
* *SP in columns 1 to 3 to blanks and comment the corresponding double
* precision statements.
*
* Revised by M. Goano, Politecnico di Torino.
* Latest modification of the revised version: April 12, 1993
************************************************************************
*   Parameters
*SP     REAL             ONE, ZERO
        DOUBLE PRECISION ONE, ZERO
*SP     PARAMETER (ONE = 1.0E+0, ZERO = 0.0E+0)
        PARAMETER (ONE = 1.0D+0, ZERO = 0.0D+0)
*   Scalar arguments
        INTEGER ITERM
*SP     REAL             RESULT, S, TERM
        DOUBLE PRECISION RESULT, S, TERM
*   Array arguments
*SP     REAL             QNUM(*), QDEN(*)
        DOUBLE PRECISION QNUM(*), QDEN(*)
*   Local scalars
        INTEGER J, K, L
*SP     REAL             C, FACTOR, FJ, FL, FTERM, RATIO
        DOUBLE PRECISION C, FACTOR, FJ, FL, FTERM, RATIO
*   Intrinsic functions
*SP     INTRINSIC REAL
        INTRINSIC DBLE
* ----------------------------------------------------------------------
        IF (ITERM.EQ.1) S = ZERO
* Get ITERM diagonal
        S = TERM + S
        L = ITERM - 1
*SP     FTERM = REAL(ITERM)
        FTERM = DBLE(ITERM)
        QDEN(ITERM) = ONE/(TERM*FTERM**2)
        QNUM(ITERM) = S*QDEN(ITERM)
        IF (ITERM.GT.1) THEN
          FACTOR = ONE
*SP       FL = REAL(L)
          FL = DBLE(L)
          RATIO = FL/FTERM
          DO 10 K = 1, L
            J = ITERM - K
*SP         FJ = REAL(J)
            FJ = DBLE(J)
            C = FACTOR*FJ/FTERM
            FACTOR = FACTOR*RATIO
            QDEN(J) = QDEN(J + 1) - C*QDEN(J)
            QNUM(J) = QNUM(J + 1) - C*QNUM(J)
   10     CONTINUE
        END IF
        RESULT = QNUM(1)/QDEN(1)
        RETURN
        END
*
* **********************************************************************
*
      SUBROUTINE MACHAR(IBETA,IT,IRND,NGRD,MACHEP,NEGEP,IEXP,MINEXP,
     1                   MAXEXP,EPS,EPSNEG,XMIN,XMAX)
C----------------------------------------------------------------------
C  This Fortran 77 subroutine is intended to determine the parameters
C   of the floating-point arithmetic system specified below.  The
C   determination of the first three uses an extension of an algorithm
C   due to M. Malcolm, CACM 15 (1972), pp. 949-951, incorporating some,
C   but not all, of the improvements suggested by M. Gentleman and S.
C   Marovich, CACM 17 (1974), pp. 276-277.  An earlier version of this
C   program was published in the book Software Manual for the
C   Elementary Functions by W. J. Cody and W. Waite, Prentice-Hall,
C   Englewood Cliffs, NJ, 1980.
C
C  If a single precision version is desired, change all occurrences of
C   CS in columns 1 and 2 to blanks and comment the corresponding double
C   precision statements.
C
C  Parameter values reported are as follows:
C
C       IBETA   - the radix for the floating-point representation
C       IT      - the number of base IBETA digits in the floating-point
C                 significand
C       IRND    - 0 if floating-point addition chops
C                 1 if floating-point addition rounds, but not in the
C                   IEEE style
C                 2 if floating-point addition rounds in the IEEE style
C                 3 if floating-point addition chops, and there is
C                   partial underflow
C                 4 if floating-point addition rounds, but not in the
C                   IEEE style, and there is partial underflow
C                 5 if floating-point addition rounds in the IEEE style,
C                   and there is partial underflow
C       NGRD    - the number of guard digits for multiplication with
C                 truncating arithmetic.  It is
C                 0 if floating-point arithmetic rounds, or if it
C                   truncates and only  IT  base  IBETA digits
C                   participate in the post-normalization shift of the
C                   floating-point significand in multiplication;
C                 1 if floating-point arithmetic truncates and more
C                   than  IT  base  IBETA  digits participate in the
C                   post-normalization shift of the floating-point
C                   significand in multiplication.
C       MACHEP  - the largest negative integer such that
C                 1.0+FLOAT(IBETA)**MACHEP .NE. 1.0, except that
C                 MACHEP is bounded below by  -(IT+3)
C       NEGEPS  - the largest negative integer such that
C                 1.0-FLOAT(IBETA)**NEGEPS .NE. 1.0, except that
C                 NEGEPS is bounded below by  -(IT+3)
C       IEXP    - the number of bits (decimal places if IBETA = 10)
C                 reserved for the representation of the exponent
C                 (including the bias or sign) of a floating-point
C                 number
C       MINEXP  - the largest in magnitude negative integer such that
C                 FLOAT(IBETA)**MINEXP is positive and normalized
C       MAXEXP  - the smallest positive power of  BETA  that overflows
C       EPS     - FLOAT(IBETA)**MACHEP.
C       EPSNEG  - FLOAT(IBETA)**NEGEPS.
C       XMIN    - the smallest non-vanishing normalized floating-point
C                 power of the radix, i.e.,  XMIN = FLOAT(IBETA)**MINEXP
C       XMAX    - the largest finite floating-point number.  In
C                 particular  XMAX = (1.0-EPSNEG)*FLOAT(IBETA)**MAXEXP
C                 Note - on some machines  XMAX  will be only the
C                 second, or perhaps third, largest number, being
C                 too small by 1 or 2 units in the last digit of
C                 the significand.
C
C  Latest modification: May 30, 1989
C
C  Author: W. J. Cody
C          Mathematics and Computer Science Division
C          Argonne National Laboratory
C          Argonne, IL 60439
C
C----------------------------------------------------------------------
      INTEGER I,IBETA,IEXP,IRND,IT,ITEMP,IZ,J,K,MACHEP,MAXEXP,
     1        MINEXP,MX,NEGEP,NGRD,NXRES
CS    REAL
      DOUBLE PRECISION
     1   A,B,BETA,BETAIN,BETAH,CONV,EPS,EPSNEG,ONE,T,TEMP,TEMPA,
     2   TEMP1,TWO,XMAX,XMIN,Y,Z,ZERO
C----------------------------------------------------------------------
CS    CONV(I) = REAL(I)
      CONV(I) = DBLE(I)
      ONE = CONV(1)
      TWO = ONE + ONE
      ZERO = ONE - ONE
C----------------------------------------------------------------------
C  Determine IBETA, BETA ala Malcolm.
C----------------------------------------------------------------------
      A = ONE
   10 A = A + A
         TEMP = A+ONE
         TEMP1 = TEMP-A
         IF (TEMP1-ONE .EQ. ZERO) GO TO 10
      B = ONE
   20 B = B + B
         TEMP = A+B
         ITEMP = INT(TEMP-A)
         IF (ITEMP .EQ. 0) GO TO 20
      IBETA = ITEMP
      BETA = CONV(IBETA)
C----------------------------------------------------------------------
C  Determine IT, IRND.
C----------------------------------------------------------------------
      IT = 0
      B = ONE
  100 IT = IT + 1
         B = B * BETA
         TEMP = B+ONE
         TEMP1 = TEMP-B
         IF (TEMP1-ONE .EQ. ZERO) GO TO 100
      IRND = 0
      BETAH = BETA / TWO
      TEMP = A+BETAH
      IF (TEMP-A .NE. ZERO) IRND = 1
      TEMPA = A + BETA
      TEMP = TEMPA+BETAH
      IF ((IRND .EQ. 0) .AND. (TEMP-TEMPA .NE. ZERO)) IRND = 2
C----------------------------------------------------------------------
C  Determine NEGEP, EPSNEG.
C----------------------------------------------------------------------
      NEGEP = IT + 3
      BETAIN = ONE / BETA
      A = ONE
      DO 200 I = 1, NEGEP
         A = A * BETAIN
  200 CONTINUE
      B = A
  210 TEMP = ONE-A
         IF (TEMP-ONE .NE. ZERO) GO TO 220
         A = A * BETA
         NEGEP = NEGEP - 1
      GO TO 210
  220 NEGEP = -NEGEP
      EPSNEG = A
C----------------------------------------------------------------------
C  Determine MACHEP, EPS.
C----------------------------------------------------------------------
      MACHEP = -IT - 3
      A = B
  300 TEMP = ONE+A
         IF (TEMP-ONE .NE. ZERO) GO TO 320
         A = A * BETA
         MACHEP = MACHEP + 1
      GO TO 300
  320 EPS = A
C----------------------------------------------------------------------
C  Determine NGRD.
C----------------------------------------------------------------------
      NGRD = 0
      TEMP = ONE+EPS
      IF ((IRND .EQ. 0) .AND. (TEMP*ONE-ONE .NE. ZERO)) NGRD = 1
C----------------------------------------------------------------------
C  Determine IEXP, MINEXP, XMIN.
C
C  Loop to determine largest I and K = 2**I such that
C         (1/BETA) ** (2**(I))
C  does not underflow.
C  Exit from loop is signaled by an underflow.
C----------------------------------------------------------------------
      I = 0
      K = 1
      Z = BETAIN
      T = ONE + EPS
      NXRES = 0
  400 Y = Z
         Z = Y * Y
C----------------------------------------------------------------------
C  Check for underflow here.
C----------------------------------------------------------------------
         A = Z * ONE
         TEMP = Z * T
         IF ((A+A .EQ. ZERO) .OR. (ABS(Z) .GE. Y)) GO TO 410
         TEMP1 = TEMP * BETAIN
         IF (TEMP1*BETA .EQ. Z) GO TO 410
         I = I + 1
         K = K + K
      GO TO 400
  410 IF (IBETA .EQ. 10) GO TO 420
      IEXP = I + 1
      MX = K + K
      GO TO 450
C----------------------------------------------------------------------
C  This segment is for decimal machines only.
C----------------------------------------------------------------------
  420 IEXP = 2
      IZ = IBETA
  430 IF (K .LT. IZ) GO TO 440
         IZ = IZ * IBETA
         IEXP = IEXP + 1
      GO TO 430
  440 MX = IZ + IZ - 1
C----------------------------------------------------------------------
C  Loop to determine MINEXP, XMIN.
C  Exit from loop is signaled by an underflow.
C----------------------------------------------------------------------
  450 XMIN = Y
         Y = Y * BETAIN
C----------------------------------------------------------------------
C  Check for underflow here.
C----------------------------------------------------------------------
         A = Y * ONE
         TEMP = Y * T
         IF (((A+A) .EQ. ZERO) .OR. (ABS(Y) .GE. XMIN)) GO TO 460
         K = K + 1
         TEMP1 = TEMP * BETAIN
         IF ((TEMP1*BETA .NE. Y) .OR. (TEMP .EQ. Y)) THEN
               GO TO 450
            ELSE
               NXRES = 3
               XMIN = Y
         END IF
  460 MINEXP = -K
C----------------------------------------------------------------------
C  Determine MAXEXP, XMAX.
C----------------------------------------------------------------------
      IF ((MX .GT. K+K-3) .OR. (IBETA .EQ. 10)) GO TO 500
      MX = MX + MX
      IEXP = IEXP + 1
  500 MAXEXP = MX + MINEXP
C----------------------------------------------------------------------
C  Adjust IRND to reflect partial underflow.
C----------------------------------------------------------------------
      IRND = IRND + NXRES
C----------------------------------------------------------------------
C  Adjust for IEEE-style machines.
C----------------------------------------------------------------------
      IF (IRND .GE. 2) MAXEXP = MAXEXP - 2
C----------------------------------------------------------------------
C  Adjust for machines with implicit leading bit in binary
C  significand, and machines with radix point at extreme
C  right of significand.
C----------------------------------------------------------------------
      I = MAXEXP + MINEXP
      IF ((IBETA .EQ. 2) .AND. (I .EQ. 0)) MAXEXP = MAXEXP - 1
      IF (I .GT. 20) MAXEXP = MAXEXP - 1
      IF (A .NE. Y) MAXEXP = MAXEXP - 2
      XMAX = ONE - EPSNEG
      IF (XMAX*ONE .NE. XMAX) XMAX = ONE - BETA * EPSNEG
      XMAX = XMAX / (BETA * BETA * BETA * XMIN)
      I = MAXEXP + MINEXP + 3
      IF (I .LE. 0) GO TO 520
      DO 510 J = 1, I
          IF (IBETA .EQ. 2) XMAX = XMAX + XMAX
          IF (IBETA .NE. 2) XMAX = XMAX * BETA
  510 CONTINUE
  520 RETURN
C---------- Last line of MACHAR ----------
      END
