C*****************************************************************************
C
C   MPFUN: A MULTIPLE PRECISION FLOATING POINT COMPUTATION PACKAGE
C
C   Standard Fortran-77 version
C   Version Date:  March 14, 1995
C
C   Author:
C
C      David H. Bailey                 Telephone:   415-604-4410
C      NASA Ames Research Center       Facsimile:   415-604-3957
C      Mail Stop T045-1                Internet:    dbailey@nas.nasa.gov
C      Moffett Field, CA 94035
C      USA
C
C   Restrictions:
C
C   This software has now been approved by NASA for unrestricted distribution.
C   However, usage of this software is subject to the following:
C
C   1. This software is offered without warranty of any kind, either expressed
C      or implied.  The author would appreciate, however, any reports of bugs
C      or other difficulties that may be encountered.
C   2. If modifications or enhancements to this software are made to this
C      software by others, NASA Ames reserves the right to obtain this enhanced
C      software at no cost and with no restrictions on its usage.
C   3. The author and NASA Ames are to be acknowledged in any published paper
C      based on computations using this software.  Accounts of practical
C      applications or other benefits resulting from this software are of
C      particular interest.  Please send a copy of such papers to the author.
C
C   Description:
C
C   The following information is a brief description of this program.  For
C   full details and instructions for usage, see the paper "A Portable High
C   Performance Multiprecision Package", available from the author.
C
C   This package of Fortran subroutines performs multiprecision floating point
C   arithmetic.  If sufficient main memory is available, the maximum precision
C   level is at least 16 million digits.  The maximum dynamic range is at
C   least 10^(+-14,000,000).  It employs advanced algorithms, including an
C   FFT-based multiplication routine and some recently discovered
C   quadratically convergent algorithms for pi, exp and log.  The package also
C   features extensive debug and self-checking facilities, so that it can be
C   used as a rigorous system integrity test.  All of the routines in this
C   package have been written to facilitate vector and parallel processing.
C
C   For users who do not wish to manually write code that calls these routines,
C   an automatic translator program is available from the author that converts
C   ordinary Fortran-77 code into code that calls these routines.  Contact the
C   author for details.
C
C   This package should run correctly on any computer with a Fortran-77
C   compiler that meets certain minimal floating point accuracy standards.
C   Any system based on the IEEE floating point standard, with a 25 bit
C   mantissa in single precision and a 53 bit mantissa in double precision,
C   easily meets these requirements.  All DEC VAX systems meet these
C   requirements.  All IBM mainframes and workstations meet these requirements.
C   Cray systems meet all of these requirements with double precision disabled
C   (i.e. by using only single precision).
C
C   Machine-specific tuning notes may be located by searching for the text
C   string C> in this program file.  It is highly recommended that these notes
C   be read before running this package on a specific system.  Also,
C   certain vectorizable DO loops that are often not recognized as such by
C   vectorizing compilers are prefaced with Cray CDIR$ IVDEP directives.  On
C   other vector systems these directives should be replaced by the
C   appropriate equivalents.
C
C   Instructions for compiling and testing this program are included in the
C   readme file that accompanies this file.
C
C*****************************************************************************
C
      BLOCK DATA
C
C   This initializes the parameters in MPCOM1 and the error codes in MPCOM2
C   with default values.
C>
C   On IEEE systems and most other 32 bit systems, set BBXC = 4096.D0,
C   NBTC = 24, NPRC = 32, and MCRC = 7.  On Cray systems, set BBXC = 2048.D0,
C   NBTC = 22, NPRC = 16, and MCRC = 8.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
C
      PARAMETER (BBXC = 4096.D0, NBTC = 24, NPRC = 32, MCRC = 7,        
     $  BDXC = BBXC ** 2, BX2C = BDXC ** 2, RBXC = 1.D0 / BBXC,         
     $  RDXC = RBXC ** 2, RX2C = RDXC ** 2, RXXC = 16.D0 * RX2C)
      DATA BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR / BBXC, BDXC,    
     $  BX2C, RBXC, RDXC, RX2C, RXXC, NBTC, NPRC/
      DATA NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS /                 
     $  16, 0, 6, 0, MCRC, 1, 1, 1, 1024/
      DATA KER /72 * 2/
      END
C
      SUBROUTINE DPADD (A, NA, B, NB, C, NC)
C
C   This adds the DPE numbers (A, NA) and (B, NB) to yield the sum (C, NC).
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION PT(64)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      SAVE PT
      DATA PT/ 64 * 0.D0/
C
      IF (IER .NE. 0) THEN
        C = 0.D0
        NC = 0
        RETURN
      ENDIF
C
C   If this is the first call to DPADD, initialize the PT table.
C
      IF (PT(1) .EQ. 0.D0) THEN
        PT(1) = 0.5D0
C
        DO 100 I = 2, 64
          PT(I) = 0.5D0 * PT(I-1)
 100    CONTINUE
C
      ENDIF
C
C   This operation reduces to five cases.
C
      IF (B .EQ. 0.D0) THEN
        C = A
        NC = NA
      ELSE IF (A .EQ. 0.D0) THEN
        C = B
        NC = NB
      ELSE IF (NA .EQ. NB) THEN
        C = A + B
        NC = NA
      ELSE IF (NA .GT. NB) THEN
        K = NA - NB
        NC = NA
        IF (K .GT. 64) THEN
          C = A
        ELSE
          C = A + B * PT(K)
        ENDIF
      ELSE
        K = NB - NA
        NC = NB
        IF (K .GT. 64) THEN
          C = B
        ELSE
          C = B + A * PT(K)
        ENDIF
      ENDIF
      IF (C .EQ. 0.D0) THEN
        NC = 0
        GOTO 130
      ENDIF
C
C   Normalize the result to a decent range if it is not.
C
 110  IF (ABS (C) .GE. BDX) THEN
        C = RDX * C
        NC = NC + NBT
        GOTO 110
      ENDIF
C
 120  IF (ABS (C) .LT. 1.D0) THEN
        C = BDX * C
        NC = NC - NBT
        GOTO 120
      ENDIF
C
 130  RETURN
      END
C
      SUBROUTINE DPDEC (A, NA, B, NB)
C
C   This converts the DPE number (A, NA) to decimal form, i.e. B * 10^NB,
C   where |B| is between 1 and 10.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (XLT = 0.3010299956639812D0)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
C
      IF (A .NE. 0.D0) THEN
        T1 = XLT * NA + LOG10 (ABS (A))
        NB = T1
        IF (T1 .LT. 0.D0) NB = NB - 1
        B = SIGN (10.D0 ** (T1 - NB), A)
      ELSE
        B = 0.D0
        NB = 0
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE DPDIV (A, NA, B, NB, C, NC)
C
C   This divides the DPE number (A, NA) by (B, NB) to yield the quotient
C   (C, NC).
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
C
      IF (IER .NE. 0) THEN
        C = 0.D0
        NC = 0
        RETURN
      ENDIF
      IF (B .EQ. 0.D0) THEN
        IF (KER(1) .NE. 0) THEN
          WRITE (LDB, 1)
 1        FORMAT ('*** DPDIV: Divisor is zero.')
          IER = 1
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
C   Divide A by B and subtract exponents, unless A is zero.
C
      IF (A .EQ. 0.D0) THEN
        C = 0.D0
        NC = 0
        GOTO 120
      ELSE
        C = A / B
        NC = NA - NB
      ENDIF
C
C   Normalize the result to a decent range if it is not.
C
 100  IF (ABS (C) .GE. BDX) THEN
        C = RDX * C
        NC = NC + NBT
        GOTO 100
      ENDIF
C
 110  IF (ABS (C) .LT. 1.D0) THEN
        C = BDX * C
        NC = NC - NBT
        GOTO 110
      ENDIF
C
 120  RETURN
      END
C
      SUBROUTINE DPMUL (A, NA, B, NB, C, NC)
C
C   This multiplies the DPE number (A, NA) by (B, NB) to yield the product
C   (C, NC).
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
C
      IF (IER .NE. 0) THEN
        C = 0.D0
        NC = 0
        RETURN
      ENDIF
C
C   Multiply A by B and add exponents, unless either is zero.
C
      IF (A .EQ. 0.D0 .OR. B .EQ. 0.D0) THEN
        C = 0.D0
        NC = 0
        GOTO 120
      ELSE
        C = A * B
        NC = NA + NB
      ENDIF
C
C   Normalize the result to a decent range if it is not.
C
 100  IF (ABS (C) .GE. BDX) THEN
        C = RDX * C
        NC = NC + NBT
        GOTO 100
      ENDIF
C
 110  IF (ABS (C) .LT. 1.D0) THEN
        C = BDX * C
        NC = NC - NBT
        GOTO 110
      ENDIF
C
 120  RETURN
      END
C
      SUBROUTINE DPPWR (A, NA, B, NB, C, NC)
C
C   This raises the DPE number (A, NA) to the (B, NB) power and places the
C   result in (C, NC).
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (CL2 = 1.4426950408889633D0)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
C
      IF (IER .NE. 0) THEN
        C = 0.D0
        NC = 0
        RETURN
      ENDIF
      IF (A .LE. 0.D0) THEN
        IF (KER(2) .NE. 0) THEN
          WRITE (LDB, 1)
 1        FORMAT ('*** DPPWR: Argument is less than or equal to zero.')
          IER = 2
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
      IF (B .EQ. 0.D0) THEN
        C = 1.D0
        NC = 0
        GOTO 120
      ENDIF
C
      IF (B .EQ. 1.D0 .AND. NB .EQ. 0) THEN
        C = A
        NC = NA
        GOTO 120
      ENDIF
C
C   Compute the base 2 logarithm of A and multiply by B.
C
      AL = CL2 * LOG (A) + NA
      CALL DPMUL (AL, 0, B, NB, T1, N1)
C
C   Check for possible overflow or underflow.
C
      IF (N1 .GT. 6) THEN
        IF (T1 .GT. 0.D0) THEN
          IF (KER(3) .NE. 0) THEN
            WRITE (LDB, 2)
 2          FORMAT ('*** DPPWR: Overflow')
            IER = 3
            IF (KER(IER) .EQ. 2) CALL MPABRT
          ENDIF
          RETURN
        ELSE
          C = 0.D0
          NC = 0
          GOTO 120
        ENDIF
      ENDIF
C
C   Compute 2 raised to the power B * Log_2 (A).
C
      T1 = T1 * 2.D0 ** N1
      NC = INT (T1)
      C = 2.D0 ** (T1 - NC)
C
C   Normalize the result to a decent range if it is not.
C
 100  IF (ABS (C) .GE. BDX) THEN
        C = RDX * C
        NC = NC + NBT
        GOTO 100
      ENDIF
C
 110  IF (ABS (C) .LT. 1.D0) THEN
        C = BDX * C
        NC = NC - NBT
        GOTO 110
      ENDIF
C
 120  RETURN
      END
C
      SUBROUTINE DPSQRT (A, NA, B, NB)
C
C   This computes the square root of the DPE number (A, NA) and places the
C   result in (B, NB).
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
C
      IF (IER .NE. 0) THEN
        B = 0.D0
        NB = 0
        RETURN
      ENDIF
      IF (A .LT. 0.D0) THEN
        IF (KER(4) .NE. 0) THEN
          WRITE (LDB, 1)
 1        FORMAT ('*** DPSQRT: Argument is negative.')
          IER = 4
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
      IF (A .EQ. 0.D0) THEN
        B = 0.D0
        NB = 0
        GOTO 120
      ENDIF
C
C   Divide the exponent of A by two and then take the square root of A.  If
C   NA is not an even number, then we have to multiply A by 10 before taking
C   the square root.
C
      NB = NA / 2
      IF (NA .EQ. 2 * NB) THEN
        B = SQRT (A)
      ELSE
        B = SQRT (2.D0 * A)
        IF (NA .LT. 0) NB = NB - 1
      ENDIF
C
C   Normalize the result to a decent range if it is not.
C
 100  IF (ABS (B) .GE. BDX) THEN
        B = RDX * B
        NB = NB + NBT
        GOTO 100
      ENDIF
C
 110  IF (ABS (B) .LT. 1.D0) THEN
        B = BDX * B
        NB = NB - NBT
        GOTO 110
      ENDIF
C
 120  RETURN
      END
C
      SUBROUTINE DPSUB (A, NA, B, NB, C, NC)
C
C   This subtracts the DPE number (B, NB) from (A, NA) to yield the difference
C   (C, NC).
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
C
      IF (IER .NE. 0) THEN
        C = 0.D0
        NC = 0
        RETURN
      ENDIF
C
      BB = -B
      CALL DPADD (A, NA, BB, NB, C, NC)
C
      RETURN
      END
C
      SUBROUTINE MPABRT
C>
C   This routine terminates execution.  Many users will want to replace the
C   default STOP with a call to a system routine that provides a traceback.
C   Examples of code that produce traceback are included here (commented out)
C   for some systems.
C
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
C
      WRITE (LDB, 1) IER
 1    FORMAT ('*** MPABRT: Execution terminated, error code =',I4)
C
C   Use this line on Cray systems.
C
C      CALL ABORT
C
C   Use this line plus the C routine TRACBK (available from author) on
C   Silicon Graphics IRIS systems.
C
C      CALL TRACBK
C
C   On other systems, merely terminate execution.
C
      STOP
      END
C
      SUBROUTINE MPADD (A, B, C)
C
C   This routine adds MP numbers A and B to yield the MP sum C.  It attempts
C   to include all significance of A and B in the result, up to the maximum
C   mantissa length NW.  Debug output starts with IDB = 9.
C
C   Max SP space for C: NW + 4 cells.  Max DP scratch space: NW + 4 cells.
C
      DOUBLE PRECISION D
      PARAMETER (NDB = 22)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM4/ D(1024)
      DIMENSION A(NW+2), B(NW+2), C(NW+4)
C
      IF (IER .NE. 0) THEN
        C(1) = 0.
        C(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 9) THEN
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 1) (A(I), I = 1, NO)
 1      FORMAT ('MPADD I'/(6F12.0))
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 1) (B(I), I = 1, NO)
      ENDIF
C
      IA = SIGN (1., A(1))
      IB = SIGN (1., B(1))
      NA = MIN (INT (ABS (A(1))), NW)
      NB = MIN (INT (ABS (B(1))), NW)
C
C   This first IF block checks for zero inputs.
C
      IF (NA .EQ. 0) THEN
C
C   A is zero -- the result is B.
C
        C(1) = SIGN (NB, IB)
C
        DO 100 I = 2, NB + 2
          C(I) = B(I)
 100    CONTINUE
C
        GOTO 420
      ELSEIF (NB .EQ. 0) THEN
C
C   B is zero -- the result is A.
C
        C(1) = SIGN (NA, IA)
C
        DO 110 I = 2, NA + 2
          C(I) = A(I)
 110    CONTINUE
C
        GOTO 420
      ENDIF
      MA = A(2)
      MB = B(2)
C
C   This IF block breaks the problem into different branches depending on
C   the relative sizes of the exponents of A and B.
C
      IF (MA .EQ. MB) THEN
C
C   A and B have the same exponent.
C
        NM = MIN (NA, NB)
        NX = MAX (NA, NB)
        IF (IA .EQ. IB) THEN
C
C   A and B have the same exponent and sign.
C
          D(1) = SIGN (NX, IA)
          D(2) = MA
          D(NX+3) = 0.D0
          D(NX+4) = 0.D0
C
          DO 120 I = 3, NM + 2
            D(I) = DBLE (A(I)) + DBLE (B(I))
 120      CONTINUE
C
          IF (NA .GT. NB) THEN
C
C   A is longer than B -- include extra words of A in C.
C
            DO 130 I = NM + 3, NA + 2
              D(I) = A(I)
 130        CONTINUE
C
          ELSEIF (NB .GT. NA) THEN
C
C   B is longer than A -- include extra words of B in C.
C
            DO 140 I = NM + 3, NB + 2
              D(I) = B(I)
 140        CONTINUE
C
          ENDIF
        ELSE
C
C   A and B have the same exponent but the opposite sign.  It is thus
C   necessary to scan through each vector until we find an unequal word.
C
          DO 150 I = 3, NM + 2
            IF (A(I) .NE. B(I)) GOTO 180
 150      CONTINUE
C
C   All words up to the common length are equal.
C
          IF (NA .EQ. NB) THEN
C
C   The length of A is the same as B -- result is zero.
C
            C(1) = 0.D0
            C(2) = 0.D0
            GOTO 420
          ELSEIF (NA .GT. NB) THEN
C
C   A is longer -- thus trailing words of A are shifted to start of C.
C
            NN = NA - NB
            D(1) = SIGN (NN, IA)
            D(2) = A(2) - NB
            D(NN+3) = 0.D0
            D(NN+4) = 0.D0
C
            DO 160 I = 3, NN + 2
              D(I) = A(I+NB)
 160        CONTINUE
C
          ELSEIF (NB .GT. NA) THEN
C
C   B is longer -- thus trailing words of B are shifted to start of C.
C
            NN = NB - NA
            D(1) = SIGN (NN, IB)
            D(2) = B(2) - NA
            D(NN+3) = 0.D0
            D(NN+4) = 0.D0
C
            DO 170 I = 3, NN + 2
              D(I) = B(I+NA)
 170        CONTINUE
C
          ENDIF
          GOTO 410
C
C   An unequal word was found.
C
 180      K = I - 3
          IF (A(K+3) .GT. B(K+3)) THEN
C
C   A is larger -- subtract B (shifted) from A.
C
            D(1) = SIGN (NX - K, IA)
            D(2) = A(2) - K
            D(NX-K+3) = 0.D0
            D(NX-K+4) = 0.D0
C
            DO 190 I = 3, NM - K + 2
              D(I) = DBLE (A(I+K)) - DBLE (B(I+K))
 190        CONTINUE
C
            DO 200 I = NB - K + 3, NA - K + 2
              D(I) = A(I+K)
 200        CONTINUE
C
            DO 210 I = NA - K + 3, NB - K + 2
              D(I) = - B(I+K)
 210        CONTINUE
C
          ELSE
C
C   B is larger -- subtract A (shifted) from B.
C
            D(1) = SIGN (NX - K, IB)
            D(2) = B(2) - K
            D(NX-K+3) = 0.D0
            D(NX-K+4) = 0.D0
C
            DO 220 I = 3, NM - K + 2
              D(I) = DBLE (B(I+K)) - DBLE (A(I+K))
 220        CONTINUE
C
            DO 230 I = NB - K + 3, NA - K + 2
              D(I) = - A(I+K)
 230        CONTINUE
C
            DO 240 I = NA - K + 3, NB - K + 2
              D(I) = B(I+K)
 240        CONTINUE
C
          ENDIF
        ENDIF
      ELSEIF (MA .GT. MB) THEN
C
C   Exponent of A is greater.  In other words, A has a larger magnitude.
C
        MC = MA - MB
        LA = MIN (MC, NA)
        LB = MIN (MC + NB, NW + 2)
        LM = MIN (NA, LB)
        LX = MIN (MAX (NA, LB), NW)
        D(1) = SIGN (LX, IA)
        D(2) = A(2)
        D(LX+3) = 0.D0
        D(LX+4) = 0.D0
C
        DO 250 I = 3, LA + 2
          D(I) = A(I)
 250    CONTINUE
C
C   If B is shifted NW + 2 or more words to the right of A then C = A.
C
        IF (MC .GE. NW + 2) THEN
          D(1) = SIGN (NA, IA)
          D(LA+3) = 0.D0
          D(LA+4) = 0.D0
          GOTO 410
        ENDIF
        IF (MC .GT. NA) THEN
C
C   There is a gap between A and the shifted B.  Fill it with zeroes.
C
          DO 260 I = NA + 3, MC + 2
            D(I) = 0.D0
 260      CONTINUE
C
          LM = MC
        ENDIF
        IF (IA .EQ. IB) THEN
C
C   A and B have the same sign -- add common words with B shifted right.
C
          DO 270 I = MC + 3, LM + 2
            D(I) = DBLE (A(I)) + DBLE (B(I-MC))
 270      CONTINUE
C
C   Include tail of A or B, whichever is longer after shift.
C
          IF (NA .GT. LB) THEN
C
            DO 280 I = LM + 3, NA + 2
              D(I) = A(I)
 280        CONTINUE
C
          ELSE
C
            DO 290 I = LM + 3, LB + 2
              D(I) = B(I-MC)
 290        CONTINUE
C
          ENDIF
        ELSE
C
C   A and B have different signs -- subtract common words with B shifted right.
C
          DO 300 I = MC + 3, LM + 2
            D(I) = DBLE (A(I)) - DBLE (B(I-MC))
 300      CONTINUE
C
C   Include tail of A or B, whichever is longer after shift.
C
          DO 310 I = LM + 3, NA + 2
            D(I) = A(I)
 310      CONTINUE
C
          DO 320 I = LM + 3, LB + 2
            D(I) = - B(I-MC)
 320      CONTINUE
C
        ENDIF
      ELSE
C
C   Exponent of B is greater.  In other words, B has a larger magnitude.
C
        MC = MB - MA
        LB = MIN (MC, NB)
        LA = MIN (MC + NA, NW + 2)
        LM = MIN (NB, LA)
        LX = MIN (MAX (NB, LA), NW)
        D(1) = SIGN (LX, IB)
        D(2) = B(2)
        D(LX+3) = 0.D0
        D(LX+4) = 0.D0
C
        DO 330 I = 3, LB + 2
          D(I) = B(I)
 330    CONTINUE
C
C   If A is shifted NW + 2 or more words to the right of B then C = B.
C
        IF (MC .GE. NW + 2) THEN
          D(1) = SIGN (NB, IB)
          D(LB+3) = 0.D0
          D(LB+4) = 0.D0
          GOTO 410
        ENDIF
        IF (MC .GT. NB) THEN
C
C   There is a gap between B and the shifted A.  Fill it with zeroes.
C
          DO 340 I = NB + 3, MC + 2
            D(I) = 0.D0
 340      CONTINUE
C
          LM = MC
        ENDIF
        IF (IB .EQ. IA) THEN
C
C   B and A have the same sign -- add common words with A shifted right.
C
          DO 350 I = MC + 3, LM + 2
            D(I) = DBLE (B(I)) + DBLE (A(I-MC))
 350      CONTINUE
C
C   Include tail of B or A, whichever is longer after shift.
C
          DO 360 I = LM + 3, NB + 2
            D(I) = B(I)
 360      CONTINUE
C
          DO 370 I = LM + 3, LA + 2
            D(I) = A(I-MC)
 370      CONTINUE
C
        ELSE
C
C   B and A have different signs -- subtract common words with A shifted right.
C
          DO 380 I = MC + 3, LM + 2
            D(I) = DBLE (B(I)) - DBLE (A(I-MC))
 380      CONTINUE
C
C   Include tail of B or A, whichever is longer after shift.
C
          DO 390 I = LM + 3, NB + 2
            D(I) = B(I)
 390      CONTINUE
C
          DO 400 I = LM + 3, LA + 2
            D(I) = - A(I-MC)
 400      CONTINUE
C
        ENDIF
      ENDIF
C
C   Fix up result, since some words may be negative or exceed BDX.
C
 410  CALL MPNORM (C)
C
 420  IF (IDB .GE. 9) THEN
        NO = MIN (INT (ABS (C(1))), NDB) + 2
        WRITE (LDB, 2) (C(I), I = 1, NO)
 2      FORMAT ('MPADD O'/(6F12.0))
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPAGMX (A, B)
C
C   This performs the arithmetic-geometric mean (AGM) iterations.  This routine
C   is called by MPLOGX.  It is not intended to be called directly by the user.
C
C   Max SP space for A and B: NW + 4 cells.  Max SP scratch space: 6.5*NW + 35
C   cells.  Max DP scratch space: 12 * NW + 6 cells.
C
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
      DIMENSION A(NW+4), B(NW+4)
C
      IF (IER .NE. 0) THEN
        A(1) = 0.
        A(2) = 0.
        B(1) = 0.
        B(2) = 0.
        RETURN
      ENDIF
      N4 = NW + 4
      NS = 2 * N4
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N4
      S(K0) = 0.
      S(K0+1) = 0.
      L1 = 0
C
 100  L1 = L1 + 1
      IF (L1 .EQ. 50) THEN
        IF (KER(5) .NE. 0) THEN
          WRITE (LDB, 1)
 1        FORMAT ('*** MPAGMX: Iteration limit exceeded.')
          IER = 5
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
      ENDIF
C
      S1 = S(K0+1)
      CALL MPADD (A, B, S(K0))
      CALL MPMULD (S(K0), 0.5D0, 0, S(K1))
      CALL MPMULX (A, B, S(K0))
      CALL MPSQRX (S(K0), B)
      CALL MPEQ (S(K1), A)
      CALL MPSUB (A, B, S(K0))
C
C   Check for convergence.
C
      IF (S(K0) .NE. 0. .AND. (S(K0+1) .LT. S1 .OR. S(K0+1) .GE. -2))   
     $  GOTO 100
C
      ICS = ISS
      IF (IDB .GE. 6) WRITE (LDB, 2) L1, S(K0+1)
 2    FORMAT ('MPAGMX: Iter., Tol. Achieved =',I5,F8.0)
      RETURN
      END
C
      SUBROUTINE MPALER
C
C   This outputs error messages when a single precision scratch space
C   allocation error is detected.
C
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
C
      IF (KER(6) .NE. 0) THEN
        WRITE (LDB, 1) ICS - 1
 1      FORMAT ('*** MPALER: Insufficient single precision scratch ',   
     $    'space.'/ 'Allocate',I10,' cells in an array in common ',     
     $    'MPCOM3 of the main '/ 'program and set IMS in common ',      
     $    'MPCOM1 to this size.')
        IER = 6
        IF (KER(IER) .EQ. 2) CALL MPABRT
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE MPANG (X, Y, PI, A)
C
C   This computes the MP angle A subtended by the MP pair (X, Y) considered as
C   a point in the x-y plane.  This is more useful than an arctan or arcsin
C   routine, since it places the result correctly in the full circle, i.e.
C   -Pi < A <= Pi.  PI is the MP value of Pi computed by a previous call to
C   MPPI.  For extra high levels of precision, use MPANGX.  The last word of
C   the result is not reliable.  Debug output starts with IDB = 5.
C
C   Max SP space for A: NW + 4 cells.  Max SP scratch space: 15 * NW + 88
C   cells.  Max DP scratch space: NW + 7 cells.
C
C   The Taylor series for Sin converges much more slowly than that of Arcsin.
C   Thus this routine does not employ Taylor series, but instead computes
C   Arccos or Arcsin by solving Cos (a) = x or Sin (a) = y using one of the
C   following Newton iterations, both of which converge to a:
C
C           z_{k+1} = z_k - [x - Cos (z_k)] / Sin (z_k)
C           z_{k+1} = z_k + [y - Sin (z_k)] / Cos (z_k)
C
C   The first is selected if Abs (x) <= Abs (y); otherwise the second is used.
C   These iterations are performed with a maximum precision level NW that
C   is dynamically changed, approximately doubling with each iteration.
C   See the comment about the parameter NIT in MPDIVX.
C
      DOUBLE PRECISION CL2, CPI, T1, T2, T3
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (CL2 = 1.4426950408889633D0, CPI = 3.141592653589793D0, 
     $  NIT = 3)
      DIMENSION A(NW+4), PI(NW+2), X(NW+2), Y(NW+2)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        A(1) = 0.
        A(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 5) THEN
        CALL MPDEB ('MPANG I', X)
        CALL MPDEB ('MPANG I', Y)
      ENDIF
C
      IX = SIGN (1., X(1))
      NX = MIN (INT (ABS (X(1))), NW)
      IY = SIGN (1., Y(1))
      NY = MIN (INT (ABS (Y(1))), NW)
C
C   Check if both X and Y are zero.
C
      IF (NX .EQ. 0 .AND. NY .EQ. 0) THEN
        IF (KER(7) .NE. 0) THEN
          WRITE (LDB, 1)
 1        FORMAT ('*** MPANG: Both arguments are zero.')
          IER = 7
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
C   Check if Pi has been precomputed.
C
      CALL MPMDC (PI, T1, N1)
      IF (N1 .NE. 0 .OR. ABS (T1 - CPI) .GT. RX2) THEN
        IF (KER(8) .NE. 0) THEN
          WRITE (LDB, 2)
 2        FORMAT ('*** MPANG: PI must be precomputed.')
          IER = 8
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
C   Check if one of X or Y is zero.
C
      IF (NX .EQ. 0) THEN
        IF (IY .GT. 0) THEN
          CALL MPMULD (PI, 0.5D0, 0, A)
        ELSE
          CALL MPMULD (PI, -0.5D0, 0, A)
        ENDIF
        GOTO 120
      ELSEIF (NY .EQ. 0) THEN
        IF (IX .GT. 0) THEN
          A(1) = 0.
          A(2) = 0.
        ELSE
          CALL MPEQ (PI, A)
        ENDIF
        GOTO 120
      ENDIF
C
      N5 = NW + 5
      NS = 5 * N5
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N5
      K2 = K1 + N5
      K3 = K2 + N5
      K4 = K3 + N5
      NWS = NW
      NW = NW + 1
C
C   Determine the least integer MQ such that 2 ^ MQ .GE. NW.
C
      T1 = NWS
      MQ = CL2 * LOG (T1) + 1.D0 - RXX
C
C   Normalize x and y so that x^2 + y^2 = 1.
C
      CALL MPMUL (X, X, S(K0))
      CALL MPMUL (Y, Y, S(K1))
      CALL MPADD (S(K0), S(K1), S(K2))
      CALL MPSQRT (S(K2), S(K3))
      CALL MPDIV (X, S(K3), S(K1))
      CALL MPDIV (Y, S(K3), S(K2))
C
C   Compute initial approximation of the angle.
C
      CALL MPMDC (S(K1), T1, N1)
      CALL MPMDC (S(K2), T2, N2)
      N1 = MAX (N1, -66)
      N2 = MAX (N2, -66)
      T1 = T1 * 2.D0 ** N1
      T2 = T2 * 2.D0 ** N2
      T3 = ATAN2 (T2, T1)
      CALL MPDMC (T3, 0, A)
C
C   The smaller of x or y will be used from now on to measure convergence.
C   This selects the Newton iteration (of the two listed above) that has the
C   largest denominator.
C
      IF (ABS (T1) .LE. ABS (T2)) THEN
        KK = 1
        CALL MPEQ (S(K1), S(K0))
      ELSE
        KK = 2
        CALL MPEQ (S(K2), S(K0))
      ENDIF
C
      NW = 3
      IQ = 0
C
C   Perform the Newton-Raphson iteration described above with a dynamically
C   changing precision level NW (one greater than powers of two).
C
      DO 110 K = 2, MQ
        NW = MIN (2 * NW - 2, NWS) + 1
 100    CONTINUE
        CALL MPCSSN (A, PI, S(K1), S(K2))
        IF (KK .EQ. 1) THEN
          CALL MPSUB (S(K0), S(K1), S(K3))
          CALL MPDIV (S(K3), S(K2), S(K4))
          CALL MPSUB (A, S(K4), S(K1))
        ELSE
          CALL MPSUB (S(K0), S(K2), S(K3))
          CALL MPDIV (S(K3), S(K1), S(K4))
          CALL MPADD (A, S(K4), S(K1))
        ENDIF
        CALL MPEQ (S(K1), A)
        IF (K .EQ. MQ - NIT .AND. IQ .EQ. 0) THEN
          IQ = 1
          GOTO 100
        ENDIF
 110  CONTINUE
C
C   Restore original precision level.
C
      NW = NWS
      ICS = ISS
      CALL MPROUN (A)
C
 120  IF (IDB .GE. 5) CALL MPDEB ('MPANG O', A)
C
      RETURN
      END
C
      SUBROUTINE MPANGX (X, Y, PI, A)
C
C   This computes the MP angle A subtended by the MP pair (X, Y) considered as
C   a point in the x-y plane.  This is more useful than an arctan or arcsin
C   routine, since it places the result correctly in the full circle, i.e.
C   -Pi < A <= Pi.  PI is the MP value of Pi computed by a previous call to
C   MPPI or MPPIX.  Before calling MPANGX, the array in MPCOM5 must be
C   initialized by calling MPINIX.  For modest levels of precision, use MPANG.
C   NW should be a power of two.  The last three words of the result are not
C   reliable.  Debug output starts with IDB = 6.
C
C   Max SP space for A: NW + 4 cells.  Max SP scratch space: 19.5 * NW + 87
C   cells.  Max DP scratch space: 12 * NW + 6 cells.
C
C   This routine employs a complex arithmetic version of the MPLOGX alogirthm.
C
      DOUBLE PRECISION CPI, T1
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (CPI = 3.141592653589793D0)
      DIMENSION A(NW+4), F0(8), F1(8), F4(8), PI(NW+2), X(NW+2), Y(NW+2)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        A(1) = 0.
        A(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 6) THEN
        CALL MPDEB ('MPANGX I', X)
        CALL MPDEB ('MPANGX I', Y)
      ENDIF
C
      IX = SIGN (1., X(1))
      NX = MIN (INT (ABS (X(1))), NW)
      IY = SIGN (1., Y(1))
      NY = MIN (INT (ABS (Y(1))), NW)
      NCR = 2 ** MCR
C
C   Check if precision level is too low to justify the advanced routine.
C
      IF (NW .LE. NCR) THEN
        CALL MPANG (X, Y, PI, A)
        GOTO 100
      ENDIF
C
C   Check if both X and Y are zero.
C
      IF (NX .EQ. 0 .AND. NY .EQ. 0) THEN
        IF (KER(9) .NE. 0) THEN
          WRITE (LDB, 1)
 1        FORMAT ('*** MPANGX: Both arguments are zero.')
          IER = 9
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
C   Check if Pi has been precomputed.
C
      CALL MPMDC (PI, T1, N1)
      IF (N1 .NE. 0 .OR. ABS (T1 - CPI) .GT. RX2) THEN
        IF (KER(10) .NE. 0) THEN
          WRITE (LDB, 2)
 2        FORMAT ('*** MPANGX: PI must be precomputed.')
          IER = 10
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
C   Check if one of X or Y is zero.
C
      IF (NX .EQ. 0) THEN
        IF (IY .GT. 0) THEN
          CALL MPMULD (PI, 0.5D0, 0, A)
        ELSE
          CALL MPMULD (PI, -0.5D0, 0, A)
        ENDIF
        GOTO 100
      ELSEIF (NY .EQ. 0) THEN
        IF (IX .GT. 0) THEN
          A(1) = 0.
          A(2) = 0.
        ELSE
          CALL MPEQ (PI, A)
        ENDIF
        GOTO 100
      ENDIF
C
C   Define scratch space.
C
      N4 = NW + 4
      N42 = 2 * N4
      NS = 4 * N42
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N42
      K2 = K1 + N42
      K3 = K2 + N42
      F0(1) = 0.
      F0(2) = 0.
      F0(3) = 0.
      F1(1) = 1.
      F1(2) = 0.
      F1(3) = 1.
      F4(1) = 1.
      F4(2) = 0.
      F4(3) = 4.
C
C   Multiply the input by a large power of two.
C
      CALL MPMDC (X, T1, N1)
      N2 = NBT * (NW / 2 + 2) - N1
      TN = N2
      CALL MPMULD (X, 1.D0, N2, S(K1))
      CALL MPMULD (Y, 1.D0, N2, S(K2))
      CALL MPMMPC (S(K1), S(K2), N4, S(K0))
C
C   Perform AGM iterations.
C
      CALL MPMMPC (F1, F0, N4, S(K1))
      CALL MPMMPC (F4, F0, N4, S(K3))
      CALL MPCDVX (N4, S(K3), S(K0), S(K2))
      CALL MPCAGX (S(K1), S(K2))
C
C   Compute A = Imag (Pi / (2 * Z)), where Z is the limit of the complex AGM.
C
      CALL MPMULD (S(K1), 2.D0, 0, S(K0))
      CALL MPMULD (S(K1+N4), 2.D0, 0, S(K0+N4))
      CALL MPMMPC (PI, F0, N4, S(K2))
      CALL MPCDVX (N4, S(K2), S(K0), S(K1))
      CALL MPEQ (S(K1+N4), A)
      ICS = ISS
C
 100  IF (IDB .GE. 6) CALL MPDEB ('MPANGX O', A)
C
      RETURN
      END
C
      SUBROUTINE MPCADD (L, A, B, C)
C
C   This computes the sum of the MPC numbers A and B and returns the MPC
C   result in C.  L is the offset between real and imaginary parts in A, B
C   and C.  L must be at least NW + 4.  Debug output starts with IDB = 9.
C
C   Max SP space for C: 2 * L cells.
C
      DIMENSION A(2*L), B(2*L), C(2*L)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
C
      IF (IER .NE. 0) THEN
        C(1) = 0.
        C(2) = 0.
        C(L+1) = 0.
        C(L+2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 9) WRITE (LDB, 1)
 1    FORMAT ('MPCADD')
C
      IF (L .LT. NW + 4) THEN
        IF (KER(11) .NE. 0) THEN
          WRITE (LDB, 2) L, NW + 4
 2        FORMAT ('*** MPCADD: Offset parameter is too small',2I8)
          IER = 11
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
      L1 = L + 1
      CALL MPADD (A, B, C)
      CALL MPADD (A(L1), B(L1), C(L1))
C
      RETURN
      END
C
      SUBROUTINE MPCAGX (A, B)
C
C   This performs the arithmetic-geometric mean (AGM) iterations.  This routine
C   is called by MPANGX.  It is not intended to be called directly by the user.
C
C   Max SP space for A and B: 2*NW + 8 cells.  Max SP scratch space:11.5*NW+55
C   cells.  Max DP scratch space: 12 * NW + 6 cells.
C
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
      DIMENSION A(2*NW+8), B(2*NW+8)
C
      IF (IER .NE. 0) THEN
        A(1) = 0.
        A(2) = 0.
        B(1) = 0.
        B(2) = 0.
        RETURN
      ENDIF
      N4 = NW + 4
      NS = 4 * N4
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + 2 * N4
      S(K0) = 0.
      S(K0+1) = 0.
      L1 = 0
C
 100  L1 = L1 + 1
      IF (L1 .EQ. 50) THEN
        IF (KER(12) .NE. 0) THEN
          WRITE (LDB, 1)
 1        FORMAT ('*** MPCAGX: Iteration limit exceeded.')
          IER = 12
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
      ENDIF
C
      S1 = S(K0+1)
      CALL MPCADD (N4, A, B, S(K0))
      CALL MPMULD (S(K0), 0.5D0, 0, S(K1))
      CALL MPMULD (S(K0+N4), 0.5D0, 0, S(K1+N4))
      CALL MPCMLX (N4, A, B, S(K0))
      CALL MPCSQX (N4, S(K0), B)
      CALL MPCEQ (N4, S(K1), A)
      CALL MPSUB (A, B, S(K0))
C
C   Check for convergence.
C
      IF (S(K0) .NE. 0. .AND. (S(K0+1) .LT. S1 .OR. S(K0+1) .GE. -2))   
     $  GOTO 100
C
      ICS = ISS
      IF (IDB .GE. 6) WRITE (LDB, 2) L1, S(K0+1)
 2    FORMAT ('MPCAGX: Iter., Tol. Achieved =',I5,F8.0)
      RETURN
      END
C
      SUBROUTINE MPCBRT (A, B)
C
C   This computes the cube root of the MP number A and returns the MP result
C   in B.  For extra high levels of precision, use MPCBRX.  Debug output
C   starts with IDB = 7.
C
C   Max SP space for B: NW + 4 cells.  Max SP scratch space: 3 * NW + 15
C   cells.  Max DP scratch space: NW + 5 cells.
C
C   This subroutine employs the following Newton-Raphson iteration, which
C   converges to A ^ (-2/3):
C
C          X_{k+1} = X_k + (1 - X_k^3 * A^2) * X_k / 3
C
C   where the muliplication () * X_k is performed with only half of the
C   normal level of precision.  These iterations are performed with a
C   maximum precision level NW that is dynamically changed, doubling with
C   each iteration.  The final iteration is performed as follows (this is
C   due to A. Karp):
C
C          Cbrt(A) = (A * X_n) + [A - (A * X_n)^3] * X_n / 3 (approx.)
C
C   where the multiplications A * X_n and [] * X_n are performed with only
C   half of the final level of precision.  See the comment about the parameter
C   NIT in MPDIVX.
C
      DOUBLE PRECISION CL2, T1, T2
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (CL2 = 1.4426950408889633D0, NDB = 22, NIT = 3)
      DIMENSION A(NW+2), B(NW+4), F(8)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        B(1) = 0.
        B(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 7) THEN
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 1) (A(I), I = 1, NO)
 1      FORMAT ('MPCBRT I'/(6F12.0))
      ENDIF
C
      IA = SIGN (1., A(1))
      NA = MIN (INT (ABS (A(1))), NW)
C
      IF (NA .EQ. 0) THEN
        B(1) = 0.
        B(2) = 0.
        GOTO 120
      ENDIF
      IF (IA .LT. 0.D0) THEN
        IF (KER(13) .NE. 0) THEN
          WRITE (LDB, 2)
 2        FORMAT ('*** MPCBRT: Argument is negative.')
          IER = 13
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
      N5 = NW + 5
      NS = 3 * N5
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N5
      K2 = K1 + N5
      NWS = NW
C
C   Determine the least integer MQ such that 2 ^ MQ .GE. NW.
C
      T1 = NW
      MQ = CL2 * LOG (T1) + 1.D0 - RXX
C
C   Compute A^2 outside of the iteration loop.
C
      NW = NWS + 1
      CALL MPMUL (A, A, S(K0))
C
C   Compute the initial approximation of A ^ (-2/3).
C
      CALL MPMDC (A, T1, N)
      N3 = - 2 * N / 3
      T2 = (T1 * 2.D0 ** (N + 3.D0 * N3 / 2.D0)) ** (-2.D0 / 3.D0)
      CALL MPDMC (T2, N3, B)
      F(1) = 1.
      F(2) = 0.
      F(3) = 1.
      NW = 3
      IQ = 0
C
C   Perform the Newton-Raphson iteration described above with a dynamically
C   changing precision level NW (one greater than powers of two).
C
      DO 110 K = 2, MQ
        NW1 = NW
        NW = MIN (2 * NW - 2, NWS) + 1
        NW2 = NW
 100    CONTINUE
        CALL MPMUL (B, B, S(K1))
        CALL MPMUL (B, S(K1), S(K2))
        CALL MPMUL (S(K0), S(K2), S(K1))
        CALL MPSUB (F, S(K1), S(K2))
        NW = NW1
        CALL MPMUL (B, S(K2), S(K1))
        CALL MPDIVD (S(K1), 3.D0, 0, S(K2))
        NW = NW2
        CALL MPADD (B, S(K2), S(K1))
        CALL MPEQ (S(K1), B)
        IF (K .EQ. MQ - NIT .AND. IQ .EQ. 0) THEN
          IQ = 1
          GOTO 100
        ENDIF
 110  CONTINUE
C
C   Perform last iteration using Karp's trick.
C
      CALL MPMUL (A, B, S(K0))
      NW1 = NW
      NW = MIN (2 * NW - 2, NWS) + 1
      NW2 = NW
      CALL MPMUL (S(K0), S(K0), S(K1))
      CALL MPMUL (S(K0), S(K1), S(K2))
      CALL MPSUB (A, S(K2), S(K1))
      NW = NW1
      CALL MPMUL (S(K1), B, S(K2))
      CALL MPDIVD (S(K2), 3.D0, 0, S(K1))
      NW = NW2
      CALL MPADD (S(K0), S(K1), S(K2))
      CALL MPEQ (S(K2), B)
C
C   Restore original precision level.
C
      NW = NWS
      ICS = ISS
      CALL MPROUN (B)
C
 120  IF (IDB .GE. 7) THEN
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 3) (B(I), I = 1, NO)
 3      FORMAT ('MPCBRT O'/(6F12.0))
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPCBRX (A, B)
C
C   This computes the cube root of the MP number A and returns the MP result
C   in B.  Before calling MPCBRX, the array in MPCOM5 must be initialized by
C   calling MPINIX.  For modest levels of precision, use MPCBRT.  NW should be
C   a power of two.  The last three words of the result are not reliable.
C   Debug output starts with IDB = 6.
C
C   Max SP space for B: NW + 4 cells.  Max SP scratch space: 4.5 * NW + 27
C   cells.  Max DP scratch space: 12 * NW + 6 cells.
C
C   This routine uses basically the same Newton iteration algorithm as MPCBRT.
C   In fact, this routine calls MPCBRT to obtain an initial approximation.
C   See the comment about the parameter NIT in MPDIVX.
C
      DOUBLE PRECISION CL2, T1
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (CL2 = 1.4426950408889633D0, NDB = 22, NIT = 1)
      DIMENSION A(NW+2), B(NW+4), F(8)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        B(1) = 0.
        B(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 6) THEN
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 1) (A(I), I = 1, NO)
 1      FORMAT ('MPCBRX I'/(6F12.0))
      ENDIF
C
      IA = SIGN (1., A(1))
      NA = MIN (INT (ABS (A(1))), NW)
      NCR = 2 ** MCR
C
      IF (NA .EQ. 0) THEN
        B(1) = 0.
        B(2) = 0.
        GOTO 120
      ENDIF
      IF (IA .LT. 0.D0) THEN
        IF (KER(14) .NE. 0) THEN
          WRITE (LDB, 2)
 2        FORMAT ('*** MPCBRX: Argument is negative.')
          IER = 14
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
C   Check if precision level is too low to justify the advanced routine.
C
      IF (NW .LE. NCR) THEN
        CALL MPCBRT (A, B)
        GOTO 120
      ENDIF
      N4 = NW + 4
      NS = 3 * N4
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N4
      K2 = K1 + N4
      NWS = NW
C
C   Determine the least integer MQ such that 2 ^ MQ .GE. NW.
C
      T1 = NW
      MQ = CL2 * LOG (T1) + 1.D0 - RXX
C
C   Compute A^2 outside of the iteration loop.
C
      CALL MPSQX (A, S(K0))
C
C   Compute the initial approximation of A ^ (-2/3).
C
      NW = NCR
      CALL MPCBRT (A, S(K1))
      CALL MPDIV (S(K1), A, B)
      F(1) = 1.
      F(2) = 0.
      F(3) = 1.
      IQ = 0
C
C   Perform the Newton-Raphson iteration described above with a dynamically
C   changing precision level NW (powers of two).
C
      DO 110 K = MCR + 1, MQ
        NW1 = NW
        NW = MIN (2 * NW, NWS)
        NW2 = NW
 100    CONTINUE
        CALL MPSQX (B, S(K1))
        CALL MPMULX (B, S(K1), S(K2))
        CALL MPMULX (S(K0), S(K2), S(K1))
        CALL MPSUB (F, S(K1), S(K2))
        NW = NW1
        CALL MPMULX (B, S(K2), S(K1))
        CALL MPDIVD (S(K1), 3.D0, 0, S(K2))
        NW = NW2
        CALL MPADD (B, S(K2), S(K1))
        CALL MPEQ (S(K1), B)
        IF (K .EQ. MQ - NIT .AND. IQ .EQ. 0) THEN
          IQ = 1
          GOTO 100
        ENDIF
 110  CONTINUE
C
C   Perform last iteration using Karp's trick.
C
      CALL MPMUL (A, B, S(K0))
      NW1 = NW
      NW = MIN (2 * NW, NWS)
      NW2 = NW
      CALL MPMUL (S(K0), S(K0), S(K1))
      CALL MPMUL (S(K0), S(K1), S(K2))
      CALL MPSUB (A, S(K2), S(K1))
      NW = NW1
      CALL MPMUL (S(K1), B, S(K2))
      CALL MPDIVD (S(K2), 3.D0, 0, S(K1))
      NW = NW2
      CALL MPADD (S(K0), S(K1), S(K2))
      CALL MPEQ (S(K2), B)
      ICS = ISS
C
 120  IF (IDB .GE. 6) THEN
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 3) (B(I), I = 1, NO)
 3      FORMAT ('MPCBRX O'/(6F12.0))
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPCDIV (L, A, B, C)
C
C   This routine divides the MP complex numbers A and B to yield the MPC
C   quotient C.  L is the offset between real and imaginary parts in A, B
C   and the result C.  L must be at least NW + 4.  For extra high levels of
C   precision, use MPCDVX.  The last word is not reliable.  Debug output
C   starts with IDB = 7
C
C   Max SP space for C: 2 * L cells.  Max SP scratch space: 5 * NW + 20
C   cells.  Max DP scratch space: NW + 4 cells.
C
C   This routine employs the formula described in MPCMUL to save multiprecision
C   multiplications.
C
      PARAMETER (NDB = 22)
      DIMENSION A(2*L), B(2*L), C(2*L), F(8)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        C(1) = 0.
        C(2) = 0.
        C(L+1) = 0.
        C(L+2) = 0.
        RETURN
      ENDIF
      L1 = L + 1
      IF (IDB .GE. 7) THEN
        WRITE (LDB, 1) L
 1      FORMAT ('MPCDIV I',I10)
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 2) (A(I), I = 1, NO)
 2      FORMAT ('MPCDIV I'/(6F12.0))
        NO = MIN (INT (ABS (A(L1))), NDB) + 2
        WRITE (LDB, 2) (A(L+I), I = 1, NO)
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 2) (B(I), I = 1, NO)
        NO = MIN (INT (ABS (B(L1))), NDB) + 2
        WRITE (LDB, 2) (B(L+I), I = 1, NO)
      ENDIF
C
      IF (L .LT. NW + 4) THEN
        IF (KER(15) .NE. 0) THEN
          WRITE (LDB, 3) L, NW + 4
 3        FORMAT ('*** MPCDIV: Offset parameter is too small',2I8)
          IER = 15
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
      IF (B(1) .EQ. 0. .AND. B(L1) .EQ. 0.) THEN
        IF (KER(16) .NE. 0) THEN
          WRITE (LDB, 4)
 4        FORMAT ('*** MPCDIV: Divisor is zero.')
          IER = 16
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
      N4 = NW + 4
      NS = 5 * N4
      ISS = ICS
      ICS = ICS + NS
      IF (ICS - 1 .GT. IMS) CALL MPALER
      IHS = MAX (ICS, IHS)
      K0 = ISS
      K1 = K0 + N4
      K2 = K1 + N4
      K3 = K2 + N4
      K4 = K3 + N4
      F(1) = 1.
      F(2) = 0.
      F(3) = 1.
C
      CALL MPMUL (A, B, S(K0))
      CALL MPMUL (A(L1), B(L1), S(K1))
      CALL MPADD (S(K0), S(K1), S(K2))
      CALL MPSUB (S(K0), S(K1), S(K3))
      CALL MPADD (A, A(L1), S(K0))
      CALL MPSUB (B, B(L1), S(K1))
      CALL MPMUL (S(K0), S(K1), S(K4))
      CALL MPSUB (S(K4), S(K3), S(K1))
      CALL MPMUL (B, B, S(K0))
      CALL MPMUL (B(L1), B(L1), S(K3))
      CALL MPADD (S(K0), S(K3), S(K4))
      CALL MPDIV (F, S(K4), S(K0))
      CALL MPMUL (S(K2), S(K0), C)
      CALL MPMUL (S(K1), S(K0), C(L1))
      ICS = ISS
C
      IF (IDB .GE. 7) THEN
        NO = MIN (INT (ABS (C(1))), NDB) + 2
        WRITE (LDB, 5) (C(I), I = 1, NO)
 5      FORMAT ('MPCDIV O'/(6F12.0))
        NO = MIN (INT (ABS (C(L1))), NDB) + 2
        WRITE (LDB, 5) (C(L+I), I = 1, NO)
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPCDVX (L, A, B, C)
C
C   This routine divides the MP complex numbers A and B to yield the MPC
C   quotient C.  L is the offset between real and imaginary parts in A, B
C   the result C.  L must be at least NW + 4.  Before calling MPCDVX, the
C   array in MPCOM5 must be initialized by calling MPINIX.  For modest levels
C   of precision, use MPCDIV.  NW should be a power of two.  The last two
C   words are not reliable.  Debug output starts with IDB = 7
C
C   Max SP space for C: 2 * L cells.  Max SP scratch space: 8 * NW + 32
C   cells.  Max DP scratch space: 12 * NW + 6 cells.
C
C   This routine employs the same scheme as MPCDIV.
C
      PARAMETER (NDB = 22)
      DIMENSION A(2*L), B(2*L), C(2*L), F(8)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        C(1) = 0.
        C(2) = 0.
        C(L+1) = 0.
        C(L+2) = 0.
        RETURN
      ENDIF
      L1 = L + 1
      IF (IDB .GE. 7) THEN
        WRITE (LDB, 1) L
 1      FORMAT ('MPCDVX I',I10)
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 2) (A(I), I = 1, NO)
 2      FORMAT ('MPCDVX I'/(6F12.0))
        NO = MIN (INT (ABS (A(L1))), NDB) + 2
        WRITE (LDB, 2) (A(L+I), I = 1, NO)
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 2) (B(I), I = 1, NO)
        NO = MIN (INT (ABS (B(L1))), NDB) + 2
        WRITE (LDB, 2) (B(L+I), I = 1, NO)
      ENDIF
C
      IF (L .LT. NW + 4) THEN
        IF (KER(17) .NE. 0) THEN
          WRITE (LDB, 3) L, NW + 4
 3        FORMAT ('*** MPCDVX: Offset parameter is too small',2I8)
          IER = 17
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
      IF (B(1) .EQ. 0. .AND. B(L1) .EQ. 0.) THEN
        IF (KER(18) .NE. 0) THEN
          WRITE (LDB, 4)
 4        FORMAT ('*** MPCDVX: Divisor is zero.')
          IER = 18
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
      N4 = NW + 4
      NS = 5 * N4
      ISS = ICS
      ICS = ICS + NS
      IF (ICS - 1 .GT. IMS) CALL MPALER
      IHS = MAX (ICS, IHS)
      K0 = ISS
      K1 = K0 + N4
      K2 = K1 + N4
      K3 = K2 + N4
      K4 = K3 + N4
      F(1) = 1.
      F(2) = 0.
      F(3) = 1.
C
      CALL MPMULX (A, B, S(K0))
      CALL MPMULX (A(L1), B(L1), S(K1))
      CALL MPADD (S(K0), S(K1), S(K2))
      CALL MPSUB (S(K0), S(K1), S(K3))
      CALL MPADD (A, A(L1), S(K0))
      CALL MPSUB (B, B(L1), S(K1))
      CALL MPMULX (S(K0), S(K1), S(K4))
      CALL MPSUB (S(K4), S(K3), S(K1))
      CALL MPSQX (B, S(K0))
      CALL MPSQX (B(L1), S(K3))
      CALL MPADD (S(K0), S(K3), S(K4))
      CALL MPDIVX (F, S(K4), S(K0))
      CALL MPMUL (S(K2), S(K0), C)
      CALL MPMUL (S(K1), S(K0), C(L1))
      ICS = ISS
C
      IF (IDB .GE. 7) THEN
        NO = MIN (INT (ABS (C(1))), NDB) + 2
        WRITE (LDB, 5) (C(I), I = 1, NO)
 5      FORMAT ('MPCDVX O'/(6F12.0))
        NO = MIN (INT (ABS (C(L1))), NDB) + 2
        WRITE (LDB, 5) (C(L+I), I = 1, NO)
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPCEQ (L, A, B)
C
C   This sets the MPC number B equal to the MPC number A.  L is the offset
C   between real and imaginary parts in A and B.  Debug output starts with
C   IDB = 10.
C
C   Max SP space for B: 2 * L cells.
C
      DIMENSION A(2*L), B(2*L)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
C
      IF (IER .NE. 0) THEN
        B(1) = 0.
        B(2) = 0.
        B(L+1) = 0.
        B(L+2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 10) WRITE (LDB, 1)
 1    FORMAT ('MPCEQ')
C
      I1 = SIGN (1., A(1))
      N1 = MIN (INT (ABS (A(1))), NW, L - 2)
      I2 = SIGN (1., A(L+1))
      N2 = MIN (INT (ABS (A(L+1))), NW, L - 2)
      B(1) = SIGN (N1, I1)
      B(L+1) = SIGN (N2, I2)
C
      DO 100 I = 2, N1 + 2
        B(I) = A(I)
 100  CONTINUE
C
      DO 110 I = 2, N2 + 2
        B(L+I) = A(L+I)
 110  CONTINUE
C
      RETURN
      END
C
      SUBROUTINE MPCFFT (IS, M, X, Y)
C
C   This routine computes the 2^M -point complex-to-complex FFT of X.  See
C   article by DHB in Intl. J. of Supercomputer Applications, Spring 1988,
C   p. 82 - 87).  X and Y are double precision.  X is both the input and the
C   output array, while Y is a scratch array.  Both X and Y must be
C   dimensioned with 2 * N cells, where N = 2^M.  The data in X are assumed
C   to have real and imaginary parts separated by N cells.  A call to MPCFFT
C   with IS = 1 (or -1) indicates a call to perform a FFT with positive (or
C   negative) exponentials.  M must be at least two.  Before calling MPCRFT,
C   the array in MPCOM5 must be initialized by calling MPINIX.
C
C   In this application, MPCFFT is called by MPRCFT and MPCRFT, which are in
C   turn called by MPMULX.  This routine is not intended to be called directly
C   by the user.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION X(*), Y(*)
C
      N = 2 ** M
C>
C   For Cray computers, it is most efficient to limit M1 to 6.  For most
C   scalar computers, it is best to limit M1 to 2.  Uncomment whichever of the
C   next two lines is appropriate.
C
C      M1 = MIN (M / 2, 6)
      M1 = MIN (M / 2, 2)
      M2 = M - M1
      N2 = 2 ** M1
      N1 = 2 ** M2
C
C   Perform one variant of the Stockham FFT.
C
      DO 100 L = 1, M1, 2
        CALL MPFFT1 (IS, L, M, X, Y)
        IF (L .EQ. M1) GOTO 120
        CALL MPFFT1 (IS, L + 1, M, Y, X)
 100  CONTINUE
C
C   Perform a transposition of X treated as a N2 x N1 x 2 matrix.
C
      CALL MPTRAN (N1, N2, X, Y)
C
C   Perform second variant of the Stockham FFT from Y to X and X to Y.
C
      DO 110 L = M1 + 1, M, 2
        CALL MPFFT2 (IS, L, M, Y, X)
        IF (L .EQ. M) GOTO 160
        CALL MPFFT2 (IS, L + 1, M, X, Y)
 110  CONTINUE
C
      GOTO 140
C
C   Perform a transposition of Y treated as a N2 x N1 x 2 matrix.
C
 120  CALL MPTRAN (N1, N2, Y, X)
C
C   Perform second variant of the Stockham FFT from X to Y and Y to X.
C
      DO 130 L = M1 + 1, M, 2
        CALL MPFFT2 (IS, L, M, X, Y)
        IF (L .EQ. M) GOTO 140
        CALL MPFFT2 (IS, L + 1, M, Y, X)
 130  CONTINUE
C
      GOTO 160
C
C   Copy Y to X.
C
 140  DO 150 I = 1, 2 * N
        X(I) = Y(I)
 150  CONTINUE
C
 160  RETURN
      END
C
      SUBROUTINE MPCMLX (L, A, B, C)
C
C   This routine multiplies the MP complex numbers A and B to yield the MPC
C   product C.  L is the offset between real and imaginary parts in A, B and
C   the result C.  L must be at least NW + 4.  Before calling MPCMLX, the
C   array in MPCOM5 must be initialized by calling MPINIX.  For modest levels
C   of precision, use MPCMUL.  NW should be a power of two.  The last word is
C   not reliable.  Debug output starts with IDB = 7.
C
C   Max SP space for C: 2 * L cells.  Max SP scratch space: 4 * NW + 16
C   cells.  Max DP scratch space: 12 * NW + 6 cells.
C
C   This routine employs the same scheme as MPCMUL.
C
      PARAMETER (NDB = 22)
      DIMENSION A(2*L), B(2*L), C(2*L)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        C(1) = 0.
        C(2) = 0.
        C(L+1) = 0.
        C(L+2) = 0.
        RETURN
      ENDIF
      L1 = L + 1
      IF (IDB .GE. 7) THEN
        WRITE (LDB, 1) L
 1      FORMAT ('MPCMLX I',I10)
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 2) (A(I), I = 1, NO)
 2      FORMAT ('MPCMLX I'/(6F12.0))
        NO = MIN (INT (ABS (A(L1))), NDB) + 2
        WRITE (LDB, 2) (A(L+I), I = 1, NO)
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 2) (B(I), I = 1, NO)
        NO = MIN (INT (ABS (B(L1))), NDB) + 2
        WRITE (LDB, 2) (B(L+I), I = 1, NO)
      ENDIF
C
      IF (L .LT. NW + 4) THEN
        IF (KER(19) .NE. 0) THEN
          WRITE (LDB, 3) L, NW + 4
 3        FORMAT ('*** MPCMLX: Offset parameter is too small',2I8)
          IER = 19
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
      N4 = NW + 4
      NS = 4 * N4
      ISS = ICS
      ICS = ICS + NS
      IF (ICS - 1 .GT. IMS) CALL MPALER
      IHS = MAX (ICS, IHS)
      K0 = ISS
      K1 = K0 + N4
      K2 = K1 + N4
      K3 = K2 + N4
C
      CALL MPMULX (A, B, S(K0))
      CALL MPMULX (A(L1), B(L1), S(K1))
      CALL MPSUB (S(K0), S(K1), C)
      CALL MPADD (S(K0), S(K1), S(K2))
      CALL MPADD (A, A(L1), S(K0))
      CALL MPADD (B, B(L1), S(K1))
      CALL MPMULX (S(K0), S(K1), S(K3))
      CALL MPSUB (S(K3), S(K2), C(L1))
      ICS = ISS
C
      IF (IDB .GE. 7) THEN
        NO = MIN (INT (ABS (C(1))), NDB) + 2
        WRITE (LDB, 4) (C(I), I = 1, NO)
 4      FORMAT ('MPCMLX O'/(6F12.0))
        NO = MIN (INT (ABS (C(L1))), NDB) + 2
        WRITE (LDB, 4) (C(L+I), I = 1, NO)
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPCMUL (L, A, B, C)
C
C   This routine multiplies the MP complex numbers A and B to yield the MPC
C   product C.  L is the offset between real and imaginary parts in A, B and
C   the result C.  L must be at least NW + 4.  For extra high levels of
C   precision, use MPCMLX.  The last word is not reliable.  Debug output
C   starts with IDB = 7.
C
C   Max SP space for C: 2 * L cells.  Max SP scratch space: 4 * NW + 16
C   cells.  Max DP scratch space: NW + 4 cells.
C
C   This routine employs the formula
C
C   (a_1 + a_2 i) (b_1 + b_2 i)  =  [a_1 b_1 - a_2 b_2]  +
C                [(a_1 + b_1) (a_2 + b_2) - (a_1 b_1 + a_2 b_2)] i
C
C   Note that this formula can be implemented with only three multiplications
C   whereas the conventional formula requires four.
C
      PARAMETER (NDB = 22)
      DIMENSION A(2*L), B(2*L), C(2*L)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        C(1) = 0.
        C(2) = 0.
        C(L+1) = 0.
        C(L+2) = 0.
        RETURN
      ENDIF
      L1 = L + 1
      IF (IDB .GE.7) THEN
        WRITE (LDB, 1) L
 1      FORMAT ('MPCMUL I',I10)
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 2) (A(I), I = 1, NO)
 2      FORMAT ('MPCMUL I'/(6F12.0))
        NO = MIN (INT (ABS (A(L1))), NDB) + 2
        WRITE (LDB, 2) (A(L+I), I = 1, NO)
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 2) (B(I), I = 1, NO)
        NO = MIN (INT (ABS (B(L1))), NDB) + 2
        WRITE (LDB, 2) (B(L+I), I = 1, NO)
      ENDIF
C
      IF (L .LT. NW + 4) THEN
        IF (KER(20) .NE. 0) THEN
          WRITE (LDB, 3) L, NW + 4
 3        FORMAT ('*** MPCMUL: Offset parameter is too small',2I8)
          IER = 20
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
      N4 = NW + 4
      NS = 4 * N4
      ISS = ICS
      ICS = ICS + NS
      IF (ICS - 1 .GT. IMS) CALL MPALER
      IHS = MAX (ICS, IHS)
      K0 = ISS
      K1 = K0 + N4
      K2 = K1 + N4
      K3 = K2 + N4
C
      CALL MPMUL (A, B, S(K0))
      CALL MPMUL (A(L1), B(L1), S(K1))
      CALL MPSUB (S(K0), S(K1), C)
      CALL MPADD (S(K0), S(K1), S(K2))
      CALL MPADD (A, A(L1), S(K0))
      CALL MPADD (B, B(L1), S(K1))
      CALL MPMUL (S(K0), S(K1), S(K3))
      CALL MPSUB (S(K3), S(K2), C(L1))
      ICS = ISS
C
      IF (IDB .GE. 7) THEN
        NO = MIN (INT (ABS (C(1))), NDB) + 2
        WRITE (LDB, 4) (C(I), I = 1, NO)
 4      FORMAT ('MPCMUL O'/(6F12.0))
        NO = MIN (INT (ABS (C(L1))), NDB) + 2
        WRITE (LDB, 4) (C(L+I), I = 1, NO)
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPCPLX (N, LA, A, X1, NX, LX, X)
C
C   This routine finds a complex root of the N-th degree polynomial whose
C   MPC coefficients are in A by Newton-Raphson iterations, beginning
C   at the complex DPE value (X1(1), NX(1)) + i (X1(2), NX(2)), and returns
C   the MPC root in X.  The N + 1 coefficients a_0, a_1, ..., a_N are
C   assumed to start in locations A(1), A(2*LA+1), A(4*LA+1), etc.  LA is the
C   offset between the real and the imaginary parts of each input coefficient.
C   Typically LA = NW + 4.  LX, also an input parameter, is the offset between
C   the real and the imaginary parts of the result to be stored in X.  LX
C   should be at least NW + 4.  Before calling MPCPLX, the array in MPCOM5
C   be initialized by calling MPINIX.  For modest levels of precision, use
C   MPCPOL.  NW should be a power of two.  The last two words of the result
C   are not reliable.  Debug output starts with IDB = 5.
C
C   Max SP space for X: 2 * LX cells.  Max SP scratch space: 18 * NW + 72
C   cells.  Max DP scratch space: 12 * NW + 6 cells.
C
C   See the note in MPPOL about repeated roots.
C
C   This routine employs the same scheme as MPCPOL.
C
      CHARACTER*8 CX
      DOUBLE PRECISION T1, X1
      DIMENSION A(2*LA,N+1), NX(2), X(2*LX), X1(2)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        X(1) = 0.
        X(2) = 0.
        X(LX+1) = 0.
        X(LX+2) = 0.
      ENDIF
      IF (IDB .GE. 5) THEN
        WRITE (LDB, 1) N, LX
 1      FORMAT ('MPCPLX I',2I6)
C
        DO 100 K = 0, N
          WRITE (CX, '(I4)') K
          CALL MPDEB (CX, A(1,K+1))
          CALL MPDEB (CX, A(LA+1,K+1))
 100    CONTINUE
C
        WRITE (LDB, 2) X1(2), NX(2)
 2      FORMAT ('MPCPLX I',F16.12,' x 10 ^',I6,F20.12,' x 10^',I6)
      ENDIF
C
C   Check if precision level is too low to justify the advanced routine.
C
      NCR = 2 ** MCR
      IF (NW .LE. NCR) THEN
        CALL MPCPOL (N, LA, A, X1, NX, LX, X)
        L1 = 0
        GOTO 150
      ENDIF
C
C   Check if the polynomial is proper.
C
      IF (A(1,1) .EQ. 0. .OR. A(1,N+1) .EQ. 0.) THEN
        IF (KER(21) .NE. 0) THEN
          WRITE (LDB, 3)
 3        FORMAT ('*** MPCPLX: Either the first or last input ',        
     $      'coefficient is zero.')
          IER = 21
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
      N4 = NW + 4
      N8 = 2 * N4
      NS = 10 * N4
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N8
      K2 = K1 + N8
      K3 = K2 + N8
      K4 = K3 + N8
      NWS = NW
C
C   Set the initial value.
C
      NW = NCR
      CALL MPCPOL (N, LA, A, X1, NX, N4, S(K0))
      TL = -4.
      L1 = 0
      LS = -10
C
C   Perform MP Newton-Raphson iterations to solve P(x) = 0.
C
 110  L1 = L1 + 1
      IF (L1 .EQ. 50) THEN
        IF (KER(22) .NE. 0) THEN
          WRITE (LDB, 4)
 4        FORMAT ('*** MPCPLX: Iteration limit exceeded.')
          IER = 22
          IF (KER(IER) .EQ. 2) CALL MPABRT
          ICS = ISS
          NW = NWS
          RETURN
        ENDIF
      ENDIF
C
C   Compute P(x).
C
      CALL MPMMPC (A(1,N+1), A(LA+1,N+1), N4, S(K1))
C
      DO 120 K = N - 1, 0, -1
        CALL MPCMLX (N4, S(K0), S(K1), S(K2))
        CALL MPADD (S(K2), A(1,K+1), S(K1))
        CALL MPADD (S(K2+N4), A(LA+1,K+1), S(K1+N4))
 120  CONTINUE
C
C   Compute P'(x).
C
      T1 = N
      CALL MPMULD (A(1,N+1), T1, 0, S(K2))
      CALL MPMULD (A(LA+1,N+1), T1, 0, S(K2+N4))
C
      DO 130 K = N - 1, 1, -1
        CALL MPCMLX (N4, S(K0), S(K2), S(K3))
        T1 = K
        CALL MPMULD (A(1,K+1), T1, 0, S(K4))
        CALL MPMULD (A(LA+1,K+1), T1, 0, S(K4+N4))
        CALL MPCADD (N4, S(K3), S(K4), S(K2))
 130  CONTINUE
C
C   Compute P(x) / P'(x) and update x.
C
      CALL MPCDVX (N4, S(K1), S(K2), S(K3))
      CALL MPCSUB (N4, S(K0), S(K3), S(K4))
C
      IF (IDB .GE. 6) THEN
        WRITE (LDB, 5) L1
 5      FORMAT ('ITERATION',I4)
        CALL MPDEB ('X', S(K0))
        CALL MPDEB (' ', S(K0+N4))
        CALL MPDEB ('P(X)', S(K1))
        CALL MPDEB (' ', S(K1+N4))
        CALL MPDEB ('P''(X)', S(K2))
        CALL MPDEB (' ', S(K2+N4))
        CALL MPDEB ('CORR', S(K3))
        CALL MPDEB (' ', S(K3+N4))
      ENDIF
      CALL MPCEQ (N4, S(K4), S(K0))
C
C   If this was the second iteration at full precision, there is no need to
C   continue (the adjusted value of x is correct); otherwise repeat.
C
      IF (L1 .EQ. LS + 1) GOTO 140
      IF (S(K3) .NE. 0. .AND. S(K3+1) .GT. TL .OR. S(K3+N4) .NE. 0.     
     $  .AND. S(K3+N4+1) .GT. TL) GOTO 110
C
C   Newton iterations have converged to current precision.  Increase precision
C   and continue.
C
      IF (NW .EQ. NWS) GOTO 140
      NW = MIN (2 * NW, NWS)
      IF (NW .EQ. NWS) LS = L1
      IF (NW .LE. 32) THEN
        TL = 2 - NW
      ELSEIF (NW .LE. 256) THEN
        TL = 3 - NW
      ELSE
        TL = 4 - NW
      ENDIF
      IF (IDB .GE. 6) THEN
        WRITE (LDB, 6) NW
 6      FORMAT (6X,'New NW =', I8)
      ENDIF
      GOTO 110
C
 140  CALL MPMMPC (S(K0), S(K0+N4), LX, X)
      ICS = ISS
C
 150  IF (IDB .GE. 5) THEN
        WRITE (LDB, 7) L1
 7      FORMAT ('Iteration count:',I5)
        CALL MPDEB ('MPCPLX O', X)
        CALL MPDEB (' ', X(LX+1))
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPCPOL (N, LA, A, X1, NX, LX, X)
C
C   This routine finds a complex root of the N-th degree polynomial whose
C   MPC coefficients are in A by Newton-Raphson iterations, beginning
C   at the complex DPE value (X1(1), NX(1)) + i (X1(2), NX(2)), and returns
C   the MPC root in X.  The N + 1 coefficients a_0, a_1, ..., a_N are
C   assumed to start in locations A(1), A(2*LA+1), A(4*LA+1), etc.  LA is the
C   offset between the real and the imaginary parts of each input coefficient.
C   Typically LA = NW + 4.  LX, also an input parameter, is the offset between
C   the real and the imaginary parts of the result to be stored in X.  LX must
C   be at least NW + 4.  For extra high levels of precision, use MPCPLX.
C   Debug output starts with IDB = 5.
C
C   Max SP space for X: 2 * LX cells.  Max SP scratch space: 15 * NW + 75
C   cells.  Max DP scratch space: NW + 5 cells.
C
C   See the note about repeated roots in MPPOL.
C
C   This routine employs the complex form of the Newton-Raphson iteration:
C
C   X_{k+1} = X_k - P(X_k) / P'(X_k)
C
C   These iterations are performed with a maximum precision level NW that is
C   dynamically changed, approximately doubling with each iteration.
C
      CHARACTER*8 CX
      DOUBLE PRECISION T1, X1
      DIMENSION A(2*LA,N+1), NX(2), X(2*LX), X1(2)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        X(1) = 0.
        X(2) = 0.
        X(LX+1) = 0.
        X(LX+2) = 0.
      ENDIF
      IF (IDB .GE. 5) THEN
        WRITE (LDB, 1) N, LX
 1      FORMAT ('MPCPOL I',2I6)
C
        DO 100 K = 0, N
          WRITE (CX, '(I4)') K
          CALL MPDEB (CX, A(1,K+1))
          CALL MPDEB (CX, A(LA+1,K+1))
 100    CONTINUE
C
        WRITE (LDB, 2) X1(1), NX(1), X1(2), NX(2)
 2      FORMAT ('MPCPOL I',F16.12,' x 10 ^',I6,F20.12,' x 10^',I6)
      ENDIF
C
C  Check if the polynomial is proper.
C
      IF (A(1,1) .EQ. 0. .OR. A(1,N+1) .EQ. 0.) THEN
        IF (KER(23) .NE. 0) THEN
          WRITE (LDB, 3)
 3        FORMAT ('*** MPCPOL: Either the first or last input ',        
     $      'coefficient is zero.')
          IER = 23
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
      N5 = NW + 5
      N10 = 2 * N5
      NS = 10 * N5
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N10
      K2 = K1 + N10
      K3 = K2 + N10
      K4 = K3 + N10
      NWS = NW
      NW = NW + 1
C
C   Set the initial value.
C
      CALL MPDMC (X1(1), NX(1), S(K0))
      CALL MPDMC (X1(2), NX(2), S(K0+N5))
      NW = 5
      TL = -4.
      L1 = 0
      LS = -10
C
C   Perform MP Newton-Raphson iterations to solve P(x) = 0.
C
 110  L1 = L1 + 1
      IF (L1 .EQ. 50) THEN
        IF (KER(24) .NE. 0) THEN
          WRITE (LDB, 4)
 4        FORMAT ('*** MPCPOL: Iteration limit exceeded.')
          IER = 24
          IF (KER(IER) .EQ. 2) CALL MPABRT
          ICS = ISS
          NW = NWS
          RETURN
        ENDIF
      ENDIF
C
C   Compute P(x).
C
      CALL MPMMPC (A(1,N+1), A(LA+1,N+1), N5, S(K1))
C
      DO 120 K = N - 1, 0, -1
        CALL MPCMUL (N5, S(K0), S(K1), S(K2))
        CALL MPADD (S(K2), A(1,K+1), S(K1))
        CALL MPADD (S(K2+N5), A(LA+1,K+1), S(K1+N5))
 120  CONTINUE
C
C   Compute P'(x).
C
      T1 = N
      CALL MPMULD (A(1,N+1), T1, 0, S(K2))
      CALL MPMULD (A(LA+1,N+1), T1, 0, S(K2+N5))
C
      DO 130 K = N - 1, 1, -1
        CALL MPCMUL (N5, S(K0), S(K2), S(K3))
        T1 = K
        CALL MPMULD (A(1,K+1), T1, 0, S(K4))
        CALL MPMULD (A(LA+1,K+1), T1, 0, S(K4+N5))
        CALL MPCADD (N5, S(K3), S(K4), S(K2))
 130  CONTINUE
C
C   Compute P(x) / P'(x) and update x.
C
      CALL MPCDIV (N5, S(K1), S(K2), S(K3))
      CALL MPCSUB (N5, S(K0), S(K3), S(K4))
C
      IF (IDB .GE. 6) THEN
        WRITE (LDB, 5) L1
 5      FORMAT ('Iteration',I4)
        CALL MPDEB ('X', S(K0))
        CALL MPDEB (' ', S(K0+N5))
        CALL MPDEB ('P(X)', S(K1))
        CALL MPDEB (' ', S(K1+N5))
        CALL MPDEB ('P''(X)', S(K2))
        CALL MPDEB (' ', S(K2+N5))
        CALL MPDEB ('CORR', S(K3))
        CALL MPDEB (' ', S(K3+N5))
      ENDIF
      CALL MPCEQ (N5, S(K4), S(K0))
C
C   If this was the second iteration at full precision, there is no need to
C   continue (the adjusted value of x is correct); otherwise repeat.
C
      IF (L1 .EQ. LS + 1) GOTO 140
      IF (S(K3) .NE. 0. .AND. S(K3+1) .GT. TL .OR. S(K3+N5) .NE. 0.     
     $  .AND. S(K3+N5+1) .GT. TL) GOTO 110
C
C   Newton iterations have converged to current precision.  Increase precision
C   and continue.
C
      IF (NW .EQ. NWS + 1) GOTO 140
      NW = MIN (2 * NW - 2, NWS) + 1
      IF (NW .EQ. NWS + 1) LS = L1
      TL = 1 - NW
      IF (IDB .GE. 6) THEN
        WRITE (LDB, 6) NW
 6      FORMAT (6X,'New NW =', I8)
      ENDIF
      GOTO 110
C
 140  CALL MPMMPC (S(K0), S(K0+N5), LX, X)
C
C   Restore original precision level.
C
      NW = NWS
      ICS = ISS
      CALL MPROUN (X)
      CALL MPROUN (X(LX+1))
C
      IF (IDB .GE. 5) THEN
        WRITE (LDB, 7) L1
 7      FORMAT ('Iteration count:',I5)
        CALL MPDEB ('MPCPOL O', X)
        CALL MPDEB (' ', X(LX+1))
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPCPR (A, B, IC)
C
C   This routine compares the MP numbers A and B and returns in IC the value
C   -1, 0, or 1 depending on whether A < B, A = B, or A > B.  It is faster
C   than merely subtracting A and B and looking at the sign of the result.
C   Debug output begins with IDB = 9.
C
      DIMENSION A(NW+4), B(NW+4)
      PARAMETER (NDB = 22)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        IC = 0
        RETURN
      ENDIF
      IF (IDB .GE. 9) THEN
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 1) (A(I), I = 1, NO)
 1      FORMAT ('MPCPR I'/(6F12.0))
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 1) (B(I), I = 1, NO)
      ENDIF
      IA = SIGN (1., A(1))
      IF (A(1) .EQ. 0.) IA = 0
      IB = SIGN (1., B(1))
      IF (B(1) .EQ. 0.) IB = 0
C
C   Compare signs.
C
      IF (IA .NE. IB) THEN
        IC = SIGN (1, IA - IB)
        GOTO 110
      ENDIF
C
C   The signs are the same.  Compare exponents.
C
      MA = A(2)
      MB = B(2)
      IF (MA .NE. MB) THEN
        IC = IA * SIGN (1, MA - MB)
        GOTO 110
      ENDIF
C
C   The signs and the exponents are the same.  Compare mantissas.
C
      NA = MIN (INT (ABS (A(1))), NW)
      NB = MIN (INT (ABS (B(1))), NW)
C
      DO 100 I = 3, MIN (NA, NB) + 2
        IF (A(I) .NE. B(I)) THEN
          IC = IA * SIGN (1., A(I) - B(I))
          GOTO 110
        ENDIF
 100  CONTINUE
C
C   The mantissas are the same to the common length.  Compare lengths.
C
      IF (NA .NE. NB) THEN
        IC = IA * SIGN (1, NA - NB)
        GOTO 110
      ENDIF
C
C   The signs, exponents, mantissas and lengths are the same.  Thus A = B.
C
      IC = 0
C
 110  IF (IDB .GE. 9) WRITE (6, 2) IC
 2    FORMAT ('MPCPR O',I4)
      RETURN
      END
C
      SUBROUTINE MPCPWR (L, A, N, B)
C
C   This computes the N-th power of the MPC number A and returns the MPC
C   result C in B.  When N is zero, 1 is returned.  When N is negative, the
C   reciprocal of A ^ |N| is returned.  L is the offset between real and
C   imaginary parts in A and B.  L should be at least NW + 4.  For extra high
C   levels of precision, use MPCPWX.  Debug output starts with IDB = 7.
C
C   Max SP space for B: 2 * L cells.  Max SP scratch space: 11 * NW + 55
C   cells.  Max DP scratch space: NW + 5 cells.
C
C   This routine employs the binary method for exponentiation.
C
      DOUBLE PRECISION CL2, T1
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (CL2 = 1.4426950408889633D0, NDB = 22)
      DIMENSION A(2*L), B(2*L), F1(8), F2(8)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        B(1) = 0.
        B(2) = 0.
        B(L+1) = 0.
        B(L+2) = 0.
        RETURN
      ENDIF
      L1 = L + 1
      IF (IDB .GE. 7) THEN
        WRITE (6, 1) L, N
 1      FORMAT ('MPCPWR I',2I10)
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 2) (A(I), I = 1, NO)
 2      FORMAT ('MPCPWR I'/(6F12.0))
        NO = MIN (INT (ABS (A(L1))), NDB) + 2
        WRITE (LDB, 2) (A(L+I), I = 1, NO)
      ENDIF
C
      NA1 = MIN (INT (ABS (A(1))), NW)
      NA2 = MIN (INT (ABS (A(L1))), NW)
      IF (NA1 .EQ. 0 .AND. NA2 .EQ. 0) THEN
        IF (N .GE. 0) THEN
          B(1) = 0.
          B(2) = 0.
          B(L1) = 0.
          B(L1+1) = 0.
          GOTO 120
        ELSE
          IF (KER(25) .NE. 0) THEN
            WRITE (LDB, 3)
 3          FORMAT ('*** MPCPWR: Argument is zero and N is negative or',
     $        ' zero.')
            IER = 25
            IF (KER(IER) .EQ. 2) CALL MPABRT
          ENDIF
          RETURN
        ENDIF
      ENDIF
C
      N5 = NW + 5
      NS = 6 * N5
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + 2 * N5
      K2 = K1 + 2 * N5
      NWS = NW
      NW = NW + 1
      NN = ABS (N)
      F1(1) = 1.
      F1(2) = 0.
      F1(3) = 1.
      F2(1) = 0.
      F2(2) = 0.
      CALL MPMMPC (A, A(L1), N5, S(K0))
      IF (NN .EQ. 0) THEN
        CALL MPMMPC (F1, F2, L, B)
        NW = NWS
        ICS = ISS
        GOTO 120
      ELSEIF (NN .EQ. 1) THEN
        CALL MPCEQ (N5, S(K0), S(K2))
        GOTO 110
      ELSEIF (NN .EQ. 2) THEN
        CALL MPCMUL (N5, S(K0), S(K0), S(K2))
        GOTO 110
      ENDIF
C
C   Determine the least integer MN such that 2 ^ MN .GT. NN.
C
      T1 = NN
      MN = CL2 * LOG (T1) + 1.D0 + RXX
      CALL MPMMPC (F1, F2, N5, S(K2))
      KN = NN
C
C   Compute B ^ N using the binary rule for exponentiation.
C
      DO 100 J = 1, MN
        KK = KN / 2
        IF (KN .NE. 2 * KK) THEN
          CALL MPCMUL (N5, S(K2), S(K0), S(K1))
          CALL MPCEQ (N5, S(K1), S(K2))
        ENDIF
        KN = KK
        IF (J .LT. MN) THEN
          CALL MPCMUL (N5, S(K0), S(K0), S(K1))
          CALL MPCEQ (N5, S(K1), S(K0))
        ENDIF
 100  CONTINUE
C
C   Compute reciprocal if N is negative.
C
 110  IF (N .LT. 0) THEN
        CALL MPMMPC (F1, F2, N5, S(K1))
        CALL MPCDIV (N5, S(K1), S(K2), S(K0))
        CALL MPCEQ (N5, S(K0), S(K2))
      ENDIF
      CALL MPMMPC (S(K2), S(N5+K2), L, B)
C
C   Restore original precision level.
C
      NW = NWS
      ICS = ISS
      CALL MPROUN (B)
      CALL MPROUN (B(L1))
C
 120  IF (IDB .GE. 7) THEN
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 4) (B(I), I = 1, NO)
 4      FORMAT ('MPCPWR O'/(6F12.0))
        NO = MIN (INT (ABS (B(L1))), NDB) + 2
        WRITE (LDB, 4) (B(L+I), I = 1, NO)
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPCPWX (L, A, N, B)
C
C   This computes the N-th power of the MPC number A and returns the MPC
C   result C in B.  When N is zero, 1 is returned.  When N is negative, the
C   reciprocal of A ^ |N| is returned.  L is the offset between real and
C   imaginary parts in A and B.  L should be at least NW + 4.  Before calling
C   MPCPWX, the array in MPCOM5 must be initialized by calling MPINIX.  For
C   modest levels of precision, use MPCPWR.  NW should be a power of two.
C   The last two words of the result are not reliable.  Debug output starts
C   with IDB = 6.
C
C   Max SP space for B: 2 * L cells.  Max SP scratch space: 14 * NW + 56
C   cells.  Max DP scratch space: 12 * NW + 6 cells.
C
C   This routine employs the binary method for exponentiation.
C
      DOUBLE PRECISION CL2, T1
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (CL2 = 1.4426950408889633D0, NDB = 22)
      DIMENSION A(2*L), B(2*L), F1(8), F2(8)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        B(1) = 0.
        B(2) = 0.
        B(L+1) = 0.
        B(L+2) = 0.
        RETURN
      ENDIF
      L1 = L + 1
      IF (IDB .GE. 6) THEN
        WRITE (6, 1) L, N
 1      FORMAT ('MPCPWX I',2I10)
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 2) (A(I), I = 1, NO)
 2      FORMAT ('MPCPWX I'/(6F12.0))
        NO = MIN (INT (ABS (A(L1))), NDB) + 2
        WRITE (LDB, 2) (A(L+I), I = 1, NO)
      ENDIF
C
      NA1 = MIN (INT (ABS (A(1))), NW)
      NA2 = MIN (INT (ABS (A(L1))), NW)
      NCR = 2 ** MCR
C
C   Check if precision level of A is too low to justify advanced routine.
C
      IF (NA1 .LE. NCR .AND. NA2 .LE. NCR) THEN
        CALL MPCPWR (L, A, N, B)
        GOTO 120
      ENDIF
      IF (NA1 .EQ. 0 .AND. NA2 .EQ. 0) THEN
        IF (N .GE. 0) THEN
          B(1) = 0.
          B(2) = 0.
          B(L1) = 0.
          B(L1+1) = 0.
          GOTO 120
        ELSE
          IF (KER(26) .NE. 0) THEN
            WRITE (LDB, 3)
 3          FORMAT ('*** MPCPWX: Argument is zero and N is negative or',
     $        ' zero.')
            IER = 26
            IF (KER(IER) .EQ. 2) CALL MPABRT
          ENDIF
          RETURN
        ENDIF
      ENDIF
C
      N4 = NW + 4
      NS = 6 * N4
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + 2 * N4
      K2 = K1 + 2 * N4
      NN = ABS (N)
      F1(1) = 1.
      F1(2) = 0.
      F1(3) = 1.
      F2(1) = 0.
      F2(2) = 0.
      CALL MPMMPC (A, A(L1), N4, S(K0))
      IF (NN .EQ. 0) THEN
        CALL MPMMPC (F1, F2, L, B)
        ICS = ISS
        GOTO 120
      ELSEIF (NN .EQ. 1) THEN
        CALL MPCEQ (N4, S(K0), S(K2))
        GOTO 110
      ELSEIF (NN .EQ. 2) THEN
        CALL MPCMLX (N4, S(K0), S(K0), S(K2))
        GOTO 110
      ENDIF
C
C   Determine the least integer MN such that 2 ^ MN .GT. NN.
C
      T1 = NN
      MN = CL2 * LOG (T1) + 1.D0 + RXX
      CALL MPMMPC (F1, F2, N4, S(K2))
      KN = NN
C
C   Compute B ^ N using the binary rule for exponentiation.
C
      DO 100 J = 1, MN
        KK = KN / 2
        IF (KN .NE. 2 * KK) THEN
          CALL MPCMLX (N4, S(K2), S(K0), S(K1))
          CALL MPCEQ (N4, S(K1), S(K2))
        ENDIF
        KN = KK
        IF (J .LT. MN) THEN
          CALL MPCMLX (N4, S(K0), S(K0), S(K1))
          CALL MPCEQ (N4, S(K1), S(K0))
        ENDIF
 100  CONTINUE
C
C   Compute reciprocal if N is negative.
C
 110  IF (N .LT. 0) THEN
        CALL MPMMPC (F1, F2, N4, S(K1))
        CALL MPCDVX (N4, S(K1), S(K2), S(K0))
        CALL MPCEQ (N4, S(K0), S(K2))
      ENDIF
      CALL MPMMPC (S(K2), S(N4+K2), L, B)
      ICS = ISS
C
 120  IF (IDB .GE. 6) THEN
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 4) (B(I), I = 1, NO)
 4      FORMAT ('MPCPWX O'/(6F12.0))
        NO = MIN (INT (ABS (B(L1))), NDB) + 2
        WRITE (LDB, 4) (B(L+I), I = 1, NO)
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPCRFT (IS, M, X, Y)
C
C   This performs an N-point complex-to-real FFT, where N = 2^M.  X and Y
C   are double precision arrays.  X is both the input and the output data
C   array, and Y is a scratch array.  N/2 + 1 complex pairs are input, with
C   real and imaginary parts separated by N/2 + 1 locations, and N real
C   values are output .  A call to MPCRFT with IS = 1 (or -1) indicates a call
C   to perform a complex-to-real FFT with positive (or negative) exponentials.
C   M must be at least three.  The arrays X and Y must be dimensioned with
C   N + 2 cells.  Before calling MPCRFT, the U array in MPCOM5 must be
C   initialized by calling MPINIX.
C
C   In this application, MPCRFT is called by MPMULX.  This routine is not
C   intended to be called directly by the user.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION X(*), Y(*)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM5/ U(1024)
C
C   Set initial parameters.
C
      K = U(1)
      MX = MOD (K, 64)
      NU = K / 64
      N = 2 ** M
      N2 = N / 2
      N21 = N2 + 1
      N4 = N / 4
      KU = N / 2
      KN = KU + NU
C
C   Check if input parameters are invalid.
C
      IF ((IS .NE. 1 .AND. IS .NE. -1) .OR. M .LT. 3 .OR. M .GT. MX)    
     $  THEN
        IF (KER(27) .NE. 0) THEN
          WRITE (LDB, 1)  IS, M, MX
 1        FORMAT ('*** MPCRFT: Either U has not been initialized'/      
     $      'or else one of the input parameters is invalid', 3I5)
          IER = 27
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
C   Construct the input to MPCFFT.
C
      Y(1) = 0.5D0 * (X(1) + X(N21))
      Y(N2+1) = 0.5D0 * (X(1) - X(N21))
      Y(N4+1) = X(N4+1)
      Y(N4+N2+1) = -IS * X(N4+N2+2)
C
CDIR$ IVDEP
      DO 100 K = 2, N4
        X11 = X(K)
        X12 = X(K+N21)
        X21 = X(N2+2-K)
        X22 = X(N+3-K)
        A1 = X11 + X21
        A2 = X11 - X21
        B1 = X12 + X22
        B2 = X12 - X22
        U1 = U(K+KU)
        U2 = IS * U(K+KN)
        T1 = U1 * B1 + U2 * A2
        T2 = U1 * A2 - U2 * B1
        Y(K) = 0.5D0 * (A1 - T1)
        Y(K+N2) = 0.5D0 * (B2 + T2)
        Y(N2+2-K) = 0.5D0 * (A1 + T1)
        Y(N+2-K) = 0.5D0 * (-B2 + T2)
 100  CONTINUE
C
C   Perform a normal N/2-point FFT on Y.
C
      CALL MPCFFT (IS, M - 1, Y, X)
C
C   Copy Y to X such that Y(k) = X(2k-1) + i X(2k).
C
CDIR$ IVDEP
      DO 110 K = 1, N2
        X(2*K-1) = Y(K)
        X(2*K) = Y(K+N2)
 110  CONTINUE
C
      RETURN
      END
C
      SUBROUTINE MPCSHX (A, PI, AL2, X, Y)
C
C   This computes the hyperbolic cosine and sine of the MP number A and
C   returns the two MP results in X and Y, respectively.  PI is the MP value
C   of Pi computed by a previous call to MPPI or MPPIX.  AL2 is the MP value
C   of Log (10) computed by a previous call to MPLOG or MPLOGX.  Before
C   calling MPCSHX, the array in MPCOM5 must be initialized by calling MPINIX.
C   For modest levels of precision, use MPCSSH.  NW should be a power of two.
C   The last four words of the result are not reliable.  Debug output starts
C   with IDB = 5.
C
C   Max SP space for X and Y: NW + 4 cells.  Max SP scratch space:16.5*NW + 75
C   cells.  Max DP scratch space: 12 * NX + 6 cells.
C
      DIMENSION A(NW+2), F(8), AL2(NW+2), PI(NW+2), X(NW+4), Y(NW+4)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        X(1) = 0.
        X(2) = 0.
        Y(1) = 0.
        Y(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 5) CALL MPDEB ('MPCSHX I', A)
C
      N4 = NW + 4
      NS = 3 * N4
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N4
      K2 = K1 + N4
      F(1) = 1.
      F(2) = 0.
      F(3) = 1.
C
      CALL MPEXPX (A, PI, AL2, S(K0))
      CALL MPDIVX (F, S(K0), S(K1))
      CALL MPADD (S(K0), S(K1), S(K2))
      CALL MPMULD (S(K2), 0.5D0, 0, X)
      CALL MPSUB (S(K0), S(K1), S(K2))
      CALL MPMULD (S(K2), 0.5D0, 0, Y)
      ICS = ISS
C
      IF (IDB .GE. 5) THEN
        CALL MPDEB ('MPCSHX O', X)
        CALL MPDEB ('MPCSHX O', Y)
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPCSQR (L, A, B)
C
C   This routine computes the complex square root of the MPC number C.  L is
C   the offset between real and imaginary parts in A and B.  L must be at
C   least NW + 4.  For extra high levels of precision, use MPCSQX.  The last
C   word is not reliable.  Debug output starts with IDB = 6.
C
C   Max SP space for B: 2 * L cells.  Max SP scratch space: 6 * NW + 27
C   cells.  Max DP scratch space: NW + 5 cells.
C
C   This routine uses the following formula, where A1 and A2 are the real and
C   imaginary parts of A, and where R = Sqrt [A1 ^ 2 + A2 ^2]:
C
C      B = Sqrt [(R + A1) / 2] + I Sqrt [(R - A1) / 2]
C
C   If the imaginary part of A is < 0, then the imaginary part of B is also
C   set to be < 0.
C
      PARAMETER (NDB = 22)
      DIMENSION A(2*L), B(2*L)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        B(1) = 0.
        B(2) = 0.
        B(L+1) = 0.
        B(L+2) = 0.
        RETURN
      ENDIF
      L1 = L + 1
      IF (IDB .GE. 6) THEN
        WRITE (LDB, 1) L
 1      FORMAT ('MPCSQR I',I10)
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 2) (A(I), I = 1, NO)
 2      FORMAT ('MPCSQR I'/(6F12.0))
        NO = MIN (INT (ABS (A(L1))), NDB) + 2
        WRITE (LDB, 2) (A(L+I), I = 1, NO)
      ENDIF
C
      IF (A(1) .EQ. 0. .AND. A(L+1) .EQ. 0.) THEN
        B(1) = 0.
        B(2) = 0.
        B(L+1) = 0.
        B(L+2) = 0.
        GOTO 100
      ENDIF
C
      N4 = NW + 4
      NS = 3 * N4
      ISS = ICS
      ICS = ICS + NS
      IF (ICS - 1 .GT. IMS) CALL MPALER
      IHS = MAX (ICS, IHS)
      K0 = ISS
      K1 = K0 + N4
      K2 = K1 + N4
C
      CALL MPMUL (A, A, S(K0))
      CALL MPMUL (A(L1), A(L1), S(K1))
      CALL MPADD (S(K0), S(K1), S(K2))
      CALL MPSQRT (S(K2), S(K0))
      CALL MPEQ (A, S(K1))
      S(K1) = ABS (S(K1))
      CALL MPADD (S(K0), S(K1), S(K2))
      CALL MPMULD (S(K2), 0.5D0, 0, S(K1))
      CALL MPSQRT (S(K1), S(K0))
      CALL MPMULD (S(K0), 2.D0, 0, S(K1))
      IF (A(1) .GE. 0.) THEN
        CALL MPEQ (S(K0), B)
        CALL MPDIV (A(L1), S(K1), B(L1))
      ELSE
        CALL MPDIV (A(L1), S(K1), B)
        B(1) = ABS (B(1))
        CALL MPEQ (S(K0), B(L1))
        B(L1) = SIGN (B(L1), A(L1))
      ENDIF
      ICS = ISS
C
 100  IF (IDB .GE. 6) THEN
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 3) (B(I), I = 1, NO)
 3      FORMAT ('MPCSQR O'/(6F12.0))
        NO = MIN (INT (ABS (B(L1))), NDB) + 2
        WRITE (LDB, 3) (B(L+I), I = 1, NO)
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPCSQX (L, A, B)
C
C   This routine computes the complex square root of the MPC number C.  L is
C   the offset between real and imaginary parts in A and B.  L must be at
C   least NW + 4.  For modest levels of precision, use MPCSQR.  The last two
C   words are not reliable.  Debug output starts with IDB = 5.
C
C   Max SP space for B: 2 * L cells.  Max SP scratch space: 7.5 * NW + 39
C   cells.  Max DP scratch space: 12 * NW + 6 cells.
C
C   This routine uses the same algorithm as MPCSQR.
C
      PARAMETER (NDB = 22)
      DIMENSION A(2*L), B(2*L)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        B(1) = 0.
        B(2) = 0.
        B(L+1) = 0.
        B(L+2) = 0.
        RETURN
      ENDIF
      L1 = L + 1
      IF (IDB .GE. 5) THEN
        WRITE (LDB, 1) L
 1      FORMAT ('MPCSQX I',I10)
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 2) (A(I), I = 1, NO)
 2      FORMAT ('MPCSQX I'/(6F12.0))
        NO = MIN (INT (ABS (A(L1))), NDB) + 2
        WRITE (LDB, 2) (A(L+I), I = 1, NO)
      ENDIF
C
      IF (A(1) .EQ. 0. .AND. A(L+1) .EQ. 0.) THEN
        B(1) = 0.
        B(2) = 0.
        B(L+1) = 0.
        B(L+2) = 0.
        GOTO 100
      ENDIF
C
      N4 = NW + 4
      NS = 3 * N4
      ISS = ICS
      ICS = ICS + NS
      IF (ICS - 1 .GT. IMS) CALL MPALER
      IHS = MAX (ICS, IHS)
      K0 = ISS
      K1 = K0 + N4
      K2 = K1 + N4
C
      CALL MPSQX (A, S(K0))
      CALL MPSQX (A(L1), S(K1))
      CALL MPADD (S(K0), S(K1), S(K2))
      CALL MPSQRX (S(K2), S(K0))
      CALL MPEQ (A, S(K1))
      S(K1) = ABS (S(K1))
      CALL MPADD (S(K0), S(K1), S(K2))
      CALL MPMULD (S(K2), 0.5D0, 0, S(K1))
      CALL MPSQRX (S(K1), S(K0))
      CALL MPMULD (S(K0), 2.D0, 0, S(K1))
      IF (A(1) .GE. 0.) THEN
        CALL MPEQ (S(K0), B)
        CALL MPDIVX (A(L1), S(K1), B(L1))
      ELSE
        CALL MPDIVX (A(L1), S(K1), B)
        B(1) = ABS (B(1))
        CALL MPEQ (S(K0), B(L1))
        B(L1) = SIGN (B(L1), A(L1))
      ENDIF
      ICS = ISS
C
 100  IF (IDB .GE. 5) THEN
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 3) (B(I), I = 1, NO)
 3      FORMAT ('MPCSQX O'/(6F12.0))
        NO = MIN (INT (ABS (B(L1))), NDB) + 2
        WRITE (LDB, 3) (B(L+I), I = 1, NO)
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPCSSH (A, AL2, X, Y)
C
C   This computes the hyperbolic cosine and sine of the MP number A and
C   returns the two MP results in X and Y, respectively.  AL2 is the MP value
C   of Log (10) computed by a previous call to MPLOG.  For extra high levels of
C   precision, use MPCSHX.  The last word of the result is not reliable.
C   Debug output starts with IDB = 5.
C
C   Max SP space for X and Y: NW + 4 cells.  Max SP scratch space: 9 * NW + 50
C   cells.  Max DP scratch space: NW + 6 cells.
C
      DIMENSION A(NW+2), F(8), AL2(NW+2), X(NW+4), Y(NW+4)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        X(1) = 0.
        X(2) = 0.
        Y(1) = 0.
        Y(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 5) CALL MPDEB ('MPCSSH I', A)
C
      N5 = NW + 5
      NS = 4 * N5
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N5
      K2 = K1 + N5
      K3 = K2 + N5
      NWS = NW
      NW = NW + 1
      F(1) = 1.
      F(2) = 0.
      F(3) = 1.
C
      CALL MPEXP (A, AL2, S(K0))
      CALL MPDIV (F, S(K0), S(K1))
      CALL MPADD (S(K0), S(K1), S(K2))
      CALL MPMULD (S(K2), 0.5D0, 0, S(K3))
      CALL MPEQ (S(K3), X)
      CALL MPSUB (S(K0), S(K1), S(K2))
      CALL MPMULD (S(K2), 0.5D0, 0, S(K3))
      CALL MPEQ (S(K3), Y)
C
C   Restore original precision level.
C
      NW = NWS
      ICS = ISS
      CALL MPROUN (X)
      CALL MPROUN (Y)
C
      IF (IDB .GE. 5) THEN
        CALL MPDEB ('MPCSSH O', X)
        CALL MPDEB ('MPCSSH O', Y)
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPCSSN (A, PI, X, Y)
C
C   This computes the cosine and sine of the MP number A and returns the two MP
C   results in X and Y, respectively.  PI is the MP value of Pi computed by a
C   previous call to MPPI.  For extra high levels of precision, use MPCSSX.
C   The last word of the result is not reliable.  Debug output starts with
C   IDB = 6.
C
C   Max SP space for X and Y: NW + 4 cells.  Max SP scratch space: 10 * NW + 53
C   cells.  Max DP scratch space: NW + 6 cells.
C
C   This routine uses the conventional Taylor's series for Sin (s):
C
C   Sin (s) =  s - s^3 / 3! + s^5 / 5! - s^7 / 7! ...
C
C   where s = t - a * pi / 2 - b * pi / 16 and the integers a and b are chosen
C   to minimize the absolute value of s.  We can then compute
C
C   Sin (t) = Sin (s + a * pi / 2 + b * pi / 16)
C   Cos (t) = Cos (s + a * pi / 2 + b * pi / 16)
C
C   by applying elementary trig identities for sums.  The sine and cosine of
C   b * pi / 16 are of the form 1/2 * Sqrt {2 +- Sqrt [2 +- Sqrt(2)]}.
C   Reducing t in this manner insures that -Pi / 32 < s <= Pi / 32, which
C   accelerates convergence in the above series.
C
      DOUBLE PRECISION CPI, T1, T2
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (CPI = 3.141592653589793D0)
      DIMENSION A(NW+2), F(8), PI(NW+2), X(NW+4), Y(NW+4)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        X(1) = 0.
        X(2) = 0.
        Y(1) = 0.
        Y(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 6) CALL MPDEB ('MPCSSN I', A)
C
      IA = SIGN (1., A(1))
      NA = MIN (INT (ABS (A(1))), NW)
      IF (NA .EQ. 0) THEN
        X(1) = 1.
        X(2) = 0.
        X(3) = 1.
        Y(1) = 0.
        Y(2) = 0.
        L1 = 0
        GOTO 120
      ENDIF
C
C   Check if Pi has been precomputed.
C
      CALL MPMDC (PI, T1, N1)
      IF (N1 .NE. 0 .OR. ABS (T1 - CPI) .GT. RX2) THEN
        IF (KER(28) .NE. 0) THEN
          WRITE (LDB, 1)
 1        FORMAT ('*** MPCSSN: PI must be precomputed.')
          IER = 28
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
      N5 = NW + 5
      NS = 7 * N5
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N5
      K2 = K1 + N5
      K3 = K2 + N5
      K4 = K3 + N5
      K5 = K4 + N5
      K6 = K5 + N5
      NWS = NW
      NW = NW + 1
      F(1) = 1.
      F(2) = 0.
      F(3) = 1.
C
C   Reduce to between - Pi and Pi.
C
      CALL MPMULD (PI, 2.D0, 0, S(K0))
      CALL MPDIV (A, S(K0), S(K1))
      CALL MPNINT (S(K1), S(K2))
      CALL MPSUB (S(K1), S(K2), S(K3))
C
C   Determine nearest multiple of Pi / 2, and within a quadrant, the nearest
C   multiple of Pi / 16.  Through most of the rest of this subroutine, KA and
C   KB are the integers a and b of the algorithm above.
C
      CALL MPMDC (S(K3), T1, N1)
      IF (N1 .GE. - NBT) THEN
        T1 = T1 * 2.D0 ** N1
        T2 = 4.D0 * T1
        KA = NINT (T2)
        KB = NINT (8.D0 * (T2 - KA))
      ELSE
        KA = 0
        KB = 0
      ENDIF
      T1 = (8 * KA + KB) / 32.D0
      CALL MPDMC (T1, 0, S(K1))
      CALL MPSUB (S(K3), S(K1), S(K2))
      CALL MPMUL (S(K0), S(K2), S(K1))
C
C   Compute cosine and sine of the reduced argument s using Taylor's series.
C
      IF (S(K1) .EQ. 0.) THEN
        S(K0) = 0.
        S(K0+1) = 0.
        L1 = 0
        GOTO 110
      ENDIF
      CALL MPEQ (S(K1), S(K0))
      CALL MPMUL (S(K0), S(K0), S(K2))
      L1 = 0
C
 100  L1 = L1 + 1
      IF (L1 .EQ. 10000) THEN
        IF (KER(29) .NE. 0) THEN
          WRITE (LDB, 2)
 2        FORMAT ('*** MPCSSN: Iteration limit exceeded.')
          IER = 29
          IF (KER(IER) .EQ. 2) CALL MPABRT
          ICS = ISS
          NW = NWS
          RETURN
        ENDIF
      ENDIF
C
      T2 = - (2.D0 * L1) * (2.D0 * L1 + 1.D0)
      CALL MPMUL (S(K2), S(K1), S(K3))
      CALL MPDIVD (S(K3), T2, 0, S(K1))
      CALL MPADD (S(K1), S(K0), S(K3))
      CALL MPEQ (S(K3), S(K0))
C
C   Check for convergence of the series.
C
      IF (S(K1) .NE. 0. .AND. S(K1+1) .GE. S(K0+1) - NW) GOTO 100
C
C   Compute Cos (s) = Sqrt [1 - Sin^2 (s)].
C
 110  CALL MPEQ (S(K0), S(K1))
      CALL MPMUL (S(K0), S(K0), S(K2))
      CALL MPSUB (F, S(K2), S(K3))
      CALL MPSQRT (S(K3), S(K0))
C
C   Compute cosine and sine of b * Pi / 16.
C
      KC = ABS (KB)
      F(3) = 2.
      IF (KC .EQ. 0) THEN
        S(K2) = 1.
        S(K2+1) = 0.
        S(K2+2) = 1.
        S(K3) = 0.
        S(K3+1) = 0.
      ELSE
        IF (KC .EQ. 1) THEN
          CALL MPSQRT (F, S(K4))
          CALL MPADD (F, S(K4), S(K5))
          CALL MPSQRT (S(K5), S(K4))
        ELSEIF (KC .EQ. 2) THEN
          CALL MPSQRT (F, S(K4))
        ELSEIF (KC .EQ. 3) THEN
          CALL MPSQRT (F, S(K4))
          CALL MPSUB (F, S(K4), S(K5))
          CALL MPSQRT (S(K5), S(K4))
        ELSEIF (KC .EQ. 4) THEN
          S(K4) = 0.
          S(K4+1) = 0.
        ENDIF
        CALL MPADD (F, S(K4), S(K5))
        CALL MPSQRT (S(K5), S(K3))
        CALL MPMULD (S(K3), 0.5D0, 0, S(K2))
        CALL MPSUB (F, S(K4), S(K5))
        CALL MPSQRT (S(K5), S(K4))
        CALL MPMULD (S(K4), 0.5D0, 0, S(K3))
      ENDIF
      IF (KB .LT. 0) S(K3) = - S(K3)
C
C   Apply the trigonometric summation identities to compute cosine and sine of
C   s + b * Pi / 16.
C
      CALL MPMUL (S(K0), S(K2), S(K4))
      CALL MPMUL (S(K1), S(K3), S(K5))
      CALL MPSUB (S(K4), S(K5), S(K6))
      CALL MPMUL (S(K1), S(K2), S(K4))
      CALL MPMUL (S(K0), S(K3), S(K5))
      CALL MPADD (S(K4), S(K5), S(K1))
      CALL MPEQ (S(K6), S(K0))
C
C   This code in effect applies the trigonometric summation identities for
C   (s + b * Pi / 16) + a * Pi / 2.
C
      IF (KA .EQ. 0) THEN
        CALL MPEQ (S(K0), X)
        CALL MPEQ (S(K1), Y)
      ELSEIF (KA .EQ. 1) THEN
        CALL MPEQ (S(K1), X)
        X(1) = - X(1)
        CALL MPEQ (S(K0), Y)
      ELSEIF (KA .EQ. -1) THEN
        CALL MPEQ (S(K1), X)
        CALL MPEQ (S(K0), Y)
        Y(1) = - Y(1)
      ELSEIF (KA .EQ. 2 .OR. KA .EQ. -2) THEN
        CALL MPEQ (S(K0), X)
        X(1) = - X(1)
        CALL MPEQ (S(K1), Y)
        Y(1) = - Y(1)
      ENDIF
C
C   Restore original precision level.
C
      NW = NWS
      ICS = ISS
      CALL MPROUN (X)
      CALL MPROUN (Y)
C
 120  IF (IDB .GE. 6) THEN
        WRITE (LDB, 3) L1
 3      FORMAT ('Iteration count:',I5)
        CALL MPDEB ('MPCSSN O', X)
        CALL MPDEB ('MPCSSN O', Y)
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPCSSX (A, PI, X, Y)
C
C   This computes the cosine and sine of the MP number A and returns the two MP
C   results in X and Y, respectively.  PI is the MP value of Pi computed by a
C   previous call to MPPI or MPPIX.  Before calling MPCSSX, the array in
C   MPCOM5 must be initialized by calling MPINIX.  For modest levels of
C   precision, use MPCSSN.  NW should be a power of two.  The last four words
C   of the result are not reliable.  Debug output starts with IDB = 5.
C
C   Max SP space for X and Y: NW + 4 cells.  Max SP scratch space:27.5*NW+119
C   cells.  Max DP scratch space: 12 * NW + 6 cells.
C
C   This routine employs a complex arithmetic version of the scheme used in
C   MPEXPX.
C
      DOUBLE PRECISION CL2, CPI, T1, T2
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (CL2 = 1.4426950408889633D0, CPI = 3.141592653589793D0, 
     $  NIT = 1)
      DIMENSION A(NW+2), F1(8), PI(NW+2), X(NW+4), Y(NW+4)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        X(1) = 0.
        X(2) = 0.
        Y(1) = 0.
        Y(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 5) CALL MPDEB ('MPCSSX I', A)
C
      IA = SIGN (1., A(1))
      NA = MIN (INT (ABS (A(1))), NW)
      NCR = 2 ** MCR
C
C   Check if precision level is too low to justify advanced routine.
C
      IF (NW .LE. NCR) THEN
        CALL MPCSSN (A, PI, X, Y)
        L1 = 0
        GOTO 120
      ENDIF
C
C   Check if input is zero.
C
      IF (NA .EQ. 0) THEN
        X(1) = 1.
        X(2) = 0.
        X(3) = 1.
        Y(1) = 0.
        Y(2) = 0.
        L1 = 0
        GOTO 120
      ENDIF
C
C   Check if Pi has been precomputed.
C
      CALL MPMDC (PI, T1, N1)
      IF (N1 .NE. 0 .OR. ABS (T1 - CPI) .GT. RX2) THEN
        IF (KER(30) .NE. 0) THEN
          WRITE (LDB, 1)
 1        FORMAT ('*** MPCSSX: PI must be precomputed.')
          IER = 30
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
      N4 = NW + 4
      N42 = 2 * N4
      NS = 4 * N42
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N42
      K2 = K1 + N42
      K3 = K2 + N42
      F1(1) = 1.
      F1(2) = 0.
      F1(3) = 1.
      NWS = NW
C
C   Reduce argument to between - Pi and Pi.
C
      CALL MPMULD (PI, 2.D0, 0, S(K0))
      CALL MPDIVX (A, S(K0), S(K1))
      CALL MPNINT (S(K1), S(K2))
      CALL MPMUL (S(K2), S(K0), S(K1))
      CALL MPSUB (A, S(K1), S(K0))
C
C   Determine the least integer MQ such that 2 ^ MQ .GE. NW.
C
      T2 = NWS
      MQ = CL2 * LOG (T2) + 1.D0 - RXX
      CALL MPEQ (F1, S(K2))
C
C   Compute initial approximation to [Cos (A), Sin (A)].
C
      NW = NCR
      CALL MPCSSN (S(K0), PI, S(K3), S(K3+N4))
      IQ = 0
C
C   Perform the Newton-Raphson iteration with a dynamically changing precision
C   level NW.
C
      DO 110 K = MCR + 1, MQ
        NW = MIN (2 * NW, NWS)
 100    CONTINUE
        CALL MPANGX (S(K3), S(K3+N4), PI, S(K1))
        CALL MPSUB (S(K0), S(K1), S(K2+N4))
        CALL MPCMLX (N4, S(K3), S(K2), S(K1))
        CALL MPCEQ (N4, S(K1), S(K3))
        IF (K .EQ. MQ - NIT .AND. IQ .EQ. 0) THEN
          IQ = 1
          GOTO 100
        ENDIF
 110  CONTINUE
C
C   The final (cos, sin) result must be normalized to have magnitude 1.
C
      CALL MPSQX (S(K3), S(K0))
      CALL MPSQX (S(K3+N4), S(K0+N4))
      CALL MPADD (S(K0), S(K0+N4), S(K1))
      CALL MPSQRX (S(K1), S(K2))
      CALL MPDIVX (S(K3), S(K2), S(K0))
      CALL MPDIVX (S(K3+N4), S(K2), S(K0+N4))
      CALL MPMPCM (N4, S(K0), X, Y)
      ICS = ISS
C
 120  IF (IDB .GE. 5) THEN
        CALL MPDEB ('MPCSSX O', X)
        CALL MPDEB ('MPCSSX O', Y)
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPCSUB (L, A, B, C)
C
C   This subracts the MPC numbers A and B and returns the MPC difference in
C   C.  L is the offset between real and imaginary parts in A, B and C.  L
C   must be at least NW + 4.  Debug output starts with IDB = 9.
C
C   Max SP space for C: 2 * L cells.
C
      DIMENSION A(2*L), B(2*L), C(2*L)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
C
      IF (IER .NE. 0) THEN
        C(1) = 0.
        C(2) = 0.
        C(L+1) = 0.
        C(L+2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 9) WRITE (LDB, 1)
 1    FORMAT ('MPCSUB')
C
      L1 = L + 1
      CALL MPSUB (A, B, C)
      CALL MPSUB (A(L1), B(L1), C(L1))
C
      RETURN
      END
C
      SUBROUTINE MPDEB (CS, A)
C
C   This outputs the character string CS, the exponent of the MP number A, and
C   the first 50 digits of A, all on one line.  CS must either be a literal
C   string not exceeding 12 characters in length or a variable of type
C   CHARACTER*n, where n does not exceed 12.
C
      CHARACTER*(*) CS
      CHARACTER*1 B(160)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      DIMENSION A(NW+2)
C
      IF (IER .NE. 0) RETURN
      IDS = IDB
      IDB = 0
      NWS = NW
      NW = MIN (NW, 10)
      CALL MPOUTC (A, B, N)
      N = MIN (N, 70)
      WRITE (LDB, 1) CS, ' ', (B(K), K = 1, 4), (B(K), K = 9, N)
 1    FORMAT (A12,67A1:/(79A1))
      IDB = IDS
      NW = NWS
      RETURN
      END
C
      SUBROUTINE MPDIV (A, B, C)
C
C   This divides the MP number A by the MP number B to yield the MP quotient C.
C   For extra high levels of precision, use MPDIVX.  Debug output starts with
C   IDB = 8.
C
C   Max SP space for C: NW + 4 cells.  Max DP scratch space: NW + 4 cells.
C
      DOUBLE PRECISION D, RB, SS, T0, T1, T2
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (NDB = 22)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM4/ D(1024)
      DIMENSION A(NW+2), B(NW+2), C(NW+4)
C
      IF (IER .NE. 0) THEN
        C(1) = 0.
        C(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 8) THEN
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 1) (A(I), I = 1, NO)
 1      FORMAT ('MPDIV I'/(6F12.0))
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 1) (B(I), I = 1, NO)
      ENDIF
C
      IA = SIGN (1., A(1))
      IB = SIGN (1., B(1))
      NA = MIN (INT (ABS (A(1))), NW)
      NB = MIN (INT (ABS (B(1))), NW)
C
C   Check if dividend is zero.
C
      IF (NA .EQ. 0) THEN
        C(1) = 0.
        C(2) = 0.
        GOTO 210
      ENDIF
      IF (NB .EQ. 1 .AND. B(3) .EQ. 1.) THEN
C
C   Divisor is 1 or -1 -- result is A or -A.
C
        C(1) = SIGN (NA, IA * IB)
        C(2) = A(2) - B(2)
C
        DO 100 I = 3, NA + 2
          C(I) = A(I)
 100    CONTINUE
C
        GOTO 210
      ENDIF
C
C   Check if divisor is zero.
C
      IF (NB .EQ. 0) THEN
        IF (KER(31) .NE. 0) THEN
          WRITE (LDB, 2)
 2        FORMAT ('*** MPDIV: Divisor is zero.')
          IER = 31
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
C   Initialize trial divisor and trial dividend.
C
      T0 = BDX * B(3)
      IF (NB .GE. 2) T0 = T0 + B(4)
      IF (NB .GE. 3) T0 = T0 + RDX * B(5)
      IF (NB .GE. 4) T0 = T0 + RX2 * B(6)
      RB = 1.D0 / T0
      MD = MIN (NA + NB, NW)
      D(1) = 0.D0
C
      DO 110 I = 2, NA + 1
        D(I) = A(I+1)
 110  CONTINUE
C
      DO 120 I = NA + 2, MD + 4
        D(I) = 0.D0
 120  CONTINUE
C
C   Perform ordinary long division algorithm.  First compute only the first
C   NA words of the quotient.
C
      DO 150 J = 2, NA + 1
        T1 = BX2 * D(J-1) + BDX * D(J) + D(J+1) + RDX * D(J+2)
        T0 = INT (RB * T1)
        J3 = J - 3
        I2 = MIN (NB, NW + 2 - J3) + 2
        IJ = I2 + J3
C
        DO 130 I = 3, I2
          I3 = I + J3
          D(I3) = D(I3) - T0 * B(I)
 130    CONTINUE
C
C   Release carries periodically to avoid overflowing the exact integer
C   capacity of double precision floating point words in D.
C
        IF (MOD (J - 1, NPR) .EQ. 0) THEN
CDIR$ IVDEP
          DO 140 I = J + 1, IJ
            T1 = D(I)
            T2 = INT (RDX * T1)
            D(I) = T1 - BDX * T2
            D(I-1) = D(I-1) + T2
 140      CONTINUE
C
        ENDIF
        D(J) = D(J) + BDX * D(J-1)
        D(J-1) = T0
 150  CONTINUE
C
C   Compute additional words of the quotient, as long as the remainder
C   is nonzero.
C
      DO 180 J = NA + 2, NW + 3
        T1 = BX2 * D(J-1) + BDX * D(J) + D(J+1)
        IF (J .LE. NW + 2) T1 = T1 + RDX * D(J+2)
        T0 = INT (RB * T1)
        J3 = J - 3
        I2 = MIN (NB, NW + 2 - J3) + 2
        IJ = I2 + J3
        SS = 0.D0
C
        DO 160 I = 3, I2
          I3 = I + J3
          D(I3) = D(I3) - T0 * B(I)
          SS = SS + ABS (D(I3))
 160    CONTINUE
C
        IF (MOD (J - 1, NPR) .EQ. 0) THEN
CDIR$ IVDEP
          DO 170 I = J + 1, IJ
            T1 = D(I)
            T2 = INT (RDX * T1)
            D(I) = T1 - BDX * T2
            D(I-1) = D(I-1) + T2
 170      CONTINUE
C
        ENDIF
        D(J) = D(J) + BDX * D(J-1)
        D(J-1) = T0
        IF (SS .EQ. 0.D0) GOTO 190
        IF (IJ .LE. NW + 1) D(IJ+3) = 0.D0
 180  CONTINUE
C
C   Set sign and exponent, and fix up result.
C
      J = NW + 3
C
 190  D(J) = 0.D0
      IF (D(1) .EQ. 0.D0) THEN
        IS = 1
      ELSE
        IS = 2
      ENDIF
      NC = MIN (J - 1, NW)
      D(NC+3) = 0.D0
      D(NC+4) = 0.D0
C
      DO 200 I = J + 1, 3, -1
        D(I) = D(I-IS)
 200  CONTINUE
C
      D(1) = SIGN (NC, IA * IB)
      D(2) = A(2) - B(2) + IS - 2
      CALL MPNORM (C)
C
 210  IF (IDB .GE. 8) THEN
        NO = MIN (INT (ABS (C(1))), NDB) + 2
        WRITE (LDB, 3) (C(I), I = 1, NO)
 3      FORMAT ('MPDIV O'/(6F12.0))
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPDIVD (A, B, N, C)
C
C   This routine divides the MP number A by the DPE number (B, N) to yield
C   the MP quotient C.   Debug output starts with IDB = 9.
C
C   Max SP space for C: NW + 4 cells.  Max DP space: NW + 4 cells.
C
      DOUBLE PRECISION B, BB, BR, D, DD, T1
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (NDB = 22)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
      COMMON /MPCOM4/ D(1024)
      DIMENSION A(NW+2), C(NW+4), F(8)
C
      IF (IER .NE. 0) THEN
        C(1) = 0.
        C(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 9) THEN
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 1) (A(I), I = 1, NO)
 1      FORMAT ('MPDIVD I'/(6F12.0))
        WRITE (LDB, 2) B, N
 2      FORMAT ('MPDIVD I',1PD25.15,I10)
      ENDIF
C
      IA = SIGN (1., A(1))
      NA = MIN (INT (ABS (A(1))), NW)
      IB = SIGN (1.D0, B)
C
C   Check if dividend is zero.
C
      IF (NA .EQ. 0) THEN
        C(1) = 0.
        C(2) = 0.
        GOTO 150
      ENDIF
C
C   Check if divisor is zero.
C
      IF (B .EQ. 0.D0) THEN
        IF (KER(32) .NE. 0) THEN
          WRITE (LDB, 3)
 3        FORMAT ('*** MPDIVD: Divisor is zero.')
          IER = 32
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
      N1 = N / NBT
      N2 = N - NBT * N1
      BB = ABS (B) * 2.D0 ** N2
C
C   Reduce BB to within 1 and BDX.
C
      IF (BB .GE. BDX) THEN
C
        DO 100 K = 1, 100
          BB = RDX * BB
          IF (BB .LT. BDX) THEN
            N1 = N1 + K
            GOTO 120
          ENDIF
 100    CONTINUE
C
      ELSEIF (BB .LT. 1.D0) THEN
C
        DO 110 K = 1, 100
          BB = BDX * BB
          IF (BB .GE. 1.D0) THEN
            N1 = N1 - K
            GOTO 120
          ENDIF
 110    CONTINUE
C
      ENDIF
C
C   If B cannot be represented exactly in a single mantissa word, use MPDIV.
C
 120  IF (BB .NE. AINT (BB)) THEN
        BB = SIGN (BB, B)
        CALL MPDMC (BB, N1 * NBT, F)
        CALL MPDIV (A, F, C)
        GOTO 150
      ENDIF
C
      BR = 1.D0 / BB
      DD = A(3)
C
C   Perform short division (not vectorizable at present).  Continue as long as
C   the remainder remains nonzero.
C
      DO 130 J = 2, NW + 3
        T1 = INT (BR * DD)
        D(J+1) = T1
        DD = BDX * (DD - T1 * BB)
        IF (J .LE. NA) THEN
          DD = DD + A(J+2)
        ELSE
          IF (DD .EQ. 0.D0) GOTO 140
        ENDIF
 130  CONTINUE
C
C   Set sign and exponent of result.
C
      J = NW + 3
C
 140  NC = MIN (J - 1, NW)
      D(1) = SIGN (NC, IA * IB)
      D(2) = A(2) - N1
      IF (J .LE. NW + 2) D(J+2) = 0.D0
      IF (J .LE. NW + 1) D(J+3) = 0.D0
      CALL MPNORM (C)
C
 150  IF (IDB .GE. 9) THEN
        NO = MIN (INT (ABS (C(1))), NDB) + 2
        WRITE (LDB, 4) (C(I), I = 1, NO)
 4      FORMAT ('MPDIVD O'/(6F12.0))
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPDIVX (A, B, C)
C
C   This divides the MP number A by the MP number B and returns the MP result
C   in C.  Before calling MPDIVX, the array in MPCOM5 must be initialized by
C   calling MPINIX.  For modest levels of precision, use MPDIV.  NW should be
C   a power of two.  The last two words of the result are not reliable.  Debug
C   output starts with IDB = 7.
C
C   Max SP space for C: NW + 4 cells.  Max SP scratch space: 3 * NW + 12
C   cells.  Max DP scratch space: 12 * NW + 6 cells.
C
C   This subroutine employs the following Newton-Raphson iteration, which
C   converges to 1 / B:
C
C          X_{k+1} = X_k + (1 - X_k * B) * X_k
C
C   where the muliplication () * X_k is performed with only half of the
C   normal level of precision.  These iterations are performed with a
C   maximum precision level NW that is dynamically changed, doubling with
C   each iteration.  The final iteration is performed as follows (this is
C   due to A. Karp):
C
C          A / B = (A * X_n) + [A - (A * X_n) * B] * X_n  (approx.)
C
C   where the multiplications A * X_n and [] * X_n are performed with only
C   half of the final level of precision.
C
C   One difficulty with this procedure is that errors often accumulate in the
C   trailing mantissa words.  This error can be controlled by repeating one of
C   the iterations.  The iteration that is repeated is controlled by setting
C   the parameter NIT below:  If NIT = 0, the last iteration is repeated (this
C   is most effective but most expensive).  If NIT = 1, then the next-to-last
C   iteration is repeated, etc.  An extra word of precision cannot be used in
C   this routine (since MPMULX prefers powers of two), so NIT = 0 or 1 is best
C   unless the user needs maximum speed.
C
      DOUBLE PRECISION CL2, T1
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (CL2 = 1.4426950408889633D0, NDB = 22, NIT = 1)
      DIMENSION A(NW+2), B(NW+2), C(NW+4), F(8)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        C(1) = 0.
        C(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 7) THEN
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 1) (A(I), I = 1, NO)
 1      FORMAT ('MPDIVX I'/(6F12.0))
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 1) (B(I), I = 1, NO)
      ENDIF
C
      IA = SIGN (1., A(1))
      IB = SIGN (1., B(1))
      NA = MIN (INT (ABS (A(1))), NW)
      NB = MIN (INT (ABS (B(1))), NW)
      NCR = 2 ** MCR
C
C   Check if dividend is zero.
C
      IF (NA .EQ. 0) THEN
        C(1) = 0.
        C(2) = 0.
        GOTO 120
      ENDIF
C
C   Check if divisor is zero.
C
      IF (NB .EQ. 0)  THEN
        IF (KER(33) .NE. 0) THEN
          WRITE (LDB, 2)
 2        FORMAT ('*** MPDIVX: Divisor is zero.')
          IER = 33
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
C   Check if precision level of divisor is too low to justify the advanced
C   routine.
C
      IF (NB .LE. NCR) THEN
        CALL MPDIV (A, B, C)
        GOTO 120
      ENDIF
      N4 = NW + 4
      NS = 3 * N4
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N4
      K2 = K1 + N4
      NWS = NW
C
C   Determine the least integer MQ such that 2 ^ MQ .GE. NW.
C
      T1 = NW
      MQ = CL2 * LOG (T1) + 1.D0 - RXX
C
C   Compute the initial approximation of 1 / B to a precision of NCR words.
C
      NW = NCR
      F(1) = 1.
      F(2) = 0.
      F(3) = 1.
      CALL MPDIV (F, B, C)
      IQ = 0
C
C   Perform the Newton-Raphson iterations described above.
C
      DO 110 K = MCR + 1, MQ - 1
        NW1 = NW
        NW = MIN (2 * NW, NWS)
        NW2 = NW
 100    CONTINUE
        CALL MPMULX (B, C, S(K0))
        CALL MPSUB (F, S(K0), S(K1))
        NW = NW1
        CALL MPMULX (C, S(K1), S(K0))
        NW = NW2
        CALL MPADD (C, S(K0), S(K1))
        CALL MPEQ (S(K1), C)
        IF (K .EQ. MQ - NIT .AND. IQ .EQ. 0) THEN
          IQ = 1
          GOTO 100
        ENDIF
 110  CONTINUE
C
C   Perform last iteration using Karp's trick.
C
      CALL MPMULX (A, C, S(K0))
      NW1 = NW
      NW = MIN (2 * NW, NWS)
      NW2 = NW
      CALL MPMULX (S(K0), B, S(K1))
      CALL MPSUB (A, S(K1), S(K2))
      NW = NW1
      CALL MPMULX (S(K2), C, S(K1))
      NW = NW2
      CALL MPADD (S(K0), S(K1), C)
      ICS = ISS
C
 120  IF (IDB .GE. 7) THEN
        NO = MIN (INT (ABS (C(1))), NDB) + 2
        WRITE (LDB, 3) (C(I), I = 1, NO)
 3      FORMAT ('MPDIVX O'/(6F12.0))
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPDMC (A, N, B)
C
C   This routine converts the DPE number (A, N) to MP form in B.  All bits of
C   A are recovered in B.  However, note for example that if A = 0.1D0 and N
C   is 0, then B will NOT be the multiprecision equivalent of 1/10.  Debug
C   output starts with IDB = 9.
C
C   Max SP space for B:  8 cells.
C
      DOUBLE PRECISION A, AA
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (NDB = 22)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      DIMENSION B(NW+4)
C
      IF (IER .NE. 0) THEN
        B(1) = 0.
        B(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 9) WRITE (LDB, 1) A, N
 1    FORMAT ('MPDMC I',1PD25.15,I10)
C
C   Check for zero.
C
      IF (A .EQ. 0.D0) THEN
        B(1) = 0.
        B(2) = 0.
        GOTO 150
      ENDIF
      N1 = N / NBT
      N2 = N - NBT * N1
      AA = ABS (A) * 2.D0 ** N2
C
C   Reduce AA to within 1 and BDX.
C
      IF (AA .GE. BDX) THEN
C
        DO 100 K = 1, 100
          AA = RDX * AA
          IF (AA .LT. BDX) THEN
            N1 = N1 + K
            GOTO 120
          ENDIF
 100    CONTINUE
C
      ELSEIF (AA .LT. 1.D0) THEN
C
        DO 110 K = 1, 100
          AA = BDX * AA
          IF (AA .GE. 1.D0) THEN
            N1 = N1 - K
            GOTO 120
          ENDIF
 110    CONTINUE
C
      ENDIF
C
C   Store successive sections of AA into B.
C
 120  B(2) = N1
      B(3) = AINT (AA)
      AA = BDX * (AA - B(3))
      B(4) = AINT (AA)
      AA = BDX * (AA - B(4))
      B(5) = AINT (AA)
      AA = BDX * (AA - B(5))
      B(6) = AINT (AA)
      B(7) = 0.
      B(8) = 0.
C
      DO 130 I = 6, 3, -1
        IF (B(I) .NE. 0.) GOTO 140
 130  CONTINUE
C
 140  AA = I - 2
      B(1) = SIGN (AA, A)
C
 150  IF (IDB .GE. 9) THEN
        NO = ABS (B(1)) + 2.
        WRITE (LDB, 2) (B(I), I = 1, NO)
 2      FORMAT ('MPDMC O'/(6F12.0))
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPEQ (A, B)
C
C   This routine sets the MP number B equal to the MP number A.  Debug output
C   starts with IDB = 10.
C
C   Max SP space for B: NW + 2 cells.
C
C   The fact that only NW + 2 cells, and not NW + 4 cells, are copied is
C   important in some routines that increase the precision level by one.
C
      PARAMETER (NDB = 22)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      DIMENSION A(NW+2), B(NW+2)
C
      IF (IER .NE. 0) THEN
        B(1) = 0.
        B(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 10) WRITE (LDB, 1)
 1    FORMAT ('MPEQ')
C
      IA = SIGN (1., A(1))
      NA = MIN (INT (ABS (A(1))), NW)
      B(1) = SIGN (NA, IA)
C
      DO 100 I = 2, NA + 2
        B(I) = A(I)
 100  CONTINUE
C
      RETURN
      END
C
      SUBROUTINE MPEXP (A, AL2, B)
C
C   This computes the exponential function of the MP number A and returns the
C   MP result in B.  AL2 is the MP value of Log(2) produced by a prior call
C   to MPLOG.  For extra high levels of precision, use MPEXPX.  The last
C   word of the result is not reliable.  Debug output starts with IDB = 7.
C
C   Max SP space for B: NW + 4 cells.  Max SP scratch space: 5 * NW + 25
C   cells.  Max DP scratch space: NW + 5 cells.
C
C   This routine uses a modification of the Taylor's series for Exp (t):
C
C   Exp (t) =  (1 + r + r^2 / 2! + r^3 / 3! + r^4 / 4! ...) ^ q * 2 ^ n
C
C   where q = 256, r = t' / q, t' = t - n Log(2) and where n is chosen so
C   that -0.5 Log(2) < t' <= 0.5 Log(2).  Reducing t mod Log(2) and
C   dividing by 256 insures that -0.001 < r <= 0.001, which accelerates
C   convergence in the above series.
C
      DOUBLE PRECISION ALT, T1, T2
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (ALT = 0.693147180559945309D0, NQ = 8)
      DIMENSION A(NW+2), B(NW+4), AL2(NW+2), F(8)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        B(1) = 0.
        B(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 7) CALL MPDEB ('MPEXP I', A)
C
      IA = SIGN (1., A(1))
      NA = MIN (INT (ABS (A(1))), NW)
      CALL MPMDC (A, T1, N1)
      T1 = T1 * 2.D0 ** N1
C
C   Unless the argument is near Log (2), Log(2) must be precomputed.  This
C   exception is necessary because MPLOG calls MPEXP to initialize Log (2).
C
      IF (ABS (T1 - ALT) .GT. RDX) THEN
        CALL MPMDC (AL2, T2, N2)
        IF (N2 .NE. - NBT .OR. ABS (T2 * 0.5D0 ** NBT - ALT) .GT. RX2)  
     $    THEN
          IF (KER(34) .NE. 0) THEN
            WRITE (LDB, 1)
 1          FORMAT ('*** MPEXP: LOG (2) must be precomputed.')
            IER = 34
            IF (KER(IER) .EQ. 2) CALL MPABRT
          ENDIF
          RETURN
        ENDIF
      ENDIF
C
C   Check for overflows and underflows.
C
      IF (T1 .GE. 1D9) THEN
        IF (T1 .GT. 0.D0) THEN
          IF (KER(35) .NE. 0) THEN
            WRITE (LDB, 2) T1, N1
 2          FORMAT ('*** MPEXP: Argument is too large',F12.6,' x 10 ^', 
     $        I8)
            IER = 35
            IF (KER(IER) .EQ. 2) CALL MPABRT
          ENDIF
          RETURN
        ELSE
          B(1) = 0.
          B(2) = 0.
          L1 = 0
          GOTO 130
        ENDIF
      ENDIF
C
      N5 = NW + 5
      NS = 4 * N5
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N5
      K2 = K1 + N5
      K3 = K2 + N5
      NWS = NW
      NW = NW + 1
      F(1) = 1.
      F(2) = 0.
      F(3) = 1.
C
C   Compute the reduced argument A' = A - Log(2) * Nint [A / Log(2)].  Save
C   NZ = Nint [A / Log(2)] for correcting the exponent of the final result.
C
      IF (ABS (T1 - ALT) .GT. RDX) THEN
        CALL MPDIV (A, AL2, S(K0))
        CALL MPNINT (S(K0), S(K1))
        CALL MPMDC (S(K1), T1, N1)
        NZ = T1 * 2.D0 ** N1 + SIGN (RXX, T1)
        CALL MPMUL (AL2, S(K1), S(K2))
        CALL MPSUB (A, S(K2), S(K0))
      ELSE
        CALL MPEQ (A, S(K0))
        NZ = 0
      ENDIF
      TL = S(K0+1) - NW
C
C   Check if the reduced argument is zero.
C
      IF (S(K0) .EQ. 0.D0) THEN
        S(K0) = 1.
        S(K0+1) = 0.
        S(K0+2) = 1.
        L1 = 0
        GOTO 120
      ENDIF
C
C   Divide the reduced argument by 2 ^ NQ.
C
      CALL MPDIVD (S(K0), 1.D0, NQ, S(K1))
C
C   Compute Exp using the usual Taylor series.
C
      CALL MPEQ (F, S(K2))
      CALL MPEQ (F, S(K3))
      L1 = 0
C
 100  L1 = L1 + 1
      IF (L1 .EQ. 10000) THEN
        IF (KER(36) .NE. 0) THEN
          WRITE (LDB, 3)
 3        FORMAT ('*** MPEXP: Iteration limit exceeded.')
          IER = 36
          IF (KER(IER) .EQ. 2) CALL MPABRT
          NW = NWS
          ICS = ISS
          RETURN
        ENDIF
      ENDIF
C
      T2 = L1
      CALL MPMUL (S(K2), S(K1), S(K0))
      CALL MPDIVD (S(K0), T2, 0, S(K2))
      CALL MPADD (S(K3), S(K2), S(K0))
      CALL MPEQ (S(K0), S(K3))
C
C   Check for convergence of the series.
C
      IF (S(K2) .NE. 0. .AND. S(K2+1) .GE. TL) GOTO 100
C
C   Raise to the (2 ^ NQ)-th power.
C
      DO 110 I = 1, NQ
        CALL MPMUL (S(K0), S(K0), S(K1))
        CALL MPEQ (S(K1), S(K0))
 110  CONTINUE
C
C  Multiply by 2 ^ NZ.
C
 120  CALL MPMULD (S(K0), 1.D0, NZ, S(K1))
      CALL MPEQ (S(K1), B)
C
C   Restore original precision level.
C
      NW = NWS
      ICS = ISS
      CALL MPROUN (B)
C
 130  IF (IDB .GE. 7) THEN
        WRITE (LDB, 4) L1
 4      FORMAT ('Iteration count:',I5)
        CALL MPDEB ('MPEXP O', B)
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPEXPX (A, PI, AL2, B)
C
C   This computes the exponential function of the MP number A and returns the
C   MP result in B.  PI is the MP value of Pi produced by a prior call to MPPI
C   or MPPIX.  AL2 is the MP value of Log(2) produced by a prior call to
C   MPLOG  or MPLOGX.  Before calling MPEXPX, the array in MPCOM5 must be
C   initialized by calling MPINIX.  NW should be a power of two.  For modest
C   levels of precision, use MPEXP.  The last four words of the result are
C   not reliable.  Debug output starts with IDB = 5.
C
C   Max SP space for B: NW + 4 cells.  Max SP scratch space: 13.5 * NW + 63
C   cells.  Max DP scratch space: 12 * NW + 6 cells.
C
C   This routine uses the Newton iteration
C
C     b_{k+1} = b_k [a + 1 - log b_k]
C
C   with a dynamically changing level of precision.  Logs are performed using
C   MPLOGX.  See the comment about the parameter NIT in MPDIVX.
C
      DOUBLE PRECISION ALT, CL2, CPI, T1, T2
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (ALT = 0.693147180559945309D0,                          
     $  CL2 = 1.4426950408889633D0, CPI = 3.141592653589793238D0,       
     $  NIT = 1)
      DIMENSION A(NW+2), AL2(NW+2), B(NW+4), F1(8), PI(NW+2)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        B(1) = 0.
        B(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 5) CALL MPDEB ('MPEXPX I', A)
C
      NCR = 2 ** MCR
      IA = SIGN (1., A(1))
      NA = MIN (INT (ABS (A(1))), NW)
      CALL MPMDC (A, T1, N1)
      T1 = T1 * 2.D0 ** N1
C
C   Check if precision level is too low to justify the advanced routine.
C
      IF (NW .LE. NCR) THEN
        CALL MPEXP (A, AL2, B)
        GOTO 120
      ENDIF
C
C   Check if Log(2) has been precomputed.
C
      CALL MPMDC (AL2, T2, N2)
      IF (N2 .NE. - NBT .OR. ABS (T2 * 0.5D0 ** NBT - ALT) .GT. RX2)    
     $  THEN
        IF (KER(37) .NE. 0) THEN
          WRITE (LDB, 1)
 1        FORMAT ('*** MPEXPX: LOG (2) must be precomputed.')
          IER = 37
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
C     Check if Pi has been precomputed.
C
      CALL MPMDC (PI, T2, N2)
      IF (N2 .NE. 0 .OR. ABS (T2 - CPI) .GT. RX2) THEN
        IF (KER(38) .NE. 0) THEN
          WRITE (LDB, 2)
 2        FORMAT ('*** MPEXPX: PI must be precomputed.')
          IER = 38
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
C   Check for overflows and underflows.
C
      IF (T1 .GE. 1D9) THEN
        IF (T1 .GT. 0.D0) THEN
          IF (KER(39) .NE. 0) THEN
            WRITE (LDB, 3) T1, N1
 3          FORMAT ('*** MPEXPX: Argument is too large',F12.6,' x 10 ^',
     $        I8)
            IER = 39
            IF (KER(IER) .EQ. 2) CALL MPABRT
          ENDIF
          RETURN
        ELSE
          B(1) = 0.
          B(2) = 0.
          GOTO 120
        ENDIF
      ENDIF
C
      N4 = NW + 4
      NS = 3 * N4
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N4
      K2 = K1 + N4
      NWS = NW
      F1(1) = 1.
      F1(2) = 0.
      F1(3) = 1.
C
C   Determine the least integer MQ such that 2 ^ MQ .GE. NW.
C
      T2 = NWS
      MQ = CL2 * LOG (T2) + 1.D0 - RXX
      CALL MPADD (A, F1, S(K0))
C
C   Compute initial approximation to Exp (A).
C
      NW = NCR
      CALL MPEXP (A, AL2, B)
      IQ = 0
C
C   Perform the Newton-Raphson iteration described above with a dynamically
C   changing precision level NW.
C
      DO 110 K = MCR + 1, MQ
        NW = MIN (2 * NW, NWS)
 100    CONTINUE
        CALL MPLOGX (B, PI, AL2, S(K1))
        CALL MPSUB (S(K0), S(K1), S(K2))
        CALL MPMULX (B, S(K2), S(K1))
        CALL MPEQ (S(K1), B)
        IF (K .EQ. MQ - NIT .AND. IQ .EQ. 0) THEN
          IQ = 1
          GOTO 100
        ENDIF
 110  CONTINUE
C
      ICS = ISS
C
 120  IF (IDB .GE. 6) CALL MPDEB ('MPEXPX O', B)
      RETURN
      END
C
      SUBROUTINE MPFFT1 (IS, L, M, X, Y)
C
C   Performs the L-th iteration of the first variant of the Stockham FFT.
C   This routine is called by MPCFFT.  It is not intended to be called directly
C   by the user.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION X(*), Y(*)
      COMMON /MPCOM5/ U(1024)
C
C   Set initial parameters.
C
      N = 2 ** M
      K = U(1)
      NU = K / 64
      N1 = N / 2
      LK = 2 ** (L - 1)
      LI = 2 ** (M - L)
      LJ = 2 * LI
      KU = LI + 1
      KN = KU + NU
C
      DO 100 K = 0, LK - 1
        I11 = K * LJ + 1
        I12 = I11 + LI
        I21 = K * LI + 1
        I22 = I21 + N1
C
CDIR$ IVDEP
        DO 100 I = 0, LI - 1
          U1 = U(KU+I)
          U2 = IS * U(KN+I)
          X11 = X(I11+I)
          X12 = X(I11+I+N)
          X21 = X(I12+I)
          X22 = X(I12+I+N)
          T1 = X11 - X21
          T2 = X12 - X22
          Y(I21+I) = X11 + X21
          Y(I21+I+N) = X12 + X22
          Y(I22+I) = U1 * T1 - U2 * T2
          Y(I22+I+N) = U1 * T2 + U2 * T1
 100  CONTINUE
C
      RETURN
      END
C
      SUBROUTINE MPFFT2 (IS, L, M, X, Y)
C
C   Performs the L-th iteration of the second variant of the Stockham FFT.
C   This routine is called by MPCFFT.  It is not intended to be called directly
C   by the user.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION X(*), Y(*)
      COMMON /MPCOM5/ U(1024)
C
C   Set initial parameters.
C
      N = 2 ** M
      K = U(1)
      NU = K / 64
      N1 = N / 2
      LK = 2 ** (L - 1)
      LI = 2 ** (M - L)
      LJ = 2 * LK
      KU = LI + 1
C
      DO 100 I = 0, LI - 1
        I11 = I * LK + 1
        I12 = I11 + N1
        I21 = I * LJ + 1
        I22 = I21 + LK
        U1 = U(KU+I)
        U2 = IS * U(KU+I+NU)
C
CDIR$ IVDEP
        DO 100 K = 0, LK - 1
          X11 = X(I11+K)
          X12 = X(I11+K+N)
          X21 = X(I12+K)
          X22 = X(I12+K+N)
          T1 = X11 - X21
          T2 = X12 - X22
          Y(I21+K) = X11 + X21
          Y(I21+K+N) = X12 + X22
          Y(I22+K) = U1 * T1 - U2 * T2
          Y(I22+K+N) = U1 * T2 + U2 * T1
 100  CONTINUE
C
      RETURN
      END
C
      SUBROUTINE MPINFR (A, B, C)
C
C   Sets B to the integer part of the MP number A and sets C equal to the
C   fractional part of A.  Note that if A = -3.3, then B = -3 and C = -0.3.
C   Debug output starts with IDB = 9.
C
C   Max SP space for B and C: NW + 4 cells.
C
      PARAMETER (NDB = 22)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      DIMENSION A(NW+2), B(NW+2), C(NW+2)
C
      IF (IER .NE. 0) THEN
        B(1) = 0.
        B(2) = 0.
        C(1) = 0.
        C(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 9) THEN
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 1) (A(I), I = 1, NO)
 1      FORMAT ('MPINFR I'/(6F12.0))
      ENDIF
C
C   Check if  A  is zero.
C
      IA = SIGN (1., A(1))
      NA = MIN (INT (ABS (A(1))), NW)
      MA = A(2)
      IF (NA .EQ. 0)  THEN
        B(1) = 0.
        B(2) = 0.
        C(1) = 0.
        C(2) = 0.
        GOTO 120
      ENDIF
C
      IF (MA .GE. NW - 1) THEN
        IF (KER(40) .NE. 0) THEN
          WRITE (LDB, 2)
 2        FORMAT ('*** MPINFR: Argument is too large.')
          IER = 40
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
C   Place integer part in  B.
C
      NB = MIN (MAX (MA + 1, 0), NA)
      IF (NB .EQ. 0) THEN
        B(1) = 0.
        B(2) = 0.
      ELSE
        B(1) = SIGN (NB, IA)
        B(2) = MA
        B(NB+3) = 0.
        B(NB+4) = 0.
C
        DO 100 I = 3, NB + 2
          B(I) = A(I)
 100    CONTINUE
C
      ENDIF
C
C   Place fractional part in C.
C
      NC = NA - NB
      IF (NC .LE. 0) THEN
        C(1) = 0.
        C(2) = 0.
      ELSE
        C(1) = SIGN (NC, IA)
        C(2) = MA - NB
        C(NC+3) = 0.
        C(NC+4) = 0.
C
        DO 110 I = 3, NC + 2
          C(I) = A(I+NB)
 110    CONTINUE
C
      ENDIF
C
C   Fix up results.  B may have trailing zeros and C may have leading zeros.
C
      CALL MPROUN (B)
      CALL MPROUN (C)
C
 120  IF (IDB .GE. 9)  THEN
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 3) (B(I), I = 1, NO)
 3      FORMAT ('MPINFR O'/(6F12.0))
        NO = MIN (INT (ABS (C(1))), NDB) + 2
        WRITE (LDB, 3) (C(I), I = 1, NO)
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPINIX (M)
C
C   This initializes the double precision array U in common MPCOM5 with roots
C   of unity required by the FFT routines, which are called by MPMULX.  Before
C   calling any of the advanced MP routines (i.e. those whose names end in X),
C   this routine must be called with M set to MX, where MX is defined as the
C   integer such that 2 ^ MX = NX, and where NX is the largest precision level
C   NW that will be used in the subsequent application.  Before calling MPINIX,
C   the user must allocate at least 2^(M + 3) double precision cells in common
C   MPCOM5, which must be placed in the user's main program.  Also, at least
C   12 * NW + 6 double precision cells must be allocated in common MPCOM4.
C   Only one call to MPINIX is required, no matter how many advanced routines
C   are called.  It is not necessary for the user to call MPINIX, to allocate
C   space in MPCOM5 or to allocate more than NW + 6 cells in MPCOM4 if the
C   advanced routines are not called.
C
      DOUBLE PRECISION PI, T1, T2, U
      PARAMETER (PI = 3.141592653589793238D0)
      COMMON /MPCOM5/ U(1024)
C
C   Initialize the U array with sines and cosines in a manner that permits
C   stride one access at each FFT iteration.
C
      MM = M + 2
      N = 2 ** MM
      NU = N
      U(1) = 64 * N + MM
      KU = 2
      KN = KU + NU
      LN = 1
C
      DO 110 J = 1, MM
        T1 = PI / LN
C
CDIR$ IVDEP
        DO 100 I = 0, LN - 1
          T2 = I * T1
          U(I+KU) = COS (T2)
          U(I+KN) = SIN (T2)
 100    CONTINUE
C
        KU = KU + LN
        KN = KU + NU
        LN = 2 * LN
 110  CONTINUE
C
      RETURN
      END
C
      SUBROUTINE MPINP (IU, A, CS)
C
C   This routine reads the MP number A from logical unit IU.  CS is a scratch
C   array of type CHARACTER*1.  CS must be dimensioned at least 7.225*NW + 100.
C   The digits of A may span more than one line.  A comma at the end of the
C   last line denotes the end of the MP number.  The input lines may not
C   exceed 120 characters in length.  Embedded blanks are allowed anywhere.
C   However, if the input number contains more than 80 embedded blanks, then
C   the dimension of CS must be increased by a corresponding amount.  The
C   exponent is optional in the input number, but if present it must appear
C   first.  Two examples:
C
C   1073741824.,
C   10 ^  -4 x  3.14159 26535 89793 23846 26433 83279
C     50288 41971 69399 37510,
C
C   Max SP space for A: NW + 4 cells.  Max SP scratch space: 5 * NW + 27 cells.
C
      CHARACTER*120 LIN
      CHARACTER*1 CS(*)
      DIMENSION A(NW+2)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        A(1) = 0.
        A(2) = 0.
        RETURN
      ENDIF
      L = 0
      ND = 7.225D0 * NW + 100.D0
C
 100  READ (IU, '(A)', END = 150) LIN
C
      DO 110 I = 120, 1, -1
        IF (LIN(I:I) .NE. ' ') GOTO 120
 110  CONTINUE
C
      GOTO 100
C
 120  K = I
      IF (L .GT. ND) GOTO 140
C
      DO 130 I = 1, K
        L = L + 1
        IF (L .GT. ND) GOTO 140
        CS(L)= LIN(I:I)
 130  CONTINUE
C
 140  IF (LIN(K:K) .NE. ',') GOTO 100
      L = L - 1
C
      CALL MPINPC (CS, L, A)
      GOTO 160
C
 150  IF (KER(72) .NE. 0) THEN
        WRITE (LDB, 1)
 1      FORMAT ('*** MPINP: End-of-file encountered.')
        IER = 72
        IF (KER(IER) .EQ. 2) CALL MPABRT
      ENDIF
C
 160  RETURN
      END
C
      SUBROUTINE MPINPC (A, N, B)
C
C   Converts the CHARACTER*1 array A of length N into the MP number B.  The
C   string A must be in the format '10^s a x tb.c' where a, b and c are digit
C   strings; s and t are '-', '+' or blank; x is either 'x' or '*'.  Blanks may
C   be embedded anywhere.  The digit string a is limited to nine digits and
C   80 total characters, including blanks.  The exponent portion (i.e. the
C   portion up to and including x) and the period may optionally be omitted.
C   Debug output starts with IDB = 7.
C
C   Max SP space for B: NW + 4 cells.  Max SP scratch space: 5 * NW + 27 cells.
C
C   The following example shows how this routine may be used to input a MP
C   number:
C
C      CHARACTER*1 CX(800)
C      READ (1, '(80A1)') (CX(I), I = 1, ND)
C      CALL MPINPC (CX, ND, B)
C
      DOUBLE PRECISION BI
      CHARACTER*1 A, AI
      CHARACTER*10 DIG
      CHARACTER*80 CA
      PARAMETER (NDB = 22, DIG = '0123456789')
      DIMENSION A(N), B(NW+4), F(8)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        B(1) = 0.
        B(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 7) THEN
        NO = MIN (N, INT (7.225 * NDB) + 20)
        WRITE (LDB, 1) (A(I), I = 1, NO)
 1      FORMAT ('MPINPC I'/(78A1))
      ENDIF
C
      N5 = NW + 5
      NS = 3 * N5
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N5
      K2 = K1 + N5
      NWS = NW
      NW = NW + 1
      I1 = 1
      NN = 0
C
C   Find the carat, period, plus or minus sign, whichever comes first.
C
      DO 100 I = 1, N
        AI = A(I)
        IF (AI .EQ. '^') GOTO 110
        IF (AI .EQ. '.' .OR. AI .EQ. '+' .OR. AI .EQ. '-') GOTO 160
 100  CONTINUE
C
      GOTO 160
C
C   Make sure number preceding the carat is 10.
C
 110  I2 = I - 1
      IF (I2 .GT. 80) GOTO 210
      CA = ' '
C
      DO 120 I = 1, I2
        AI = A(I)
        IF (AI .EQ. ' ') THEN
          GOTO 120
        ELSEIF (INDEX (DIG, AI) .EQ. 0) THEN
          GOTO 210
        ENDIF
        CA(I:I) = AI
 120  CONTINUE
C
      READ (CA, '(BN,I80)') NN
      IF (NN .NE. 10) GOTO 210
      I1 = I2 + 2
C
C   Find the x or *.
C
      DO 130 I = I1, N
        AI = A(I)
        IF (AI .EQ. 'x' .OR. AI .EQ. '*') GOTO 140
 130  CONTINUE
C
      GOTO 210
C
C   Convert the exponent.
C
 140  I2 = I - 1
      L1 = I2 - I1 + 1
      IF (L1 .GT. 80) GOTO 210
      CA = ' '
      ID = 0
      IS = 1
C
      DO 150 I = 1, L1
        AI = A(I+I1-1)
        IF (AI .EQ. ' ' .OR. AI .EQ. '+') THEN
          GOTO 150
        ELSEIF (AI .EQ. '-' .AND. ID .EQ. 0) THEN
          ID = 1
          IS = -1
          CA(I:I) = ' '
        ELSE
          IF (INDEX (DIG, AI) .EQ. 0) GOTO 210
          ID = 1
          CA(I:I) = AI
        ENDIF
 150  CONTINUE
C
      READ (CA, '(BN,I80)') NN
      NN = IS * NN
      I1 = I2 + 2
C
C   Find the next nonblank character.
C
 160  DO 170 I = I1, N
        IF (A(I) .NE. ' ') GOTO 180
 170  CONTINUE
C
      GOTO 210
C
C   Check if the nonblank character is a plus or minus sign.
C
 180  I1 = I
      IF (A(I1) .EQ. '+') THEN
        I1 = I1 + 1
        IS = 1
      ELSEIF (A(I1) .EQ. '-') THEN
        I1 = I1 + 1
        IS = -1
      ELSE
        IS = 1
      ENDIF
      NB = 0
      IB = 0
      ID = 0
      IP = 0
      S(K2) = 0.
      S(K2+1) = 0.
      F(1) = 1.
      F(2) = 0.
      IT = 0
C
 190  IP = 0
      CA(1:6) = '000000'
C
C   Scan for digits, looking for the period also.  On the first pass we just
C   count, so that on the second pass it will come out right.
C
      DO 200 I = I1, N
        AI = A(I)
        IF (AI .EQ. ' ') THEN
        ELSEIF (AI .EQ. '.') THEN
          IF (IP .NE. 0) GOTO 210
          IP = ID
        ELSEIF (INDEX (DIG, AI) .EQ. 0) THEN
          GOTO 210
        ELSE
          IB = IB + 1
          ID = ID + 1
          CA(IB:IB) = AI
        ENDIF
        IF (IB .EQ. 6 .OR. I .EQ. N .AND. IB .NE. 0) THEN
          IF (IT .NE. 0) THEN
            NB = NB + 1
            READ (CA(1:6), '(F6.0)') BI
            CALL MPMULD (S(K2), 1.D6, 0, S(K0))
            IF (BI .NE. 0) THEN
              F(1) = 1.
              F(3) = BI
            ELSE
              F(1) = 0.
            ENDIF
            CALL MPADD (S(K0), F, S(K2))
            CA(1:6) = '000000'
          ENDIF
          IF (I .NE. N) IB = 0
        ENDIF
 200  CONTINUE
C
      IF (IT .EQ. 0) THEN
        IB = 6 - IB
        IF (IB .EQ. 6) IB = 0
        IT = 1
        GOTO 190
      ENDIF
      IF (IS .EQ. -1) S(K2) = - S(K2)
      IF (IP .EQ. 0) IP = ID
      NN = NN + IP - ID
      F(1) = 1.
      F(3) = 10.
      CALL MPNPWR (F, NN, S(K0))
      CALL MPMUL (S(K2), S(K0), S(K1))
      CALL MPEQ (S(K1), B)
      NW = NWS
      CALL MPROUN (B)
      ICS = ISS
C
      IF (IDB .GE. 7) THEN
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 2) (B(I), I = 1, NO)
 2      FORMAT ('MPINPC O'/(6F12.0))
      ENDIF
      GOTO 220
C
 210  IF (KER(41) .NE. 0) THEN
        WRITE (LDB, 3)
 3      FORMAT ('*** MPINPC: Syntax error in literal string.')
        IER = 41
        IF (KER(IER) .EQ. 2) CALL MPABRT
        NW = NWS
        ICS = ISS
      ENDIF
C
 220  RETURN
      END
C
      SUBROUTINE MPINQP (IA, IB)
C
C   This routine returns the value of the parameter whose name is IA in common
C   MPCOM1.  By using this routine instead of merely including the MPCOM1
C   block in the code, a user may eliminate the possibility of confusion with
C   a variable name in his or her program.  IA is of type CHARACTER and IB
C   is the value.
C
      CHARACTER*(*) IA
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
C
      IF (IA .EQ. 'NW' .OR. IA .EQ. 'nw') THEN
        IB = NW
      ELSEIF (IA .EQ. 'IDB' .OR. IA .EQ. 'idb') THEN
        IB = IDB
      ELSEIF (IA .EQ. 'LDB' .OR. IA .EQ. 'ldb') THEN
        IB = LDB
      ELSEIF (IA .EQ. 'IER' .OR. IA .EQ. 'ier') THEN
        IB = IER
      ELSEIF (IA .EQ. 'MCR' .OR. IA .EQ. 'mcr') THEN
        IB = MCR
      ELSEIF (IA .EQ. 'IRD' .OR. IA .EQ. 'ird') THEN
        IB = IRD
      ELSEIF (IA .EQ. 'ICS' .OR. IA .EQ. 'ics') THEN
        IB = ICS
      ELSEIF (IA .EQ. 'IHS' .OR. IA .EQ. 'ihs') THEN
        IB = IHS
      ELSEIF (IA .EQ. 'IMS' .OR. IA .EQ. 'ims') THEN
        IB = IMS
      ELSE
        IB = 0
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE MPINRL (N, LX, X, MN, MT, LR, R, IQ)
C
C   This routine searches for integer relations among the entries of the
C   N-long MP vector X.  An integer relation is an n-long vector r such that
C   r_1 x_1 + r_2 x_2 + ... + r_n x_n = 0.  The entries of x are assumed to
C   start at X(1), X(LX+1), X(2*LX+1), etc.  MN is the Log_10 of the maximum
C   Euclidean norm of an acceptable relation.  IQ is set to 1 if the routine
C   succeeds in recovering a relation that (1) produces zero to within the
C   relative tolerance 10^MT and (2) has Euclidean norm less than 10^MN.  If
C   no relation is found that meets these standards, IQ is set to 0.  When a
C   valid relation vector is recovered, it is placed in R, beginning at R(1),
C   R(LR+1), R(2*LR+1), etc., where LR, like LX, is an input parameter.  LR
C   should be at least MN/6 + 3.  For extra-high levels of precision, call
C   MPINRX.  Debug output starts with IDB = 4.  When IDB = 5, norm bounds are
C   output within which no relation can exist.
C
C   Max SP space for R: LR * N cells.  Max SP scratch space:
C   (4 * N^2 + 5 * N + 12) * (NW + 4) cells.  Max DP scratch space: NW + 4
C   cells.
C
C   A typical application of this routine is to determine if a given computed
C   real number r is the root of any algebraic equation of degree n - 1 with
C   integer coefficients.  One merely sets x_k = r^(k-1) for k = 1 to n and
C   calls MPINRL.  If an integer relation is found, this relation is the vector
C   of coefficients of a polynomial satisfied by r.  If MPINRL outputs a norm
C   bound of B, then r is not the root of any polynomial of degree n or less
C   with integer coefficients, where the Euclidean norm of the vector of
C   coefficients is less than B.
C
C   It sometimes happens that the "precision exhausted" message is output
C   before finding a relation that is known to exist.  If this happens,
C   increase NW, the working precision level, as well as scratch space
C   allocations if necessary, and try again.  Typically MT is set to roughly
C   10 - 6 * NX, where NX is the precision level used to compute X.  Repeating
C   a run with somewhat higher precision is highly recommended to certify that
C   bounds results are valid.
C
C   This routine allocates the scratch space array S for arrays.  Otherwise the
C   indexing in MPINRQ is too complicated.
C
      CHARACTER*8 CX
      PARAMETER (IB = 6)
      DIMENSION R(LR,N), X(LX,N)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        IQ = 0
        RETURN
      ENDIF
      IF (IDB .GE. 5) THEN
        WRITE (LDB, 1) N, LX, MN, LR
 1      FORMAT ('MPINRL I',4I6)
C
        DO 100 K = 1, N
          WRITE (CX, '(I4)') K
          CALL MPDEB (CX, X(1,K))
 100    CONTINUE
C
      ENDIF
C
C   Check if enough space is allowed for R.
C
      IF (LR .LT. MN / IB + 3) THEN
        IF (KER(42) .NE. 0) THEN
          WRITE (LDB, 2)
 2        FORMAT ('*** MPINRL: Argument LR must be larger to match MN.')
          IER = 42
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
      N4 = NW + 4
      NS = (4 * N ** 2 + 5 * N + 7) * N4
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      KBN = N * (N + 1)
      KBS = N + 1 + KBN
      KC = N * (N + 1) + KBS
      KU = N * N + KC
      CALL MPINRQ (N, LX, X, MN, MT, LR, R, IQ, S(K0), S(KBN*N4+K0),    
     $  S(KBS*N4+K0), S(KC*N4+K0), S(KU*N4+K0))
      ICS = ISS
C
      IF (IDB .GE. 5) THEN
        WRITE (LDB, 3) IQ
 3      FORMAT ('MPINRL O',I2)
        IF (IQ .EQ. 1) THEN
C
          DO 110 K = 1, N
            WRITE (CX, '(I4)') K
            CALL MPDEB (CX, R(1,K))
 110      CONTINUE
C
        ENDIF
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPINRQ (N, LX, X, MN, MT, LR, R, IQ, B, BN, BS, C, U)
C
C   This routine implements the "Small Integer Relation Algorithm" described
C   in Hastad, Just, Lagarias, and Schnorr, "Polynomial Time Algorithms for
C   Finding Integer Relations Among Real Numbers", to appear in SIAM J. on
C   Computing.  This routine is called by MPINRL.  It is not intended to be
C   called directly by the user.
C
C   IMX = Number of iterations after which run is terminated.
C   ITP = Print interval.  Also the interval at which norm bounds are computed.
C   LB  = Reduction in log_10 (BN(N)) from previous iteration.  Used to detect
C         that a tentative relation has been found.
C
      DOUBLE PRECISION AB, BNN, BNS, BNZ, BX, BY, T1, T2, T3, T4, TB
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (IMX = 10000, ITP = 10, ITZ = 100, LB = 20)
      DIMENSION B(NW+4,N,0:N), BN(NW+4,0:N), BS(NW+4,N,0:N),            
     $  C(NW+4,N,N), R(LR,N), U(NW+4,0:N,0:N), X(LX,N)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
C   Step 1: Initialization.
C
      N4 = NW + 4
      NS = 5 * N4
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N4
      K2 = K1 + N4
      K3 = K2 + N4
      K4 = K3 + N4
      NWS = NW
      TL = 2 - NW
      BNS = 0.D0
      BNZ = 0.D0
      IBS = 0
      IBZ = 0
      II = 0
      IQ = 0
C
      DO 100 I = 1, N
        CALL MPEQ (X(1,I), B(1,I,0))
 100  CONTINUE
C
      DO 120 J = 1, N
C
        DO 110 I = 1, N
          B(1,I,J) = 0.
          B(2,I,J) = 0.
          C(1,I,J) = 0.
          C(2,I,J) = 0.
 110    CONTINUE
C
        B(1,J,J) = 1.
        B(2,J,J) = 0.
        B(3,J,J) = 1.
        C(1,J,J) = 1.
        C(2,J,J) = 0.
        C(3,J,J) = 1.
 120  CONTINUE
C
      DO 180 I = 0, N
C
        DO 130 K = 1, N
          CALL MPEQ (B(1,K,I), BS(1,K,I))
 130    CONTINUE
C
        DO 160 J = 0, I - 1
          S(K0) = 0.
          S(K0+1) = 0.
C
          DO 140 K = 1, N
            CALL MPMUL (B(1,K,I), BS(1,K,J), S(K1))
            CALL MPADD (S(K0), S(K1), S(K2))
            CALL MPEQ (S(K2), S(K0))
 140      CONTINUE
C
          IF (BN(1,J) .EQ. 0. .OR. BN(2,J) .LT. TL) THEN
            U(1,I,J) = 0.
            U(2,I,J) = 0.
          ELSE
            CALL MPDIV (S(K0), BN(1,J), U(1,I,J))
          ENDIF
          U(1,J,I) = 0.
          U(2,J,I) = 0.
C
          DO 150 K = 1, N
            CALL MPMUL (U(1,I,J), BS(1,K,J), S(K0))
            CALL MPSUB (BS(1,K,I), S(K0), S(K1))
            CALL MPEQ (S(K1), BS(1,K,I))
 150      CONTINUE
 160    CONTINUE
C
        S(K0) = 0.
        S(K0+1) = 0.
C
        DO 170 K = 1, N
          CALL MPMUL (BS(1,K,I), BS(1,K,I), S(K1))
          CALL MPADD (S(K0), S(K1), S(K2))
          CALL MPEQ (S(K2), S(K0))
 170    CONTINUE
C
        CALL MPEQ (S(K0), BN(1,I))
        U(1,I,I) = 1.
        U(2,I,I) = 0.
        U(3,I,I) = 1.
 180  CONTINUE
C
C   Step 2: Termination test.
C
 190  II = II + 1
      IF (IER .NE. 0) RETURN
      IF (II .GT. IMX) THEN
        IF (KER(43) .NE. 0) THEN
          WRITE (LDB, 1) II
 1        FORMAT ('*** MPINRQ: Iteration limit exceeded',I6)
          IER = 43
          IF (KER(IER) .EQ. 2) CALL MPABRT
          ICS = ISS
          RETURN
        ENDIF
      ENDIF
      BX = 0.D0
      BY = 0.D0
      IX = -10000
      IY = -10000
C
      DO 200 I = 1, N - 1
        CALL MPMDC (BN(1,I), AB, IB)
        IF ((AB .GT. BX .AND. IB .EQ. IX) .OR. (AB .NE. 0.D0            
     $    .AND. IB .GT. IX)) THEN
          BX = AB
          IX = IB
        ENDIF
        CALL DPMUL (2.D0 ** I, 0, AB, IB, T1, N1)
        IF ((T1 .GT. BY .AND. N1 .EQ. IY) .OR. (T1 .NE. 0.D0            
     $    .AND. N1 .GT. IY)) THEN
          BY = T1
          IY = N1
          I1 = I
        ENDIF
 200  CONTINUE
C
      CALL DPSQRT (BX, IX, T1, N1)
      CALL DPDIV (1.D0, 0, T1, N1, T2, N2)
      CALL DPDEC (T2, N2, TB, NB)
      CALL MPMDC (BN(1,N), T2, N2)
      CALL DPDEC (T2, N2, BNN, IBN)
      IF ((IDB .GE. 5 .AND. MOD (II, ITP) .EQ. 0) .OR. IDB .GE. 6) THEN
        WRITE (LDB, 2) II, TB, NB, BNN, IBN
 2      FORMAT ('Iteration', I6/ 'Norm bound =', F10.6, ' x 10^', I6,   
     $    4X, 'BN(N) =', F10.6, ' x 10^', I6)
        IF (IDB .GE. 6) THEN
          WRITE (LDB, 3)
 3        FORMAT ('BSTAR square norms:')
          CALL MPMOUT (1, N, BN(1,1))
        ENDIF
        IF (IDB .GE. 7) THEN
          WRITE (LDB, 4)
 4        FORMAT ('B Matrix')
          CALL MPMOUT (N, N + 1, B)
          WRITE (LDB, 5)
 5        FORMAT ('U Matrix')
          CALL MPMOUT (N + 1, N + 1, U)
        ENDIF
      ENDIF
      IF (NB .GT. MN) GOTO 280
C
C   Test if current BN(N) is 10^LB times the previous BN(N).
C
      IF (BNN .NE. 0.D0 .AND. IBN .GT. IBS + LB) THEN
        IF (IDB .GE. 5) WRITE (LDB, 6) II, BNN, IBN
 6      FORMAT (/'Tentative relation, iteration', I6, 4X, 'BN(N) =',    
     $    F10.6, ' x 10^', I6)
C
C   Compute residual and norm of tentative relation.
C
        DO 220 K = N, 1, -1
          T2 = 0.D0
          N2 = 0
          S(K0) = 0.
          S(K0+1) = 0.
C
          DO 210 J = 1, N
            NW = LR - 2
            CALL MPEQ (C(1,J,K), R(1,J))
            NW = NWS
            CALL MPMDC (R(1,J), T1, N1)
            CALL DPMUL (T1, N1, T1, N1, T3, N3)
            CALL DPADD (T2, N2, T3, N3, T4, N4)
            T2 = T4
            N2 = N4
            CALL MPMUL (R(1,J), X(1,J), S(K1))
            CALL MPADD (S(K0), S(K1), S(K2))
            CALL MPEQ (S(K2), S(K0))
 210      CONTINUE
C
C   If the residual is zero or within tolerance 10^MT of zero, it is a real
C   relation.  Otherwise it was a false alarm.
C
          CALL MPMDC (S(K0), T3, N3)
          CALL DPDEC (T3, N3, T1, N1)
          IF (T1 .EQ. 0.D0 .OR. N1 .LE. MT) THEN
            IF (IDB .GE. 4) THEN
              CALL DPSQRT (T2, N2, T3, N3)
              CALL DPDEC (T3, N3, T1, N1)
              CALL MPMDC (S(K0), T4, N4)
              CALL DPDEC (T4, N4, T2, N2)
              WRITE (LDB, 7) K, T1, N1, T2, N2
 7            FORMAT ('Relation in column',I4,3X,'Norm =',F10.6,        
     $          ' x 10^',I6/'Residual =',F10.6,' x 10^',I6)
            ENDIF
            IQ = 1
            GOTO 280
          ENDIF
 220    CONTINUE
C
      ENDIF
C
C   Test if BN(N) is the same as ITZ iterations ago.
C
      IF (MOD (II, ITZ) .EQ. 0) THEN
        IF (BNN .EQ. BNZ .AND. IBN .EQ. IBZ) THEN
          IF (KER(44) .NE. 0) THEN
            WRITE (LDB, 8) INT (LOG10 (BDX) * (NW + 3))
 8          FORMAT ('*** MPINRQ: Numeric overflow has occurred.  Call ',
     $        'MPINRL with at least',I8/'digits precision.')
            IER = 44
            IF (KER(IER) .EQ. 2) CALL MPABRT
          ENDIF
          ICS = ISS
          RETURN
        ENDIF
        BNZ = BNN
        IBZ = IBN
      ENDIF
      BNS = BNN
      IBS = IBN
C
C   Step 3: Update B and C for transformation and then exchange B and C.
C
      I2 = I1 + 1
C
C   Check if U(i2,i1) can be converted exactly to an integer.  The error
C   number and message are the same as the previous one.
C
      IF (ABS (U(2,I2,I1)) .GE. NW - 1) THEN
        IF (KER(45) .NE. 0) THEN
          WRITE (LDB, 8) INT (LOG10 (BDX) * ABS (U(2,I2,I1)))
          IER = 45
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        ICS = ISS
        RETURN
      ENDIF
      CALL MPNINT (U(1,I2,I1), S(K0))
C
      DO 230 K = 1, N
        CALL MPMUL (S(K0), B(1,K,I1), S(K1))
        CALL MPSUB (B(1,K,I2), S(K1), S(K2))
        CALL MPEQ (S(K2), B(1,K,I2))
        CALL MPMUL (S(K0), C(1,K,I2), S(K1))
        CALL MPADD (C(1,K,I1), S(K1), S(K2))
        CALL MPEQ (S(K2), C(1,K,I1))
 230  CONTINUE
C
      DO 240 K = 1, N
        CALL MPEQ (B(1,K,I1), S(K1))
        CALL MPEQ (B(1,K,I2), B(1,K,I1))
        CALL MPEQ (S(K1), B(1,K,I2))
        CALL MPEQ (C(1,K,I1), S(K1))
        CALL MPEQ (C(1,K,I2), C(1,K,I1))
        CALL MPEQ (S(K1), C(1,K,I2))
 240  CONTINUE
C
C   Update U for transformation.
C
      DO 250 J = 0, I1
        CALL MPMUL (S(K0), U(1,I1,J), S(K1))
        CALL MPSUB (U(1,I2,J), S(K1), S(K2))
        CALL MPEQ (S(K2), U(1,I2,J))
 250  CONTINUE
C
C   Update BN and U for exchange.
C
      CALL MPEQ (U(1,I2,I1), S(K0))
      CALL MPMUL (S(K0), S(K0), S(K1))
      CALL MPMUL (S(K1), BN(1,I1), S(K2))
      CALL MPADD (BN(1,I2), S(K2), S(K1))
      IF (S(K1) .NE. 0. .AND. S(K1+1) .GT. TL) THEN
        CALL MPDIV (BN(1,I1), S(K1), S(K2))
        CALL MPMUL (BN(1,I2), S(K2), S(K3))
        CALL MPEQ (S(K3), BN(1,I2))
        CALL MPMUL (S(K0), S(K2), S(K3))
        CALL MPEQ (S(K3), U(1,I2,I1))
      ELSE
        CALL MPEQ (BN(1,I1), BN(1,I2))
        U(1,I2,I1) = 0.
        U(2,I2,I1) = 0.
      ENDIF
      CALL MPEQ (S(K1), BN(1,I1))
C
      DO 260 J = 1, I1 - 1
        CALL MPEQ (U(1,I1,J), S(K1))
        CALL MPEQ (U(1,I2,J), U(1,I1,J))
        CALL MPEQ (S(K1), U(1,I2,J))
 260  CONTINUE
C
      S(K1) = 1.
      S(K1+1) = 0.
      S(K1+2) = 1.
C
      DO 270 J = I1 + 2, N
        CALL MPMUL (U(1,J,I1), U(1,I2,I1), S(K2))
        CALL MPMUL (S(K0), U(1,I2,I1), S(K3))
        CALL MPSUB (S(K1), S(K3), S(K4))
        CALL MPMUL (U(1,J,I2), S(K4), S(K3))
        CALL MPADD (S(K2), S(K3), S(K4))
        CALL MPMUL (S(K0), U(1,J,I2), S(K2))
        CALL MPSUB (U(1,J,I1), S(K2), U(1,J,I2))
        CALL MPEQ (S(K4), U(1,J,I1))
 270  CONTINUE
C
      GOTO 190
C
 280  IF (IDB .GE. 4) WRITE (6, 9) II, TB, NB
 9    FORMAT ('No. iterations =',I6/'Max. bound =',1PD15.6,             
     $  ' x 10^',I5)
      ICS = ISS
      RETURN
      END
C
      SUBROUTINE MPINRX (N, LX, X, MN, MT, LR, R, IQ)
C
C   This routine searches for integer relations among the entries of the
C   N-long MP vector X.  An integer relation is an n-long vector r such that
C   r_1 x_1 + r_2 x_2 + ... + r_n x_n = 0.  The entries of x are assumed to
C   start at X(1), X(LX+1), X(2*LX+1), etc.  MN is the Log_10 of the maximum
C   Euclidean norm of an acceptable relation.  IQ is set to 1 if the routine
C   succeeds in recovering a relation that (1) produces zero to within the
C   relative tolerance 10^MT and (2) has Euclidean norm less than 10^MN.  If
C   no relation is found that meets these standards, IQ is set to 0.  When a
C   valid relation vector is recovered, it is placed in R, beginning at R(1),
C   R(LR+1), R(2*LR+1), etc., where LR, like LX, is an input parameter.  LR
C   should be at least MN/6 + 3.  Before calling MPINRX, the array in MPCOM5
C   must be initialized by calling MPINIX.  For modest levels of precision,
C   call MPINRL.  Debug output starts with IDB = 4.  When IDB = 5, norm bounds
C   are output within which no relation can exist.
C
C   Max SP space for R: LR * N cells.  Max SP scratch space:
C   (4 * N^2 + 5 * N + 15) * (NW + 4) cells.  Max DP scratch space: 12 * NW + 6
C   cells.
C
C   See the comments in MPINRL about applying this routine.
C
C   This allocates the scratch space array S for arrays.  Otherwise the
C   indexing in MPINRZ is too complicated.
C
      CHARACTER*8 CX
      PARAMETER (IB = 6)
      DIMENSION R(LR,N), X(LX,N)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        IQ = 0
        RETURN
      ENDIF
      IF (IDB .GE. 5) THEN
        WRITE (LDB, 1) N, LX, MN, LR
 1      FORMAT ('MPINRX I',4I6)
C
        DO 100 K = 1, N
          WRITE (CX, '(I4)') K
          CALL MPDEB (CX, X(1,K))
 100    CONTINUE
C
      ENDIF
C
C   Check if enough space is allowed for R.
C
      IF (LR .LE. MN / IB) THEN
        IF (KER(46) .NE. 0) THEN
          WRITE (LDB, 2)
 2        FORMAT ('*** MPINRX: Argument LR must be larger to match MN.')
          IER = 46
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
C   Check if the precision level is too low to justify the advanced routine.
C
      NCR = 2 ** MCR
      IF (NW .LE. NCR) THEN
        CALL MPINRL (N, LX, X, MN, MT, LR, R, IQ)
        GOTO 110
      ENDIF
C
C   Compute pointers for arrays.
C
      N4 = NW + 4
      NS = (4 * N ** 2 + 5 * N + 7) * N4
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      KBN = N * (N + 1)
      KBS = N + 1 + KBN
      KC = N * (N + 1) + KBS
      KU = N * N + KC
      CALL MPINRZ (N, LX, X, MN, MT, LR, R, IQ, S(K0), S(KBN*N4+K0),    
     $  S(KBS*N4+K0), S(KC*N4+K0), S(KU*N4+K0))
      ICS = ISS
C
 110  IF (IDB .GE. 5) THEN
        WRITE (LDB, 3) IQ
 3      FORMAT ('MPINRX O',I2)
        IF (IQ .EQ. 1) THEN
C
          DO 120 K = 1, N
            WRITE (CX, '(I4)') K
            CALL MPDEB (CX, R(1,K))
 120      CONTINUE
C
        ENDIF
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPINRZ (N, LX, X, MN, MT, LR, R, IQ, B, BN, BS, C, U)
C
C   This is the extra-high precision version of MPINRQ.  See the comments
C   there for details.
C
      DOUBLE PRECISION AB, BNN, BNS, BNZ, BX, BY, T1, T2, T3, T4, TB
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (IMX = 10000, ITP = 25, ITZ = 250, LB = 20)
      DIMENSION B(NW+4,N,0:N), BN(NW+4,0:N), BS(NW+4,N,0:N),            
     $  C(NW+4,N,N), R(LR,N), U(NW+4,0:N,0:N), X(LX,N)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
C   Step 1: Initialization.
C
      N4 = NW + 4
      NS = 5 * N4
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N4
      K2 = K1 + N4
      K3 = K2 + N4
      K4 = K3 + N4
      NWS = NW
      IF (NW .LE. 32) THEN
        TL = 2 - NW
      ELSEIF (NW .LE. 256) THEN
        TL = 3 - NW
      ELSE
        TL = 4 - NW
      ENDIF
      BNS = 0.D0
      BNZ = 0.D0
      IBS = 0
      IBZ = 0
      II = 0
      IQ = 0
C
      DO 100 I = 1, N
        CALL MPEQ (X(1,I), B(1,I,0))
 100  CONTINUE
C
      DO 120 J = 1, N
C
        DO 110 I = 1, N
          B(1,I,J) = 0.
          B(2,I,J) = 0.
          C(1,I,J) = 0.
          C(2,I,J) = 0.
 110    CONTINUE
C
        B(1,J,J) = 1.
        B(2,J,J) = 0.
        B(3,J,J) = 1.
        C(1,J,J) = 1.
        C(2,J,J) = 0.
        C(3,J,J) = 1.
 120  CONTINUE
C
      DO 180 I = 0, N
C
        DO 130 K = 1, N
          CALL MPEQ (B(1,K,I), BS(1,K,I))
 130    CONTINUE
C
        DO 160 J = 0, I - 1
          S(K0) = 0.
          S(K0+1) = 0.
C
          DO 140 K = 1, N
            CALL MPMULX (B(1,K,I), BS(1,K,J), S(K1))
            CALL MPADD (S(K0), S(K1), S(K2))
            CALL MPEQ (S(K2), S(K0))
 140      CONTINUE
C
          IF (BN(1,J) .EQ. 0. .OR. BN(2,J) .LT. TL) THEN
            U(1,I,J) = 0.
            U(2,I,J) = 0.
          ELSE
            CALL MPDIVX (S(K0), BN(1,J), U(1,I,J))
          ENDIF
          U(1,J,I) = 0.
          U(2,J,I) = 0.
C
          DO 150 K = 1, N
            CALL MPMULX (U(1,I,J), BS(1,K,J), S(K0))
            CALL MPSUB (BS(1,K,I), S(K0), S(K1))
            CALL MPEQ (S(K1), BS(1,K,I))
 150      CONTINUE
 160    CONTINUE
C
        S(K0) = 0.
        S(K0+1) = 0.
C
        DO 170 K = 1, N
          CALL MPSQX (BS(1,K,I), S(K1))
          CALL MPADD (S(K0), S(K1), S(K2))
          CALL MPEQ (S(K2), S(K0))
 170    CONTINUE
C
        CALL MPEQ (S(K0), BN(1,I))
        U(1,I,I) = 1.
        U(2,I,I) = 0.
        U(3,I,I) = 1.
 180  CONTINUE
C
C   Step 2: Termination test.
C
 190  II = II + 1
      IF (IER .NE. 0) RETURN
      IF (II .GT. IMX) THEN
        IF (KER(47) .NE. 0) THEN
          WRITE (LDB, 1) II
 1        FORMAT ('*** MPINRZ: Iteration limit exceeded',I6)
          IER = 47
          IF (KER(IER) .EQ. 2) CALL MPABRT
          ICS = ISS
          RETURN
        ENDIF
      ENDIF
      BX = 0.D0
      BY = 0.D0
      IX = -10000
      IY = -10000
C
      DO 200 I = 1, N - 1
        CALL MPMDC (BN(1,I), AB, IB)
        IF ((AB .GT. BX .AND. IB .EQ. IX) .OR. (AB .NE. 0.D0            
     $    .AND. IB .GT. IX)) THEN
          BX = AB
          IX = IB
        ENDIF
        CALL DPMUL (2.D0 ** I, 0, AB, IB, T1, N1)
        IF ((T1 .GT. BY .AND. N1 .EQ. IY) .OR. (T1 .NE. 0.D0            
     $    .AND. N1 .GT. IY)) THEN
          BY = T1
          IY = N1
          I1 = I
        ENDIF
 200  CONTINUE
C
      CALL DPSQRT (BX, IX, T1, N1)
      CALL DPDIV (1.D0, 0, T1, N1, T2, N2)
      CALL DPDEC (T2, N2, TB, NB)
      CALL MPMDC (BN(1,N), T2, N2)
      CALL DPDEC (T2, N2, BNN, IBN)
      IF ((IDB .GE. 5 .AND. MOD (II, ITP) .EQ. 0) .OR. IDB .GE. 6) THEN
        WRITE (LDB, 2) II, TB, NB, BNN, IBN
 2      FORMAT ('Iteration', I6/ 'Norm bound =', F10.6, ' x 10^', I6,   
     $    4X, 'BN(N) =', F10.6, ' x 10^', I6)
        IF (IDB .GE. 6) THEN
          WRITE (LDB, 3)
 3        FORMAT ('BSTAR square norms:')
          CALL MPMOUT (1, N, BN(1,1))
        ENDIF
        IF (IDB .GE. 7) THEN
          WRITE (LDB, 4)
 4        FORMAT ('B Matrix')
          CALL MPMOUT (N, N + 1, B)
          WRITE (LDB, 5)
 5        FORMAT ('U Matrix')
          CALL MPMOUT (N + 1, N + 1, U)
        ENDIF
      ENDIF
      IF (NB .GT. MN) GOTO 280
C
C   Test if current BN(N) is 10^LB times the previous BN(N).
C
      IF (BNN .NE. 0.D0 .AND. IBN .GT. IBS + LB) THEN
        IF (IDB .GE. 5) WRITE (LDB, 6) II, BNN, IBN
 6      FORMAT (/'Tentative relation, iteration', I6, 4X, 'BN(N) =',    
     $    F10.6, ' x 10^', I6)
C
C   Compute residual and norm of tentative relation.
C
        DO 220 K = N, 1, -1
          T2 = 0.D0
          N2 = 0
          S(K0) = 0.
          S(K0+1) = 0.
C
          DO 210 J = 1, N
            NW = LR - 2
            CALL MPEQ (C(1,J,K), R(1,J))
            NW = NWS
            CALL MPMDC (R(1,J), T1, N1)
            CALL DPMUL (T1, N1, T1, N1, T3, N3)
            CALL DPADD (T2, N2, T3, N3, T4, N4)
            T2 = T4
            N2 = N4
            CALL MPMULX (R(1,J), X(1,J), S(K1))
            CALL MPADD (S(K0), S(K1), S(K2))
            CALL MPEQ (S(K2), S(K0))
 210      CONTINUE
C
C   If the residual is zero or within tolerance 10^MT of zero, it is a real
C   relation.  Otherwise it was a false alarm.
C
          CALL MPMDC (S(K0), T3, N3)
          CALL DPDEC (T3, N3, T1, N1)
          IF (T1 .EQ. 0.D0 .OR. N1 .LE. MT) THEN
            IF (IDB .GE. 4) THEN
              CALL DPSQRT (T2, N2, T3, N3)
              CALL DPDEC (T3, N3, T1, N1)
              CALL MPMDC (S(K0), T4, N4)
              CALL DPDEC (T4, N4, T2, N2)
              WRITE (LDB, 7) K, T1, N1, T2, N2
 7            FORMAT ('Relation in column',I4,3X,'Norm =',F10.6,        
     $          ' x 10^',I6/'Residual =',F10.6,' x 10^',I6)
            ENDIF
            IQ = 1
            GOTO 280
          ENDIF
 220    CONTINUE
C
      ENDIF
C
C   Test if BN(N) is the same as ITZ iterations ago.
C
      IF (MOD (II, ITZ) .EQ. 0) THEN
        IF (BNN .EQ. BNZ .AND. IBN .EQ. IBZ) THEN
          IF (KER(48) .NE. 0) THEN
            WRITE (LDB, 8) INT (LOG10 (BDX) * (NW + 3))
 8          FORMAT ('*** MPINRZ: Numeric overflow has occurred.  Call ',
     $        'MPINRX with at least',I8/'digits precision.')
            IER = 48
            IF (KER(IER) .EQ. 2) CALL MPABRT
          ENDIF
          ICS = ISS
          RETURN
        ENDIF
        BNZ = BNN
        IBZ = IBN
      ENDIF
      BNS = BNN
      IBS = IBN
C
C   Step 3: Update B and C for transformation and then exchange B and C.
C
      I2 = I1 + 1
C
C   Check if U(i2,i1) can be converted exactly to an integer.  The error
C   number and message are the same as the previous one.
C
      IF (ABS (U(2,I2,I1)) .GE. NW - 1) THEN
        IF (KER(49) .NE. 0) THEN
          WRITE (LDB, 8) INT (LOG10 (BDX) * ABS (U(2,I2,I1)))
          IER = 49
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        ICS = ISS
        RETURN
      ENDIF
      CALL MPNINT (U(1,I2,I1), S(K0))
C
      DO 230 K = 1, N
        CALL MPMULX (S(K0), B(1,K,I1), S(K1))
        CALL MPSUB (B(1,K,I2), S(K1), S(K2))
        CALL MPEQ (S(K2), B(1,K,I2))
        CALL MPMULX (S(K0), C(1,K,I2), S(K1))
        CALL MPADD (C(1,K,I1), S(K1), S(K2))
        CALL MPEQ (S(K2), C(1,K,I1))
 230  CONTINUE
C
      DO 240 K = 1, N
        CALL MPEQ (B(1,K,I1), S(K1))
        CALL MPEQ (B(1,K,I2), B(1,K,I1))
        CALL MPEQ (S(K1), B(1,K,I2))
        CALL MPEQ (C(1,K,I1), S(K1))
        CALL MPEQ (C(1,K,I2), C(1,K,I1))
        CALL MPEQ (S(K1), C(1,K,I2))
 240  CONTINUE
C
C   Update U for transformation.
C
      DO 250 J = 0, I1
        CALL MPMULX (S(K0), U(1,I1,J), S(K1))
        CALL MPSUB (U(1,I2,J), S(K1), S(K2))
        CALL MPEQ (S(K2), U(1,I2,J))
 250  CONTINUE
C
C   Update BN and U for exchange.
C
      CALL MPEQ (U(1,I2,I1), S(K0))
      CALL MPSQX (S(K0), S(K1))
      CALL MPMULX (S(K1), BN(1,I1), S(K2))
      CALL MPADD (BN(1,I2), S(K2), S(K1))
      IF (S(K1) .NE. 0. .AND. S(K1+1) .GT. TL) THEN
        CALL MPDIVX (BN(1,I1), S(K1), S(K2))
        CALL MPMULX (BN(1,I2), S(K2), S(K3))
        CALL MPEQ (S(K3), BN(1,I2))
        CALL MPMULX (S(K0), S(K2), S(K3))
        CALL MPEQ (S(K3), U(1,I2,I1))
      ELSE
        CALL MPEQ (BN(1,I1), BN(1,I2))
        U(1,I2,I1) = 0.
        U(2,I2,I1) = 0.
      ENDIF
      CALL MPEQ (S(K1), BN(1,I1))
C
      DO 260 J = 1, I1 - 1
        CALL MPEQ (U(1,I1,J), S(K1))
        CALL MPEQ (U(1,I2,J), U(1,I1,J))
        CALL MPEQ (S(K1), U(1,I2,J))
 260  CONTINUE
C
      S(K1) = 1.
      S(K1+1) = 0.
      S(K1+2) = 1.
C
      DO 270 J = I1 + 2, N
        CALL MPMULX (U(1,J,I1), U(1,I2,I1), S(K2))
        CALL MPMULX (S(K0), U(1,I2,I1), S(K3))
        CALL MPSUB (S(K1), S(K3), S(K4))
        CALL MPMULX (U(1,J,I2), S(K4), S(K3))
        CALL MPADD (S(K2), S(K3), S(K4))
        CALL MPMULX (S(K0), U(1,J,I2), S(K2))
        CALL MPSUB (U(1,J,I1), S(K2), U(1,J,I2))
        CALL MPEQ (S(K4), U(1,J,I1))
 270  CONTINUE
C
      GOTO 190
C
 280  IF (IDB .GE. 4) WRITE (6, 9) II, TB, NB
 9    FORMAT ('No. iterations =',I6/'Max. bound =',1PD15.6,             
     $  ' x 10^',I5)
      ICS = ISS
      RETURN
      END
C
      SUBROUTINE MPLOG (A, AL2, B)
C
C   This computes the natural logarithm of the MP number A and returns the MP
C   result in B.  AL2 is the MP value of Log(2) produced by a prior call to
C   MPLOG.  For extra high levels of precision, use MPLOGX.  The last word of
C   the result is not reliable.  Debug output starts with IDB = 6.
C
C   Max SP space for B: NW + 4 cells.  Max SP scratch space: 8 * NW + 45
C   cells.  Max DP scratch space: NW + 6 cells.
C
C   The Taylor series for Log converges much more slowly than that of Exp.
C   Thus this routine does not employ Taylor series, but instead computes
C   logarithms by solving Exp (b) = a using the following Newton iteration,
C   which converges to b:
C
C           x_{k+1} = x_k + [a - Exp (x_k)] / Exp (x_k)
C
C   These iterations are performed with a maximum precision level NW that
C   is dynamically changed, approximately doubling with each iteration.
C   See the comment about the parameter NIT in MPDIVX.
C
      DOUBLE PRECISION ALT, CL2, T1, T2
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (ALT = 0.693147180559945309D0,                          
     $  CL2 = 1.4426950408889633D0, NIT = 3)
      DIMENSION A(NW+2), AL2(NW+2), B(NW+4)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        B(1) = 0.
        B(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 6) CALL MPDEB ('MPLOG I', A)
C
      IA = SIGN (1., A(1))
      NA = MIN (INT (ABS (A(1))), NW)
C
      IF (IA .LT. 0 .OR. NA .EQ. 0) THEN
        IF (KER(50) .NE. 0) THEN
          WRITE (LDB, 1)
 1        FORMAT ('*** MPLOG: Argument is less than or equal to zero.')
          IER = 50
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
C  Unless the input is close to 2, Log (2) must have been precomputed.
C
      CALL MPMDC (A, T1, N1)
      IF (ABS (T1 - 2.D0) .GT. 1D-3 .OR. N1 .NE. 0) THEN
        CALL MPMDC (AL2, T2, N2)
        IF (N2 .NE. - NBT .OR. ABS (T2 * 0.5D0 ** NBT - ALT) .GT. RX2)  
     $    THEN
          IF (KER(51) .NE. 0) THEN
            WRITE (LDB, 2)
 2          FORMAT ('*** MPLOG: LOG (2) must be precomputed.')
            IER = 51
            IF (KER(IER) .EQ. 2) CALL MPABRT
          ENDIF
          RETURN
        ENDIF
      ENDIF
C
C   Check if input is exactly one.
C
      IF (A(1) .EQ. 1. .AND. A(2) .EQ. 0. .AND. A(3) .EQ. 1.) THEN
        B(1) = 0.
        B(2) = 0.
        GOTO 120
      ENDIF
C
      N5 = NW + 5
      NS = 3 * N5
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N5
      K2 = K1 + N5
      NWS = NW
C
C   Determine the least integer MQ such that 2 ^ MQ .GE. NW.
C
      T2 = NWS
      MQ = CL2 * LOG (T2) + 1.D0 - RXX
C
C   Compute initial approximation of Log (A).
C
      T1 = LOG (T1) + N1 * ALT
      CALL MPDMC (T1, 0, B)
      NW = 3
      IQ = 0
C
C   Perform the Newton-Raphson iteration described above with a dynamically
C   changing precision level NW (one greater than powers of two).
C
      DO 110 K = 2, MQ
        NW = MIN (2 * NW - 2, NWS) + 1
 100    CONTINUE
        CALL MPEXP (B, AL2, S(K0))
        CALL MPSUB (A, S(K0), S(K1))
        CALL MPDIV (S(K1), S(K0), S(K2))
        CALL MPADD (B, S(K2), S(K1))
        CALL MPEQ (S(K1), B)
        IF (K .EQ. MQ - NIT .AND. IQ .EQ. 0) THEN
          IQ = 1
          GOTO 100
        ENDIF
 110  CONTINUE
C
C   Restore original precision level.
C
      NW = NWS
      ICS = ISS
      CALL MPROUN (B)
C
 120  IF (IDB .GE. 6) CALL MPDEB ('MPLOG O', B)
C
      RETURN
      END
C
      SUBROUTINE MPLOGX (A, PI, AL2, B)
C
C   This computes the natural logarithm of the MP number A and returns the MP
C   result in B.  PI is the MP value of Pi produced by a prior call to MPPI or
C   MPPIX.  AL2 is the MP value of Log(2) produced by a prior call to MPLOG
C   or MPLOGX.  Before calling MPLOGX, the array in MPCOM5 must be
C   initialized by calling MPINIX.  For modest levels of precision, use MPLOG.
C   NW should be a power of two.  The last three words of the result are not
C   reliable.  Debug output starts with IDB = 6.
C
C   Max SP space for B: NW + 4 cells.  Max SP scratch space: 10.5 * NW + 51
C   cells.  Max DP scratch space: 12 * NW + 6 cells.
C
C   This uses the following algorithm, which is due to Salamin.  If a is
C   extremely close to 1, use a Taylor series.  Otherwise select n such that
C   z = x 2^n is at least 2^m, where m is the number of bits of desired
C   precision in the result.  Then
C
C   Log(x) = Pi / [2 AGM (1, 4/x)]
C
      DOUBLE PRECISION ALT, CPI, ST, T1, T2, TN
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (MZL = -5, ALT = 0.693147180559945309D0,                
     $  CPI = 3.141592653589793D0)
      DIMENSION AL2(NW+2), F1(8), F4(8), PI(NW+2), A(NW+4), B(NW+4)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        B(1) = 0.
        B(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 6) CALL MPDEB ('MPLOGX I', A)
C
      IA = SIGN (1., A(1))
      NA = MIN (INT (ABS (A(1))), NW)
      NCR = 2 ** MCR
C
C   Check if precision level is too low to justify the advanced routine.
C
      IF (NW .LE. NCR) THEN
        CALL MPLOG (A, AL2, B)
        GOTO 120
      ENDIF
C
      IF (IA .LT. 0 .OR. NA .EQ. 0) THEN
C
C   Input is less than or equal to zero.
C
        IF (KER(52) .NE. 0) THEN
          WRITE (LDB, 1)
 1        FORMAT ('*** MPLOGX: Argument is less than or equal to zero.')
          IER = 52
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
C   Check if Pi has been precomputed.
C
      CALL MPMDC (PI, T1, N1)
      IF (N1 .NE. 0 .OR. ABS (T1 - CPI) .GT. RX2) THEN
        IF (KER(53) .NE. 0) THEN
          WRITE (LDB, 2)
 2        FORMAT ('*** MPLOGX: PI must be precomputed.')
          IER = 53
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
C  Unless the input is 2, Log (2) must have been precomputed.
C
      IF (A(1) .NE. 1. .OR. A(2) .NE. 0. .OR. A(3) .NE. 2.) THEN
        IT2 = 0
        CALL MPMDC (AL2, T2, N2)
        IF (N2 .NE. - NBT .OR. ABS (T2 * 0.5D0 ** NBT - ALT) .GT. RX2)  
     $    THEN
          IF (KER(54) .NE. 0) THEN
            WRITE (LDB, 3)
 3          FORMAT ('*** MPLOGX: Log (2) must be precomputed.')
            IER = 54
            IF (KER(IER) .EQ. 2) CALL MPABRT
          ENDIF
          RETURN
        ENDIF
      ELSE
        IT2 = 1
      ENDIF
C
C   Define sections of the scratch array.
C
      N4 = NW + 4
      NS = 4 * N4
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N4
      K2 = K1 + N4
      K3 = K2 + N4
      F1(1) = 1.
      F1(2) = 0.
      F1(3) = 1.
      F4(1) = 1.
      F4(2) = 0.
      F4(3) = 4.
C
C   If argument is 1, the result is zero.  If the argument is extremely close
C   to 1.  If so, employ a Taylor's series instead.
C
      CALL MPSUB (A, F1, S(K0))
      IF (S(K0) .EQ. 0.) THEN
        B(1) = 0.
        B(2) = 0.
        GOTO 110
      ELSEIF (S(K0+1) .LE. MZL) THEN
        CALL MPEQ (S(K0), S(K1))
        CALL MPEQ (S(K1), S(K2))
        I1 = 1
        IS = 1
        TL = S(K0+1) - NW - 1
C
 100    I1 = I1 + 1
        IS = - IS
        ST = IS * I1
        CALL MPMULX (S(K1), S(K2), S(K3))
        CALL MPDIVD (S(K3), ST, 0, S(K2))
        CALL MPADD (S(K0), S(K2), S(K3))
        CALL MPEQ (S(K3), S(K0))
        IF (S(K2+1) .GE. TL) GOTO 100
C
        CALL MPEQ (S(K0), B)
        GOTO 110
      ENDIF
C
C   If input is exactly 2, set the exponent to a large value.  Otherwise
C   multiply the input by a large power of two.
C
      CALL MPMDC (A, T1, N1)
      N2 = NBT * (NW / 2 + 2) - N1
      TN = N2
      IF (IT2 .EQ. 1) THEN
        CALL MPDMC (1.D0, N2, S(K0))
      ELSE
        CALL MPMULD (A, 1.D0, N2, S(K0))
      ENDIF
C
C   Perform AGM iterations.
C
      CALL MPEQ (F1, S(K1))
      CALL MPDIVX (F4, S(K0), S(K2))
      CALL MPAGMX (S(K1), S(K2))
C
C   Compute B = Pi / (2 * A), where A is the limit of the AGM iterations.
C
      CALL MPMULD (S(K1), 2.D0, 0, S(K0))
      CALL MPDIVX (PI, S(K0), S(K1))
C
C  If the input was exactly 2, divide by TN.  Otherwise subtract TN * Log(2).
C
      IF (IT2 .EQ. 1) THEN
        CALL MPDIVD (S(K1), TN, 0, S(K0))
      ELSE
        CALL MPMULD (AL2, TN, 0, S(K2))
        CALL MPSUB (S(K1), S(K2), S(K0))
      ENDIF
      CALL MPEQ (S(K0), B)
C
 110  ICS = ISS
 120  IF (IDB .GE. 6) CALL MPDEB ('MPLOGX O', B)
      RETURN
      END
C
      SUBROUTINE MPMDC (A, B, N)
C
C   This converts the MP number A to the DPE form (B, N), accurate to between
C   14 and 17 digits, depending on system.  B will be between 1 and BDX.
C   Debug output starts with IDB = 9.
C
      DOUBLE PRECISION AA, B
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (NDB = 22)
      DIMENSION A(NW+2)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
C
      IF (IER .NE. 0) THEN
        B = 0.D0
        N = 0
        RETURN
      ENDIF
      IF (IDB .GE. 9) THEN
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 1) (A(I), I = 1, NO)
 1      FORMAT ('MPMDC I'/(6F12.0))
      ENDIF
C
      IF (A(1) .EQ. 0.)  THEN
        B = 0.D0
        N = 0
        GOTO 100
      ENDIF
C
      NA = ABS (A(1))
      AA = A(3)
      IF (NA .GE. 2) AA = AA + RDX * A(4)
      IF (NA .GE. 3) AA = AA + RX2 * A(5)
      IF (NA .GE. 4) AA = AA + RDX * RX2 * A(6)
C
      N = NBT * A(2)
      B = SIGN (AA, DBLE (A(1)))
C
 100  IF (IDB .GE. 9) WRITE (LDB, 2) B, N
 2    FORMAT ('MPMDC O',F10.0,I10)
      RETURN
      END
C
      SUBROUTINE MPMMPC (A, B, L, C)
C
C   This converts MP numbers A and B to MPC form in C, i.e. C = A + B i.
C   L (an input parameter) is the offset between real and imaginary parts in
C   C.  Debug output starts with IDB = 10.
C
C   Max SP space for C: 2 * L cells.
C
      DIMENSION A(NW+2), B(NW+2), C(2*L)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
C
      IF (IER .NE. 0) THEN
        C(1) = 0.
        C(2) = 0.
        C(L+1) = 0.
        C(L+2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 10) WRITE (LDB, 1)
 1    FORMAT ('MPMMPC')
C
      I1 = SIGN (1., A(1))
      N1 = MIN (INT (ABS (A(1))), NW, L - 2)
      I2 = SIGN (1., B(1))
      N2 = MIN (INT (ABS (B(1))), NW, L - 2)
      C(1) = SIGN (N1, I1)
      C(L+1) = SIGN (N2, I2)
C
      DO 100 I = 2, N1 + 2
        C(I) = A(I)
 100  CONTINUE
C
      DO 110 I = 2, N2 + 2
        C(L+I) = B(I)
 110  CONTINUE
C
      RETURN
      END
C
      SUBROUTINE MPMOUT (N1, N2, A)
C
C   This produces a compact printout of the N1 x N1 MP array A.  It is called
C   MPINRQ and MPINRZ.  It is not indended to be called directly by the user.
C
      DOUBLE PRECISION T1, T2
      DIMENSION A(NW+4,N1,N2)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      DIMENSION T1(100), I1(100)
C
      IF (IER .NE. 0) RETURN
C
      DO 110 J = 1, N1
        WRITE (LDB, 1) J
 1      FORMAT ('Row', I3)
C
        DO 100 K = 1, N2
          CALL MPMDC (A(1,J,K), T2, M2)
          CALL DPDEC (T2, M2, T1(K), I1(K))
 100    CONTINUE
C
        WRITE (LDB, 2) (T1(K), I1(K), K = 1, N2)
 2      FORMAT (4(F10.6,I6))
 110  CONTINUE
C
      RETURN
      END
C
      SUBROUTINE MPMPCM (L, A, B, C)
C
C   This converts the MPC number A to its MP real and imaginary parts, i.e.
C   B = Real (A) and C = Imag (A).  L is the offset between real and
C   imaginary parts in A.  Debug output starts with IDB = 10.
C
C   Max SP space for B and C: NW + 2 cells.
C
      DIMENSION A(2*L), B(NW+2), C(NW+2)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
C
      IF (IER .NE. 0) THEN
        B(1) = 0.
        B(2) = 0.
        C(1) = 0.
        C(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 10) WRITE (LDB, 1)
 1    FORMAT ('MPMPCM')
C
      I1 = SIGN (1., A(1))
      N1 = MIN (INT (ABS (A(1))), NW, L - 2)
      I2 = SIGN (1., A(L+1))
      N2 = MIN (INT (ABS (A(L+1))), NW, L - 2)
      B(1) = SIGN (N1, I1)
      C(1) = SIGN (N2, I2)
C
      DO 100 I = 2, N1 + 2
        B(I) = A(I)
 100  CONTINUE
C
      DO 110 I = 2, N2 + 2
        C(I) = A(L+I)
 110  CONTINUE
C
      RETURN
      END
C
      SUBROUTINE MPMUL (A, B, C)
C
C   This routine multiplies MP numbers A and B to yield the MP product C.
C   When one of the arguments has a much higher level of precision than the
C   other, this routine is slightly more efficient if A has the lower level of
C   precision.  For extra high levels of precision, use MPMULX.  Debug output
C   starts with IDB = 8.
C
C   Max SP space for C: NW + 4 cells.  Max DP scratch space: NW + 4 cells.
C
C   This routine returns up to NW mantissa words of the product.  If the
C   complete double-long product of A and B is desired (for example in large
C   integer applications), then NW must be at least as large as the sum of the
C   mantissa lengths of A and B.  In other words, if the precision levels of A
C   and B are both 64 words, then NW must be at least 128 words to obtain the
C   complete double-long product in C.
C
      DOUBLE PRECISION D, T1, T2
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (NDB = 22)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM4/ D(1024)
      DIMENSION A(NW+2), B(NW+2), C(NW+4)
C
      IF (IER .NE. 0) THEN
        C(1) = 0.
        C(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 8) THEN
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 1) (A(I), I = 1, NO)
 1      FORMAT ('MPMUL I'/(6F12.0))
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 1) (B(I), I = 1, NO)
      ENDIF
C
      IA = SIGN (1., A(1))
      IB = SIGN (1., B(1))
      NA = MIN (INT (ABS (A(1))), NW)
      NB = MIN (INT (ABS (B(1))), NW)
      IF (NA .EQ. 0 .OR. NB .EQ. 0) THEN
C
C   One of the inputs is zero -- result is zero.
C
        C(1) = 0.
        C(2) = 0.
        GOTO 170
      ENDIF
      IF (NA .EQ. 1 .AND. A(3) .EQ. 1.) THEN
C
C   A is 1 or -1 -- result is B or -B.
C
        C(1) = SIGN (NB, IA * IB)
        C(2) = A(2) + B(2)
C
        DO 100 I = 3, NB + 2
          C(I) = B(I)
 100    CONTINUE
C
        GOTO 170
      ELSEIF (NB .EQ. 1 .AND. B(3) .EQ. 1.) THEN
C
C   B is 1 or -1 -- result is A or -A.
C
        C(1) = SIGN (NA, IA * IB)
        C(2) = A(2) + B(2)
C
        DO 110 I = 3, NA + 2
          C(I) = A(I)
 110    CONTINUE
C
        GOTO 170
      ENDIF
C
      NC = MIN (NA + NB, NW)
      D2 = A(2) + B(2)
C
      DO 120 I = 1, NC + 4
        D(I) = 0.D0
 120  CONTINUE
C
C   Perform ordinary long multiplication algorithm.  Accumulate at most NW + 4
C   mantissa words of the product.
C
      DO 150 J = 3, NA + 2
        T1 = A(J)
        J3 = J - 3
        N2 = MIN (NB + 2, NW + 4 - J3)
C
        DO 130 I = 3, N2
          D(I+J3) = D(I+J3) + T1 * B(I)
 130    CONTINUE
C
C   Release carries periodically to avoid overflowing the exact integer
C   capacity of double precision floating point words in D.
C
        IF (MOD (J - 2, NPR) .EQ. 0) THEN
          I1 = MAX (3, J - NPR)
          I2 = N2 + J3
C
CDIR$ IVDEP
          DO 140 I = I1, I2
            T1 = D(I)
            T2 = INT (RDX * T1)
            D(I) = T1 - BDX * T2
            D(I-1) = D(I-1) + T2
 140      CONTINUE
C
        ENDIF
 150  CONTINUE
C
C   If D(2) is nonzero, shift the result one cell right.
C
      IF (D(2) .NE. 0.D0) THEN
        D2 = D2 + 1.
C
CDIR$ IVDEP
        DO 160 I = NC + 4, 3, -1
          D(I) = D(I-1)
 160    CONTINUE
C
      ENDIF
      D(1) = SIGN (NC, IA * IB)
      D(2) = D2
C
C   Fix up result, since some words may be negative or exceed BDX.
C
      CALL MPNORM (C)
C
 170  IF (IDB .GE. 8) THEN
        NO = MIN (INT (ABS (C(1))), NDB) + 2
        WRITE (LDB, 2) (C(I), I = 1, NO)
 2      FORMAT ('MPMUL O'/(6F12.0))
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPMULD (A, B, N, C)
C
C   This routine multiplies the MP number A by the DPE number (B, N) to yield
C   the MP product C.  Debug output starts with IDB = 9.
C
C   Max SP space for C: NW + 4 cells.  Max DP space: NW + 4 cells.
C
      DOUBLE PRECISION B, BB, D
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (NDB = 22)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM3/ S(1024)
      COMMON /MPCOM4/ D(1024)
      DIMENSION A(NW+2), C(NW+4), F(8)
C
      IF (IER .NE. 0) THEN
        C(1) = 0.
        C(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 9) THEN
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 1) (A(I), I = 1, NO)
 1      FORMAT ('MPMULD I'/(6F12.0))
        WRITE (LDB, 2) B, N
 2      FORMAT ('MPMULD I',1PD25.15,I10)
      ENDIF
C
C   Check for zero inputs.
C
      IA = SIGN (1., A(1))
      NA = MIN (INT (ABS (A(1))), NW)
      IB = SIGN (1.D0, B)
      IF (NA .EQ. 0 .OR. B .EQ. 0.D0) THEN
        C(1) = 0.
        C(2) = 0.
        GOTO 140
      ENDIF
      N1 = N / NBT
      N2 = N - NBT * N1
      BB = ABS (B) * 2.D0 ** N2
C
C   Reduce BB to within 1 and BDX.
C
      IF (BB .GE. BDX) THEN
C
        DO 100 K = 1, 100
          BB = RDX * BB
          IF (BB .LT. BDX) THEN
            N1 = N1 + K
            GOTO 120
          ENDIF
 100    CONTINUE
C
      ELSEIF (BB .LT. 1.D0) THEN
C
        DO 110 K = 1, 100
          BB = BDX * BB
          IF (BB .GE. 1.D0) THEN
            N1 = N1 - K
            GOTO 120
          ENDIF
 110    CONTINUE
C
      ENDIF
C
C   If B cannot be represented exactly in a single mantissa word, use MPMUL.
C
 120  IF (BB .NE. AINT (BB)) THEN
        BB = SIGN (BB, B)
        CALL MPDMC (BB, N1 * NBT, F)
        CALL MPMUL (F, A, C)
        GOTO 140
      ENDIF
C
C   Perform short multiply operation.
C
CDIR$ IVDEP
      DO 130 I = 3, NA + 2
        D(I) = BB * A(I)
 130  CONTINUE
C
C   Set the exponent and fix up the result.
C
      D(1) = SIGN (NA, IA * IB)
      D(2) = A(2) + N1
      D(NA+3) = 0.D0
      D(NA+4) = 0.D0
      CALL MPNORM (C)
C
 140  IF (IDB .GE. 9) THEN
        NO = MIN (INT (ABS (C(1))), NDB) + 2
        WRITE (LDB, 3) (C(I), I = 1, NO)
 3      FORMAT ('MPMULD O'/(6F12.0))
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPMULX (A, B, C)
C
C   This routine multiplies MP numbers A and B to yield the MP product C.
C   Before calling MPMULX, the array in MPCOM5 must be initialized by calling
C   MPINIX.  For modest levels of precision, use MPMUL.  NW should be a power
C   of two.  Debug output starts with IDB = 8.
C
C   Max SP space for C: NW + 4 cells.  Max DP scratch space: 12 * NW + 6 cells.
C   The fact that all advanced routines require this amount of DP scratch
C   space derives from the requirement in this routine, which all of them call.
C
C   This routine returns up to NW mantissa words of the product.  If the
C   complete double-long product of A and B is desired (for example in large
C   integer applications), then NW must be at least as large as the sum of the
C   mantissa lengths of A and B.  In other words, if the precision levels of A
C   and B are both 256 words, then NW must be at least 512 words to obtain the
C   complete double-long product in C.
C
C   This subroutine uses an advanced technique involving the fast Fourier
C   transform (FFT).  For high precision it is significantly faster than the
C   conventional scheme used in MPMUL.
C>
C   Two machine-dependent parameters are set in this routine:
C
C     ERM   Maximum tolerated FFT roundoff error.  On IEEE systems ERM =
C           0.438D0.  It is not necessary to specify ERM for modest levels of
C           precision -- see comments below.
C     MBT   Number of mantissa bits in double precision data.  MBT = 53 on
C           IEEE systems, and MBT = 48 (i.e. single precision) on Crays.
C           It is not necessary to specify MBT for modest levels of precision.
C
      DOUBLE PRECISION AN, CL2, D, ERM, T1, T2, T3, T4
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (CL2 = 1.4426950408889633D0, ERM = 0.438D0, MBT = 53,   
     $  NDB = 22)
      DIMENSION A(NW+2), B(NW+2), C(NW+4)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM4/ D(1024)
C
      IF (IER .NE. 0) THEN
        C(1) = 0.
        C(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 8)  THEN
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 1) (A(I), I = 1, NO)
 1      FORMAT ('MPMULX I'/(6F12.0))
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 1) (B(I), I = 1, NO)
      ENDIF
C
      IA = SIGN (1., A(1))
      IB = SIGN (1., B(1))
      NA = MIN (INT (ABS (A(1))), NW)
      NB = MIN (INT (ABS (B(1))), NW)
      NCR = 2 ** MCR
C
      IF (NA .EQ. 0 .OR. NB .EQ. 0) THEN
C
C   One of the inputs is zero -- result is zero.
C
        C(1) = 0.
        C(2) = 0.
        GOTO 190
      ENDIF
C
C   Check if precision level of one of the arguments is too low to justify the
C   advanced routine.
C
      IF (NA .LE. NCR .OR. NB .LE. NCR) THEN
        CALL MPMUL (A, B, C)
        GOTO 190
      ENDIF
C
C   Determine N1, the smallest power of two at least as large as NA and NB.
C
      T1 = NA
      T2 = NB
      M1 = CL2 * LOG (T1) + 1.D0 - RXX
      M2 = CL2 * LOG (T2) + 1.D0 - RXX
      M1 = MAX (M1, M2)
      N1 = 2 ** M1
      M2 = M1 + 2
      N2 = 2 * N1
      N4 = 2 * N2
      N6 = 3 * N2
      N8 = 4 * N2
      N21 = N2 + 1
      N42 = N4 + 2
      N63 = N6 + 3
      N84 = N8 + 4
C
C   Place the input data in A and B into separate sections of the scratch
C   array D.  This code also splits the input data into half-sized words.
C
CDIR$ IVDEP
      DO 100 I = 1, NA
        T1 = A(I+2)
        T2 = INT (RBX * T1)
        D(2*I-1) = T2
        D(2*I) = T1 - BBX * T2
 100  CONTINUE
C
      DO 110 I = 2 * NA + 1, N2
        D(I) = 0.D0
 110  CONTINUE
C
CDIR$ IVDEP
      DO 120 I = 1, NB
        T1 = B(I+2)
        T2 = INT (RBX * T1)
        D(2*I-1+N42) = T2
        D(2*I+N42) = T1 - BBX * T2
 120  CONTINUE
C
      DO 130 I = 2 * NB + 1, N2
        D(I+N42) = 0.D0
 130  CONTINUE
C
C   Set the second half of each input vector in D to zero.
C
CDIR$ IVDEP
      DO 140 I = N2 + 1, N4
        D(I) = 0.D0
        D(I+N42) = 0.D0
 140  CONTINUE
C
C   Perform forward real-to-complex FFTs on the two vectors in D.  The complex
C   results are placed in (D(I), I = 1, N4+2) and (D(I), I = N4 + 3, N8 + 4).
C
      CALL MPRCFT (1, M2, D, D(N84+1))
      CALL MPRCFT (1, M2, D(N42+1), D(N84+1))
C
C   Multiply the resulting complex vectors.
C
CDIR$ IVDEP
      DO 150 I = 1, N21
        T1 = D(I)
        T2 = D(I+N21)
        T3 = D(I+N42)
        T4 = D(I+N63)
        D(I+N42) = T1 * T3 - T2 * T4
        D(I+N63) = T1 * T4 + T2 * T3
 150  CONTINUE
C
C   Perform an inverse complex-to-real FFT on the resulting data.
C
      CALL MPCRFT (-1, M2, D(N42+1), D(N84+1))
C
C   Divide by N8, recombine words and release carries.
C
      NC = MIN (NA + NB, NW)
      NC1 = MIN (NW + 1, NA + NB - 1)
      D(1) = SIGN (NC, IA * IB)
      D(2) = A(2) + B(2) + 1
      AN = 1.D0 / N8
      T1 = AN * D(N42+1)
      D(3) = AINT (T1 + 0.5D0)
      D(NC+3) = 0.D0
      D(NC+4) = 0.D0
      D(N42+1) = 0.D0
C
CDIR$ IVDEP
      DO 160 I = 1, NC1
        T1 = AN * D(N42+2*I)
        T2 = AN * D(N42+2*I+1)
        T3 = AINT (T1 + 0.5D0)
        T4 = AINT (T2 + 0.5D0)
C        D(N42+2*I) = ABS (T3 - T1)
C        D(N42+2*I+1) = ABS (T4 - T2)
        T1 = INT (RDX * T3)
        T2 = T3 - BDX * T1
        T3 = INT (RDX * T4)
        T4 = T4 - BDX * T3
        D(I+3) = BBX * T2 + T4
        D(I+2) = D(I+2) + BBX * T1 + T3
 160  CONTINUE
C
C   Find the largest FFT roundoff error.  Roundoff error is minimal unless
C   exceedingly high precision (i.e. over one million digits) is used.  Thus
C   this test may be disabled in normal use.  To disable this test, uncomment
C   the next line of code and comment out the two lines of the previous loop
C   that begin D(N42..
C
C   This code can be used as a rigorous system integrity test.  First set
C   MBT according to the system being used, and then set ERM to be fairly
C   small, say 0.001 or whatever is somewhat larger than the largest FFT
C   roundoff error typically encountered for a given precision level on the
C   computer being used.  Enable this test as explained in the previous
C   paragraph.  Then if an anomalously large roundoff error is detected, a
C   hardware or compiler error has likely occurred.
C
      GOTO 180
      T1 = 0.D0
C
      DO 170 I = 1, 2 * NC1 + 1
        IF (D(N42+I) .GT. T1) THEN
          I1 = I
          T1 = D(N42+I)
        ENDIF
 170  CONTINUE
C
C   Check if maximum roundoff error exceeds the limit ERM, which is set above.
C   Also determine the number of fractional bits and how large the error is in
C   terms of units in the last place (ulp).
C
      IF (T1 .GT. ERM)  THEN
        IF (KER(55) .NE. 0) THEN
          T2 = AN * D(I1)
          I2 = CL2 * LOG (T1) + 1.D0 + RXX
          I3 = CL2 * LOG (T2) + 1.D0 + RXX
          I4 = MBT + I2 - I3
          I5 = T1 * 2 ** I4 + RXX
          WRITE (LDB, 2) I1, T1, I4, I5
 2        FORMAT ('*** MPMULX: Excessive FFT roundoff error',I10,F10.6, 
     $      2I6)
          IER = 55
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
      ENDIF
C
C   Fix up the result.
C
 180  CALL MPNORM (C)
C
 190  IF (IDB .GE. 8) THEN
        NO = MIN (INT (ABS (C(1))), NDB) + 2
        WRITE (LDB, 3) (C(I), I = 1, NO)
 3      FORMAT ('MPMULX O'/(6F12.0))
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPNINT (A, B)
C
C   This sets B equal to the integer nearest to the MP number A.  Debug output
C   starts with IDB = 8.
C
C   Max SP space for B: NW + 4 cells.  Max SP scratch space: NW + 4 cells.
C
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (NDB = 22)
      DIMENSION A(NW+2), B(NW+2), F(8)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        B(1) = 0.
        B(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 8) THEN
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 1) (A(I), I = 1, NO)
 1      FORMAT ('MPNINT I'/(6F12.0))
      ENDIF
C
      IA = SIGN (1., A(1))
      NA = MIN (INT (ABS (A(1))), NW)
      MA = A(2)
      IF (NA .EQ. 0)  THEN
C
C   A is zero -- result is zero.
C
        B(1) = 0.
        B(2) = 0.
        GOTO 110
      ENDIF
      IF (MA .GE. NW) THEN
C
C   A cannot be represented exactly as an integer.
C
        IF (KER(56) .NE. 0) THEN
          WRITE (LDB, 2)
 2        FORMAT ('*** MPNINT: Argument is too large.')
          IER = 56
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
      NS = NW + 4
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      F(1) = 1.
      F(2) = -1.
      F(3) = 0.5D0 * BDX
C
C   Add or subtract 1/2 from the input, depending on its sign.
C
      IF (IA .EQ. 1) THEN
        CALL MPADD (A, F, S(K0))
      ELSE
        CALL MPSUB (A, F, S(K0))
      ENDIF
      IC = SIGN (1., S(K0))
      NC = ABS (S(K0))
      MC = S(K0+1)
C
C   Place integer part of S in B.
C
      NB = MIN (MAX (MC + 1, 0), NC)
      IF (NB .EQ. 0) THEN
        B(1) = 0.
        B(2) = 0.
      ELSE
        B(1) = SIGN (NB, IC)
        B(2) = MC
        B(NB+3) = 0.
        B(NB+4) = 0.
C
        DO 100 I = 3, NB + 2
          B(I) = S(I+K0-1)
 100    CONTINUE
C
      ENDIF
      ICS = ISS
C
 110  IF (IDB .GE. 8) THEN
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 3) (B(I), I = 1, NO)
 3      FORMAT ('MPNINT O'/(6F12.0))
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPNORM (A)
C
C   This converts the MP number in array D of MPCOM4 to the standard
C   normalized form in A.  The MP routines often leave negative numbers or
C   values exceeding the radix BDX in result arrays, and this fixes them.
C   MPNORM assumes that two extra mantissa words are input at the end of D.
C   This reduces precision loss when it is necessary to shift the result to
C   the left.  This routine is not intended to be called directly by the user.
C   The output is placed in the SP array A.  Debug output starts with IDB = 10.
C
C   Max SP space for A: NW + 4 cells.
C
      DOUBLE PRECISION D, R1, S1, T1, T2
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (NDB = 22)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM4/ D(1024)
      DIMENSION A(NW+4)
C
      IF (IER .NE. 0) THEN
        A(1) = 0.
        A(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 10) THEN
        NO = MIN (INT (ABS (D(1))), NDB) + 4
        WRITE (LDB, 1) (D(I), I = 1, NO)
 1      FORMAT ('MPNORM I'/(4F18.0))
      ENDIF
C
      IA = SIGN (1.D0, D(1))
      NA = MIN (INT (ABS (D(1))), NW)
      IF (NA .EQ. 0)  GOTO 170
      N4 = NA + 4
      A2 = D(2)
      D(2) = 0.D0
      R1 = 2.D0 + 0.125D0 * RDX
C>
C   Try a vectorized fixup loop three times, unless A is very short.  This
C   should handle 99% of the inputs.  On scalar computers, it is more
C   efficient to completely bypass this loop, by uncommenting the next line.
C
      GOTO 120
      IF (NA .LE. 8) GOTO 120
C
      DO 110 K = 1, 3
        S1 = 0.D0
C
CDIR$ IVDEP
        DO 100 I = 3, N4
          T1 = INT (D(I) * RDX + R1) - 2.D0
          D(I) = D(I) - T1 * BDX
          D(I-1) = D(I-1) + T1
          S1 = S1 + ABS (T1)
 100    CONTINUE
C
        IF (S1 .EQ. 0.D0) GOTO 140
 110  CONTINUE
C
C   Still not fixed - use recursive loop.  This loop is not vectorizable,
C   but it is guaranteed to complete the job in one pass.
C
 120  T1 = 0.D0
C
      DO 130 I = N4, 3, -1
        T2 = T1 + D(I)
        T1 = INT (T2 * RDX + R1) - 2.D0
        D(I) = T2 - T1 * BDX
 130  CONTINUE
C
      D(2) = D(2) + T1
C
 140  IF (D(2) .NE. 0.) THEN
C
C   The fixup loops above "spilled" a nonzero number into D(2).  Shift the
C   entire number right one cell.  The exponent and length of the result
C   are increased by one.
C
        DO 150 I = N4, 3, -1
          A(I) = D(I-1)
 150    CONTINUE
C
        NA = MIN (NA + 1, NW)
        A2 = A2 + 1.
      ELSE
C
        DO 160 I = 3, N4
          A(I) = D(I)
 160    CONTINUE
C
      ENDIF
C
C   Perform rounding and truncation.
C
      A(1) = SIGN (NA, IA)
      A(2) = A2
      CALL MPROUN (A)
C
 170  IF (IDB .GE. 10) THEN
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 2) (A(I), I = 1, NO)
 2      FORMAT ('MPNORM O'/(6F12.0))
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPNPWR (A, N, B)
C
C   This computes the N-th power of the MP number A and returns the MP result
C   in B.  When N is zero, 1 is returned.  When N is negative, the reciprocal
C   of A ^ |N| is returned.  For extra high levels of precision, use MPNPWX.
C   Debug output starts with IDB = 7.
C
C   Max SP space for B: NW + 4 cells.  Max SP scratch space: 2 * NW + 10
C   cells.  Max DP scratch space: NW + 5 cells.
C
C   This routine employs the binary method for exponentiation.
C
      DOUBLE PRECISION CL2, T1
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (CL2 = 1.4426950408889633D0, NDB = 22)
      DIMENSION A(NW+2), B(NW+4), F1(8)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        B(1) = 0.
        B(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 7) THEN
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 1) N, (A(I), I = 1, NO)
 1      FORMAT ('MPNPWR I',I5/(6F12.0))
      ENDIF
C
      NA = MIN (INT (ABS (A(1))), NW)
      IF (NA .EQ. 0) THEN
        IF (N .GE. 0) THEN
          B(1) = 0.
          B(2) = 0.
          GOTO 120
        ELSE
          IF (KER(57) .NE. 0) THEN
            WRITE (LDB, 2)
 2          FORMAT ('*** MPNPWR: Argument is zero and N is negative or',
     $        ' zero.')
            IER = 57
            IF (KER(IER) .EQ. 2) CALL MPABRT
          ENDIF
          RETURN
        ENDIF
      ENDIF
C
      N5 = NW + 5
      NS = 2 * N5
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N5
      NWS = NW
      NW = NW + 1
      NN = ABS (N)
      F1(1) = 1.
      F1(2) = 0.
      F1(3) = 1.
      IF (NN .EQ. 0) THEN
        CALL MPEQ (F1, B)
        NW = NWS
        ICS = ISS
        GOTO 120
      ELSEIF (NN .EQ. 1) THEN
        CALL MPEQ (A, B)
        GOTO 110
      ELSEIF (NN .EQ. 2) THEN
        CALL MPMUL (A, A, S(K0))
        CALL MPEQ (S(K0), B)
        GOTO 110
      ENDIF
C
C   Determine the least integer MN such that 2 ^ MN .GT. NN.
C
      T1 = NN
      MN = CL2 * LOG (T1) + 1.D0 + RXX
      CALL MPEQ (F1, B)
      CALL MPEQ (A, S(K0))
      KN = NN
C
C   Compute B ^ N using the binary rule for exponentiation.
C
      DO 100 J = 1, MN
        KK = KN / 2
        IF (KN .NE. 2 * KK) THEN
          CALL MPMUL (B, S(K0), S(K1))
          CALL MPEQ (S(K1), B)
        ENDIF
        KN = KK
        IF (J .LT. MN) THEN
          CALL MPMUL (S(K0), S(K0), S(K1))
          CALL MPEQ (S(K1), S(K0))
        ENDIF
 100  CONTINUE
C
C   Compute reciprocal if N is negative.
C
 110  IF (N .LT. 0) THEN
        CALL MPDIV (F1, B, S(K0))
        CALL MPEQ (S(K0), B)
      ENDIF
C
C   Restore original precision level.
C
      NW = NWS
      ICS = ISS
      CALL MPROUN (B)
C
 120  IF (IDB .GE. 7) THEN
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 3) (B(I), I = 1, NO)
 3      FORMAT ('MPNPWR O'/(6F12.0))
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPNPWX (A, N, B)
C
C   This computes the N-th power of the MP number A and returns the MP result
C   in B.  When N is zero, 1 is returned.  When N is negative, the reciprocal
C   of A ^ |N| is returned.  Before calling MPNPWX, the array in MPCOM5 must
C   be initialized by calling MPINIX.  For modest levels of precision, use
C   MPNPWR.  NW should be a power of two.  The last two words of the result
C   are not reliable.  Debug output starts with IDB = 6.
C
C   Max SP space for B: NW + 4 cells.  Max SP scratch space: 5 * NW + 20
C   cells.  Max DP scratch space: 12 * NW + 6 cells.
C
C   This routine employs the binary method for exponentiation.
C
      DOUBLE PRECISION CL2, T1
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (CL2 = 1.4426950408889633D0, NDB = 22)
      DIMENSION A(NW+2), B(NW+4), F1(8)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        B(1) = 0.
        B(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 6) THEN
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 1) N, (A(I), I = 1, NO)
 1      FORMAT ('MPNPWX I',I5/(6F12.0))
      ENDIF
C
      NCR = 2 ** MCR
      NA = MIN (INT (ABS (A(1))), NW)
C
C   Check if precision level of A is too low to justify the advanced routine.
C
      IF (NA .LE. NCR) THEN
        CALL MPNPWR (A, N, B)
        GOTO 120
      ENDIF
      IF (NA .EQ. 0) THEN
        IF (N .GE. 0) THEN
          B(1) = 0.
          B(2) = 0.
          GOTO 120
        ELSE
          IF (KER(58) .NE. 0) THEN
            WRITE (LDB, 2)
 2          FORMAT ('*** MPNPWX: argument is zero and N is negative or',
     $        ' zero.')
            IER = 58
            IF (KER(IER) .EQ. 2) CALL MPABRT
          ENDIF
          RETURN
        ENDIF
      ENDIF
C
      N4 = NW + 4
      NS = 2 * N4
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N4
      NN = ABS (N)
      F1(1) = 1.
      F1(2) = 0.
      F1(3) = 1.
      IF (NN .EQ. 0) THEN
        CALL MPEQ (F1, B)
        ICS = ISS
        GOTO 120
      ELSEIF (NN .EQ. 1) THEN
        CALL MPEQ (A, B)
        GOTO 110
      ELSEIF (NN .EQ. 2) THEN
        CALL MPSQX (A, B)
        GOTO 110
      ENDIF
C
C   Determine the least integer MN such that 2 ^ MN .GT. NN.
C
      T1 = NN
      MN = CL2 * LOG (T1) + 1.D0 + RXX
      CALL MPEQ (F1, B)
      CALL MPEQ (A, S(K0))
      KN = NN
C
C   Compute B ^ N using the binary rule for exponentiation.
C
      DO 100 J = 1, MN
        KK = KN / 2
        IF (KN .NE. 2 * KK) THEN
          CALL MPMULX (B, S(K0), S(K1))
          CALL MPEQ (S(K1), B)
        ENDIF
        KN = KK
        IF (J .LT. MN) THEN
          CALL MPSQX (S(K0), S(K1))
          CALL MPEQ (S(K1), S(K0))
        ENDIF
 100  CONTINUE
C
C   Compute reciprocal if N is negative.
C
 110  IF (N .LT. 0) THEN
        CALL MPDIVX (F1, B, S(K0))
        CALL MPEQ (S(K0), B)
      ENDIF
      ICS = ISS
C
 120  IF (IDB .GE. 6) THEN
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 3) (B(I), I = 1, NO)
 3      FORMAT ('MPNPWX O'/(6F12.0))
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPNRT (A, N, B)
C
C   This computes the N-th root of the MP number A and returns the MP result
C   in B.  N must be at least one and must not exceed 2 ^ 30.  For extra high
C   levels of precision, use MPNRTX.  Debug output starts with IDB = 7.
C
C   Max SP space for B: NW + 4 cells.  Max SP scratch space: 6 * NW + 32
C   cells.  Max DP scratch space: NW + 6 cells.
C
C   This subroutine employs the following Newton-Raphson iteration, which
C   converges to A ^ (-1/N):
C
C          X_{k+1} = X_k + (X_k / N) * (1 - A * X_k^N)
C
C   The reciprocal of the final approximation to A ^ (-1/N) is the N-th root.
C   These iterations are performed with a maximum precision level NW that
C   is dynamically changed, approximately doubling with each iteration.
C   See the comment about the parameter NIT in MPDIVX.
C
C   When N is large and A is very near one, the following binomial series is
C   employed instead of the Newton scheme:
C
C       (1 + x)^(1/N)  =  1  +  x / N  +  x^2 * (1 - N) / (2! N^2)  +  ...
C
      DOUBLE PRECISION ALT, CL2, T1, T2, TN
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (ALT = 0.693147180559945309D0,                          
     $  CL2 = 1.4426950408889633D0, NDB = 22, NIT = 3, N30 = 2 ** 30)
      DIMENSION A(NW+2), B(NW+4), F1(8), F2(8)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        B(1) = 0.
        B(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 7) THEN
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 1) N, (A(I), I = 1, NO)
 1      FORMAT ('MPNRT I',I5/(6F12.0))
      ENDIF
C
      IA = SIGN (1., A(1))
      NA = MIN (INT (ABS (A(1))), NW)
C
      IF (NA .EQ. 0) THEN
        B(1) = 0.
        B(2) = 0.
        GOTO 140
      ENDIF
      IF (IA .LT. 0) THEN
        IF (KER(59) .NE. 0) THEN
          WRITE (LDB, 2)
 2        FORMAT ('*** MPNRT: Argument is negative.')
          IER = 59
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
      IF (N .LE. 0 .OR. N .GT. N30) THEN
        IF (KER(60) .NE. 0) THEN
          WRITE (LDB, 3) N
 3        FORMAT ('*** MPNRT: Improper value of N',I10)
          IER = 60
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
C   If N = 1, 2 or 3, call MPEQ, MPSQRT or MPCBRT.  These are faster.
C
      IF (N .EQ. 1) THEN
        CALL MPEQ (A, B)
        GOTO 140
      ELSEIF (N .EQ. 2) THEN
        CALL MPSQRT (A, B)
        GOTO 140
      ELSEIF (N .EQ. 3) THEN
        CALL MPCBRT (A, B)
        GOTO 140
      ENDIF
C
      N5 = NW + 5
      NS = 4 * N5
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N5
      K2 = K1 + N5
      K3 = K2 + N5
      NWS = NW
      F1(1) = 1.
      F1(2) = 0.
      F1(3) = 1.
C
C   Determine the least integer MQ such that 2 ^ MQ .GE. NW.
C
      T1 = NW
      MQ = CL2 * LOG (T1) + 1.D0 - RXX
C
C   Check how close A is to 1.
C
      CALL MPSUB (A, F1, S(K0))
      IF (S(K0) .EQ. 0.) THEN
        CALL MPEQ (F1, B)
        ICS = ISS
        GOTO 140
      ENDIF
      CALL MPMDC (S(K0), T1, N1)
      N2 = CL2 * LOG (ABS (T1))
      T1 = T1 * 0.5D0 ** N2
      N1 = N1 + N2
      IF (N1 .LE. -30) THEN
        T2 = N
        N2 = CL2 * LOG (T2) + 1.D0 + RXX
        N3 = - NBT * NW / N1
        IF (N3 .LT. 1.25 * N2) THEN
C
C   A is so close to 1 that it is cheaper to use the binomial series.
C
          NW = NW + 1
          CALL MPDIVD (S(K0), T2, 0, S(K1))
          CALL MPADD (F1, S(K1), S(K2))
          K = 0
C
 100      K = K + 1
          T1 = 1 - K * N
          T2 = (K + 1) * N
          CALL MPMULD (S(K1), T1, 0, S(K3))
          CALL MPDIVD (S(K3), T2, 0, S(K1))
          CALL MPMUL (S(K0), S(K1), S(K3))
          CALL MPEQ (S(K3), S(K1))
          CALL MPADD (S(K1), S(K2), S(K3))
          CALL MPEQ (S(K3), S(K2))
          IF (S(K1) .NE. 0. .AND. S(K1+1) .GE. - NW) GOTO 100
C
          CALL MPEQ (S(K2), B)
          CALL MPDIV (F1, S(K2), S(K0))
          GOTO 130
        ENDIF
      ENDIF
C
C   Compute the initial approximation of A ^ (-1/N).
C
      TN = N
      CALL MPMDC (A, T1, N1)
      N2 = - N1 / TN
      T2 = EXP (-1.D0 / TN * (LOG (T1) + (N1 + TN * N2) * ALT))
      T2 = (T1 * 2.D0 ** (N1 + TN * N2)) ** (- 1.D0 / TN)
      CALL MPDMC (T2, N2, B)
      CALL MPDMC (TN, 0, F2)
      NW = 3
      IQ = 0
C
C   Perform the Newton-Raphson iteration described above with a dynamically
C   changing precision level NW (one greater than powers of two).
C
      DO 120 K = 2, MQ
        NW = MIN (2 * NW - 2, NWS) + 1
 110    CONTINUE
        CALL MPNPWR (B, N, S(K0))
        CALL MPMUL (A, S(K0), S(K1))
        CALL MPSUB (F1, S(K1), S(K0))
        CALL MPMUL (B, S(K0), S(K1))
        CALL MPDIVD (S(K1), TN, 0, S(K0))
        CALL MPADD (B, S(K0), S(K1))
        CALL MPEQ (S(K1), B)
        IF (K .EQ. MQ - NIT .AND. IQ .EQ. 0) THEN
          IQ = 1
          GOTO 110
        ENDIF
 120  CONTINUE
C
C   Take the reciprocal to give final result.
C
      CALL MPDIV (F1, B, S(K1))
      CALL MPEQ (S(K1), B)
C
C   Restore original precision level.
C
 130  NW = NWS
      ICS = ISS
      CALL MPROUN (B)
C
 140  IF (IDB .GE. 7) THEN
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 4) (B(I), I = 1, NO)
 4      FORMAT ('MPNRT O'/(6F12.0))
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPNRTX (A, N, B)
C
C   This computes the N-th root of the MP number A and returns the MP result
C   in B.  N must be at least one and must not exceed 2 ^ 30.  Before calling
C   MPNRTX, the array in MPCOM5 must be initialized by calling MPINIX.  For
C   modest levels of precision, use MPNRT.  NW should be a power of two.  The
C   last three words of the result are not reliable.  Debug output starts with
C   IDB = 6.
C
C   Max SP space for B: NW + 4 cells.  Max SP scratch space: 9 * NW + 36
C   cells.  Max DP scratch space: 12 * NW + 6 cells.
C
C   This routine uses basically the same Newton iteration algorithm as MPNRT.
C   In fact, this routine calls MPNRT to obtain an initial approximation.
C   See the comment about the parameter NIT in MPDIVX.
C
      DOUBLE PRECISION CL2, T1, T2, TN
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (CL2 = 1.4426950408889633D0, NDB = 22, NIT = 1,         
     $  N30 = 2 ** 30)
      DIMENSION A(NW+2), B(NW+4), F1(8), F2(8)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        B(1) = 0.
        B(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 6) THEN
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 1) N, (A(I), I = 1, NO)
 1      FORMAT ('MPNRTX I',I5/(6F12.0))
      ENDIF
C
      NCR = 2 ** MCR
      IA = SIGN (1., A(1))
      NA = MIN (INT (ABS (A(1))), NW)
C
      IF (NA .EQ. 0) THEN
        B(1) = 0.
        B(2) = 0.
        GOTO 140
      ENDIF
      IF (IA .LT. 0) THEN
        IF (KER(61) .NE. 0) THEN
          WRITE (LDB, 2)
 2        FORMAT ('*** MPNRTX: Argument is negative.')
          IER = 61
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
      IF (N .LE. 0 .OR. N .GT. N30) THEN
        IF (KER(62) .NE. 0) THEN
          WRITE (LDB, 3) N
 3        FORMAT ('*** MPNRTX: Improper value of N',I10)
          IER = 62
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
C   Check if precision level is too low to justify the advanced routine.
C
      IF (NW .LE. NCR) THEN
        CALL MPNRT (A, N, B)
        GOTO 140
      ENDIF
C
C   If N = 1, 2 or 3, call MPEQ, MPSQRX or MPCBRX.  These are faster.
C
      IF (N .EQ. 1) THEN
        CALL MPEQ (A, B)
        GOTO 140
      ELSEIF (N .EQ. 2) THEN
        CALL MPSQRX (A, B)
        GOTO 140
      ELSEIF (N .EQ. 3) THEN
        CALL MPCBRX (A, B)
        GOTO 140
      ENDIF
C
      N4 = NW + 4
      NS = 4 * N4
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N4
      K2 = K1 + N4
      K3 = K2 + N4
      NWS = NW
      F1(1) = 1.
      F1(2) = 0.
      F1(3) = 1.
C
C   Determine the least integer MQ such that 2 ^ MQ .GE. NW.
C
      T1 = NW
      MQ = CL2 * LOG (T1) + 1.D0 - RXX
C
C   Check how close A is to 1.
C
      CALL MPSUB (A, F1, S(K0))
      IF (S(K0) .EQ. 0.) THEN
        CALL MPEQ (F1, B)
        GOTO 130
      ENDIF
      CALL MPMDC (S(K0), T1, N1)
      N2 = CL2 * LOG (ABS (T1))
      T1 = T1 * 0.5D0 ** N2
      N1 = N1 + N2
      IF (N1 .LE. -30) THEN
        T2 = N
        N2 = CL2 * LOG (T2) + 1.D0 + RXX
        N3 = - NBT * NW / N1
        IF (N3 .LT. 1.25 * N2) THEN
C
C   A is so close to 1 that it is cheaper to use the binomial series.
C
          CALL MPDIVD (S(K0), T2, 0, S(K1))
          CALL MPADD (F1, S(K1), S(K2))
          K = 0
C
 100      K = K + 1
          T1 = 1 - K * N
          T2 = (K + 1) * N
          CALL MPMULD (S(K1), T1, 0, S(K3))
          CALL MPDIVD (S(K3), T2, 0, S(K1))
          CALL MPMULX (S(K0), S(K1), S(K3))
          CALL MPEQ (S(K3), S(K1))
          CALL MPADD (S(K1), S(K2), S(K3))
          CALL MPEQ (S(K3), S(K2))
          IF (S(K1) .NE. 0. .AND. S(K1+1) .GE. - NW) GOTO 100
C
          CALL MPEQ (S(K2), B)
          GOTO 130
        ENDIF
      ENDIF
C
C   Compute the initial approximation of A ^ (-1/N).
C
      NW = NCR
      CALL MPNRT (A, N, S(K0))
      CALL MPDIV (F1, S(K0), B)
      TN = N
      CALL MPDMC (TN, 0, F2)
      IQ = 0
C
C   Perform the Newton-Raphson iteration described above with a dynamically
C   changing precision level NW (powers of two).
C
      DO 120 K = MCR + 1, MQ
        AN = NW
        NW = MIN (2 * NW, NWS)
 110    CONTINUE
        CALL MPNPWX (B, N, S(K0))
        CALL MPMULX (A, S(K0), S(K1))
        CALL MPSUB (F1, S(K1), S(K0))
        S(K1) = MIN (S(K1), AN)
        CALL MPMULX (B, S(K0), S(K1))
        CALL MPDIVD (S(K1), TN, 0, S(K0))
        CALL MPADD (B, S(K0), S(K1))
        CALL MPEQ (S(K1), B)
        IF (K .EQ. MQ - NIT .AND. IQ .EQ. 0) THEN
          IQ = 1
          GOTO 110
        ENDIF
 120  CONTINUE
C
C   Take the reciprocal to give final result.
C
      CALL MPDIVX (F1, B, S(K0))
      CALL MPEQ (S(K0), B)
C
 130  ICS = ISS
C
 140  IF (IDB .GE. 6) THEN
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 4) (B(I), I = 1, NO)
 4      FORMAT ('MPNRTX O'/(6F12.0))
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPOUT (IU, A, LA, CS)
C
C   This routine writes the exponent plus LA mantissa digits of the MP number
C   A to logical unit IU.  CS is a scratch array of type CHARACTER*1.  CS must
C   be dimensioned at least LA + 25.  The digits of A may span more than one
C   line.  A comma is placed at the end of the last line to denote the end of
C   the MP number.  Here is an example of the output:
C
C   10 ^        -4 x  3.14159265358979323846264338327950288419716939937510,
C
C   Max SP scratch space: 4 * NW + 22 cells.
C
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      CHARACTER*1 CS(LA+25)
      DIMENSION A(NW+2)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
C
      IF (IER .NE. 0) RETURN
C
      NWS = NW
      LL = LA / LOG10 (BDX) + 2.D0
      NW = MIN (NW, LL)
      CALL MPOUTC (A, CS, L)
      NW = NWS
      L = MIN (L, LA + 20) + 1
      CS(L) = ','
      WRITE (IU, '(78A1)') (CS(I), I = 1, L)
C
      RETURN
      END
C
      SUBROUTINE MPOUTC (A, B, N)
C
C   Converts the MP number A into character form in the CHARACTER*1 array B.
C   N (an output parameter) is the length of the output.  In other words, B is
C   contained in B(1), ..., B(N).  The format is analogous to the Fortran
C   exponential format (E format), except that the exponent is placed first.
C   Debug output starts with IDB = 7.
C
C   Max CHARACTER*1 space for B: 7.225 * NW + 30 cells.  Max SP scratch space:
C   4 * NW + 22 cells.
C
C   This routine is called by MPOUT, but it may be directly called by the user
C   if desired for custom output.  Example:
C
C      CHARACTER*1 CX(800)
C      CALL MPOUTC (A, CX, ND)
C      WRITE (1, '(20A1/(72A1))') (CX(I), I = 1, ND)
C
      DOUBLE PRECISION AA, AL2, T1
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      CHARACTER*1 B
      CHARACTER*16 CA
      PARAMETER (AL2 = 0.301029995663981195D0, CON = 0.8304820235D0,    
     $  NDB = 22)
      DIMENSION A(NW+2), B(*), F(8)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        B(1) = ' '
        N = 0
        RETURN
      ENDIF
      IF (IDB .GE. 7) THEN
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 1) (A(I), I = 1, NO)
 1      FORMAT ('MPOUTC I'/(6F12.0))
      ENDIF
C
      IA = SIGN (1., A(1))
      NA = MIN (INT (ABS (A(1))), NW)
      N5 = NW + 5
      NS = 2 * N5
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N5
      NWS = NW
      NW = NW + 1
      F(1) = 1.
      F(2) = 0.
      F(3) = 10.
C
C   Determine exact power of ten for exponent.
C
      IF (NA .NE. 0) THEN
        AA = A(3)
        IF (NA .GE. 2) AA = AA + RDX * A(4)
        IF (NA .GE. 3) AA = AA + RX2 * A(5)
        IF (NA .GE. 4) AA = AA + RDX * RX2 * A(6)
        T1 = AL2 * NBT * A(2) + LOG10 (AA)
        IF (T1 .GE. 0.D0) THEN
          NX = T1
        ELSE
          NX = T1 - 1.D0
        ENDIF
        CALL MPNPWR (F, NX, S(K0))
        CALL MPDIV (A, S(K0), S(K1))
C
C   If we didn't quite get it exactly right, multiply or divide by 10 to fix.
C
 100    IF (S(K1+1) .LT. 0.) THEN
          NX = NX - 1
          CALL MPMULD (S(K1), 10.D0, 0, S(K0))
          CALL MPEQ (S(K0), S(K1))
          GOTO 100
        ELSEIF (S(K1+2) .GE. 10.) THEN
          NX = NX + 1
          CALL MPDIVD (S(K1), 10.D0, 0, S(K0))
          CALL MPEQ (S(K0), S(K1))
          GOTO 100
        ENDIF
        S(K1) = ABS (S(K1))
      ELSE
        NX = 0
      ENDIF
C
C   Place exponent first instead of at the very end as in Fortran.
C
      B(1) = '1'
      B(2) = '0'
      B(3) = ' '
      B(4) = '^'
      WRITE (CA, '(I10)') NX
C
      DO 110 I = 1, 10
        B(I+4) = CA(I:I)
 110  CONTINUE
C
      B(15) = ' '
      B(16) = 'x'
      B(17) = ' '
C
C   Insert sign and first digit.
C
      IF (IA .EQ. -1) THEN
        B(18) = '-'
      ELSE
        B(18) = ' '
      ENDIF
      IF (NA .NE. 0) THEN
        NN = S(K1+2)
      ELSE
        NN = 0
      ENDIF
      WRITE (CA, '(I1)') NN
      B(19) = CA(1:1)
      B(20) = '.'
      IX = 20
      IF (NA .EQ. 0) GOTO 190
      F(3) = NN
      CALL MPSUB (S(K1), F, S(K0))
      IF (S(K0) .EQ. 0) GOTO 190
      CALL MPMULD (S(K0), 1.D6, 0, S(K1))
      NL = MAX (NW * LOG10 (BDX) / 6.D0 - 1.D0, 1.D0)
C
C   Insert the digits of the remaining words.
C
      DO 130 J = 1, NL
        IF (S(K1+1) .EQ. 0.) THEN
          NN = S(K1+2)
          F(1) = 1.
          F(3) = NN
        ELSE
          F(1) = 0.
          NN = 0
        ENDIF
        WRITE (CA, '(I6.6)') NN
C
        DO 120 I = 1, 6
          B(I+IX) = CA(I:I)
 120    CONTINUE
C
        IX = IX + 6
        CALL MPSUB (S(K1), F, S(K0))
        CALL MPMULD (S(K0), 1.D6, 0, S(K1))
        IF (S(K1) .EQ. 0.) GOTO 140
 130  CONTINUE
C
C   Check if trailing zeroes should be trimmed.
C
      J = NL + 1
C
 140  L = IX
      IF (B(L) .EQ. '0' .OR. (J .GT. NL .AND. B(L-1) .EQ. '0' .AND.     
     $  B(L-2) .EQ. '0' .AND. B(L-3) .EQ. '0')) THEN
        B(L) = ' '
C
        DO 150 I = L - 1, 21, -1
          IF (B(I) .NE. '0') THEN
            IX = I
            GOTO 190
          ENDIF
          B(I) = ' '
 150    CONTINUE
C
        IX = 20
C
C   Check if trailing nines should be rounded up.
C
      ELSEIF (J .GT. NL .AND. B(L-1) .EQ. '9' .AND. B(L-2) .EQ. '9'     
     $    .AND. B(L-3) .EQ. '9') THEN
        B(L) = ' '
C
        DO 160 I = L - 1, 21, -1
          IF (B(I) .NE. '9') GOTO 180
          B(I) = ' '
 160    CONTINUE
C
C   We have rounded away all digits to the right of the decimal point, and the
C   digit to the left of the digit is a 9.  Set the digit to 1 and increase
C   the exponent by one.
C
        IX = 20
        IF (B(19) .EQ. '9') THEN
          B(19) = '1'
          WRITE (CA, '(I10)') NX + 1
C
          DO 170 I = 1, 10
            B(I+4) = CA(I:I)
 170      CONTINUE
C
        ELSE
          CA = B(19)
          READ (CA, '(I1)') NN
          WRITE (CA, '(I1)') NN + 1
          B(19) = CA(1:1)
        ENDIF
        GOTO 190
C
 180    CA = B(I)
        READ (CA, '(I1)') NN
        WRITE (CA, '(I1)') NN + 1
        B(I) = CA(1:1)
        IX = I
      ENDIF
C
 190  N = IX
      NW = NWS
      ICS = ISS
      IF (IDB .GE. 7) THEN
        NO = MIN (N, 6 * NDB + 20)
        WRITE (LDB, 2) (B(I), I = 1, NO)
 2      FORMAT ('MPOUTC O'/(78A1))
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPPI (PI)
C
C   This computes Pi to available precision (NW mantissa words).  For extra
C   high levels of precision, use MPPIX.  Debug output starts with IDB = 7.
C
C   Max SP space for PI: NW + 4 cells.  Max SP scratch space: 8 * NW + 43
C   cells.  Max DP scratch space: NW + 6 cells.
C
C   The algorithm that is used for computing Pi, which is due to Salamin
C   and Brent, is as follows:
C
C   Set  A_0 = 1,  B_0 = 1/Sqrt(2)  and  D_0 = Sqrt(2) - 1/2.
C
C   Then from k = 1 iterate the following operations:
C
C      A_k = 0.5 * (A_{k-1} + B_{k-1})
C      B_k = Sqrt (A_{k-1} * B_{k-1})
C      D_k = D_{k-1} - 2^k * (A_k - B_k) ^ 2
C
C   Then  P_k = (A_k + B_k) ^ 2 / D_k  converges quadratically to Pi.
C   In other words, each iteration approximately doubles the number of correct
C   digits, providing all iterations are done with the maximum precision.
C
      DOUBLE PRECISION CL2, T1
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (CL2 = 1.4426950408889633D0)
      DIMENSION F(8), PI(NW+4)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        PI(1) = 0.
        PI(2) = 0.
        RETURN
      ENDIF
C
C   Perform calculations to one extra word accuracy.
C
      N5 = NW + 5
      NS = 5 * N5
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N5
      K2 = K1 + N5
      K3 = K2 + N5
      K4 = K3 + N5
      NWS = NW
      NW = NW + 1
C
C   Determine the number of iterations required for the given precision level.
C   This formula is good only for this Pi algorithm.
C
      T1 = NWS * LOG10 (BDX)
      MQ = CL2 * (LOG (T1) - 1.D0) + 1.D0
C
C  Initialize as above.
C
      S(K0) = 1.
      S(K0+1) = 0.
      S(K0+2) = 1.
      F(1) = 1.
      F(2) = 0.
      F(3) = 2.
      CALL MPSQRT (F, S(K2))
      CALL MPMULD (S(K2), 0.5D0, 0, S(K1))
      F(2) = -1.
      F(3) = 0.5D0 * BDX
      CALL MPSUB (S(K2), F, S(K4))
C
C   Perform iterations as described above.
C
      DO 100 K = 1, MQ
        CALL MPADD (S(K0), S(K1), S(K2))
        CALL MPMUL (S(K0), S(K1), S(K3))
        CALL MPSQRT (S(K3), S(K1))
        CALL MPMULD (S(K2), 0.5D0, 0, S(K0))
        CALL MPSUB (S(K0), S(K1), S(K2))
        CALL MPMUL (S(K2), S(K2), S(K3))
        T1 = 2.D0 ** K
        CALL MPMULD (S(K3), T1, 0, S(K2))
        CALL MPSUB (S(K4), S(K2), S(K3))
        CALL MPEQ (S(K3), S(K4))
 100  CONTINUE
C
C   Complete computation.
C
      CALL MPADD (S(K0), S(K1), S(K2))
      CALL MPMUL (S(K2), S(K2), S(K3))
      CALL MPDIV (S(K3), S(K4), S(K2))
      CALL MPEQ (S(K2), PI)
C
C   Restore original precision level.
C
      NW = NWS
      ICS = ISS
      CALL MPROUN (PI)
C
      IF (IDB .GE. 7) CALL MPDEB ('MPPI O', PI)
      RETURN
      END
C
      SUBROUTINE MPPIX (PI)
C
C   This computes Pi to available precision (NW mantissa words).  Before
C   calling MPPIX, the array in MPCOM5 must be initialized by calling MPINIX.
C   For modest levels of precision, use MPPI.  NW should be a power of two.
C   The last three words of the result are not reliable.  Debug output starts
C   with IDB = 7.
C
C   Max SP space for PI: NW + 4 cells.  Max SP scratch space: 9.5 * NW + 47
C   cells.  Max DP scratch space: 12 * NW + 6 cells.
C
C   This routine uses basically the same algorithm as MPPI.
C
      DOUBLE PRECISION CL2, T1
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (CL2 = 1.4426950408889633D0)
      DIMENSION F(8), PI(NW+4)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        PI(1) = 0.
        PI(2) = 0.
        RETURN
      ENDIF
      NCR = 2 ** MCR
C
C   Check if precision level is too low to justify the advanced routine.
C
      IF (NW .LE. NCR) THEN
        CALL MPPI (PI)
        GOTO 110
      ENDIF
      N4 = NW + 4
      NS = 5 * N4
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N4
      K2 = K1 + N4
      K3 = K2 + N4
      K4 = K3 + N4
C
C   Determine the number of iterations required for the given precision level.
C   This formula is good only for this Pi algorithm.
C
      T1 = NW * LOG10 (BDX)
      MQ = CL2 * (LOG (T1) - 1.D0) + 1.D0
C
C  Initialize as above.
C
      S(K0) = 1.
      S(K0+1) = 0.
      S(K0+2) = 1.
      F(1) = 1.
      F(2) = 0.
      F(3) = 2.
      CALL MPSQRX (F, S(K2))
      CALL MPMULD (S(K2), 0.5D0, 0, S(K1))
      F(2) = -1.
      F(3) = 0.5D0 * BDX
      CALL MPSUB (S(K2), F, S(K4))
C
C   Perform iterations as described above.
C
      DO 100 K = 1, MQ
        CALL MPADD (S(K0), S(K1), S(K2))
        CALL MPMULX (S(K0), S(K1), S(K3))
        CALL MPSQRX (S(K3), S(K1))
        CALL MPMULD (S(K2), 0.5D0, 0, S(K0))
        CALL MPSUB (S(K0), S(K1), S(K2))
        CALL MPSQX (S(K2), S(K3))
        T1 = 2.D0 ** K
        CALL MPMULD (S(K3), T1, 0, S(K2))
        CALL MPSUB (S(K4), S(K2), S(K3))
        CALL MPEQ (S(K3), S(K4))
 100  CONTINUE
C
C   Complete computation.
C
      CALL MPADD (S(K0), S(K1), S(K2))
      CALL MPSQX (S(K2), S(K3))
      CALL MPDIVX (S(K3), S(K4), S(K2))
      CALL MPEQ (S(K2), PI)
      ICS = ISS
C
 110  IF (IDB .GE. 7) CALL MPDEB ('MPPIX O', PI)
      RETURN
      END
C
      SUBROUTINE MPPOL (N, L, A, X1, NX, X)
C
C   This finds a real root of the N-th degree polynomial whose MP coefficients
C   are in A by Newton-Raphson iterations, beginning at the DPE value (X1, NX)
C   and returns the MP root in X.  The N + 1 coefficients a_0, a_1, ..., a_N
C   are assumed to start in locations A(1), A(L+1), A(2*L+1), etc.  For extra
C   high levels of precision, use MPPOLX.  Debug output starts with IDB = 6.
C
C   Max SP space for X: NW + 4 cells.  Max SP scratch space: 5 * NW + 25
C   cells.  Max DP scratch space: NW + 5 cells.
C
C   One requirement for this routine to work is that the desired root is not
C   a repeated root.  If one wishes to apply this routine to find a repeated
C   root, it is first necessary to reduce the polynomial to one that has only
C   simple roots.  This can be done by performing the Euclidean algorithm in
C   the ring of polynomials to determine the greatest common divisor Q(t) of
C   P(t) and P'(t).  Here P(t) is the polynomial a_0 + a_1 t + a_2 t^2 +
C   ... + a_n t^n, and P'(t) is the derivative of P(t).  Then R(t) = P(t)/Q(t)
C   is a polynomial that has only simple roots.
C
C   This routine employs the standard form of the Newton-Raphson iteration:
C
C   X_{k+1} = X_k - P(X_k) / P'(X_k)
C
C   These iterations are performed with a maximum precision level NW that is
C   dynamically changed, approximately doubling with each iteration.
C
      CHARACTER*8 CX
      DOUBLE PRECISION T1, X1
      DIMENSION A(L,N+1), X(NW+4)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        X(1) = 0.
        X(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 6) THEN
        WRITE (LDB, 1) N
 1      FORMAT ('MPPOL I',I4)
C
        DO 100 K = 0, N
          WRITE (CX, '(I4)') K
          CALL MPDEB (CX, A(1,K+1))
 100    CONTINUE
C
        WRITE (LDB, 2) X1, NX
 2      FORMAT ('MPPOL I',F16.12,' x 10 ^',I6)
      ENDIF
C
C  Check if the polynomial is proper.
C
      IF (A(1,1) .EQ. 0. .OR. A(1,N+1) .EQ. 0.) THEN
        IF (KER(63) .NE. 0) THEN
          WRITE (LDB, 3)
 3        FORMAT ('*** MPPOL: Either the first or last input ',         
     $      'coefficient is zero.')
          IER = 63
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
      N5 = NW + 5
      NS = 5 * N5
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N5
      K2 = K1 + N5
      K3 = K2 + N5
      K4 = K3 + N5
      NWS = NW
      NW = NW + 1
C
C   Set the initial value.
C
      CALL MPDMC (X1, NX, S(K0))
      NW = 5
      TL = -4.
      L1 = 0
      LS = -10
C
C   Perform MP Newton-Raphson iterations to solve P(x) = 0.
C
 110  L1 = L1 + 1
      IF (L1 .EQ. 50) THEN
        IF (KER(64) .NE. 0) THEN
          WRITE (LDB, 4)
 4        FORMAT ('*** MPPOL: Iteration limit exceeded.')
          IER = 64
          IF (KER(IER) .EQ. 2) CALL MPABRT
          NW = NWS
          ICS = ISS
          RETURN
        ENDIF
      ENDIF
C
C   Compute P(x).
C
      CALL MPEQ (A(1,N+1), S(K1))
C
      DO 120 K = N - 1, 0, -1
        CALL MPMUL (S(K0), S(K1), S(K2))
        CALL MPADD (S(K2), A(1,K+1), S(K1))
 120  CONTINUE
C
C   Compute P'(x).
C
      T1 = N
      CALL MPMULD (A(1,N+1), T1, 0, S(K2))
C
      DO 130 K = N - 1, 1, -1
        CALL MPMUL (S(K0), S(K2), S(K3))
        T1 = K
        CALL MPMULD (A(1,K+1), T1, 0, S(K4))
        CALL MPADD (S(K3), S(K4), S(K2))
 130  CONTINUE
C
C   Compute P(x) / P'(x) and update x.
C
      CALL MPDIV (S(K1), S(K2), S(K3))
      CALL MPSUB (S(K0), S(K3), S(K4))
C
      IF (IDB .GE. 7) THEN
        WRITE (LDB, 5) L1
 5      FORMAT ('Iteration',I4)
        CALL MPDEB ('X', S(K0))
        CALL MPDEB ('P(X)', S(K1))
        CALL MPDEB ('P''(X)', S(K2))
        CALL MPDEB ('CORR', S(K3))
      ENDIF
      CALL MPEQ (S(K4), S(K0))
C
C   If this was the second iteration at full precision, there is no need to
C   continue (the adjusted value of x is correct); otherwise repeat.
C
      IF (L1 .EQ. LS + 1) GOTO 140
      IF (S(K3) .NE. 0. .AND. S(K3+1) .GT. TL) GOTO 110
C
C   Newton iterations have converged to current precision.  Increase precision
C   and continue.
C
      IF (NW .EQ. NWS + 1) GOTO 140
      NW = MIN (2 * NW - 2, NWS) + 1
      IF (NW .EQ. NWS + 1) LS = L1
      TL = 1 - NW
      IF (IDB .GE. 7) THEN
        WRITE (LDB, 6) NW
 6      FORMAT (6X,'New NW =', I8)
      ENDIF
      GOTO 110
C
 140  CALL MPEQ (S(K0), X)
C
C   Restore original precision level.
C
      NW = NWS
      ICS = ISS
      CALL MPROUN (X)
C
      IF (IDB .GE. 6) THEN
        WRITE (LDB, 7) L1
 7      FORMAT ('Iteration count:',I5)
        CALL MPDEB ('MPPOL O', X)
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPPOLX (N, L, A, X1, NX, X)
C
C   This finds a real root of the N-th degree polynomial whose MP coefficients
C   are in A by Newton-Raphson iterations, beginning at the DP value (X1, NX)
C   and returns the MP root in X.  The N + 1 coefficients a_0, a_1, ..., a_N
C   are assumed to start in locations A(1), A(L+1), A(2*L+1), etc.  Before
C   calling MPPOLX, the array in MPCOM5 must be initialized by calling MPINIX.
C   For modest levels of precision, use MPPOL.  NW should be a power of two.
C   The last three words of the result are not reliable.  Debug output starts
C   with IDB = 5.
C
C   Max SP space for X: NW + 4 cells.  Max SP scratch space: 8 * NW + 32
C   cells.  Max DP scratch space: 12 * NW + 6 cells.
C
C   For a discussion of the algorithm and usage, see MPPOL.  This routine uses
C   basically the same Newton iteration algorithm as MPPOL.  In fact, this
C   routine calls MPPOL to obtain an initial approximation.
C
      CHARACTER*8 CX
      DOUBLE PRECISION T1, X1
      DIMENSION A(L,N+1), X(NW+4)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        X(1) = 0.
        X(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 5) THEN
        WRITE (LDB, 1) N
 1      FORMAT ('MPPOLX I',I4)
C
        DO 100 K = 0, N
          WRITE (CX, '(I4)') K
          CALL MPDEB (CX, A(1,K+1))
 100    CONTINUE
C
        WRITE (LDB, 2) X1, NX
 2      FORMAT ('MPPOLX I',F16.12,' x 10 ^',I6)
      ENDIF
C
C   Check if precision level is too low to justify the advanced routine.
C
      NCR = 2 ** MCR
      IF (NW .LE. NCR) THEN
        CALL MPPOL (N, L, A, X1, NX, X)
        L1 = 0
        GOTO 150
      ENDIF
C
C  Check if the polynomial is proper.
C
      IF (A(1,1) .EQ. 0. .OR. A(1,N+1) .EQ. 0.) THEN
        IF (KER(65) .NE. 0) THEN
          WRITE (LDB, 3)
 3        FORMAT ('*** MPPOLX: Either the first or last input ',        
     $      'coefficient is zero.')
          IER = 65
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
      N4 = NW + 4
      NS = 5 * N4
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N4
      K2 = K1 + N4
      K3 = K2 + N4
      K4 = K3 + N4
      NWS = NW
C
C   Compute the initial approximation.
C
      NW = NCR
      CALL MPPOL (N, L, A, X1, NX, X)
      CALL MPEQ (X, S(K0))
      TL = 2 - NW
      L1 = 0
      LS = -10
C
C   Perform MP Newton-Raphson iterations to solve P(x) = 0.
C
 110  L1 = L1 + 1
      IF (L1 .EQ. 50) THEN
        IF (KER(66) .NE. 0) THEN
          WRITE (LDB, 4)
 4        FORMAT ('*** MPPOLX: Iteration limit exceeded.')
          IER = 66
          IF (KER(IER) .EQ. 2) CALL MPABRT
          NW = NWS
          ICS = ISS
          RETURN
        ENDIF
      ENDIF
C
C   Compute P(x).
C
      CALL MPEQ (A(1,N+1), S(K1))
C
      DO 120 K = N - 1, 0, -1
        CALL MPMULX (S(K0), S(K1), S(K2))
        CALL MPADD (S(K2), A(1,K+1), S(K1))
 120  CONTINUE
C
C   Compute P'(x).
C
      T1 = N
      CALL MPMULD (A(1,N+1), T1, 0, S(K2))
C
      DO 130 K = N - 1, 1, -1
        CALL MPMULX (S(K0), S(K2), S(K3))
        T1 = K
        CALL MPMULD (A(1,K+1), T1, 0, S(K4))
        CALL MPADD (S(K3), S(K4), S(K2))
 130  CONTINUE
C
C   Compute P(x) / P'(x) and update x.
C
      CALL MPDIVX (S(K1), S(K2), S(K3))
      CALL MPSUB (S(K0), S(K3), S(K4))
C
      IF (IDB .GE. 6) THEN
        WRITE (LDB, 5) L1
 5      FORMAT ('Iteration',I4)
        CALL MPDEB ('X', S(K0))
        CALL MPDEB ('P(X)', S(K1))
        CALL MPDEB ('P''(X)', S(K2))
        CALL MPDEB ('CORR', S(K3))
      ENDIF
      CALL MPEQ (S(K4), S(K0))
C
C   If this was the second iteration at full precision, there is no need to
C   continue (the adjusted value of x is correct); otherwise repeat.
C
      IF (L1 .EQ. LS + 1) GOTO 140
      IF (S(K3) .NE. 0. .AND. S(K3+1) .GT. TL) GOTO 110
C
C   Newton iterations have converged to current precision.  Increase precision
C   and continue.
C
      IF (NW .EQ. NWS) GOTO 140
      NW = MIN (2 * NW, NWS)
      IF (NW .EQ. NWS) LS = L1
      IF (NW .LE. 32) THEN
        TL = 2 - NW
      ELSEIF (NW .LE. 256) THEN
        TL = 3 - NW
      ELSE
        TL = 4 - NW
      ENDIF
      IF (IDB .GE. 6) THEN
        WRITE (LDB, 6) NW
 6      FORMAT (6X,'New NW =', I8)
      ENDIF
      GOTO 110
C
 140  CALL MPEQ (S(K0), X)
      ICS = ISS
C
 150  IF (IDB .GE. 5) THEN
        WRITE (LDB, 7) L1
 7      FORMAT ('Iteration count:',I5)
        CALL MPDEB ('MPPOLX O', X)
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPRAND (A)
C>
C   This returns a pseudo-random MP number A between 0 and 1.  This routine
C   calls the pseudo-random number generator routine MPRANQ in the file below.
C   Better routines than MPRANQ are available for this purpose on some
C   computer systems.  If so, it is suggested that the call to MPRANQ here be
C   replaced by a call to its equivalent on the host system.  Note, however,
C   that test no. 55 of the TESTMP test suite will fail if another generator
C   is used.  Debug output starts with IDB = 9.
C
C   Max SP space for A: NW + 4 cells.
C
      DOUBLE PRECISION MPRANQ, SD, S0
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (NDB = 22, S0 = 314159265.D0)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      DIMENSION A(NW+4)
      SAVE SD
      DATA SD/S0/
C
      IF (IER .NE. 0) THEN
        A(1) = 0.
        A(2) = 0.
        RETURN
      ENDIF
      A(1) = NW
      A(2) = -1.
C
      DO 100 I = 3, NW + 4
        A(I) = AINT (BDX * MPRANQ (SD))
 100  CONTINUE
C
      CALL MPROUN (A)
C
      IF (IDB .GE. 9) THEN
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 1) (A(I), I = 1, NO)
 1      FORMAT ('MPRAND O'/(6F12.0))
      ENDIF
      RETURN
      END
C
      FUNCTION MPRANQ (SD)
C
C   This routine returns a pseudorandom DP floating number uniformly
C   distributed between 0 and 1, computed from the seed SD, which is updated
C   after each reference.  The initial value of SD should be an odd whole
C   number in the range (1, 2 ^ 30).  2 ^ 28 pseudorandom numbers with 30 bits
C   each are returned before repeating.  The same sequence is generated on any
C   computer system.
C
      DOUBLE PRECISION F7, R30, SD, T1, T2, T30, MPRANQ
      PARAMETER (F7 = 78125.D0)
      SAVE R30, T30
      DATA R30/0.D0/
C
C   If this is the first time MPRANQ has been called, compute R30 = 2^(-30)
C   and T30 = 2^30.  This must be done in a loop rather than by merely using
C   ** in order to insure the results are exact on all systems.
C
      IF (R30 .EQ. 0.D0) THEN
        R30 = 1.D0
        T30 = 1.D0
C
        DO 100 I = 1, 30
          R30 = 0.5D0 * R30
          T30 = 2.D0 * T30
 100    CONTINUE
C
      ENDIF
C
C   Generate a pseudorandom number using a linear congruential scheme.
C
      T1 = F7 * SD
      T2 = AINT (R30 * T1)
      SD = T1 - T30 * T2
      MPRANQ = R30 * SD
C
      RETURN
      END
C
      SUBROUTINE MPRCFT (IS, M, X, Y)
C
C   This performs an N-point real-to-complex FFT, where N = 2^M.  X and Y
C   are double precision arrays.  X is both the input and the output data
C   array, and Y is a scratch array.  N real values are input and N/2 + 1
C   complex pairs are output, with real and imaginary parts separated by
C   N/2 + 1 locations.  A call to MPRCFT with IS = 1 (or -1) indicates a call
C   to perform a complex-to-real FFT with positive (or negative) exponentials.
C   M must be at least three.  The arrays X and Y must be dimensioned with
C   N + 2 cells.  Before calling MPRCFT, the U array in MPCOM5 must be
C   initialized by calling MPINIX.
C
C   In this application, MPRCFT is called by MPMULX.  This routine is not
C   intended to be called directly by the user.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION X(*), Y(*)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM5/ U(1024)
C
C   Set initial parameters.
C
      K = U(1)
      MX = MOD (K, 64)
      NU = K / 64
      N = 2 ** M
      N2 = N / 2
      N21 = N2 + 1
      N4 = N / 4
      KU = N / 2
      KN = KU + NU
C
C   Check if input parameters are invalid.
C
      IF ((IS .NE. 1 .AND. IS .NE. -1) .OR. M .LT. 3 .OR. M .GT. MX)    
     $  THEN
        IF (KER(67) .NE. 0) THEN
          WRITE (LDB, 1)  IS, M, MX
 1        FORMAT ('*** MPRCFT: either U has not been initialized'/      
     $      'or else one of the input parameters is invalid',3I5)
          IER = 67
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
C   Copy X to Y such that Y(k) = X(2k-1) + i X(2k).
C
CDIR$ IVDEP
      DO 100 K = 1, N2
        Y(K) = X(2*K-1)
        Y(K+N2) = X(2*K)
 100  CONTINUE
C
C   Perform a normal N/2-point FFT on Y.
C
      CALL MPCFFT (IS, M - 1, Y, X)
C
C   Reconstruct the FFT of X.
C
      X(1) = 2.D0 * (Y(1) + Y(N21))
      X(N21+1) = 0.D0
      X(N4+1) = 2.D0 * Y(N4+1)
      X(N4+1+N21) = 2.D0 * IS * Y(N4+N2+1)
      X(N21) = 2.D0 * (Y(1) - Y(N21))
      X(N+2) = 0.D0
C
CDIR$ IVDEP
      DO 110 K = 2, N4
        Y11 = Y(K)
        Y12 = Y(K+N2)
        Y21 = Y(N2+2-K)
        Y22 = Y(N+2-K)
        A1 = Y11 + Y21
        A2 = Y11 - Y21
        B1 = Y12 + Y22
        B2 = Y12 - Y22
        U1 = U(K+KU)
        U2 = IS * U(K+KN)
        T1 = U1 * B1 + U2 * A2
        T2 = - U1 * A2 + U2 * B1
        X(K) = A1 + T1
        X(K+N21) = B2 + T2
        X(N2+2-K) = A1 - T1
        X(N+3-K) = -B2 + T2
 110  CONTINUE
C
      RETURN
      END
C
      SUBROUTINE MPROUN (A)
C
C   This performs rounding and truncation of the MP number A.  It is called
C   by MPNORM, and also by other subroutines when the precision level is
C   reduced by one.  It is not intended to be directly called by the user.
C
C   Maximum SP space for A:  NW + 4 cells.
C
C   The parameter AMX is the absolute value of the largest exponent word
C   allowed for MP numbers.
C
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (AMX = 2.E6)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      DIMENSION A(NW+4)
C
      IF (IER .NE. 0) THEN
        A(1) = 0.
        A(2) = 0.
        RETURN
      ENDIF
C
C   Check for initial zeroes.
C
      A2 = A(2)
      A(2) = 0.
      IA = SIGN (1., A(1))
      NA = MIN (INT (ABS (A(1))), NW)
      N4 = NA + 4
      IF (A(3) .EQ. 0.) THEN
C
C   Find the first nonzero word and shift the entire number left.  The length
C   of the result is reduced by the length of the shift.
C
        DO 100 I = 4, N4
          IF (A(I) .NE. 0.)  GOTO 110
 100    CONTINUE
C
        A(1) = 0.
        A(2) = 0.
        GOTO 170
C
 110    K = I - 3
C
CDIR$ IVDEP
        DO 120 I = 3, N4 - K
          A(I) = A(I+K)
 120    CONTINUE
C
        A2 = A2 - K
        NA = NA - MAX (K - 2, 0)
      ENDIF
C
C   Perform rounding depending on IRD.
C
      IF (NA .EQ. NW .AND. IRD .GE. 1) THEN
        IF (IRD .EQ. 1 .AND. A(NA+3) .GE. 0.5D0 * BDX .OR. IRD .EQ. 2   
     $    .AND. A(NA+3) .GE. 1.) A(NA+2) = A(NA+2) + 1.
C
C   Release carries as far as necessary due to rounding.
C
        DO 130 I = NA + 2, 3, -1
          IF (A(I) .LT. BDX) GOTO 140
          A(I) = A(I) - BDX
          A(I-1) = A(I-1) + 1.
 130    CONTINUE
C
C   Release of carries due to rounding continued all the way to the start --
C   i.e. number was entirely 9's.
C
        A(3) = A(2)
        NA = 1
        A2 = A2 + 1.
      ENDIF
C
 140  IF (A(NA+2) .EQ. 0.) THEN
C
C   At least the last mantissa word is zero.  Find the last nonzero word
C   and adjust the length of the result accordingly.
C
        DO 150 I = NA + 2, 3, -1
          IF (A(I) .NE. 0.) GOTO 160
 150    CONTINUE
C
        A(1) = 0.
        A(2) = 0.
        GOTO 170
C
 160    NA = I - 2
      ENDIF
C
C   Check for overflow and underflow.
C
      IF (A2 .LT. - AMX) THEN
        IF (KER(68) .NE. 0) THEN
          WRITE (LDB, 1)
 1        FORMAT ('*** MPROUN: Exponent underflow.')
          IER = 68
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
      ELSEIF (A2 .GT. AMX) THEN
        IF (KER(69) .NE. 0) THEN
          WRITE (LDB, 2)
 2        FORMAT ('*** MPROUN: Exponent overflow.')
          IER = 69
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
      ENDIF
C
C   Check for zero.
C
      IF (A(3) .EQ. 0.) THEN
        A(1) = 0.
        A(2) = 0.
      ELSE
        A(1) = SIGN (NA, IA)
        A(2) = A2
        A(NA+3) = 0.
        A(NA+4) = 0.
      ENDIF
C
 170  RETURN
      END
C
      SUBROUTINE MPSETP (IA, IB)
C
C   This routine sets the parameter whose name is IA in common MPCOM1 to the
C   value IB.  By using this routine instead of merely including the MPCOM1
C   block in the code, a user may eliminate the possibility of confusion with
C   a variable name in his or her program.  IA is of type CHARACTER and IB
C   is the integer value.
C
      CHARACTER*(*) IA
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
C
      IF (IA .EQ. 'NW' .OR. IA .EQ. 'nw') THEN
        NW = IB
      ELSEIF (IA .EQ. 'IDB' .OR. IA .EQ. 'idb') THEN
        IDB = IB
      ELSEIF (IA .EQ. 'LDB' .OR. IA .EQ. 'ldb') THEN
        LDB = IB
      ELSEIF (IA .EQ. 'IER' .OR. IA .EQ. 'ier') THEN
        IER = IB
      ELSEIF (IA .EQ. 'MCR' .OR. IA .EQ. 'mcr') THEN
        MCR = IB
      ELSEIF (IA .EQ. 'IRD' .OR. IA .EQ. 'ird') THEN
        IRD = IB
      ELSEIF (IA .EQ. 'ICS' .OR. IA .EQ. 'ics') THEN
        ICS = IB
      ELSEIF (IA .EQ. 'IHS' .OR. IA .EQ. 'ihs') THEN
        IHS = IB
      ELSEIF (IA .EQ. 'IMS' .OR. IA .EQ. 'ims') THEN
        IMS = IB
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE MPSORT (N, LA, A, IP)
C
C   This routine sorts the entries of the N-long MP vector A into ascending
C   order using the quicksort algorithm.  The entries of A are assumed to
C   start at A(1), A(LA+1), A(2*LA+1), etc. The permutation vector that would
C   sort the vector is returned in IP.  Debug output starts with IDB = 7.
C
C   Max integer space for IP: N cells.  Max SP scratch space: 2 * NW + 8 cells.
C
      CHARACTER*8 CX
      DIMENSION A(LA,N), IP(N), IK(50), JK(50)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
C
        DO 100 I = 1, N
          IP(I) = I
 100    CONTINUE
C
        RETURN
      ENDIF
      IF (IDB .GE. 7) THEN
        WRITE (LDB, 1) N, LA
 1      FORMAT ('MPSORT I',2I6)
C
        DO 110 K = 1, N
          WRITE (CX, '(I4)') K
          CALL MPDEB (CX, A(1,K))
 110    CONTINUE
C
      ENDIF
C
      N4 = NW + 4
      NS = 2 * N4
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N4
C
      DO 120 I = 1, N
        IP(I) = I
 120  CONTINUE
C
      K = 1
      IK(1) = 1
      JK(1) = N
C
 130  I = IK(K)
      J = JK(K)
      IQ = I
      JQ = J
      IT = (I + J + 1) / 2
      L = IP(J)
      IP(J) = IP(IT)
      IP(IT) = L
      CALL MPEQ (A(1,IP(J)), S(K0))
      J = J - 1
C
 140  DO 150 L = I, J
        CALL MPCPR (S(K0), A(1,IP(L)), IC)
        IF (IC .LT. 0) GOTO 160
 150  CONTINUE
C
      I = J
      GOTO 190
C
 160  I = L
C
      DO 170 L = J, I, -1
        CALL MPCPR (S(K0), A(1,IP(L)), IC)
        IF (IC .GT. 0) GOTO 180
 170  CONTINUE
C
      J = I
      GOTO 190
C
 180  J = L
      IF (I .GE. J)  GOTO 190
      L = IP(I)
      IP(I) = IP(J)
      IP(J) = L
      GOTO 140
C
 190  CALL MPCPR (S(K0), A(1,IP(I)), IC)
      IF (IC .GE. 0) GOTO 200
      L = IP(JQ)
      IP(JQ) = IP(I)
      IP(I) = L
C
 200  K = K - 1
      JZ = 0
      IF (J .EQ. IQ)  GOTO 210
      K = K + 1
      JK(K) = J
      JZ = 1
C
 210  I = I + 1
      IF (I .EQ. JQ)  GOTO 220
      K = K + 1
      IK(K) = I
      JK(K) = JQ
      IF (JZ .EQ. 0)  GOTO 220
      IF (J - IQ .GE. JQ - I)  GOTO 220
      IK(K-1) = I
      JK(K-1) = JQ
      IK(K) = IQ
      JK(K) = J
C
 220  IF (K .GT. 0)  GOTO 130
C
      ICS = ISS
      IF (IDB .GE. 7) WRITE (6, 2) IP
 2    FORMAT ('MPSORT O'/(8I9))
      RETURN
      END
C
      SUBROUTINE MPSQRT (A, B)
C
C   This computes the square root of the MP number A and returns the MP result
C   in B.  For extra high levels of precision, use MPSQRX.  Debug output
C   starts with IDB = 7.
C
C   Max SP space for B: NW + 4 cells.  Max SP scratch space: 3 * NW + 15
C   cells.  Max DP scratch space: NW + 5 cells.
C
C   This subroutine employs the following Newton-Raphson iteration, which
C   converges to 1 / Sqrt(A):
C
C          X_{k+1} = X_k + 0.5 * (1 - X_k^2 * A) * X_k
C
C   where the muliplication () * X_k is performed with only half of the
C   normal level of precision.  These iterations are performed with a
C   maximum precision level NW that is dynamically changed, doubling with
C   each iteration.  The final iteration is performed as follows (this is
C   due to A. Karp):
C
C          Sqrt(A) = (A * X_n) + 0.5 * [A - (A * X_n)^2] * X_n  (approx.)
C
C   where the multiplications A * X_n and [] * X_n are performed with only
C   half of the final level of precision.  See the comment about the parameter
C   NIT is MPDIVX.
C
      DOUBLE PRECISION CL2, T1, T2
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (CL2 = 1.4426950408889633D0, NDB = 22, NIT = 3)
      DIMENSION A(NW+2), B(NW+4), F(8)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        B(1) = 0.
        B(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 7) THEN
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 1) (A(I), I = 1, NO)
 1      FORMAT ('MPSQRT I'/(6F12.0))
      ENDIF
C
      IA = SIGN (1., A(1))
      NA = MIN (INT (ABS (A(1))), NW)
C
      IF (NA .EQ. 0) THEN
        B(1) = 0.
        B(2) = 0.
        GOTO 120
      ENDIF
      IF (IA .LT. 0.D0) THEN
        IF (KER(70) .NE. 0) THEN
          WRITE (LDB, 2)
 2        FORMAT ('*** MPSQRT: Argument is negative.')
          IER = 70
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
      N5 = NW + 5
      NS = 3 * N5
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N5
      K2 = K1 + N5
      NWS = NW
C
C   Determine the least integer MQ such that 2 ^ MQ .GE. NW.
C
      T1 = NW
      MQ = CL2 * LOG (T1) + 1.D0 - RXX
C
C   Compute the initial approximation of 1 / Sqrt(A).
C
      CALL MPMDC (A, T1, N)
      N2 = - N / 2
      T2 = SQRT (T1 * 2.D0 ** (N + 2 * N2))
      T1 = 1.D0 / T2
      CALL MPDMC (T1, N2, B)
      F(1) = 1.
      F(2) = 0.
      F(3) = 1.
      NW = 3
      IQ = 0
C
C   Perform the Newton-Raphson iteration described above with a dynamically
C   changing precision level NW (one greater than powers of two).
C
      DO 110 K = 2, MQ - 1
        NW1 = NW
        NW = MIN (2 * NW - 2, NWS) + 1
        NW2 = NW
 100    CONTINUE
        CALL MPMUL (B, B, S(K0))
        CALL MPMUL (A, S(K0), S(K1))
        CALL MPSUB (F, S(K1), S(K0))
        NW = NW1
        CALL MPMUL (B, S(K0), S(K1))
        CALL MPMULD (S(K1), 0.5D0, 0, S(K0))
        NW = NW2
        CALL MPADD (B, S(K0), S(K1))
        CALL MPEQ (S(K1), B)
        IF (K .EQ. MQ - NIT .AND. IQ .EQ. 0) THEN
          IQ = 1
          GOTO 100
        ENDIF
 110  CONTINUE
C
C   Perform last iteration using Karp's trick.
C
      CALL MPMUL (A, B, S(K0))
      NW1 = NW
      NW = MIN (2 * NW - 2, NWS) + 1
      NW2 = NW
      CALL MPMUL (S(K0), S(K0), S(K1))
      CALL MPSUB (A, S(K1), S(K2))
      NW = NW1
      CALL MPMUL (S(K2), B, S(K1))
      CALL MPMULD (S(K1), 0.5D0, 0, S(K2))
      NW = NW2
      CALL MPADD (S(K0), S(K2), S(K1))
      CALL MPEQ (S(K1), B)
C
C   Restore original precision level.
C
      NW = NWS
      ICS = ISS
      CALL MPROUN (B)
C
 120  IF (IDB .GE. 7) THEN
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 3) (B(I), I = 1, NO)
 3      FORMAT ('MPSQRT O'/(6F12.0))
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPSQRX (A, B)
C
C   This computes the square root of the MP number A and returns the MP result
C   in B.  Before calling MPSQRX, the array in MPCOM5 must be initialized by
C   calling MPINIX.  For modest levels of precision, use MPSQRT.  NW should be
C   a power of two.  The last three words of the result are not reliable.
C   Debug output starts with IDB = 6.
C
C   Max SP space for B: NW + 4 cells.  Max SP scratch space: 4.5 * NW + 27
C   cells.  Max DP scratch space: 12 * NW + 6 cells.
C
C   This routine uses basically the same Newton iteration algorithm as MPSQRT.
C   In fact, this routine calls MPSQRT to obtain an initial approximation.
C   See the comment about the parameter NIT in MPDIVX.
C
      DOUBLE PRECISION CL2, T1
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (CL2 = 1.4426950408889633D0, NDB = 22, NIT = 1)
      DIMENSION A(NW+2), B(NW+4), F(8)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM3/ S(1024)
C
      IF (IER .NE. 0) THEN
        B(1) = 0.
        B(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 6) THEN
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 1) (A(I), I = 1, NO)
 1      FORMAT ('MPSQRX I'/(6F12.0))
      ENDIF
C
      IA = SIGN (1., A(1))
      NA = MIN (INT (ABS (A(1))), NW)
      NCR = 2 ** MCR
C
      IF (NA .EQ. 0) THEN
        B(1) = 0.
        B(2) = 0.
        GOTO 120
      ENDIF
      IF (IA .LT. 0.D0) THEN
        IF (KER(71) .NE. 0) THEN
          WRITE (LDB, 2)
 2        FORMAT ('*** MPSQRX: Argument is negative.')
          IER = 71
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
        RETURN
      ENDIF
C
C   Check if precision level is too low to justify the advanced routine.
C
      IF (NW .LE. NCR) THEN
        CALL MPSQRT (A, B)
        GOTO 120
      ENDIF
      N4 = NW + 4
      NS = 3 * N4
      ISS = ICS
      ICS = ICS + NS
      IHS = MAX (ICS, IHS)
      IF (ICS - 1 .GT. IMS) CALL MPALER
      K0 = ISS
      K1 = K0 + N4
      K2 = K1 + N4
      NWS = NW
C
C   Determine the least integer MQ such that 2 ^ MQ .GE. NW.
C
      T1 = NW
      MQ = CL2 * LOG (T1) + 1.D0 - RXX
C
C   Compute the initial approximation of 1 / Sqrt(A).
C
      NW = NCR
      CALL MPSQRT (A, S(K0))
      CALL MPDIV (S(K0), A, B)
      F(1) = 1.
      F(2) = 0.
      F(3) = 1.
      IQ = 0
C
C   Perform the Newton-Raphson iteration described above with a dynamically
C   changing precision level NW (powers of two).
C
      DO 110 K = MCR + 1, MQ - 1
        NW1 = NW
        NW = MIN (2 * NW, NWS)
        NW2 = NW
 100    CONTINUE
        CALL MPSQX (B, S(K0))
        CALL MPMULX (A, S(K0), S(K1))
        CALL MPSUB (F, S(K1), S(K0))
        NW = NW1
        CALL MPMULX (B, S(K0), S(K1))
        CALL MPMULD (S(K1), 0.5D0, 0, S(K0))
        NW = NW2
        CALL MPADD (B, S(K0), S(K1))
        CALL MPEQ (S(K1), B)
        IF (K .EQ. MQ - NIT .AND. IQ .EQ. 0) THEN
          IQ = 1
          GOTO 100
        ENDIF
 110  CONTINUE
C
C   Perform last iteration using Karp's trick.
C
      CALL MPMULX (A, B, S(K0))
      NW1 = NW
      NW = MIN (2 * NW, NWS)
      NW2 = NW
      CALL MPSQX (S(K0), S(K1))
      CALL MPSUB (A, S(K1), S(K2))
      NW = NW1
      CALL MPMULX (S(K2), B, S(K1))
      CALL MPMULD (S(K1), 0.5D0, 0, S(K2))
      NW = NW2
      CALL MPADD (S(K0), S(K2), S(K1))
      CALL MPEQ (S(K1), B)
      ICS = ISS
C
 120  IF (IDB .GE. 6) THEN
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 3) (B(I), I = 1, NO)
 3      FORMAT ('MPSQRX O'/(6F12.0))
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPSQX (A, B)
C
C   This routine squares the MP number A to yield the MP product B.
C   Before calling MPSQX, the array in MPCOM5 must be initialized by calling
C   MPINIX.  For modest levels of precision, use MPMUL.  NW should be a power
C   of two.  Debug output starts with IDB = 8.
C
C   Max SP space for B: NW + 4 cells.  Max DP scratch space: 8 * NW + 4 cells.
C
C   This subroutine uses the same FFT technique as MPMULX.  It is faster
C   because only one forward FFT has to be computed.  See the comments in
C   MPMULX about obtaining the complete double-long result.
C>
C   See comments in MPMULX about the machine-dependent parameters ERM and MBT.
C
      DOUBLE PRECISION AN, CL2, D, ERM, T1, T2, T3, T4
      DOUBLE PRECISION BBX, BDX, BX2, RBX, RDX, RX2, RXX
      PARAMETER (CL2 = 1.4426950408889633D0, ERM = 0.438D0, MBT = 53,   
     $  NDB = 22)
      DIMENSION A(NW+2), B(NW+4)
      COMMON /MPCOM0/ BBX, BDX, BX2, RBX, RDX, RX2, RXX, NBT, NPR
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
      COMMON /MPCOM2/ KER(72)
      COMMON /MPCOM4/ D(1024)
C
      IF (IER .NE. 0) THEN
        B(1) = 0.
        B(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 8)  THEN
        NO = MIN (INT (ABS (A(1))), NDB) + 2
        WRITE (LDB, 1) (A(I), I = 1, NO)
 1      FORMAT ('MPSQX I'/(6F12.0))
      ENDIF
C
      IA = SIGN (1., A(1))
      NA = MIN (INT (ABS (A(1))), NW)
      NCR = 2 ** MCR
C
      IF (NA .EQ. 0) THEN
        B(1) = 0.
        B(2) = 0.
        GOTO 170
      ENDIF
C
C   Check if precision level of the argument is too low to justify this
C   advanced routine.
C
      IF (NA .LE. NCR) THEN
        CALL MPMUL (A, A, B)
        GOTO 170
      ENDIF
C
C   Determine N1, the smallest power of two at least as large as NA.
C
      T1 = NA
      M1 = CL2 * LOG (T1) + 1.D0 - RXX
      N1 = 2 ** M1
      M2 = M1 + 2
      N2 = 2 * N1
      N4 = 2 * N2
      N6 = 3 * N2
      N8 = 4 * N2
      N21 = N2 + 1
      N42 = N4 + 2
      N63 = N6 + 3
C
C   Place the input data in A into the scratch array D.
C   This code also splits the input data into half-sized words.
C
CDIR$ IVDEP
      DO 100 I = 1, NA
        T1 = A(I+2)
        T2 = INT (RBX * T1)
        D(2*I-1) = T2
        D(2*I) = T1 - BBX * T2
 100  CONTINUE
C
      DO 110 I = 2 * NA + 1, N2
        D(I) = 0.D0
 110  CONTINUE
C
C   Set the second half of the input vector in D to zero.
C
      DO 120 I = N2 + 1, N4
        D(I) = 0.D0
 120  CONTINUE
C
C   Perform a forward real-to-complex FFT on the vector in D.
C
      CALL MPRCFT (1, M2, D, D(N42+1))
C
C   Square the resulting complex vector.
C
CDIR$ IVDEP
      DO 130 I = 1, N21
        T1 = D(I)
        T2 = D(I+N21)
        D(I+N42) = T1 * T1 - T2 * T2
        D(I+N63) = 2.D0 * T1 * T2
 130  CONTINUE
C
C   Perform an inverse complex-to-real FFT on the resulting data.
C
      CALL MPCRFT (-1, M2, D(N42+1), D)
C
C   Divide by N8, recombine words and release carries.
C
      NB = MIN (2 * NA, NW)
      NB1 = MIN (NW + 1, 2 * NA - 1)
      D(1) = ABS (NB)
      D(2) = 2 * A(2) + 1
      AN = 1.D0 / N8
      T1 = AN * D(N42+1)
      D(3) = AINT (T1 + 0.5D0)
      D(NB+3) = 0.D0
      D(NB+4) = 0.D0
      D(N42+1) = 0.D0
C
CDIR$ IVDEP
      DO 140 I = 1, NB1
        T1 = AN * D(N42+2*I)
        T2 = AN * D(N42+2*I+1)
        T3 = AINT (T1 + 0.5D0)
        T4 = AINT (T2 + 0.5D0)
C        D(N42+2*I) = ABS (T3 - T1)
C        D(N42+2*I+1) = ABS (T4 - T2)
        T1 = INT (RDX * T3)
        T2 = T3 - BDX * T1
        T3 = INT (RDX * T4)
        T4 = T4 - BDX * T3
        D(I+3) = BBX * T2 + T4
        D(I+2) = D(I+2) + BBX * T1 + T3
 140  CONTINUE
C
C   Find the largest FFT roundoff error.  See comments in MPMULX about this
C   test.  To disable this test, uncomment the next line of code and comment
C   out the two lines of the previous loop that begin D(N42..
C
      GOTO 160
      T1 = 0.D0
C
      DO 150 I = 1, 2 * NB1 + 1
        IF (D(N42+I) .GT. T1) THEN
          I1 = I
          T1 = D(N42+I)
        ENDIF
 150  CONTINUE
C
C   Check if maximum roundoff error exceeds the limit ERM, which is set above.
C
      IF (T1 .GT. ERM)  THEN
        IF (KER(55) .NE. 0) THEN
          T2 = AN * D(I1)
          I2 = CL2 * LOG (T1) + 1.D0 + RXX
          I3 = CL2 * LOG (T2) + 1.D0 + RXX
          I4 = MBT + I2 - I3
          I5 = T1 * 2 ** I4 + RXX
          WRITE (LDB, 2) I1, T1, I4, I5
 2        FORMAT ('*** MPSQX: Excessive FFT roundoff error',I10,F10.6,  
     $      2I6)
          IER = 55
          IF (KER(IER) .EQ. 2) CALL MPABRT
        ENDIF
      ENDIF
C
C   Fix up the result.
C
 160  CALL MPNORM (B)
C
 170  IF (IDB .GE. 8) THEN
        NO = MIN (INT (ABS (B(1))), NDB) + 2
        WRITE (LDB, 3) (B(I), I = 1, NO)
 3      FORMAT ('MPSQX O'/(6F12.0))
      ENDIF
      RETURN
      END
C
      SUBROUTINE MPSUB (A, B, C)
C
C   This routine subtracts MP numbers A and B to yield the MP difference C,
C   by negating B and adding.  Debug output starts with IDB = 9.
C
C   Max SP space for C: NW + 4 cells.
C
      DIMENSION A(NW+2), B(NW+2), C(NW+2)
      COMMON /MPCOM1/ NW, IDB, LDB, IER, MCR, IRD, ICS, IHS, IMS
C
      IF (IER .NE. 0) THEN
        C(1) = 0.
        C(2) = 0.
        RETURN
      ENDIF
      IF (IDB .GE. 9) WRITE (LDB, 1)
 1    FORMAT ('MPSUB')
C
C   Check if A = B.  This is necessary because A and B might be same array,
C   in which case negating B below won't work.
C
      IF (A(1) .NE. B(1)) GOTO 110
C
      DO 100 I = 2, INT (ABS (A(1))) + 2
        IF (A(I) .NE. B(I)) GOTO 110
 100  CONTINUE
C
C   A = B.  Result is zero.
C
      C(1) = 0.
      C(2) = 0.
      IF (IDB .GE. 9) WRITE (LDB, 2) (C(I), I = 1, 2)
 2    FORMAT ('MPSUB O'/2F9.0)
      GOTO 120
C
C   Save the sign of B, and then negate B.
C
 110  B1 = B(1)
      B(1) = - B1
C
C   Perform addition and restore the sign of B.
C
      CALL MPADD (A, B, C)
      B(1) = B1
C
 120  RETURN
      END
C
      SUBROUTINE MPTRAN (N1, N2, X, Y)
C
C   Performs a transpose of the vector X, returning the result in Y.  X is
C   treated as a N1 x N2 complex matrix, and Y is treated as a N2 x N1 complex
C   matrix.  The complex data is assumed stored with real and imaginary parts
C   separated by N1 x N2 locations.
C
C   This routine is called by MPCFFT.  It is not intended to be called directly
C   by the user.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NA = 32, NC = 32)
      DIMENSION X(2*N1*N2), Y(2*N1*N2), Z(NC,2*NC)
C
      N = N1 * N2
C>
C   Use different techniques, depending on the system, N1 and N2.  For Cray
C   systems, uncomment the next line.
C
C      GOTO 100
C
C   This strategy is good for many scalar cache memory computers.  The
C   value of NC (i.e. the size of Z) may have to be changed depending on
C   how large the cache is.
C
      IF (N1 .LE. NC .OR. N2 .LE. NC) THEN
        IF (N1 .GE. N2) THEN
          GOTO 110
        ELSE
          GOTO 130
        ENDIF
      ELSE
        GOTO 150
      ENDIF
C
C   This strategy is best for Cray systems.
C
 100  IF (N1 .LT. NA .OR. N2 .LT. NA) THEN
        IF (N1 .GE. N2) THEN
          GOTO 110
        ELSE
          GOTO 130
        ENDIF
      ELSE
        GOTO 220
      ENDIF
C
C   Scheme 1:  Perform a simple transpose in the usual way.
C
 110  DO 120 J = 0, N2 - 1
        J1 = J + 1
        J2 = J * N1 + 1
C
CDIR$ IVDEP
        DO 120 I = 0, N1 - 1
          Y(I*N2+J1) = X(I+J2)
          Y(I*N2+J1+N) = X(I+J2+N)
 120  CONTINUE
C
      GOTO 260
C
C   Scheme 2:  Perform a simple transpose with the loops reversed.
C
 130  DO 140 I = 0, N1 - 1
        I1 = I * N2 + 1
        I2 = I + 1
C
CDIR$ IVDEP
        DO 140 J = 0, N2 - 1
          Y(J+I1) = X(J*N1+I2)
          Y(J+I1+N) = X(J*N1+I2+N)
 140  CONTINUE
C
      GOTO 260
C
C   Scheme 3:  Perform a transpose using the intermediate array Z.  This gives
C   better performance than schemes 1 and 2 on certain cache memory systems.
C   The size of the array Z (i.e. the parameter NC above) may have to be
C   adjusted for optimal performance.
C
 150  DO 210 JJ = 0, N2 - 1, NC
        DO 200 II = 0, N1 - 1, NC
C
          DO 170 J = 1, NC
            J1 = II + (J - 1 + JJ) * N1
C
CDIR$ IVDEP
            DO 160 I = 1, NC
              Z(J,I) = X(I+J1)
              Z(J,I+NC) = X(I+J1+N)
 160        CONTINUE
C
 170      CONTINUE
C
          DO 190 I = 1, NC
            I1 = JJ + (I - 1 + II) * N2
C
CDIR$ IVDEP
            DO 180 J = 1, NC
              Y(J+I1) = Z(J,I)
              Y(J+I1+N) = Z(J,I+NC)
 180        CONTINUE
C
 190      CONTINUE
C
 200    CONTINUE
 210  CONTINUE
C
      GOTO 260
C
C   Scheme 4:  Perform the transpose along diagonals to insure odd strides.
C   This works well on moderate vector, variable stride computers, when both
C   N1 and N2 are divisible by reasonably large powers of two (32 or larger on
C   Cray computers).
C
 220  N11 = N1 + 1
      N21 = N2 + 1
      IF (N1 .GE. N2) THEN
        K1 = N1
        K2 = N2
        I11 = N1
        I12 = 1
        I21 = 1
        I22 = N2
      ELSE
        K1 = N2
        K2 = N1
        I11 = 1
        I12 = N2
        I21 = N1
        I22 = 1
      ENDIF
C
      DO 230 J = 0, K2 - 1
        J1 = J * I11 + 1
        J2 = J * I12 + 1
C
CDIR$ IVDEP
        DO 230 I = 0, K2 - 1 - J
          Y(N21*I+J2) = X(N11*I+J1)
          Y(N21*I+J2+N) = X(N11*I+J1+N)
 230  CONTINUE
C
      DO 240 J = 1, K1 - K2 - 1
        J1 = J * I21 + 1
        J2 = J * I22 + 1
C
CDIR$ IVDEP
        DO 240 I = 0, K2 - 1
          Y(N21*I+J2) = X(N11*I+J1)
          Y(N21*I+J2+N) = X(N11*I+J1+N)
 240  CONTINUE
C
      DO 250 J = K1 - K2, K1 - 1
        J1 = J * I21 + 1
        J2 = J * I22 + 1
C
CDIR$ IVDEP
        DO 250 I = 0, K1 - 1 - J
          Y(N21*I+J2) = X(N11*I+J1)
          Y(N21*I+J2+N) = X(N11*I+J1+N)
 250  CONTINUE
C
 260  RETURN
      END
