*DECK DXPQNU
      SUBROUTINE DXPQNU (NU1, NU2, MU, THETA, ID, PQA, IPQA, IERROR)
C***BEGIN PROLOGUE  DXPQNU
C***SUBSIDIARY
C***PURPOSE  To compute the values of Legendre functions for DXLEGF.
C            This subroutine calculates initial values of P or Q using
C            power series, then performs forward nu-wise recurrence to
C            obtain P(-MU,NU,X), Q(0,NU,X), or Q(1,NU,X). The nu-wise
C            recurrence is stable for P for all mu and for Q for mu=0,1.
C***LIBRARY   SLATEC
C***CATEGORY  C3A2, C9
C***TYPE      DOUBLE PRECISION (XPQNU-S, DXPQNU-D)
C***KEYWORDS  LEGENDRE FUNCTIONS
C***AUTHOR  Smith, John M., (NBS and George Mason University)
C***ROUTINES CALLED  DXADD, DXADJ, DXPSI
C***COMMON BLOCKS    DXBLK1
C***REVISION HISTORY  (YYMMDD)
C   820728  DATE WRITTEN
C   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
C   901019  Revisions to prologue.  (DWL and WRB)
C   901106  Changed all specific intrinsics to generic.  (WRB)
C           Corrected order of sections in prologue and added TYPE
C           section.  (WRB)
C   920127  Revised PURPOSE section of prologue.  (DWL)
C***END PROLOGUE  DXPQNU
      DOUBLE PRECISION A,NU,NU1,NU2,PQ,PQA,DXPSI,R,THETA,W,X,X1,X2,XS,
     1 Y,Z
      DOUBLE PRECISION DI,DMU,PQ1,PQ2,FACTMU,FLOK
      DIMENSION PQA(*),IPQA(*)
      COMMON /DXBLK1/ NBITSF
      SAVE /DXBLK1/
!$OMP THREADPRIVATE(/DXBLK1/)
C
C        J0, IPSIK, AND IPSIX ARE INITIALIZED IN THIS SUBROUTINE.
C        J0 IS THE NUMBER OF TERMS USED IN SERIES EXPANSION
C        IN SUBROUTINE DXPQNU.
C        IPSIK, IPSIX ARE VALUES OF K AND X RESPECTIVELY
C        USED IN THE CALCULATION OF THE DXPSI FUNCTION.
C
C***FIRST EXECUTABLE STATEMENT  DXPQNU
      IERROR=0
      J0=NBITSF
      IPSIK=1+(NBITSF/10)
      IPSIX=5*IPSIK
      IPQ=0
C        FIND NU IN INTERVAL [-.5,.5) IF ID=2  ( CALCULATION OF Q )
      NU=MOD(NU1,1.D0)
      IF(NU.GE..5D0) NU=NU-1.D0
C        FIND NU IN INTERVAL (-1.5,-.5] IF ID=1,3, OR 4  ( CALC. OF P )
      IF(ID.NE.2.AND.NU.GT.-.5D0) NU=NU-1.D0
C        CALCULATE MU FACTORIAL
      K=MU
      DMU=MU
      IF(MU.LE.0) GO TO 60
      FACTMU=1.D0
      IF=0
      DO 50 I=1,K
      FACTMU=FACTMU*I
   50 CALL DXADJ(FACTMU,IF,IERROR)
      IF (IERROR.NE.0) RETURN
   60 IF(K.EQ.0) FACTMU=1.D0
      IF(K.EQ.0) IF=0
C
C        X=COS(THETA)
C        Y=SIN(THETA/2)**2=(1-X)/2=.5-.5*X
C        R=TAN(THETA/2)=SQRT((1-X)/(1+X)
C
      X=COS(THETA)
      Y=SIN(THETA/2.D0)**2
      R=TAN(THETA/2.D0)
C
C        USE ASCENDING SERIES TO CALCULATE TWO VALUES OF P OR Q
C        FOR USE AS STARTING VALUES IN RECURRENCE RELATION.
C
      PQ2=0.0D0
      DO 100 J=1,2
      IPQ1=0
      IF(ID.EQ.2) GO TO 80
C
C        SERIES FOR P ( ID = 1, 3, OR 4 )
C        P(-MU,NU,X)=1./FACTORIAL(MU)*SQRT(((1.-X)/(1.+X))**MU)
C                *SUM(FROM 0 TO J0-1)A(J)*(.5-.5*X)**J
C
      IPQ=0
      PQ=1.D0
      A=1.D0
      IA=0
      DO 65 I=2,J0
      DI=I
      A=A*Y*(DI-2.D0-NU)*(DI-1.D0+NU)/((DI-1.D0+DMU)*(DI-1.D0))
      CALL DXADJ(A,IA,IERROR)
      IF (IERROR.NE.0) RETURN
      IF(A.EQ.0.D0) GO TO 66
      CALL DXADD(PQ,IPQ,A,IA,PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
   65 CONTINUE
   66 CONTINUE
      IF(MU.LE.0) GO TO 90
      X2=R
      X1=PQ
      K=MU
      DO 77 I=1,K
      X1=X1*X2
   77 CALL DXADJ(X1,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
      PQ=X1/FACTMU
      IPQ=IPQ-IF
      CALL DXADJ(PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
      GO TO 90
C
C        Z=-LN(R)=.5*LN((1+X)/(1-X))
C
   80 Z=-LOG(R)
      W=DXPSI(NU+1.D0,IPSIK,IPSIX)
      XS=1.D0/SIN(THETA)
C
C        SERIES SUMMATION FOR Q ( ID = 2 )
C        Q(0,NU,X)=SUM(FROM 0 TO J0-1)((.5*LN((1+X)/(1-X))
C    +DXPSI(J+1,IPSIK,IPSIX)-DXPSI(NU+1,IPSIK,IPSIX)))*A(J)*(.5-.5*X)**J
C
C        Q(1,NU,X)=-SQRT(1./(1.-X**2))+SQRT((1-X)/(1+X))
C             *SUM(FROM 0 T0 J0-1)(-NU*(NU+1)/2*LN((1+X)/(1-X))
C                 +(J-NU)*(J+NU+1)/(2*(J+1))+NU*(NU+1)*
C     (DXPSI(NU+1,IPSIK,IPSIX)-DXPSI(J+1,IPSIK,IPSIX))*A(J)*(.5-.5*X)**J
C
C        NOTE, IN THIS LOOP K=J+1
C
      PQ=0.D0
      IPQ=0
      IA=0
      A=1.D0
      DO 85 K=1,J0
      FLOK=K
      IF(K.EQ.1) GO TO 81
      A=A*Y*(FLOK-2.D0-NU)*(FLOK-1.D0+NU)/((FLOK-1.D0+DMU)*(FLOK-1.D0))
      CALL DXADJ(A,IA,IERROR)
      IF (IERROR.NE.0) RETURN
   81 CONTINUE
      IF(MU.GE.1) GO TO 83
      X1=(DXPSI(FLOK,IPSIK,IPSIX)-W+Z)*A
      IX1=IA
      CALL DXADD(PQ,IPQ,X1,IX1,PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
      GO TO 85
   83 X1=(NU*(NU+1.D0)*(Z-W+DXPSI(FLOK,IPSIK,IPSIX))+(NU-FLOK+1.D0)
     1  *(NU+FLOK)/(2.D0*FLOK))*A
      IX1=IA
      CALL DXADD(PQ,IPQ,X1,IX1,PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
   85 CONTINUE
      IF(MU.GE.1) PQ=-R*PQ
      IXS=0
      IF(MU.GE.1) CALL DXADD(PQ,IPQ,-XS,IXS,PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
      IF(J.EQ.2) MU=-MU
      IF(J.EQ.2) DMU=-DMU
   90 IF(J.EQ.1) PQ2=PQ
      IF(J.EQ.1) IPQ2=IPQ
      NU=NU+1.D0
  100 CONTINUE
      K=0
      IF(NU-1.5D0.LT.NU1) GO TO 120
      K=K+1
      PQA(K)=PQ2
      IPQA(K)=IPQ2
      IF(NU.GT.NU2+.5D0) RETURN
  120 PQ1=PQ
      IPQ1=IPQ
      IF(NU.LT.NU1+.5D0) GO TO 130
      K=K+1
      PQA(K)=PQ
      IPQA(K)=IPQ
      IF(NU.GT.NU2+.5D0) RETURN
C
C        FORWARD NU-WISE RECURRENCE FOR F(MU,NU,X) FOR FIXED MU
C        USING
C        (NU+MU+1)*F(MU,NU,X)=(2.*NU+1)*F(MU,NU,X)-(NU-MU)*F(MU,NU-1,X)
C        WHERE F(MU,NU,X) MAY BE P(-MU,NU,X) OR IF MU IS REPLACED
C        BY -MU THEN F(MU,NU,X) MAY BE Q(MU,NU,X).
C        NOTE, IN THIS LOOP, NU=NU+1
C
  130 X1=(2.D0*NU-1.D0)/(NU+DMU)*X*PQ1
      X2=(NU-1.D0-DMU)/(NU+DMU)*PQ2
      CALL DXADD(X1,IPQ1,-X2,IPQ2,PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
      CALL DXADJ(PQ,IPQ,IERROR)
      IF (IERROR.NE.0) RETURN
      NU=NU+1.D0
      PQ2=PQ1
      IPQ2=IPQ1
      GO TO 120
C
      END
