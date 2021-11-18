      FUNCTION CLOGAM(Z)
C
C     this routine computes the logarithm of the gamma function gamma(z)
C     for any complex argument 'Z' to any accuracy preset by CALL LOGAM
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 Z,U,V,H,R,CLOGAM,CDIGAM,LOGAM,SER
      DIMENSION B(15),BN(15),BD(15)
C
      DATA LERR /6/, NX0 /6/, NB /15/
      DATA ZERO,ONE,TWO /0D+0,1D+0,2D+0/
      DATA FOUR,HALF,QUART /4D+0,.5D+0,.25D+0/
      DATA B(1),B(2) /8.3333333333333329D-002,-2.7777777777777779D-003/
      DATA B(3),B(4) /7.9365079365079365E-004,-5.9523809523809529E-004/
      DATA B(5),B(6) /8.4175084175084171E-004,-1.9175269175269176E-003/
      DATA B(7),B(8) /6.4102564102564100E-003,-2.9550653594771242E-002/
      DATA B(9),B(10) /0.17964437236883057,-1.3924322169059011/
      DATA B(11),B(12) /13.402864044168393,-156.84828462600203/
      DATA PI /3.1415926535897931/, ALPI /1.1447298858494002/
      DATA HL2P /0.91893853320467267/, ACCUR /2D-16/
      DATA NT /12/, NX0 /6/
      DATA BN(1),BD(1)    / +1D+0,   6D+0 /,
     X     BN(2),BD(2)    / -1D+0,  30D+0 /,
     X     BN(3),BD(3)    / +1D+0,  42D+0 /,
     X     BN(4),BD(4)    / -1D+0,  30D+0 /,
     X     BN(5),BD(5)    / +5D+0,  66D+0 /,
     X     BN(6),BD(6)    /          -691D+0,  2730D+0/,
     X     BN(7),BD(7)    /          +  7D+0,     6D+0/,
     X     BN(8),BD(8)    /         -3617D+0,   510D+0/,
     X     BN(9),BD(9)    /         43867D+0,   798D+0/,
     X     BN(10),BD(10)  /       -174611D+0,   330D+0/,
     X     BN(11),BD(11)  /        854513D+0,   138D+0/,
     X     BN(12),BD(12)  /    -236364091D+0,  2730D+0/,
     X     BN(13),BD(13)  /     + 8553103D+0,     6D+0/,
     X     BN(14),BD(14)  /  -23749461029D+0,   870D+0/,
     X     BN(15),BD(15)  / 8615841276005D+0, 14322D+0/
      DATA FPLMIN / -140D+0 /

C
      X=DBLE(Z)
      T=IMAG(Z)
      MX = INT(DBLE(ACCUR*100 - X))
      IF(ABS(ABS(X)-MX) + ABS(T).LT.ACCUR*50) GO TO 60
      F=ABS(T)
      V=DCMPLX(X,F)
      IF(X .LT. ZERO) V=ONE-V
      H=ZERO
      C=DBLE(V)
      N=NX0-INT(C)
      IF(N .LT. 0) GO TO 30
      H=V
      D=IMAG(V)
      A=ATAN2(D,C)
      IF(N .EQ. 0) GO TO 20
      DO 10 I = 1,N
      C=C+ONE
      V=DCMPLX(C,D)
      H=H*V
   10 A=A+ATAN2(D,C)
   20 H=DCMPLX(HALF*LOG(DBLE(H)**2+IMAG(H)**2),A)
      V=V+ONE
   30 R=ONE/V**2
      SER = B(NT)
      DO 40 J=2,NT
        K = NT+1 - J
   40 SER = B(K) + R*SER
      CLOGAM = HL2P+(V-HALF)*LOG(V)-V + SER/V - H
      IF(X .GE. ZERO) GO TO 50
C
      A= INT(X)-ONE
      C=PI*(X-A)
      D=PI*F
C     E=EXP(-TWO*D)
        E = ZERO
        F = -TWO*D
        IF(F.GT.FPLMIN) E = EXP(F)
      F=SIN(C)
      E= D + HALF*LOG(E*F**2+QUART*(ONE-E)**2)
      F=ATAN2(COS(C)*TANH(D),F)-A*PI
      CLOGAM=ALPI-DCMPLX(E,F)-CLOGAM
C
   50 IF(SIGN(ONE,T) .LT. -HALF) CLOGAM=CONJG(CLOGAM)

      RETURN
C
   60 WRITE(LERR,1000) 'CLOGAM',X
 1000 FORMAT(1X,A6,' ... ARGUMENT IS NON POSITIVE INTEGER = ',F20.2)
      CLOGAM = ZERO
      RETURN
C
      ENTRY CDIGAM(Z)
C
C     this routine computes the logarithmic derivative of the gamma
C     function  psi(Z) = digamma(Z) = d (ln gamma(Z))/dZ  for any
C     complex argument Z, to any accuracy preset by CALL LOGAM(ACC)
C
      U=Z
      X=DBLE(U)
      A=ABS(X)
      IF(ABS(IMAG(U)) + ABS(A + INT(X)) .LT. ACCUR) GO TO 110
      IF(X .LT. ZERO) U=-U
      V=U
      H=ZERO
      N=NX0-INT(A)
      IF(N .LT. 0) GO TO 90
      H=ONE/V
      IF(N .EQ. 0) GO TO 80
      DO 70 I = 1,N
      V=V+ONE
   70 H=H+ONE/V
   80 V=V+ONE
   90 R=ONE/V**2
      SER = B(NT) * (2*NT-1)
      DO 100 J=2,NT
        K = NT+1 - J
  100 SER = B(K)*(2*K-1) + R*SER
      CDIGAM = LOG(V) - HALF/V - R*SER - H
      IF(X .GE. ZERO) RETURN
      H=PI*U
      CDIGAM = CDIGAM + ONE/U + PI*COS(H)/SIN(H)
      RETURN
C
  110 WRITE(LERR,1000) 'CDIGAM',X
      CDIGAM=ZERO
      RETURN
C
      ENTRY LOGAM(ACC)
C
C      initialisation call for calculations to accuracy 'ACC'
C
      LOGAM=ZERO
      RETURN
      END
