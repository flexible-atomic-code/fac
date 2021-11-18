
*
************************************************************************
*
      SUBROUTINE ACOFZ1(Z,AM,NU,NL,AC,ANC,NDL,IOPT)
*
*  Subroutine to calculate electric dipole radial integrals and
*  transition probabilities for hydrogenic ions.  Results are obtained
*  for all transitions NU,LU - NL,LL for a given NU, NL (NU.ne.NL).
*  Direct recurrence on the radial matrix elements is used.
*  Note that an initialisation call is required.  See IOPT=2 below.
*
*  Reference: P.J. Storey and D.G. Hummer, 1991. CPC.
*
*** Non-standard FORTRAN statement follows
      IMPLICIT NONE
*
*         Subroutine arguments are:
*             Z     =  nuclear charge
*             AM    =  nuclear mass in atomic mass units
*             NU    =  upper state principal quantum number
*             NL    =  lower state principal quantum number
*             AC    =  array of transition probabilities (or radial
*                      matrix elements) for all allowable angular
*                      momentum quantum numbers
*             ANC   =  total transition probability from level NU to NL
*             NDL   =  first dimension of AC array in calling program
*             IOPT  =  option switch:
*                   =  0,  Electric dipole transition probabilities
*                          in sec**(-1) are stored in AC and ANC
*                   =  1,  Electric dipole radial matrix elements
*                          in atomic units are stored in AC
*                   =  2,  Initialisation of exponentials and factorials
*                          Input maximum principal quantum number in NU
*
*         Storage of transition probabilities (or radial integrals)
*         is such that
*           AC(L,1)  is the probability for transition NU,L-1 to NL,L
*     and   AC(L,2)  is the probability for transition NU,L   to NL,L-1
*
*  Import
      INTEGER NU,NL,IOPT,NDL
      DOUBLE PRECISION Z,AM
*  Export
      DOUBLE PRECISION AC,ANC
*  Local
      INTEGER I,L,MAX
      DOUBLE PRECISION FAL,ALO,CON,ANU,ANL,
     :                 AI,RP,RM,FP,FM,AL,CU,CL,BU,BL,TWOL,REM,
     :                 ZERO,ONE,TWO,FOUR,TEN,FAL1, ALO1
*
      INTEGER NMAX
      PARAMETER (NMAX=1024)
      DIMENSION FAL(NMAX),ALO(NMAX),AC(NDL,2)
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, FOUR=4.D0, TEN=10.D0)
*
*              Save arrays of factorials and exponentials
*
      SAVE FAL,ALO
!$OMP THREADPRIVATE(FAL,ALO)
*
*              Initialise arrays of factorials and exponentials
*
      IF(IOPT.EQ.2) THEN
           MAX=2*NU
           IF(MAX.GT.NMAX) THEN
             WRITE(6,101)
             IOPT=-1
             RETURN
           ENDIF
           FAL(1) = ZERO
           DO 2 I =2,MAX
              AI = REAL(I-1)
              ALO(I-1) = LOG10(AI)
              FAL(I) = ALO(I-1) + FAL(I-1)
    2      CONTINUE
           RETURN
*
      ELSE
*
*              Reduced mass
*
           REM=ONE/(ONE+ONE/(1822.889D0*AM))
*
           DO 3 I=1,NL+1
              AC(I,1)=ZERO
              AC(I,2)=ZERO
    3      CONTINUE
           ANC=ZERO
*
*              Calculation of radial matrix elements, R
*              ****************************************
*
*              Using the notation R(ll,lu) = R(nl,ll,nu,lu)
*
*              Start with direct calculation of R(nl-1,nl)
*              In addition, use                 R(nl,nl-1) = 0
*              Also initialise the factors cu, cl
*
           ANU=REAL(NU)
           ANL=REAL(NL)
           CON=2.6775015D9*REM**3*Z**6*
     :               (ONE/(ANL*ANL) - ONE/(ANU*ANU))**3
*
           IF (NU .GT. NL) THEN
              FAL1 = FAL(NU-NL)
              ALO1 = ALO(NU-NL)
           ELSE
              FAL1 = 0.0
              ALO1 = 0.0
           ENDIF
           FM = FAL(NU+NL+1) - FAL1 - FAL(NL+NL)
           FM = FM/TWO
     :         + (ANL+TWO)*(ALO(NU)+ALO(NL)+0.602059991327962D0)
     :         + (ANU-ANL-TWO)*ALO1 - (ANU+ANL+TWO)*ALO(NU+NL)
           FM = TEN**FM/Z/REM/FOUR
           FP = ZERO
           CU=SQRT((ANU-ANL)*(ANU+ANL))/(ANU*ANL)
           CL=ONE
*
*               Recur on R :
*               Use            R(l,l+1) and R(l+1,l)      (rm, rp)
*               to construct   R(l-1,l) and R(l,l-1)      (fm, fp)
*
           DO 5 L=NL,1,-1
               AL=REAL(L)
               TWOL=TWO*AL
               IF(L.EQ.NL) GO TO 6
               CU=SQRT((ANU-AL)*(ANU+AL))/(ANU*AL)
               CL=SQRT((ANL-AL)*(ANL+AL))/(ANL*AL)
               FP = ( (TWOL+ONE)*BL*RP + BU*RM) / (TWOL*CU)
               FM = ( BL*RP + (TWOL+ONE)*BU*RM) / (TWOL*CL)
*
*               Calculate transition probabilities
*
    6          AC(L,1) = AL/(TWOL-ONE)*FP*FP*CON
               AC(L,2) = AL/(TWOL+ONE)*FM*FM*CON
               IF(IOPT.EQ.1) THEN
                   AC(L,1)=FP
                   AC(L,2)=FM
               ENDIF
               RP=FP
               RM=FM
               BU=CU
               BL=CL
    5      CONTINUE
*
*               Calculate transition probability for nu to nl
*
           IF(IOPT.EQ.0) THEN
               DO 7 L=1,NL
                   AL=REAL(L)
                   ANC = ANC + (TWO*AL-ONE)*AC(L,1)
     :                       + (TWO*AL+ONE)*AC(L,2)
    7          CONTINUE
               ANC=ANC/(ANU*ANU)
               RETURN
           ENDIF
      ENDIF
*
  101 FORMAT(' INSUFFICIENT WORK SPACE IN ACOFZ1'/
     :       ' INCREASE DIMENSIONS OF FAL AND ALO TO AT LEAST',
     :       ' 2*(MAXIMUM PRINCIPAL QUANTUM NUMBER)'/)
*
      END
*
************************************************************************
*
      SUBROUTINE PIXZ1(Z,AM,NE,NL,PHE,PC,PCP,PNC,NDE,NDL,IOPT,IPCP)
*
*  Subroutine to calculate electric dipole radial integrals and  photo-
*  ionization cross-sections for hydrogenic ions.  Results are obtained
*  for all transitions NL,LL - EU,LU for a given NL and free-electron
*  energy EU.  Direct recurrence on the radial matrix elements is used.
*  Note that an initialisation call is required.   See IOPT=2 below.
*
*  Reference: P.J. Storey and D.G. Hummer, 1991. CPC.
*

********************************************************************
*  Added the option IPCP to indicate whether the partial cross sections
*  are needed. M. F. Gu, 10/25/2001.
*********************************************************************

*** Non-standard FORTRAN statement follows
      IMPLICIT NONE
*
*        Subroutine arguments are:
*            Z     =  nuclear charge
*            AM    =  nuclear mass in atomic mass units
*            NE    =  number of free electron energies
*                     (or photon energies) at which the cross-
*                     -section is to be calculated
*            NL    =  bound state principal quantum number
*            PHE   =  array of free electron energies (in a.u.) or
*                     photon frequencies in Hertz
*            PC    =  array of total photoionization cross-sections
*                     as a function of orbital angular momentum
*                     and frequency
*            PCP   =  array of partial photoionization cross-sections
*                     as a function of orbital angular momentum
*                     and frequency
*            PNC   =  array of photoionization cross-sections
*                     for level NL as a function of frequncy
*            NDE   =  dimension of PHE and PNC arrays in calling program
*            NDL   =  first dimension of PC array in calling program
*            IOPT  =  option switch:
*                  =  0,  array PHE contains free electron energies
*                  =  1,  array PHE contains photon frequencies
*                  =  2,  initialisation of exponentials and factorials,
*                         input maximum principal quantum number in NE
*            IPCP  =  0,  no PCP needed.
*                  =  1,  PCP needed.

*
*      Storage of photoionization cross-sections for states NL,L
*      is such that
*      PC(L+1,IE) is the total cross-section from NL,L summed over
*                    final states at energy IE
*      PCP(L,1,IE) is the cross-section from NL,L   to L-1, at energy IE
*      PCP(L,2,IE) is the cross-section from NL,L-1 to L  , at energy IE
********
* NOTE * If the partial cross-sections PCP are not required, storage can
******** be reduced by removing the array PCP from the calling sequence
*        of PIXZ1,and from the DIMENSION and DOUBLE PRECISION statements
*        Statements using PCP in PIXZ1 must also be deleted; they are
*        marked by **.  The total cross-section for each bound state
*        NL,L will still be returned in PC.
*
*  Import
      INTEGER NL,NE,IOPT,IPCP,NDE,NDL
      DOUBLE PRECISION Z,AM,PHE
*  Export
      DOUBLE PRECISION PC,PCP,PNC
*  Local
      INTEGER I,IE,L,MAX,MUL,POW
      DOUBLE PRECISION FAL,ALO,AW,RYD,CON,KAP,ANL,EU,PCP1,PCP2,
     :                 AI,RP,RM,FP,FM,AL,CU,CL,BU,BL,TWOL,REM,
     :                 TM,TP,FMUL,R0,ZERO,ONE,TWO,TEN,D10,DM10
*
      DIMENSION FAL(1000),ALO(1000),PHE(NDE),PC(NDL,NDE),PCP(NDL,2,NDE),
     :          PNC(NDE)
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0, TEN=10.D0,
     :           D10=1.D10, DM10=1.D-10)
*
*            Save factorials and exponentials
*
      SAVE FAL,ALO
!$OMP THREADPRIVATE(FAL,ALO)
*
*            Initialise arrays of factorials and exponentials
*
      IF(IOPT.EQ.2) THEN
           MAX=2*NE
           IF(MAX.GT.1000) THEN
             WRITE(6,101)
             IOPT=-1
             RETURN
           ENDIF
           FAL(1) = ZERO
           DO 2 I =2,MAX
              AI = REAL(I-1)
              ALO(I-1) = LOG10(AI)
              FAL(I) = ALO(I-1) + FAL(I-1)
    2      CONTINUE
*
      ELSE
*
*
*          Reduced mass
*
         REM=ONE/(ONE+ONE/(1822.889D0*AM))
         DO 3 IE=1,NE
            PNC(IE)=ZERO
            DO 3 I=1,NL+1
               PC(I,IE)=ZERO
**
               IF (IPCP .EQ. 1) THEN 
                 PCP(I,1,IE)=ZERO
                 PCP(I,2,IE)=ZERO
               ENDIF
**
    3    CONTINUE
*
*              Calculation of radial matrix elements, R
*              ****************************************
*
*              Using the notation R(eu,ll,lu) = R(nl,ll,eu,lu)
*              where eu is the free-electron energy
*
*              Start with direct calculation of R(eu,nl-1,nl)
*              In addition, use                 R(eu,nl,nl-1) = 0
*
         ANL=REAL(NL)
* Need bound-free OS in Atomic Units. change the conversion factor.
* M. F. Gu, 10/25/2001.
*        CON=8.5596557D-19*REM*REM*Z*Z/(ANL*ANL)
         CON = 0.6666666667*REM*REM*Z*Z/(ANL*ANL)

*
*              Evaluate log10 ( R(0,nl-1,nl) ) -
*              the value of the radial integral at threshold
*
         R0=(ANL+TWO)*(0.6020599913279624D0+ALO(NL))-FAL(NL+NL)/TWO
     :          -0.8685889638065036D0*ANL-0.504000052812886D0

*
*              Obtain R(eu,nl-1,nl) from R(0,nl-1,nl)
*
         DO 8 IE=1,NE
*
              IF(IOPT.EQ.0) EU=TWO*PHE(IE)/(REM*Z*Z)
              IF(IOPT.EQ.1) THEN
                  EU=3.0396597D-16*PHE(IE)/(REM*Z*Z)-ONE/(ANL*ANL)
              ENDIF
              IF(EU.LT.ZERO) GO TO 8
              KAP=SQRT(EU)
              FM=ONE
              MUL=0
              DO 4 I=1,NL
                  AI=REAL(I)
                  AI=ONE+AI*AI*EU
                  FM=FM*AI
                  IF(FM.GT.D10) THEN
                     POW=LOG10(FM)
                     FM=FM/TEN**POW
                     MUL=MUL+POW
                  ENDIF
    4         CONTINUE
              FP=ZERO
              IF((ANL*KAP).GT.1.D-20) THEN
                  FP=(ANL-ATAN(ANL*KAP)/KAP)
              ENDIF
              FM = R0 + LOG10(FM)/TWO - (ANL+TWO)*LOG10(AI)
     :                + 0.8685889638065036D0*FP + REAL(MUL)/TWO
              FP=ONE
              IF(KAP.GE.0.1D0) THEN
                  FP=6.283185307179586D0/KAP
                  FP=ONE-EXP(-FP)
                  FP=SQRT(ONE/FP)
              ENDIF
              MUL=FM
              FM=FM-MUL
              FM=FP*TEN**FM/(REM*REM*Z*Z)
              FP = 0.D0
              CU=SQRT(ONE+ANL*ANL*KAP*KAP)/(ANL)
              CL=ONE
*
*               Recur on R :
*               Use            R(eu,l,l+1) and R(eu,l+1,l)      (rm, rp)
*               to construct   R(eu,l-1,l) and R(eu,l,l-1)      (fm, fp)
*
              DO 5 L=NL,1,-1
                  AL=REAL(L)
                  TWOL=TWO*AL
                  IF(L.EQ.NL) GO TO 6
                  CU=SQRT(ONE+AL*AL*KAP*KAP)/(AL)
                  CL=SQRT((ANL-AL)*(ANL+AL))/(ANL*AL)
                  FP = ( (TWOL+ONE)*BL*RP + BU*RM) / (TWOL*CU)
                  FM = ( BL*RP + (TWOL+ONE)*BU*RM) / (TWOL*CL)
*
*               Calculate photoionization cross-sections
*
    6             TP = TEN**MUL
                  TM = FM*TP
                  TP = FP*TP
                  PCP1 = AL/(TWOL+ONE)*TP*AI*TP*CON
                  PCP2 = AL/(TWOL-ONE)*TM*AI*TM*CON
**
                  IF (IPCP .EQ. 1) THEN 
                     PCP(L,1,IE) = PCP1
                     PCP(L,2,IE) = PCP2
                  ENDIF
**
                  IF(L.LT.NL)  PC(L+1,IE) = PC(L+1,IE) + PCP1
                  PC(L,IE)   = PC(L,IE) + PCP2
                  FMUL=ONE
                  IF(FM.GT.D10) THEN
                      FMUL=DM10
                      MUL=MUL+10
                  ENDIF
                  RP=FP*FMUL
                  RM=FM*FMUL
                  BU=CU
                  BL=CL
    5         CONTINUE
*
*               Calculate total cross-section from nl.
*
              PNC(IE) = 0.D0
              DO 7 L=1,NL
                  AL=REAL(L)
                  PNC(IE) = PNC(IE) + (TWO*AL-ONE)*PC(L,IE)
    7         CONTINUE
              PNC(IE) = PNC(IE)/(ANL*ANL)
    8    CONTINUE
         RETURN
      ENDIF
*
  101 FORMAT(' INSUFFICIENT WORK SPACE IN PIXZ1'/
     :       ' INCREASE DIMENSIONS OF FAL AND ALO TO AT LEAST'/
     :       ' 2*(MAXIMUM PRINCIPAL QUANTUM NUMBER)'/)
*
      END
*
************************************************************************
*
      SUBROUTINE GFFZ1(Z,T,XLF,G,IOPT,NE,NDE,IFLAG)
*
*          Generates thermally averaged free-free non-relativistic Gaunt
*          factor for a hydrogenic ion of charge Z, with a maximum
*          relative error of .007, (rms fitting error = .001) for
*          temperature and frequency in intervals:
*                     10**-4 le u le 10**1.5,
*                     10**-3 le gams le 10**3,
*          where u = h*nu/k*t  and  gams = z**2*ryd/k*t.  To obtain the
*          stated accuracy, the full number of significant figures in
*          the coefficients must be retained.
*
*          This subroutine uses a two-dimensional Chebyshev expansion
*          computed from expressions given by Karzas and Latter (Ap.J.
*          Suppl., v.6, p.167, 1961) augmented by various limiting forms
*          of energy-specific gaunt factors.
*          Reference: D. G. Hummer, Ap. J. 327, p477, 1988.
*
*          Subroutine arguments are:
*              Z    =  nuclear charge
*              T    =  temperature in degrees kelvin
*              XLF  =  array containing log10 of input frequencies
*                      specified according to value of IOPT
*              G    =  array containing output values of gff
*              IOPT =  option switch for specification of frequencies:
*                   =  0,  array XLF contains log10(nu)
*                   =  1,  array XLF contains log10(h*nu/Ryd)
*                   =  2,  array XLF contains log10(h*nu/k*T) = log10(u)
*              NE   =  length of XLF and G arrays
*              NDE  =  dimension of XLF and G arrays in calling program
*              IFLAG=  1 if gams is out of range; 2 if u is out of range
*
*  Non-standard FORTRAN statement follows
      IMPLICIT NONE
*
*  Import
      INTEGER NE,NDE,IFLAG,IOPT
      DOUBLE PRECISION Z,T,XLF
*  Export
      DOUBLE PRECISION G
*  Local
      INTEGER I,IR,J
      DOUBLE PRECISION XLRKT,TXU,TXG,CON,B,D,C
*
      DIMENSION XLF(NDE),G(NDE),D(8,11),B(11),C(8)
      DATA D/
     .  8.986940175D+00, -4.009515855D+00,  8.808871266D-01,
     .  2.640245111D-02, -4.580645915D-02, -3.568055702D-03,
     .  2.827798067D-03,  3.365860195D-04, -8.006936989D-01,
     .  9.466021705D-01,  9.043402532D-02, -9.608451450D-02,
     . -1.885629865D-02,  1.050313890D-02,  2.800889961D-03,
     . -1.078209202D-03, -3.781305103D-01,  1.102726332D-01,
     . -1.543619180D-02,  8.310561114D-03,  2.179620525D-02,
     .  4.259726289D-03, -4.181588794D-03, -1.770208330D-03,
     .  1.877213132D-02, -1.004885705D-01, -5.483366378D-02,
     . -4.520154409D-03,  8.366530426D-03,  3.700273930D-03,
     .  6.889320423D-04,  9.460313195D-05,  7.300158392D-02,
     .  3.576785497D-03, -4.545307025D-03, -1.017965604D-02,
     . -9.530211924D-03, -3.450186162D-03,  1.040482914D-03,
     .  1.407073544D-03, -1.744671550D-03,  2.864013856D-02,
     .  1.903394837D-02,  7.091074494D-03, -9.668371391D-04,
     . -2.999107465D-03, -1.820642230D-03, -3.874082085D-04,
     . -1.707268366D-02, -4.694254776D-03,  1.311691517D-03,
     .  5.316703136D-03,  5.178193095D-03,  2.451228935D-03,
     . -2.277321615D-05, -8.182359057D-04,  2.567331664D-04,
     . -9.155339970D-03, -6.997479192D-03, -3.571518641D-03,
     . -2.096101038D-04,  1.553822487D-03,  1.509584686D-03,
     .  6.212627837D-04,  4.098322531D-03,  1.635218463D-03,
     . -5.918883504D-04, -2.333091048D-03, -2.484138313D-03,
     . -1.359996060D-03, -5.371426147D-05,  5.553549563D-04,
     .  3.837562402D-05,  2.938325230D-03,  2.393747064D-03,
     .  1.328839809D-03,  9.135013312D-05, -7.137252303D-04,
     . -7.656848158D-04, -3.504683798D-04, -8.491991820D-04,
     . -3.615327726D-04,  3.148015257D-04,  8.909207650D-04,
     .  9.869737522D-04,  6.134671184D-04,  1.068883394D-04,
     . -2.046080100D-04/
*
*          Compute temperature-dependent coefficients for u series
*
      IFLAG=0
*           xlrkt is log(ryd/kt)
      XLRKT=5.1983649D0-DLOG10(T)
      TXG=0.66666667D0*(2.0D0*DLOG10(Z)+XLRKT)
      IF(DABS(TXG).GT.2.0001D0) THEN
           IFLAG=1
           RETURN
      ENDIF
      DO 2 J=1,8
           B(11)=D(J,11)
           B(10)=TXG*B(11)+D(J,10)
           DO 1 IR=9,1,-1
                B(IR)=TXG*B(IR+1)-B(IR+2)+D(J,IR)
    1      CONTINUE
           C(J)=0.25D0*(B(1)-B(3))
    2 CONTINUE
*
*          Sum u series
*
      IF(IOPT.EQ.0) THEN
           CON=0.72727273D0*XLRKT-10.376127D0
      ELSEIF(IOPT.EQ.1) THEN
           CON=0.72727273D0*XLRKT+0.90909091D0
      ELSEIF(IOPT.EQ.2) THEN
           CON=0.90909091D0
      ELSE
           WRITE(6,5) IOPT
           STOP
      ENDIF
      DO 4 I=1,NE
           TXU=0.72727273D0*XLF(I)+CON
           IF(DABS(TXU).GT.2.0001D0) THEN
                IFLAG=2
                G(I)=0.0D0
                GOTO 4
           ENDIF
           B(8)=C(8)
           B(7)=TXU*B(8)+C(7)
           DO 3 IR=6,1,-1
                B(IR)=TXU*B(IR+1)-B(IR+2)+C(IR)
    3      CONTINUE
           G(I)=B(1)-B(3)
    4 CONTINUE
      RETURN
    5 FORMAT(' IOPT=',I2,' IS NOT PERMITTED, GFFZ1 ABORTED')
      END
