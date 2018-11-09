      DOUBLE PRECISION FUNCTION EIONE(IOPT,ARG)                    
C  FOR DOUBLE PRECISION ARG, THIS FUNCTION CALCULATES DOUBLE            
C  PRECISION APPROXIMATIONS TO INTEGRALS RELATED TO THE EXPONENTIAL     
C  INTEGRAL FUNCTION, THE PRECISE INTEGRAL EVALUATED DEPENDING ON       
C  THE VALUE OF IOPT, AS FOLLOWS:                                       
C                                                                       
C     FOR IOPT = 1, THE INTEGRAL (FROM -INFINITY TO                     
C       ARG) OF EXP(T)/T DT WILL BE EVALUATED IF ARG                    
C       IS GREATER THAN 0. IF ARG IS LESS THAN 0.0,                     
C       (-1)*THE INTEGRAL (FROM -ARG TO INFINITY)                       
C       OF EXP(-T)/T DT WILL BE EVALUATED.                              
C     FOR IOPT = 2, THE INTEGRAL (FROM ARG TO                           
C       INFINITY) OF EXP(-T)/T DT WILL BE EVALUATED.                    
C       ARG MUST BE GREATER THAN 0.                                     
C     FOR IOPT = 3, EXP(-ARG)*THE INTEGRAL (FROM                        
C       -INFINITY TO ARG) OF EXP(T)/T DT WILL BE                        
C       EVALUATED IF ARG IS GREATER THAN 0. IF ARG                      
C       IS LESS THAN 0.0, EXP(-ARG)*(-1)*THE                            
C       INTEGRAL (FROM -ARG TO INFINITY) OF                             
C       EXP(-T)/T DT WILL BE EVALUATED.                                 
C                                                                       
C  IER IS AN ERROR PARAMETER USED INTERNALLY, WITH THE                  
C      FOLLOWING MEANING:                                               
C       IER = 0 INDICATES A NORMAL RETURN.  IF IER .NE. 0,              
C         AN ERROR MESSAGE AND THE VALUE OF IER ARE WRITTEN             
C         ON LOGICAL UNIT NUMBER IOUT. SHOULD AN UNDERFLOW              
C         OCCUR, EIONE IS SET TO 0.0 IF IOPT = 1 OR 2.                  
C       IER = 4 INDICATES THAT ARG WAS NEGATIVE                         
C         FOR IOPT = 2. CALCULATION CONTINUES USING                     
C         ABS(ARG).                                                     
C       IER = 5 INDICATES THAT IOPT WAS LESS THAN                       
C         1 OR GREATER THAN 3. EIONE IS SET TO                          
C         MACHINE INFINITY.                                             
C       IER = 6 INDICATES THAT ARG WAS EQUAL TO                         
C         0.0. EIONE IS SET TO MACHINE INFINITY IF                      
C         IOPT = 2 AND NEGATIVE MACHINE INFINITY IF                     
C         IOPT = 1 OR 3.                                                
C       IER = 7 INDICATES THAT AN OVERFLOW WOULD                        
C         HAVE OCCURRED. IF IOPT = 1, EIONE IS SET                      
C         TO MACHINE INFINITY.                                          
C                                                                       
C  IOUT IS THE LOGICAL UNIT NUMBER FOR WRITING.                         
C                                                                       
C  REFERENCES:                                                          
C                                                                       
C       1) W.J. CODY AND HENRY C. THACHER, JR., MATH. COMP. 22, 641     
C               (1968) AND  MATH. COMP. 23, 289 (1969).                 
C       2) NATS-FUNPACK, ARGONNE NATIONAL LABORATORY, ARGONNE CODE      
C               CENTER, ARGONNNE, ILLINOIS, 1974.                       
C                                                                       
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION A(6),B(6),C(8),D(8),E(8),F(8),P1(9),Q1(9),              
     1 P2(9),Q2(8),P3(10),Q3(9),P4(10),Q4(9),P0(6),                     
     2 Q0(6),PX(9),QX(9)                                                
      DATA  SIX/6.D0/,TWELVE/12.D0/,THREE/3.D0/,TWO/2.D0/,              
     1 ONE/1.D0/,HALF/.5D0/,TWENT4/24.D0/,FOUR/4.D0/,                   
     2 FORTY/40.D0/,ZERO/0.D0/                                          
C                                                                       
C  MACHINE-DEPENDENT CONSTANTS                                          
C  XINF = MACHINE INFINITY                                              
C  XMAX = LOGE(MACHINE INFINITY)                                        
C  XMIN = LOGE(SMALLEST POSITIVE REAL CONSTANT)                         
C  XOVR = 83.D0 (SET TO BE APPROX. XMAX - 5.D0 IN ORDER TO              
C                 AVOID PREMATURE OVERFLOWS)                            
C                                                                       
C N.B. IN ALL LIKELIHOOD, THESE CONSTANTS WILL NOT NEED TO              
C      BE RESET WHEN USING EIONE FOR LANDAU-ZENER CALCULATIONS.         
C                                                                       
      DATA  XINF/1.700000000D+38/                                       
      DATA  XMAX/88.02969193111305D0/                                   
      DATA  XMIN/-88.72240904541719D0/                                  
      DATA  XOVR/83.0D0/                                                
      DATA  DEXP40/0.2353852668370200D18/                               
      DATA  X0/0.3725074107813666D00/                                   
      DATA  X01/0.3725074107805994D00/                                  
      DATA  X02/0.7671772501993940D-12/                                 
      DATA  A(1)/ -0.5772156649015328D 00/                              
      DATA  A(2)/  0.7541643136630163D 00/                              
      DATA  A(3)/  0.1298492329273731D 00/                              
      DATA  A(4)/  0.2406813556839774D-01/                              
      DATA  A(5)/  0.1320843092096093D-02/                              
      DATA  A(6)/  0.6577393997532639D-04/                              
      DATA  B(1)/  0.1000000000000000D 01/                              
      DATA  B(2)/  0.4258991938115897D 00/                              
      DATA  B(3)/  0.7977947184102281D-01/                              
      DATA  B(4)/  0.8302084760987714D-02/                              
      DATA  B(5)/  0.4864271383930161D-03/                              
      DATA  B(6)/  0.1306551958228487D-04/                              
      DATA  C(1)/  0.8677459548384432D-07/                              
      DATA  C(2)/  0.9999955193013902D 00/                              
      DATA  C(3)/  0.1184831055549458D 02/                              
      DATA  C(4)/  0.4559306442533897D 02/                              
      DATA  C(5)/  0.6992794512910029D 02/                              
      DATA  C(6)/  0.4252020347688406D 02/                              
      DATA  C(7)/  0.8836718088038437D 01/                              
      DATA  C(8)/  0.4013776649406646D 00/                              
      DATA  D(1)/  0.1000000000000000D 01/                              
      DATA  D(2)/  0.1284819353791566D 02/                              
      DATA  D(3)/  0.5644335695618032D 02/                              
      DATA  D(4)/  0.1066451837699138D 03/                              
      DATA  D(5)/  0.8973110971252896D 02/                              
      DATA  D(6)/  0.3149718491704406D 02/                              
      DATA  D(7)/  0.3795590037621223D 01/                              
      DATA  D(8)/  0.9088045691888691D-01/                              
      DATA  E(1)/ -0.9999999999999731D 00/                              
      DATA  E(2)/ -0.3440619950066848D 02/                              
      DATA  E(3)/ -0.4275326712019885D 03/                              
      DATA  E(4)/ -0.2396019432474904D 04/                              
      DATA  E(5)/ -0.6168852100554763D 04/                              
      DATA  E(6)/ -0.6576096987480218D 04/                              
      DATA  E(7)/ -0.2106077371426331D 04/                              
      DATA  E(8)/ -0.1489908499729481D 02/                              
      DATA  F(1)/  0.1000000000000000D 01/                              
      DATA  F(2)/  0.3640619950064596D 02/                              
      DATA  F(3)/  0.4943450702099040D 03/                              
      DATA  F(4)/  0.3190272374895431D 04/                              
      DATA  F(5)/  0.1033707530858408D 05/                              
      DATA  F(6)/  0.1632414535577834D 05/                              
      DATA  F(7)/  0.1114977528710966D 05/                              
      DATA  F(8)/  0.2378138991021601D 04/                              
      DATA  P0(1)/  0.1525388359511120D 03/                             
      DATA  P0(2)/  0.3402682862739600D 03/                             
      DATA  P0(3)/  0.2597386446160079D 03/                             
      DATA  P0(4)/  0.7787096586760712D 02/                             
      DATA  P0(5)/  0.7460510544921461D 01/                             
      DATA  P0(6)/  0.5574718225325585D-01/                             
      DATA  Q0(1)/  0.1525388359511120D 03/                             
      DATA  Q0(2)/  0.4165377042495159D 03/                             
      DATA  Q0(3)/  0.4171612180903951D 03/                             
      DATA  Q0(4)/  0.1857403824840772D 03/                             
      DATA  Q0(5)/  0.3490362129565328D 02/                             
      DATA  Q0(6)/  0.2000000000000000D 01/                             
      DATA  P1(1)/  0.5531977362081977D 01/                             
      DATA  P1(2)/  0.2063133336244559D 03/                             
      DATA  P1(3)/  0.1427234409068234D 05/                             
      DATA  P1(4)/  0.3685738950923286D 05/                             
      DATA  P1(5)/  0.4493416458218790D 07/                             
      DATA  P1(6)/ -0.1684497821007958D 07/                             
      DATA  P1(7)/  0.3529608047950282D 09/                             
      DATA  P1(8)/ -0.1251949974431755D 09/                             
      DATA  P1(9)/  0.2997849734461850D 10/                             
      DATA  Q1(1)/  0.2562890625000000D 02/                             
      DATA  Q1(2)/ -0.1512618411191135D 04/                             
      DATA  Q1(3)/  0.4268855000903744D 05/                             
      DATA  Q1(4)/ -0.7478772860127960D 06/                             
      DATA  Q1(5)/  0.8857915400539992D 07/                             
      DATA  Q1(6)/ -0.7249035719651191D 08/                             
      DATA  Q1(7)/  0.4014138914734781D 09/                             
      DATA  Q1(8)/ -0.1398366755614922D 10/                             
      DATA  Q1(9)/  0.1279632488038080D 10/                             
      DATA  P2(1)/ -0.2469409834483613D 01/                             
      DATA  P2(2)/ -0.3677831134783113D 02/                             
      DATA  P2(3)/  0.2327302338390390D 02/                             
      DATA  P2(4)/  0.7894722092944569D 01/                             
      DATA  P2(5)/ -0.1941329675144305D 02/                             
      DATA  P2(6)/  0.5886582407532809D 01/                             
      DATA  P2(7)/  0.4181024225628565D 01/                             
      DATA  P2(8)/  0.5731167057445080D 01/                             
      DATA  P2(9)/  0.9989576665165515D 00/                             
      DATA  Q2(1)/  0.2639830073180245D 01/                             
      DATA  Q2(2)/  0.9654052174292799D 03/                             
      DATA  Q2(3)/ -0.8387670841896405D 01/                             
      DATA  Q2(4)/  0.3172794892543692D 03/                             
      DATA  Q2(5)/  0.5231655687345586D 02/                             
      DATA  Q2(6)/  0.3413652125243753D 03/                             
      DATA  Q2(7)/ -0.1991496002312352D 03/                             
      DATA  Q2(8)/  0.1146252532490162D 01/                             
      DATA  P3(1)/ -0.1647721172463462D 01/                             
      DATA  P3(2)/ -0.1860092121726437D 02/                             
      DATA  P3(3)/ -0.1000641913989284D 02/                             
      DATA  P3(4)/ -0.2105740799548040D 02/                             
      DATA  P3(5)/ -0.9134835699998741D 00/                             
      DATA  P3(6)/ -0.3323612579343960D 02/                             
      DATA  P3(7)/  0.2495487730402057D 02/                             
      DATA  P3(8)/  0.2652575818452800D 02/                             
      DATA  P3(9)/ -0.1845086232391277D 01/                             
      DATA  P3(10)/ 0.9999933106160568D 00/                             
      DATA  Q3(1)/  0.9792403599217288D 02/                             
      DATA  Q3(2)/  0.6403800405352414D 02/                             
      DATA  Q3(3)/  0.5994932325667409D 02/                             
      DATA  Q3(4)/  0.2538819315630708D 03/                             
      DATA  Q3(5)/  0.4429413178337928D 02/                             
      DATA  Q3(6)/  0.1192832423968600D 04/                             
      DATA  Q3(7)/  0.1991004470817741D 03/                             
      DATA  Q3(8)/ -0.1093556195391090D 02/                             
      DATA  Q3(9)/  0.1001533852045341D 01/                             
      DATA  P4(1)/  0.1753388012654660D 03/                             
      DATA  P4(2)/ -0.2231276707776324D 03/                             
      DATA  P4(3)/ -0.1819496649298688D 02/                             
      DATA  P4(4)/ -0.2797985286243052D 02/                             
      DATA  P4(5)/ -0.7631477016202536D 01/                             
      DATA  P4(6)/ -0.1528566236369296D 02/                             
      DATA  P4(7)/ -0.7068109778950293D 01/                             
      DATA  P4(8)/ -0.5000066404131309D 01/                             
      DATA  P4(9)/ -0.3000000003209813D 01/                             
      DATA  P4(10)/  0.1000000000000010D 01/                            
      DATA  Q4(1)/  0.3978459771674147D 05/                             
      DATA  Q4(2)/  0.3972771091004144D 01/                             
      DATA  Q4(3)/  0.1377903902357480D 03/                             
      DATA  Q4(4)/  0.1171792205020864D 03/                             
      DATA  Q4(5)/  0.7048318471804246D 02/                             
      DATA  Q4(6)/ -0.1201877635471546D 02/                             
      DATA  Q4(7)/ -0.7992435957763395D 01/                             
      DATA  Q4(8)/ -0.2999998940403249D 01/                             
      DATA  Q4(9)/  0.1999999999990480D 01/                             
C  FIRST EXECUTABLE STATEMENT                                           
      IEND = 8                                                          
      JEND = 9                                                          
      KEND = 6                                                          
      IENDM1 = IEND-1                                                   
      IENDP1 = IEND+1                                                   
      JENDP1 = JEND+1                                                   
      KENDM1 = KEND-1                                                   
      X = ARG                                                           
      IER = 0                                                           
      IF (IOPT.LT.1.OR.IOPT.GT.3) GOTO 250                              
      GOTO (10,200,10), IOPT                                            
C  IOPT = 1 OR 3                                                        
   10 IF (X) 120,230,20                                                 
   20 IF (X.GE.TWELVE) GOTO 70                                          
      IF (X.GE.SIX) GOTO 50                                             
C  X GREATER THAN OR EQUAL TO 0 AND                                     
C  LESS THAN 6 - RATIONAL APPROXIMATION                                 
C  USED IS EXPRESSED IN TERMS OF                                        
C  CHEBYSHEV POLYNOMIALS TO IMPROVE                                     
C  CONDITIONING                                                         
      T = X+X                                                           
      T = T/THREE-TWO                                                   
      PX(1) = ZERO                                                      
      QX(1) = ZERO                                                      
      PX(2) = P1(1)                                                     
      QX(2) = Q1(1)                                                     
      DO 30 I=2,IEND                                                    
      PX(I+1) = T*PX(I)-PX(I-1)+P1(I)                                   
      QX(I+1) = T*QX(I)-QX(I-1)+Q1(I)                                   
   30 CONTINUE                                                          
      SUMP = HALF*T*PX(IENDP1)-PX(IEND)+P1(IENDP1)                      
      SUMQ = HALF*T*QX(IENDP1)-QX(IEND)+Q1(IENDP1)                      
      FRAC = SUMP/SUMQ                                                  
      XMX0 = (X-X01)-X02                                                
      IF (DABS(XMX0).LT.0.037D0) GOTO 40                                
      XX0 = X/X0                                                        
      EIONE = DLOG(XX0)+XMX0*FRAC                                       
      IF (IOPT.EQ.3) EIONE = DEXP(-X)*EIONE                             
      RETURN                                                            
C  EVALUATE APPROXIMATION FOR LN(X/X0)                                  
C  FOR X CLOSE TO X0                                                    
   40 Y = XMX0/X0                                                       
      SUMP = ((((P0(6)*Y+P0(5))*Y+P0(4))*Y+P0(3))*Y+P0(2))*Y+P0(1)      
      SUMQ = ((((Q0(6)*Y+Q0(5))*Y+Q0(4))*Y+Q0(3))*Y+Q0(2))*Y+Q0(1)      
      EIONE = (SUMP/(SUMQ*X0)+FRAC)*XMX0                                
      IF (IOPT.EQ.3) EIONE = DEXP(-X)*EIONE                             
      RETURN                                                            
C  X GREATER THAN OR EQUAL TO 6 AND                                     
C  LESS THAN 12                                                         
   50 FRAC = ZERO                                                       
      DO 60 I=1,IEND                                                    
      FRAC = Q2(I)/(P2(I)+X+FRAC)                                       
   60 CONTINUE                                                          
      EIONE = (P2(IENDP1)+FRAC)/X                                       
      IF (IOPT.NE.3) EIONE = EIONE*DEXP(X)                              
      RETURN                                                            
C  X GREATER THAN OR EQUAL TO 12 AND                                    
C  LESS THAN 24                                                         
   70 IF (X.GE.TWENT4) GOTO 90                                          
      FRAC = ZERO                                                       
      DO 80 I=1,JEND                                                    
      FRAC = Q3(I)/(P3(I)+X+FRAC)                                       
   80 CONTINUE                                                          
      EIONE = (P3(JENDP1)+FRAC)/X                                       
      IF (IOPT.NE.3) EIONE = EIONE*DEXP(X)                              
      RETURN                                                            
C  X GREATER THAN OR EQUAL TO 24                                        
   90 IF ((X.GE.XMAX).AND.(IOPT.LT.3)) GOTO 220                         
      Y = ONE/X                                                         
      FRAC = ZERO                                                       
      DO 100 I=1,JEND                                                   
      FRAC = Q4(I)/(P4(I)+X+FRAC)                                       
  100 CONTINUE                                                          
      FRAC = P4(JENDP1)+FRAC                                            
      EIONE = Y+Y*Y*FRAC                                                
      IF (IOPT.EQ.3) RETURN                                             
      IF (X.GT.XOVR) GOTO 110                                           
      EIONE = EIONE*DEXP(X)                                             
      RETURN                                                            
C  CALCULATION REFORMULATED TO AVOID                                    
C  PREMATURE OVERFLOW                                                   
  110 EIONE = (EIONE*DEXP(X-FORTY))*DEXP40                              
      RETURN                                                            
C  ORIGINAL X WAS NEGATIVE.                                             
  120 Y = -X                                                            
  130 W = ONE/Y                                                         
      IF (Y.GT.FOUR) GOTO 170                                           
      IF (Y.GT.ONE) GOTO 150                                            
C  -X GREATER THAN 0 AND LESS THAN OR                                   
C  EQUAL TO 1                                                           
      SUMB = B(KEND)                                                    
      DO 140 I=1,KENDM1                                                 
      J = KEND-I                                                        
      SUMB = (SUMB*Y)+B(J)                                              
  140 CONTINUE                                                          
      EIONE = DLOG(Y)-(((((A(6)*Y+A(5))*Y+A(4))*Y+A(3))*Y+A(2))*Y+A(1)) 
     1/SUMB                                                             
      IF (IOPT.EQ.3) EIONE = EIONE*DEXP(Y)                              
      GOTO 190                                                          
C  -X GREATER THAN -1 AND LESS THAN OR                                  
C  EQUAL TO 4                                                           
  150 SUMC = C(IEND)                                                    
      SUMD = D(IEND)                                                    
      DO 160 I=1,IENDM1                                                 
      J = IEND-I                                                        
      SUMC = (SUMC*W)+C(J)                                              
      SUMD = (SUMD*W)+D(J)                                              
  160 CONTINUE                                                          
      EIONE = -SUMC/SUMD                                                
      IF (IOPT.EQ.3) RETURN                                             
      EIONE = EIONE*DEXP(-Y)                                            
      GOTO 190                                                          
C  -X GREATER THAN 4                                                    
  170 IF ((-DABS(X).LT.XMIN).AND.(IOPT.LT.3)) GOTO 210                  
      SUME = E(IEND)                                                    
      SUMF = F(IEND)                                                    
      DO 180 I=1,IENDM1                                                 
      J = IEND-I                                                        
      SUME = (SUME*W)+E(J)                                              
      SUMF = (SUMF*W)+F(J)                                              
  180 CONTINUE                                                          
      EIONE = -W*(1.0D0+W*SUME/SUMF)                                    
      IF (IOPT.EQ.3) RETURN                                             
      EIONE = EIONE*DEXP(-Y)                                            
  190 IF (IOPT.EQ.2) EIONE = -EIONE                                     
      IF (IER.EQ.4) GOTO 260                                            
      RETURN                                                            
  200 Y = X                                                             
      IF (Y) 240,230,130                                                
C  ARG IS LESS THAN XMIN CAUSING UNDERFLOW                              
C  FOR IOPT = 1 OR 2. RETURN WITH EIONE = 0.0                           
  210 EIONE = ZERO                                                      
      RETURN                                                            
C  FATAL ERROR - X IS GREATER THAN                                      
C  XMAX CAUSING OVERFLOW                                                
  220 EIONE = XINF                                                      
      IER = 7                                                           
      GOTO 260                                                          
C  FATAL ERROR - ARG = 0                                                
  230 EIONE = -XINF                                                     
      IF (IOPT.EQ.2) EIONE = -EIONE                                     
      IER = 6                                                           
      GOTO 260                                                          
C  FATAL ERROR - ARG IS LESS THAN 0.0 FOR IOPT = 2                      
C  CONTINUE CURRENT CALCULATION WITH ABS(ARG)                           
C  AND STOP WITH IER = 4                                                
  240 IER = 4                                                           
      GOTO 120                                                          
C  FATAL ERROR - IOPT IS OUT OF RANGE                                   
  250 EIONE = XINF                                                      
      IER = 5                                                           
  260 CONTINUE                                                          
C  STOP IF IER. GT. 3   (I.E. UNDERFLOWS ARE ALLOWED                    
C  WITH EIONE BEING SET TO 0.0 FOR IOPT = 1 OR 2.)                      
      WRITE(*,1000) IER                                              
 1000 FORMAT(5X,'ERROR IN FUNCTION EIONE: IER =',I2)                    
      STOP                                                              
      END                                
