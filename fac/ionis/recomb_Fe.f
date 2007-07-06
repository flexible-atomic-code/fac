C
      SUBROUTINE RECOMBFe(Z,S,TE,RR,RDIEL)
C 
C=========================================================
C THIS SUBROUTINE COMPUTES THE RECOMBINATION RATES
C (RADIATIVE AND DIELECTRONIQUE) 
C !!!!!!!!!!!!!!!!IRON ONLY!!!!!!!!!!!!!!!!!!!!!!!!!!!
C 
C S IS THE IONICITY OF THE RECOMBINING ION ex: S=2 for FeII
C Z IS THE NUCLEAR CHARGE here Z=26
C TE IS THE TEMPERATURE IN KELVIN
C RR IS THE RADIATIVE RECOMBINATION RATE in CM+3/S
C RDIEL IS THE DIELECTRONIC RECOMBINATION RATE in CM+3/S
C
C !!!!!!!!!! WORK IN DOUBLE PRECISION !!!!!!!!!!!!!!!
C=========================================================
C 
      IMPLICIT REAL*8 (A-H,O-Z) 
      INTEGER Z,S
      DIMENSION AR(26),ETA(26),BETA(26)
      DATA AR/1.42,10.2,33.2,78.,151.,262.,412.,605.,813.,1090.,
     ]1330.,1640.,2000.,1460.,1680.,1910.,2240.,2590.,2960.,3160.,
     ]3490.,3910.,4330.,4770.,5840.,0./ 
      DATA ETA/.891,.843,.746,.682,.699,.728,.759,
     $.790,.810,.829,.828,.834,.836,.597,.602,.601,.579,.567,.557,.534,
     $.521,.523,.531,.537,.546,.0/ 
      DATA BETA/13*0.,5.22,5.07,5.10,5.49,5.65,5.79,6.02,6.22,6.15
     $,5.77,5.49,4.02,0./  
C
      RR=0. 
      RDIEL=0.
      IZ=S-1
C IONS HYDROGENOIDES ( SEATON)
C------------------------------ 
C 
      IF(S.GT.Z.OR.Z.EQ.1) THEN     
      ALAM=157890.*Z*Z/TE
      RR=5.197D-14*Z*DSQRT(ALAM)*(.4288+.5*DLOG(ALAM)+.469*
     &(ALAM**(-1./3.))) 
      RETURN 
      ENDIF
C 
C AUTRES IONS 
C---------------
C 
      TE4=TE/10000.
      RR=1.E-13*AR(IZ)/(TE4**(ETA(IZ)+BETA(IZ)*1.D-2*DLOG10(TE4))) 

C RECOMBINAISON DIELECTRONIQUE
C-----------------------------
      RDIEL=alphadi(Z,S,TE)
C
      RETURN
      END


      function alphadi(n,j,t)

c  returns value of dielectronic recombination rate coefficient.

c  for now, does only iron (n=26)
c  j is the ionicity : j=1 for Fe I, 2 for Fe II ...
c  t is in Kelvin
c  result is cm**3 s**-1

c  data arrays contain coefficients for a sum of terms
c        t**-1.5 * a * exp(-e1*11590./t)

c  Fe XXV - Karim&Bhalla : XXV - Chen : XXIV -  McLH & Romanik
c  XXIII-XXII - Badnell : XXI-XIX - Roszman
c  XVIII - Dasgupta&Whitney : XVII Smith etal scaled to Chen
c  XVI - II : Sarazin, Shull&Woods ?

      IMPLICIT REAL*8 (A-H,O-Z) 
      dimension a(4,25),e1(4,25)

c  data for Shull and van Steenberg for Fe XVI-II

c     data a/0.43,3*0., .256,.452,2*0., .011,.0488,.0801,.529,
c    1  .131,0.,.613,.0849,  .129,.00092,.912,.192,
c    2  .0185,.0953,.079,1.23, .016,.0717,.0906,.739,
c    3  .00567,.0782,.0318,1.263, 
c    4  .0336,.00253,1.92,.181, 1.23,3*0., 
c    5  .12,.12,2*0., .19,.09,2*0., .26,.16,2*0., .24,.17,2*0.,
c    6  .20,.30,2*0., .10,.53,2*0., .14,.46,2*0., .15,.63,2*0.,
c    7  .29,.067,2*0., .25,.025,2*0., .12,.043,2*0., .038,.016,2*0.,
c    8  .015,.0047,2*0., .0084,.0027,2*0., .0016,.000736,2*0./


c  e1 for Fe XIII looks odd

c     data e1/5300.,3*0., 4625.,6000.,2*0., 0.1,36.2,306.,928.,
c    1  73.2,0.,877.,316.,  80.3,39.1,919.,392.,
c    2  13.2,66.6,297.,714., 23.7,85.1,329.,787.,
c    3  16.2,96.0,330.,729.,
c    4  117.,22.5,683.,341., 560.,3*0., 
c    5  24.6,248.,2*0.,39.4,198.,2*0.,36.3,193.,2*0.,75.0,205.,2*0.,
c    6  24.5,155.,2*0.,22.2,144.,2*0.,21.6,136.,2*0.,22.6,136.,2*0.,
c    7  66.7,123.,2*0.,64.7,105.,2*0.,54.2,100.,2*0.,37.3,67.4,2*0.,
c    8  28.6,52.1,2*0.,16.7,31.4,2*0.,5.12,12.9,2*0./

c  SvS scaled for Fe XIII- XVI, with Zhdanov inner shell added :
c  Hahn for Fe XII and Fe VII
c  SvS high T scaled to Hahn for Fe X,XI 
c  Burgess formula for Fe IX
c  SvS scaled to Hahn for Fe VIII, Fe II-VI

      data a/0.43,3*0., .256,.452,2*0., .011,.0488,.0801,.529,
     1  .131,0.,.613,.0849,  .129,.00092,.912,.192,
     2  .0185,.0953,.079,1.23, .016,.0717,.0906,.739,
     3  .00567,.0782,.0318,1.263,
     4  .0336,.00253,1.92,.181, 1.23,3*0.,
     5  .12,.12,.6,0., .19,.09,2*0., .26,.16,2*0., .24,.17,2*0.,
     6  .225,.231,2*0., .10,.28,2*0., .14,.26,2*0., .18,.07,2*0.,
     7  .16,.036,2*0., .092,.041,2*0., .08,.024,2*0., .038,.016,2*0.,
     8  .015,.0047,2*0., .0023,.0027,2*0., .00022,.0001,2*0./

      data e1/5300.,3*0., 4625.,6000.,2*0., 0.1,36.2,306.,928.,
     1  73.2,0.,877.,316.,  80.3,39.1,919.,392.,
     2  13.2,66.6,297.,714., 23.7,85.1,329.,787.,
     3  16.2,96.0,330.,729.,
     4  117.,22.5,683.,341., 560.,3*0.,
     5  24.6,248.,560.,0.,39.4,198.,2*0.,36.3,193.,2*0.,75.0,205.,2*0.,
     6  59.6,362.,2*0.,22.2,144.,2*0.,21.6,136.,2*0.,66.1,129.,2*0.,
     7  66.7,123.,2*0.,45.5,360.,2*0.,54.2,100.,2*0.,37.3,67.4,2*0.,
     8  28.6,52.1,2*0.,16.7,31.4,2*0.,5.12,12.9,2*0./
 



      ad = 0.
      do 10 i = 1,4
   10 ad = ad + a(i,n-j+1) * dexp(-e1(i,n-j+1)*11590./t)

      alphadi = ad * t ** (-1.5)

      return
      end
