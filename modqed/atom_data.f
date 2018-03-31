c       =================================================
        subroutine atom_data
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine stores the data for the nuclear 
c       charge distribution and some other constants
c
c       knucl ..... Nuclear model:
c          Knucl=1  Fermi model
c          Knucl=2  Uniform sphere model
c          Knucl=3  Gaussian model
c       z     ..... Nuclear charge
c       cl    ..... Speed of light
c       R2    ..... The size of the box for the evaluation of wave functions
c       rnucl ..... Nuclear radius
c       rms   ..... Nuclear root-mean-square radius
c      
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /z/z/r2/r2/cl/cl
        common /name/name
        common /knucl/knucl/rnucl/rnucl
        character*1 let(11)
        common /Rrms/Rrms(120)/name_at/name_at(120)/nprms/rms0
        character*2 name,name_at
        data
     1  let /'s','p','d','f','g','h','i','k','l','m','n'/
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call init_atom_data
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        iz=z+0.5d0
        rms=Rrms(iz)
        if (rms0 .gt. 0) rms = rms0
        name=name_at(iz)
c        write( *,15) z
c        write(11,15) z
c15      format (2x,'Z=',f10.4)
c        write( *,25) rms
c        write(11,25) rms
c25      format (2x,'Rms=',f8.4,' fm')
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Defult value of the lihgt velocity.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        cl=137.03599911d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Default value of Atomic box.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        r2=150.0/z
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c        write( *,'(2x,a,f16.9)') 'C=',cl
c        write(11,'(2x,a,f16.9)') 'C=',cl
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (knucl.eq.1.and.iz.le.15) knucl=2
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (dabs(rms).lt.1.d-10) then
          knucl=0
          rnucl=0.d0
c          write( *,'(/2x,a)') 'Point nuuclei.'
c          write(11,'(/2x,a)') 'Point nuuclei.'
        else
          fermi2at=1.d-5/0.529177249d0
          rfermi=rms*dsqrt(5.d0/3.d0)
          rnucl=rfermi*fermi2at
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (knucl.lt.0.or.knucl.gt.3) then
          write( *,'(/2x,a,3x,a,i4)') 'Wrong Nuclear model !!!',
     &    'Nucl=',knucl
          write(11,'(/2x,a,3x,a,i4)') 'Wrong Nuclear model !!!',
     &    'Nucl=',knucl
         call exit1
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
 1000   return
        call exit1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        end
c       =================================================
        subroutine init_atom_data
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /Rrms/Rrms(120)/name_at/name_at(120)
        common /Am/Am(120)
        character*2 name_at
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       The subroutine contains pre-stored data for the nuclear
c       charge radii of different elements
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Nuclear radii are taken from
c       I. Angeli, K. P. Marinova, Atomic Data and Nuclear
c       Data Tables {\bf 99}, 69-95 (2013).
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       If no data for the element are available, we use
c       the empirical formula for Z < 91 from:
c       W.R. Johnson, G.Soff, Atomic Data and Nudear Data
c                     Tables {\bf 33}, p.405 (1985)
c       and for Z> 90 from:
c       I. Goidenko, Private communication.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        name_at( 1)=' H'
        name_at( 2)='He'
        name_at( 3)='Li'
        name_at( 4)='Be'
        name_at( 5)=' B'
        name_at( 6)=' C'
        name_at( 7)=' N'
        name_at( 8)=' O'
        name_at( 9)=' F'
        name_at(10)='Ne'
        name_at(11)='Na'
        name_at(12)='Mg'
        name_at(13)='Al'
        name_at(14)='Al'
        name_at(15)=' P'
        name_at(16)=' S'
        name_at(17)='Cl'
        name_at(18)='Ar'
        name_at(19)='K '
        name_at(20)='Ca'
        name_at(21)='Sc'
        name_at(22)='Ti'
        name_at(23)=' V'
        name_at(24)='Cr'
        name_at(25)='Mn'
        name_at(26)='Fe'
        name_at(27)='Co'
        name_at(28)='Ni'
        name_at(29)='Cu'
        name_at(30)='Zn'
        name_at(31)='Ga'
        name_at(32)='Ge'
        name_at(33)='As'
        name_at(34)='Se'
        name_at(35)='Br'
        name_at(36)='Kr'
        name_at(37)='Rb'
        name_at(38)='Sr'
        name_at(39)=' Y'
        name_at(40)='Zr'
        name_at(41)='Nb'
        name_at(42)='Mo'
        name_at(43)='Tc'
        name_at(44)='Ru'
        name_at(45)='Rh'
        name_at(46)='Pd'
        name_at(47)='Ag'
        name_at(48)='Cd'
        name_at(49)='In'
        name_at(50)='Sn'
        name_at(51)='Sb'
        name_at(52)='Te'
        name_at(53)=' I'
        name_at(54)='Xe'
        name_at(55)='Cs'
        name_at(56)='Ba'
        name_at(57)='La'
        name_at(58)='Ce'
        name_at(59)='Pr'
        name_at(60)='Nd'
        name_at(61)='Pm'
        name_at(62)='Sm'
        name_at(63)='Eu'
        name_at(64)='Gd'
        name_at(65)='Tb'
        name_at(66)='Dy'
        name_at(67)='Ho'
        name_at(68)='Er'
        name_at(69)='Tm'
        name_at(70)='Yb'
        name_at(71)='Lu'
        name_at(72)='Hf'
        name_at(73)='Ta'
        name_at(74)=' W'
        name_at(75)='Re'
        name_at(76)='Os'
        name_at(77)='Ir'
        name_at(78)='Pt'
        name_at(79)='Au'
        name_at(80)='Hg'
        name_at(81)='Tl'
        name_at(82)='Pb'
        name_at(83)='Bi'
        name_at(84)='Po'
        name_at(85)='At'
        name_at(86)='Rn'
        name_at(87)='Fr'
        name_at(88)='Ra'
        name_at(89)='Ac'
        name_at(90)='Th'
        name_at(91)='Pa'
        name_at(92)=' U'
        name_at(93)='Np'
        name_at(94)='Pu'
        name_at(95)='Am'
        name_at(96)='Cm'
        name_at(97)='Bk'
        name_at(98)='Cf'
        name_at(99)='Es'
        name_at(100)='Fm'
        name_at(101)='Md'
        name_at(102)='No'
        name_at(103)='Lr'
c       Rutherfordium
        name_at(104)='Rf'
        name_at(105)='Db'
        name_at(106)='Sg'
        name_at(107)='Bh'
        name_at(108)='Hs'
        name_at(109)='Mt'
        name_at(110)='Ds'
        name_at(111)='Rg'
        name_at(112)='Cn' 
        name_at(113)='Ut'
        name_at(114)='Fl'
        name_at(115)='Up'
        name_at(116)='Lv'
        name_at(117)='Us'
        name_at(118)='Uo'
        name_at(119)='Ue'
        name_at(120)='Ub'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       From NIST: http://www.nist.gov/pml/data/comp.cfm
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        Am(  1) =  1.0
        Am(  2) =  4.0
        Am(  3) =  6.0
        Am(  4) =  9.0
        Am(  5) =  10.0
        Am(  6) =  12.0
        Am(  7) =  14.0
        Am(  8) =  16.0
        Am(  9) =  19.0
        Am( 10) =  20.0
        Am( 11) =  23.0
        Am( 12) =  24.0
        Am( 13) =  27.0
        Am( 14) =  28.0
        Am( 15) =  31.0
        Am( 16) =  32.0
        Am( 17) =  35.0
        Am( 18) =  38.0
        Am( 19) =  39.0
        Am( 20) =  40.0
        Am( 21) =  45.0
        Am( 22) =  48.0
        Am( 23) =  51.0
        Am( 24) =  52.0
        Am( 25) =  55.0
        Am( 26) =  56.0
        Am( 27) =  59.0
        Am( 28) =  60.0
        Am( 29) =  63.0
        Am( 30) =  66.0
        Am( 31) =  69.0
        Am( 32) =  74.0
        Am( 33) =  75.0
        Am( 34) =  80.0
        Am( 35) =  79.0
        Am( 36) =  86.0
        Am( 37) =  87.0
        Am( 38) =  88.0
        Am( 39) =  89.0
        Am( 40) =  90.0
        Am( 41) =  93.0
        Am( 42) =  92.0
        Am( 43) =  97.0
        Am( 44) = 104.0
        Am( 45) = 103.0
        Am( 46) = 108.0
        Am( 47) = 109.0
        Am( 48) = 114.0
        Am( 49) = 115.0
        Am( 50) = 120.0
        Am( 51) = 121.0
        Am( 52) = 130.0
        Am( 53) = 127.0
        Am( 54) = 136.0
        Am( 55) = 133.0
        Am( 56) = 138.0
        Am( 57) = 139.0
        Am( 58) = 140.0
        Am( 59) = 141.0
        Am( 60) = 142.0
        Am( 61) = 145.0
        Am( 62) = 144.0
        Am( 63) = 145.0
        Am( 64) = 160.0
        Am( 65) = 159.0
        Am( 66) = 148.0
        Am( 67) = 165.0
        Am( 68) = 170.0
        Am( 69) = 169.0
        Am( 70) = 176.0
        Am( 71) = 175.0
        Am( 72) = 178.0
        Am( 73) = 181.0
        Am( 74) = 184.0
        Am( 75) = 185.0
        Am( 76) = 192.0
        Am( 77) = 191.0
        Am( 78) = 194.0
        Am( 79) = 197.0
        Am( 80) = 198.0
        Am( 81) = 205.0
        Am( 82) = 208.0
        Am( 83) = 209.0
        Am( 84) = 208.0
        Am( 85) = 210.0
        Am( 86) = 212.0
        Am( 87) = 212.0
        Am( 88) = 214.0
        Am( 89) = 227.0
        Am( 90) = 232.0
        Am( 91) = 232.0
        Am( 92) = 238.0
        Am( 93) = 237.0
        Am( 94) = 239.0
        Am( 95) = 243.0
        Am( 96) = 244.0
        Am( 97) = 247.0
        Am( 98) = 251.0
        Am( 99) = 252.0
        Am(100) = 257.0
        Am(101) = 258.0
        Am(102) = 259.0
        Am(103) = 262.0
        Am(104) = 265.0
        Am(105) = 268.0
        Am(106) = 271.0
        Am(107) = 272.0
        Am(108) = 270.0
        Am(109) = 276.0
        Am(110) = 281.0
        Am(111) = 280.0
        Am(112) = 285.0
        Am(113) = 284.0
        Am(114) = 289.0
        Am(115) = 288.0
        Am(116) = 293.0
        Am(117) = 292.0
        Am(118) = 294.0
        Am(119) = 295.0
        Am(120) = 295.0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do iz=1,120
          if (iz.le.90) then
            rms=(0.836d0*Am(iz)**(1.d0/3.d0)+0.570d0)
          else
            rms=(0.77d0*Am(iz)**(1.d0/3.d0)+0.980d0)
          endif
          Rrms(iz)=rms
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Nuclear R_rms in fm.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        Rrms  (  1) =         0.8783 d0 !#       A=1
        Rrms  (  2) =         1.6755 d0 !#       A=4        
        Rrms  (  3) =         2.5890 d0 !#       A=6
        Rrms  (  4) =         2.5190 d0 !#       A=9
        Rrms  (  5) =         2.4277 d0 !#       A=10
        Rrms  (  6) =         2.4702 d0 !#       A=12
        Rrms  (  7) =         2.5582 d0 !#       A=14
        Rrms  (  8) =         2.6991 d0 !#       A=16
        Rrms  (  9) =         2.8976 d0 !#       A=19
        Rrms  ( 10) =         3.0055 d0 !#       A=20
        Rrms  ( 11) =         2.9936 d0 !#       A=23
        Rrms  ( 12) =         3.0570 d0 !#       A=24
        Rrms  ( 13) =         3.0610 d0 !#       A=27
        Rrms  ( 14) =         3.1224 d0 !#       A=28
        Rrms  ( 15) =         3.1889 d0 !#       A=31
        Rrms  ( 16) =         3.3468 d0 !#       A=32
        Rrms  ( 17) =         3.3654 d0 !#       A=35
        Rrms  ( 18) =         3.4028 d0 !#       A=38
        Rrms  ( 19) =         3.4349 d0 !#       A=39
        Rrms  ( 20) =         3.4776 d0 !#       A=40
        Rrms  ( 21) =         3.5459 d0 !#       A=45
        Rrms  ( 22) =         3.5921 d0 !#       A=48
        Rrms  ( 23) =         3.6002 d0 !#       A=51
        Rrms  ( 24) =         3.6452 d0 !#       A=52
        Rrms  ( 25) =         3.7057 d0 !#       A=55
        Rrms  ( 26) =         3.7377 d0 !#       A=56
        Rrms  ( 27) =         3.7875 d0 !#       A=59
        Rrms  ( 28) =         3.8118 d0 !#       A=60
        Rrms  ( 29) =         3.8823 d0 !#       A=63
        Rrms  ( 30) =         3.9491 d0 !#       A=66
        Rrms  ( 31) =         3.9973 d0 !#       A=69
        Rrms  ( 32) =         4.0742 d0 !#       A=74
        Rrms  ( 33) =         4.0968 d0 !#       A=75
        Rrms  ( 34) =         4.1400 d0 !#       A=80
        Rrms  ( 35) =         4.1629 d0 !#       A=79
        Rrms  ( 36) =         4.1835 d0 !#       A=86
        Rrms  ( 37) =         4.1989 d0 !#       A=87
        Rrms  ( 38) =         4.2240 d0 !#       A=88
        Rrms  ( 39) =         4.2430 d0 !#       A=89
        Rrms  ( 40) =         4.2694 d0 !#       A=90
        Rrms  ( 41) =         4.3240 d0 !#       A=93
        Rrms  ( 42) =         4.3151 d0 !#       A=92
        Rrms  ( 43) =         4.4240 d0 !#       A=97
        Rrms  ( 44) =         4.5098 d0 !#       A=104
        Rrms  ( 45) =         4.4945 d0 !#       A=103
        Rrms  ( 46) =         4.5563 d0 !#       A=108
        Rrms  ( 47) =         4.5638 d0 !#       A=109
        Rrms  ( 48) =         4.6087 d0 !#       A=114
        Rrms  ( 49) =         4.6156 d0 !#       A=115
        Rrms  ( 50) =         4.6519 d0 !#       A=120
        Rrms  ( 51) =         4.6802 d0 !#       A=121
        Rrms  ( 52) =         4.7423 d0 !#       A=130
        Rrms  ( 53) =         4.7964 d0 !#       A=127
        Rrms  ( 54) =         4.5098 d0 !#       A=136
        Rrms  ( 55) =         4.8041 d0 !#       A=133
        Rrms  ( 56) =         4.8378 d0 !#       A=138
        Rrms  ( 57) =         4.8550 d0 !#       A=139
        Rrms  ( 58) =         4.8771 d0 !#       A=140
        Rrms  ( 59) =         4.8919 d0 !#       A=141
        Rrms  ( 60) =         4.9123 d0 !#       A=142
        Rrms  ( 61) =         4.9620 d0 !#       A=145
        Rrms  ( 62) =         4.9524 d0 !#       A=144
        Rrms  ( 63) =         4.9663 d0 !#       A=145
        Rrms  ( 64) =         5.1734 d0 !#       A=160
        Rrms  ( 65) =         5.0600 d0 !#       A=159
        Rrms  ( 66) =         5.0455 d0 !#       A=148
        Rrms  ( 67) =         5.2022 d0 !#       A=165
        Rrms  ( 68) =         5.2789 d0 !#       A=170
        Rrms  ( 69) =         5.2256 d0 !#       A=169
        Rrms  ( 70) =         5.3215 d0 !#       A=176
        Rrms  ( 71) =         5.3770 d0 !#       A=175
        Rrms  ( 72) =         5.3371 d0 !#       A=178
        Rrms  ( 73) =         5.3507 d0 !#       A=181
        Rrms  ( 74) =         5.3658 d0 !#       A=184
        Rrms  ( 75) =         5.3596 d0 !#       A=185
        Rrms  ( 76) =         5.4126 d0 !#       A=192
        Rrms  ( 77) =         5.3968 d0 !#       A=191
        Rrms  ( 78) =         5.4236 d0 !#       A=194
        Rrms  ( 79) =         5.4371 d0 !#       A=197
        Rrms  ( 80) =         5.4463 d0 !#       A=198
        Rrms  ( 81) =         5.4759 d0 !#       A=205
        Rrms  ( 82) =         5.5012 d0 !#       A=208
        Rrms  ( 83) =         5.5211 d0 !#       A=209
        Rrms  ( 84) =         5.5584 d0 !#       A=208
        Rrms  ( 85) =         5.5390 d0 !#       A=210
        Rrms  ( 86) =         5.5915 d0 !#       A=212
        Rrms  ( 87) =         5.5915 d0 !#       A=212
        Rrms  ( 88) =         5.6079 d0 !#       A=214
        Rrms  ( 89) =         5.6700 d0 !#       A=227
        Rrms  ( 90) =         5.7848 d0 !#       A=232
        Rrms  ( 91) =         5.7000 d0 !#
        Rrms  ( 92) =         5.8571 d0 !#       A=238
        Rrms  ( 93) =         5.7440 d0 !#
        Rrms  ( 94) =         5.8601 d0 !#       A=239
        Rrms  ( 95) =         5.9048 d0 !#       A=243
        Rrms  ( 96) =         5.8420 d0 !#       A=244
        Rrms  (100) =         5.8570 d0 !#  fm 
        Rrms  (105) =         5.9190 d0 !#  fm 
        Rrms  (110) =         5.9930 d0 !#  fm 
        Rrms  (115) =         6.0880 d0 !#  fm 
        Rrms  (120) =         6.1750 d0 !#  fm 
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
