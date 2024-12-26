

      subroutine fdm (xnu, alpha, fd)
      implicit double precision (a-h,o-z)

c w. fullerton, imsl houston, apr 1983 edition.
c
c fermi-dirac function of order xnu and argument alpha, where
c xnu is an integer multiple of 0.5 between -0.5 and 4.0
c inclusive.  this routine is accurate everywhere to approx-
c imately 10 digits.
c
      dimension fdscs(21,10), fdmcs(33,10), fdbcs(26,10),
     1  nscs(10), nmcs(10), nbcs(10), alfmax(10)
c
c  ixer, jxer, kxer flags for xerr returns
      data ixer,jxer,kxer/-1,2,2/
c
c series for fm5s    on the interval  0.00000e-01 to  2.71828e+00
c                                        with weighted error   5.75e-11
c                                         log weighted error  10.24
c                               significant figures required  10.32
c                                    decimal places required  10.90
c
      data fdscs(  1, 1) /   2.1468096512e0 /
      data fdscs(  2, 1) /   -.5099961735 5e0 /
      data fdscs(  3, 1) /    .1357694225 8e0 /
      data fdscs(  4, 1) /   -.0379066840 1e0 /
      data fdscs(  5, 1) /    .0108700231 1e0 /
      data fdscs(  6, 1) /   -.0031712624 6e0 /
      data fdscs(  7, 1) /    .0009364562 1e0 /
      data fdscs(  8, 1) /   -.0002790240 1e0 /
      data fdscs(  9, 1) /    .0000837159 4e0 /
      data fdscs( 10, 1) /   -.0000252565 8e0 /
      data fdscs( 11, 1) /    .0000076541 8e0 /
      data fdscs( 12, 1) /   -.0000023283 6e0 /
      data fdscs( 13, 1) /    .0000007105 2e0 /
      data fdscs( 14, 1) /   -.0000002174 1e0 /
      data fdscs( 15, 1) /    .0000000666 8e0 /
      data fdscs( 16, 1) /   -.0000000204 9e0 /
      data fdscs( 17, 1) /    .0000000063 1e0 /
      data fdscs( 18, 1) /   -.0000000019 4e0 /
      data fdscs( 19, 1) /    .0000000006 0e0 /
      data fdscs( 20, 1) /   -.0000000001 8e0 /
      data fdscs( 21, 1) /    .0000000000 57e0 /
c
c series for fm5m    on the interval  1.00000e+00 to  4.00000e+00
c                                        with weighted error   5.27e-11
c                                         log weighted error  10.28
c                               significant figures required  10.55
c                                    decimal places required  10.94
c
      data fdmcs(  1, 1) /   3.7374179482 1e0 /
      data fdmcs(  2, 1) /    .0702840192 74e0 /
      data fdmcs(  3, 1) /    .0049804498 35e0 /
      data fdmcs(  4, 1) /   -.0102460511 86e0 /
      data fdmcs(  5, 1) /    .0046134524 88e0 /
      data fdmcs(  6, 1) /   -.0015084055 45e0 /
      data fdmcs(  7, 1) /    .0004453628 58e0 /
      data fdmcs(  8, 1) /   -.0001323433 48e0 /
      data fdmcs(  9, 1) /    .0000410236 76e0 /
      data fdmcs( 10, 1) /   -.0000130973 69e0 /
      data fdmcs( 11, 1) /    .0000042172 38e0 /
      data fdmcs( 12, 1) /   -.0000013563 51e0 /
      data fdmcs( 13, 1) /    .0000004356 20e0 /
      data fdmcs( 14, 1) /   -.0000001400 58e0 /
      data fdmcs( 15, 1) /    .0000000451 47e0 /
      data fdmcs( 16, 1) /   -.0000000145 91e0 /
      data fdmcs( 17, 1) /    .0000000047 25e0 /
      data fdmcs( 18, 1) /   -.0000000015 33e0 /
      data fdmcs( 19, 1) /    .0000000004 97e0 /
      data fdmcs( 20, 1) /   -.0000000001 61e0 /
      data fdmcs( 21, 1) /    .0000000000 527e0 /
c
c series for fm5b    on the interval  0.00000e-01 to  2.50000e-01
c                                        with weighted error   2.77e-11
c                                         log weighted error  10.56
c                               significant figures required  10.85
c                                    decimal places required  11.26
c
      data fdbcs(  1, 1) /   3.9538986547 2e0 /
      data fdbcs(  2, 1) /   -.0314341746 06e0 /
      data fdbcs(  3, 1) /   -.0086537614 36e0 /
      data fdbcs(  4, 1) /   -.0000437949 13e0 /
      data fdbcs(  5, 1) /    .0003173634 01e0 /
      data fdbcs(  6, 1) /    .0000743852 22e0 /
      data fdbcs(  7, 1) /   -.0000283890 50e0 /
      data fdbcs(  8, 1) /   -.0000089391 79e0 /
      data fdbcs(  9, 1) /    .0000045170 94e0 /
      data fdbcs( 10, 1) /    .0000007855 81e0 /
      data fdbcs( 11, 1) /   -.0000009131 13e0 /
      data fdbcs( 12, 1) /    .0000000758 35e0 /
      data fdbcs( 13, 1) /    .0000001572 30e0 /
      data fdbcs( 14, 1) /   -.0000000671 69e0 /
      data fdbcs( 15, 1) /   -.0000000100 31e0 /
      data fdbcs( 16, 1) /    .0000000194 28e0 /
      data fdbcs( 17, 1) /   -.0000000059 47e0 /
      data fdbcs( 18, 1) /   -.0000000020 95e0 /
      data fdbcs( 19, 1) /    .0000000025 45e0 /
      data fdbcs( 20, 1) /   -.0000000007 47e0 /
      data fdbcs( 21, 1) /   -.0000000002 91e0 /
      data fdbcs( 22, 1) /    .0000000003 71e0 /
      data fdbcs( 23, 1) /   -.0000000001 31e0 /
      data fdbcs( 24, 1) /   -.0000000000 290e0 /
      data fdbcs( 25, 1) /    .0000000000 573e0 /
      data fdbcs( 26, 1) /   -.0000000000 277e0 /
c
c series for f00s    on the interval  0.00000e-01 to  2.71828e+00
c                                        with weighted error   3.18e-11
c                                         log weighted error  10.50
c                               significant figures required  10.36
c                                    decimal places required  11.15
c
      data fdscs(  1, 2) /   1.3659874054e0 /
      data fdscs(  2, 2) /   -.2438973101 9e0 /
      data fdscs(  3, 2) /    .0547680108 5e0 /
      data fdscs(  4, 2) /   -.0135159352 3e0 /
      data fdscs(  5, 2) /    .0035158670 3e0 /
      data fdscs(  6, 2) /   -.0009461111 9e0 /
      data fdscs(  7, 2) /    .0002607200 1e0 /
      data fdscs(  8, 2) /   -.0000731250 4e0 /
      data fdscs(  9, 2) /    .0000207911 1e0 /
      data fdscs( 10, 2) /   -.0000059759 5e0 /
      data fdscs( 11, 2) /    .0000017329 5e0 /
      data fdscs( 12, 2) /   -.0000005062 5e0 /
      data fdscs( 13, 2) /    .0000001488 2e0 /
      data fdscs( 14, 2) /   -.0000000439 8e0 /
      data fdscs( 15, 2) /    .0000000130 5e0 /
      data fdscs( 16, 2) /   -.0000000038 9e0 /
      data fdscs( 17, 2) /    .0000000011 6e0 /
      data fdscs( 18, 2) /   -.0000000003 4e0 /
      data fdscs( 19, 2) /    .0000000001 0e0 /
      data fdscs( 20, 2) /   -.0000000000 31e0 /
c
c series for f00m    on the interval  1.00000e+00 to  4.00000e+00
c                                        with weighted error   6.63e-11
c                                         log weighted error  10.18
c                               significant figures required  10.22
c                                    decimal places required  10.85
c
      data fdmcs(  1, 2) /   2.1728598421e0 /
      data fdmcs(  2, 2) /   -.1263072405 5e0 /
      data fdmcs(  3, 2) /    .0627040952 7e0 /
      data fdmcs(  4, 2) /   -.0248015923 3e0 /
      data fdmcs(  5, 2) /    .0086910218 0e0 /
      data fdmcs(  6, 2) /   -.0028969250 7e0 /
      data fdmcs(  7, 2) /    .0009558659 0e0 /
      data fdmcs(  8, 2) /   -.0003167553 1e0 /
      data fdmcs(  9, 2) /    .0001054756 2e0 /
      data fdmcs( 10, 2) /   -.0000351869 0e0 /
      data fdmcs( 11, 2) /    .0000117376 0e0 /
      data fdmcs( 12, 2) /   -.0000039134 7e0 /
      data fdmcs( 13, 2) /    .0000013044 3e0 /
      data fdmcs( 14, 2) /   -.0000004347 7e0 /
      data fdmcs( 15, 2) /    .0000001449 1e0 /
      data fdmcs( 16, 2) /   -.0000000483 0e0 /
      data fdmcs( 17, 2) /    .0000000161 0e0 /
      data fdmcs( 18, 2) /   -.0000000053 6e0 /
      data fdmcs( 19, 2) /    .0000000017 8e0 /
      data fdmcs( 20, 2) /   -.0000000005 9e0 /
      data fdmcs( 21, 2) /    .0000000001 9e0 /
      data fdmcs( 22, 2) /   -.0000000000 66e0 /
c
c series for f00b    on the interval  0.00000e-01 to  2.50000e-01
c                                        with weighted error   2.19e-11
c                                         log weighted error  10.66
c                               significant figures required  10.66
c                                    decimal places required  11.33
c
      data fdbcs(  1, 2) /   2.0021325464 9e0 /
      data fdbcs(  2, 2) /    .0018258034 26e0 /
      data fdbcs(  3, 2) /    .0011225105 54e0 /
      data fdbcs(  4, 2) /    .0004555439 28e0 /
      data fdbcs(  5, 2) /    .0000880185 91e0 /
      data fdbcs(  6, 2) /   -.0000137060 73e0 /
      data fdbcs(  7, 2) /   -.0000090015 63e0 /
      data fdbcs(  8, 2) /    .0000013512 84e0 /
      data fdbcs(  9, 2) /    .0000010449 37e0 /
      data fdbcs( 10, 2) /   -.0000003164 95e0 /
      data fdbcs( 11, 2) /   -.0000001071 28e0 /
      data fdbcs( 12, 2) /    .0000000783 47e0 /
      data fdbcs( 13, 2) /   -.0000000017 22e0 /
      data fdbcs( 14, 2) /   -.0000000148 37e0 /
      data fdbcs( 15, 2) /    .0000000055 82e0 /
      data fdbcs( 16, 2) /    .0000000010 85e0 /
      data fdbcs( 17, 2) /   -.0000000017 29e0 /
      data fdbcs( 18, 2) /    .0000000005 16e0 /
      data fdbcs( 19, 2) /    .0000000001 80e0 /
      data fdbcs( 20, 2) /   -.0000000002 22e0 /
      data fdbcs( 21, 2) /    .0000000000 691e0 /
      data fdbcs( 22, 2) /    .0000000000 219e0 /
c
c series for f05s    on the interval  0.00000e-01 to  2.71828e+00
c                                        with weighted error   2.87e-11
c                                         log weighted error  10.54
c                               significant figures required  10.38
c                                    decimal places required  11.18
c
      data fdscs(  1, 3) /   1.3284029124e0 /
      data fdscs(  2, 3) /   -.1783849741 0e0 /
      data fdscs(  3, 3) /    .0338871501 3e0 /
      data fdscs(  4, 3) /   -.0074116985 5e0 /
      data fdscs(  5, 3) /    .0017528135 4e0 /
      data fdscs(  6, 3) /   -.0004358562 3e0 /
      data fdscs(  7, 3) /    .0001122556 1e0 /
      data fdscs(  8, 3) /   -.0000296748 6e0 /
      data fdscs(  9, 3) /    .0000080041 3e0 /
      data fdscs( 10, 3) /   -.0000021938 6e0 /
      data fdscs( 11, 3) /    .0000006092 5e0 /
      data fdscs( 12, 3) /   -.0000001710 4e0 /
      data fdscs( 13, 3) /    .0000000484 6e0 /
      data fdscs( 14, 3) /   -.0000000138 4e0 /
      data fdscs( 15, 3) /    .0000000039 8e0 /
      data fdscs( 16, 3) /   -.0000000011 5e0 /
      data fdscs( 17, 3) /    .0000000003 3e0 /
      data fdscs( 18, 3) /   -.0000000000 97e0 /
      data fdscs( 19, 3) /    .0000000000 28e0 /
c
c series for f05m    on the interval  1.00000e+00 to  4.00000e+00
c                                        with weighted error   2.90e-11
c                                         log weighted error  10.54
c                               significant figures required  10.52
c                                    decimal places required  11.23
c
      data fdmcs(  1, 3) /   1.8324909204 1e0 /
      data fdmcs(  2, 3) /   -.2831869649 90e0 /
      data fdmcs(  3, 3) /    .1228778370 71e0 /
      data fdmcs(  4, 3) /   -.0473020885 11e0 /
      data fdmcs(  5, 3) /    .0172527816 16e0 /
      data fdmcs(  6, 3) /   -.0061561464 96e0 /
      data fdmcs(  7, 3) /    .0021769305 63e0 /
      data fdmcs(  8, 3) /   -.0007657618 44e0 /
      data fdmcs(  9, 3) /    .0002681487 29e0 /
      data fdmcs( 10, 3) /   -.0000935040 12e0 /
      data fdmcs( 11, 3) /    .0000324836 07e0 /
      data fdmcs( 12, 3) /   -.0000112489 79e0 /
      data fdmcs( 13, 3) /    .0000038849 42e0 /
      data fdmcs( 14, 3) /   -.0000013385 73e0 /
      data fdmcs( 15, 3) /    .0000004602 70e0 /
      data fdmcs( 16, 3) /   -.0000001579 79e0 /
      data fdmcs( 17, 3) /    .0000000541 36e0 /
      data fdmcs( 18, 3) /   -.0000000185 24e0 /
      data fdmcs( 19, 3) /    .0000000063 30e0 /
      data fdmcs( 20, 3) /   -.0000000021 60e0 /
      data fdmcs( 21, 3) /    .0000000007 36e0 /
      data fdmcs( 22, 3) /   -.0000000002 50e0 /
      data fdmcs( 23, 3) /    .0000000000 854e0 /
      data fdmcs( 24, 3) /   -.0000000000 290e0 /
c
c series for f05b    on the interval  0.00000e-01 to  2.50000e-01
c                                        with weighted error   2.38e-11
c                                         log weighted error  10.62
c                               significant figures required  10.46
c                                    decimal places required  11.27
c
      data fdbcs(  1, 3) /   1.3737552319e0 /
      data fdbcs(  2, 3) /    .0271865732 7e0 /
      data fdbcs(  3, 3) /    .0071381409 2e0 /
      data fdbcs(  4, 3) /    .0001635361 7e0 /
      data fdbcs(  5, 3) /   -.0000117682 0e0 /
      data fdbcs(  6, 3) /   -.0000141643 9e0 /
      data fdbcs(  7, 3) /   -.0000002415 6e0 /
      data fdbcs(  8, 3) /    .0000012719 9e0 /
      data fdbcs(  9, 3) /   -.0000000281 9e0 /
      data fdbcs( 10, 3) /   -.0000001602 3e0 /
      data fdbcs( 11, 3) /    .0000000306 4e0 /
      data fdbcs( 12, 3) /    .0000000187 1e0 /
      data fdbcs( 13, 3) /   -.0000000101 5e0 /
      data fdbcs( 14, 3) /   -.0000000003 8e0 /
      data fdbcs( 15, 3) /    .0000000021 2e0 /
      data fdbcs( 16, 3) /   -.0000000007 1e0 /
      data fdbcs( 17, 3) /   -.0000000001 7e0 /
      data fdbcs( 18, 3) /    .0000000002 3e0 /
      data fdbcs( 19, 3) /   -.0000000000 70e0 /
      data fdbcs( 20, 3) /   -.0000000000 23e0 /
c
c series for f10s    on the interval  0.00000e-01 to  2.71828e+00
c                                        with weighted error   3.46e-11
c                                         log weighted error  10.46
c                               significant figures required  10.38
c                                    decimal places required  11.09
c
      data fdscs(  1, 4) /   1.6098847156e0 /
      data fdscs(  2, 4) /   -.1624060718 4e0 /
      data fdscs(  3, 4) /    .0261468226 3e0 /
      data fdscs(  4, 4) /   -.0050782609 1e0 /
      data fdscs(  5, 4) /    .0010937471 2e0 /
      data fdscs(  6, 4) /   -.0002516893 5e0 /
      data fdscs(  7, 4) /    .0000606609 8e0 /
      data fdscs(  8, 4) /   -.0000151303 0e0 /
      data fdscs(  9, 4) /    .0000038751 7e0 /
      data fdscs( 10, 4) /   -.0000010136 9e0 /
      data fdscs( 11, 4) /    .0000002697 7e0 /
      data fdscs( 12, 4) /   -.0000000728 3e0 /
      data fdscs( 13, 4) /    .0000000199 0e0 /
      data fdscs( 14, 4) /   -.0000000054 9e0 /
      data fdscs( 15, 4) /    .0000000015 3e0 /
      data fdscs( 16, 4) /   -.0000000004 3e0 /
      data fdscs( 17, 4) /    .0000000001 2e0 /
      data fdscs( 18, 4) /   -.0000000000 34e0 /
c
c series for f10m    on the interval  1.00000e+00 to  4.00000e+00
c                                        with weighted error   3.92e-11
c                                         log weighted error  10.41
c                               significant figures required  10.44
c                                    decimal places required  11.11
c
      data fdmcs(  1, 4) /   1.8823916245e0 /
      data fdmcs(  2, 4) /   -.4963882506 6e0 /
      data fdmcs(  3, 4) /    .2216981727 0e0 /
      data fdmcs(  4, 4) /   -.0903287142 0e0 /
      data fdmcs(  5, 4) /    .0352537623 9e0 /
      data fdmcs(  6, 4) /   -.0134370871 3e0 /
      data fdmcs(  7, 4) /    .0050410084 6e0 /
      data fdmcs(  8, 4) /   -.0018681458 3e0 /
      data fdmcs(  9, 4) /    .0006853986 5e0 /
      data fdmcs( 10, 4) /   -.0002493647 5e0 /
      data fdmcs( 11, 4) /    .0000900867 6e0 /
      data fdmcs( 12, 4) /   -.0000323503 7e0 /
      data fdmcs( 13, 4) /    .0000115572 5e0 /
      data fdmcs( 14, 4) /   -.0000041103 4e0 /
      data fdmcs( 15, 4) /    .0000014560 9e0 /
      data fdmcs( 16, 4) /   -.0000005140 2e0 /
      data fdmcs( 17, 4) /    .0000001808 9e0 /
      data fdmcs( 18, 4) /   -.0000000634 8e0 /
      data fdmcs( 19, 4) /    .0000000222 2e0 /
      data fdmcs( 20, 4) /   -.0000000077 6e0 /
      data fdmcs( 21, 4) /    .0000000027 0e0 /
      data fdmcs( 22, 4) /   -.0000000009 4e0 /
      data fdmcs( 23, 4) /    .0000000003 2e0 /
      data fdmcs( 24, 4) /   -.0000000001 1e0 /
      data fdmcs( 25, 4) /    .0000000000 39e0 /
c
c series for f10b    on the interval  0.00000e-01 to  2.50000e-01
c                                        with weighted error   3.24e-11
c                                         log weighted error  10.49
c                               significant figures required  10.22
c                                    decimal places required  11.12
c
      data fdbcs(  1, 4) /   1.0766097168e0 /
      data fdbcs(  2, 4) /    .0509708968 7e0 /
      data fdbcs(  3, 4) /    .0125669010 0e0 /
      data fdbcs(  4, 4) /   -.0001333933 2e0 /
      data fdbcs(  5, 4) /   -.0000390220 0e0 /
      data fdbcs(  6, 4) /   -.0000033829 5e0 /
      data fdbcs(  7, 4) /    .0000018566 7e0 /
      data fdbcs(  8, 4) /    .0000003251 2e0 /
      data fdbcs(  9, 4) /   -.0000001931 9e0 /
      data fdbcs( 10, 4) /   -.0000000183 8e0 /
      data fdbcs( 11, 4) /    .0000000281 7e0 /
      data fdbcs( 12, 4) /   -.0000000030 6e0 /
      data fdbcs( 13, 4) /   -.0000000037 4e0 /
      data fdbcs( 14, 4) /    .0000000016 2e0 /
      data fdbcs( 15, 4) /    .0000000001 6e0 /
      data fdbcs( 16, 4) /   -.0000000003 7e0 /
      data fdbcs( 17, 4) /    .0000000001 1e0 /
      data fdbcs( 18, 4) /    .0000000000 32e0 /
c
c series for f15s    on the interval  0.00000e-01 to  2.71828e+00
c                                        with weighted error   5.17e-11
c                                         log weighted error  10.29
c                               significant figures required  10.34
c                                    decimal places required  10.90
c
      data fdscs(  1, 5) /   2.2601818297e0 /
      data fdscs(  2, 5) /   -.1708722765 0e0 /
      data fdscs(  3, 5) /    .0233363666 4e0 /
      data fdscs(  4, 5) /   -.0040304134 4e0 /
      data fdscs(  5, 5) /    .0007916285 5e0 /
      data fdscs(  6, 5) /   -.0001687845 2e0 /
      data fdscs(  7, 5) /    .0000381078 6e0 /
      data fdscs(  8, 5) /   -.0000089765 6e0 /
      data fdscs(  9, 5) /    .0000021848 5e0 /
      data fdscs( 10, 5) /   -.0000005458 3e0 /
      data fdscs( 11, 5) /    .0000001393 0e0 /
      data fdscs( 12, 5) /   -.0000000361 8e0 /
      data fdscs( 13, 5) /    .0000000095 4e0 /
      data fdscs( 14, 5) /   -.0000000025 4e0 /
      data fdscs( 15, 5) /    .0000000006 8e0 /
      data fdscs( 16, 5) /   -.0000000001 8e0 /
      data fdscs( 17, 5) /    .0000000000 51e0 /
c
c series for f15m    on the interval  1.00000e+00 to  4.00000e+00
c                                        with weighted error   5.45e-11
c                                         log weighted error  10.26
c                               significant figures required  10.43
c                                    decimal places required  10.97
c
      data fdmcs(  1, 5) /   2.2310484506e0 /
      data fdmcs(  2, 5) /   -.8435086752 5e0 /
      data fdmcs(  3, 5) /    .4030118348 2e0 /
      data fdmcs(  4, 5) /   -.1765850337 0e0 /
      data fdmcs(  5, 5) /    .0738231089 9e0 /
      data fdmcs(  6, 5) /   -.0299237390 3e0 /
      data fdmcs(  7, 5) /    .0118550148 7e0 /
      data fdmcs(  8, 5) /   -.0046128514 7e0 /
      data fdmcs(  9, 5) /    .0017688631 8e0 /
      data fdmcs( 10, 5) /   -.0006701511 6e0 /
      data fdmcs( 11, 5) /    .0002513292 1e0 /
      data fdmcs( 12, 5) /   -.0000934452 7e0 /
      data fdmcs( 13, 5) /    .0000344851 9e0 /
      data fdmcs( 14, 5) /   -.0000126440 2e0 /
      data fdmcs( 15, 5) /    .0000046095 2e0 /
      data fdmcs( 16, 5) /   -.0000016719 6e0 /
      data fdmcs( 17, 5) /    .0000006037 1e0 /
      data fdmcs( 18, 5) /   -.0000002171 0e0 /
      data fdmcs( 19, 5) /    .0000000777 9e0 /
      data fdmcs( 20, 5) /   -.0000000277 8e0 /
      data fdmcs( 21, 5) /    .0000000098 9e0 /
      data fdmcs( 22, 5) /   -.0000000035 1e0 /
      data fdmcs( 23, 5) /    .0000000012 4e0 /
      data fdmcs( 24, 5) /   -.0000000004 3e0 /
      data fdmcs( 25, 5) /    .0000000001 5e0 /
      data fdmcs( 26, 5) /   -.0000000000 54e0 /
c
c series for f15b    on the interval  0.00000e-01 to  2.50000e-01
c                                        with weighted error   2.21e-11
c                                         log weighted error  10.66
c                               significant figures required  10.32
c                                    decimal places required  11.28
c
      data fdbcs(  1, 5) /    .9138458031 3e0 /
      data fdbcs(  2, 5) /    .0756461485 3e0 /
      data fdbcs(  3, 5) /    .0185325720 6e0 /
      data fdbcs(  4, 5) /   -.0002173856 5e0 /
      data fdbcs(  5, 5) /   -.0000237328 8e0 /
      data fdbcs(  6, 5) /    .0000042673 3e0 /
      data fdbcs(  7, 5) /    .0000012018 7e0 /
      data fdbcs(  8, 5) /   -.0000002038 1e0 /
      data fdbcs(  9, 5) /   -.0000000983 9e0 /
      data fdbcs( 10, 5) /    .0000000298 3e0 /
      data fdbcs( 11, 5) /    .0000000073 9e0 /
      data fdbcs( 12, 5) /   -.0000000054 3e0 /
      data fdbcs( 13, 5) /    .0000000001 9e0 /
      data fdbcs( 14, 5) /    .0000000008 2e0 /
      data fdbcs( 15, 5) /   -.0000000002 9e0 /
      data fdbcs( 16, 5) /   -.0000000000 48e0 /
      data fdbcs( 17, 5) /    .0000000000 77e0 /
      data fdbcs( 18, 5) /   -.0000000000 22e0 /
c
c series for f20s    on the interval  0.00000e-01 to  2.71828e+00
c                                        with weighted error   2.47e-11
c                                         log weighted error  10.61
c                               significant figures required  10.86
c                                    decimal places required  11.22
c
      data fdscs(  1, 6) /   3.5445815749 6e0 /
      data fdscs(  2, 6) /   -.2001497509 36e0 /
      data fdscs(  3, 6) /    .0231937129 11e0 /
      data fdscs(  4, 6) /   -.0035654858 18e0 /
      data fdscs(  5, 6) /    .0006393090 63e0 /
      data fdscs(  6, 6) /   -.0001264180 87e0 /
      data fdscs(  7, 6) /    .0000267615 70e0 /
      data fdscs(  8, 6) /   -.0000059580 71e0 /
      data fdscs(  9, 6) /    .0000013790 87e0 /
      data fdscs( 10, 6) /   -.0000003292 55e0 /
      data fdscs( 11, 6) /    .0000000806 22e0 /
      data fdscs( 12, 6) /   -.0000000201 61e0 /
      data fdscs( 13, 6) /    .0000000051 32e0 /
      data fdscs( 14, 6) /   -.0000000013 26e0 /
      data fdscs( 15, 6) /    .0000000003 47e0 /
      data fdscs( 16, 6) /   -.0000000000 921e0 /
      data fdscs( 17, 6) /    .0000000000 246e0 /
c
c series for f20m    on the interval  1.00000e+00 to  4.00000e+00
c                                        with weighted error   2.78e-11
c                                         log weighted error  10.56
c                               significant figures required  10.91
c                                    decimal places required  11.28
c
      data fdmcs(  1, 6) /   2.9777839001 0e0 /
      data fdmcs(  2, 6) /  -1.4577413536 2e0 /
      data fdmcs(  3, 6) /    .7528269512 35e0 /
      data fdmcs(  4, 6) /   -.3549647428 05e0 /
      data fdmcs(  5, 6) /    .1584014924 88e0 /
      data fdmcs(  6, 6) /   -.0680073485 74e0 /
      data fdmcs(  7, 6) /    .0283569566 67e0 /
      data fdmcs(  8, 6) /   -.0115545568 43e0 /
      data fdmcs(  9, 6) /    .0046209871 95e0 /
      data fdmcs( 10, 6) /   -.0018197198 20e0 /
      data fdmcs( 11, 6) /    .0007073370 31e0 /
      data fdmcs( 12, 6) /   -.0002719114 95e0 /
      data fdmcs( 13, 6) /    .0001035295 30e0 /
      data fdmcs( 14, 6) /   -.0000390900 34e0 /
      data fdmcs( 15, 6) /    .0000146509 87e0 /
      data fdmcs( 16, 6) /   -.0000054554 02e0 /
      data fdmcs( 17, 6) /    .0000020195 19e0 /
      data fdmcs( 18, 6) /   -.0000007436 80e0 /
      data fdmcs( 19, 6) /    .0000002725 59e0 /
      data fdmcs( 20, 6) /   -.0000000994 63e0 /
      data fdmcs( 21, 6) /    .0000000361 53e0 /
      data fdmcs( 22, 6) /   -.0000000130 94e0 /
      data fdmcs( 23, 6) /    .0000000047 26e0 /
      data fdmcs( 24, 6) /   -.0000000017 01e0 /
      data fdmcs( 25, 6) /    .0000000006 10e0 /
      data fdmcs( 26, 6) /   -.0000000002 18e0 /
      data fdmcs( 27, 6) /    .0000000000 780e0 /
      data fdmcs( 28, 6) /   -.0000000000 277e0 /
c
c series for f20b    on the interval  0.00000e-01 to  2.50000e-01
c                                        with weighted error   1.41e-11
c                                         log weighted error  10.85
c                               significant figures required  10.48
c                                    decimal places required  11.47
c
      data fdbcs(  1, 6) /    .8211121274 8e0 /
      data fdbcs(  2, 6) /    .1030146856 2e0 /
      data fdbcs(  3, 6) /    .0258442758 9e0 /
      data fdbcs(  4, 6) /    .0000739478 7e0 /
      data fdbcs(  5, 6) /    .0000269632 3e0 /
      data fdbcs(  6, 6) /    .0000055393 5e0 /
      data fdbcs(  7, 6) /   -.0000000666 0e0 /
      data fdbcs(  8, 6) /   -.0000002863 2e0 /
      data fdbcs(  9, 6) /    .0000000098 7e0 /
      data fdbcs( 10, 6) /    .0000000250 2e0 /
      data fdbcs( 11, 6) /   -.0000000043 8e0 /
      data fdbcs( 12, 6) /   -.0000000022 7e0 /
      data fdbcs( 13, 6) /    .0000000011 2e0 /
      data fdbcs( 14, 6) /    .0000000000 40e0 /
      data fdbcs( 15, 6) /   -.0000000001 9e0 /
      data fdbcs( 16, 6) /    .0000000000 60e0 /
      data fdbcs( 17, 6) /    .0000000000 14e0 /
c
c series for f25s    on the interval  0.00000e-01 to  2.71828e+00
c                                        with weighted error   4.95e-11
c                                         log weighted error  10.31
c                               significant figures required  10.79
c                                    decimal places required  10.91
c
      data fdscs(  1, 7) /   6.0776352655 5e0 /
      data fdscs(  2, 7) /   -.2553089335 90e0 /
      data fdscs(  3, 7) /    .0250962593 28e0 /
      data fdscs(  4, 7) /   -.0034359138 78e0 /
      data fdscs(  5, 7) /    .0005628501 70e0 /
      data fdscs(  6, 7) /   -.0001033045 44e0 /
      data fdscs(  7, 7) /    .0000205192 58e0 /
      data fdscs(  8, 7) /   -.0000043206 22e0 /
      data fdscs(  9, 7) /    .0000009516 32e0 /
      data fdscs( 10, 7) /   -.0000002172 43e0 /
      data fdscs( 11, 7) /    .0000000510 64e0 /
      data fdscs( 12, 7) /   -.0000000122 98e0 /
      data fdscs( 13, 7) /    .0000000030 23e0 /
      data fdscs( 14, 7) /   -.0000000007 56e0 /
      data fdscs( 15, 7) /    .0000000001 92e0 /
      data fdscs( 16, 7) /   -.0000000000 495e0 /
c
c series for f25m    on the interval  1.00000e+00 to  4.00000e+00
c                                        with weighted error   4.15e-11
c                                         log weighted error  10.38
c                               significant figures required  10.96
c                                    decimal places required  11.11
c
      data fdmcs(  1, 7) /   4.4112295532 9e0 /
      data fdmcs(  2, 7) /  -2.6064119236 9e0 /
      data fdmcs(  3, 7) /   1.4541147675 7e0 /
      data fdmcs(  4, 7) /   -.7348922045 75e0 /
      data fdmcs(  5, 7) /    .3485252773 58e0 /
      data fdmcs(  6, 7) /   -.1579050700 24e0 /
      data fdmcs(  7, 7) /    .0690927288 72e0 /
      data fdmcs(  8, 7) /   -.0294110574 69e0 /
      data fdmcs(  9, 7) /    .0122427694 58e0 /
      data fdmcs( 10, 7) /   -.0050026362 88e0 /
      data fdmcs( 11, 7) /    .0020124730 44e0 /
      data fdmcs( 12, 7) /   -.0007988327 90e0 /
      data fdmcs( 13, 7) /    .0003134423 09e0 /
      data fdmcs( 14, 7) /   -.0001217494 74e0 /
      data fdmcs( 15, 7) /    .0000468708 54e0 /
      data fdmcs( 16, 7) /   -.0000179017 70e0 /
      data fdmcs( 17, 7) /    .0000067890 45e0 /
      data fdmcs( 18, 7) /   -.0000025582 83e0 /
      data fdmcs( 19, 7) /    .0000009584 71e0 /
      data fdmcs( 20, 7) /   -.0000003572 13e0 /
      data fdmcs( 21, 7) /    .0000001324 92e0 /
      data fdmcs( 22, 7) /   -.0000000489 26e0 /
      data fdmcs( 23, 7) /    .0000000179 94e0 /
      data fdmcs( 24, 7) /   -.0000000065 93e0 /
      data fdmcs( 25, 7) /    .0000000024 07e0 /
      data fdmcs( 26, 7) /   -.0000000008 76e0 /
      data fdmcs( 27, 7) /    .0000000003 17e0 /
      data fdmcs( 28, 7) /   -.0000000001 15e0 /
      data fdmcs( 29, 7) /    .0000000000 415e0 /
c
c series for f25b    on the interval  0.00000e-01 to  2.50000e-01
c                                        with weighted error   4.96e-11
c                                         log weighted error  10.30
c                               significant figures required   9.92
c                                    decimal places required  10.91
c
      data fdbcs(  1, 7) /    .7721541990 3e0 /
      data fdbcs(  2, 7) /    .1348979022 5e0 /
      data fdbcs(  3, 7) /    .0353565117 1e0 /
      data fdbcs(  4, 7) /    .0009476728 1e0 /
      data fdbcs(  5, 7) /    .0001277291 1e0 /
      data fdbcs(  6, 7) /    .0000007218 8e0 /
      data fdbcs(  7, 7) /   -.0000009206 8e0 /
      data fdbcs(  8, 7) /   -.0000001154 7e0 /
      data fdbcs(  9, 7) /    .0000000580 6e0 /
      data fdbcs( 10, 7) /    .0000000047 4e0 /
      data fdbcs( 11, 7) /   -.0000000060 8e0 /
      data fdbcs( 12, 7) /    .0000000005 1e0 /
      data fdbcs( 13, 7) /    .0000000006 5e0 /
      data fdbcs( 14, 7) /   -.0000000002 4e0 /
      data fdbcs( 15, 7) /   -.0000000000 28e0 /
      data fdbcs( 16, 7) /    .0000000000 49e0 /
c
c series for f30s    on the interval  0.00000e-01 to  2.71828e+00
c                                        with weighted error   2.87e-11
c                                         log weighted error  10.54
c                               significant figures required  11.29
c                                    decimal places required  11.14
c
      data fdscs(  1, 8) /  11.2341939777e0 /
      data fdscs(  2, 8) /   -.3495789894 02e0 /
      data fdscs(  3, 8) /    .0291275872 60e0 /
      data fdscs(  4, 8) /   -.0035525827 96e0 /
      data fdscs(  5, 8) /    .0005319821 79e0 /
      data fdscs(  6, 8) /   -.0000906823 61e0 /
      data fdscs(  7, 8) /    .0000169110 38e0 /
      data fdscs(  8, 8) /   -.0000033697 23e0 /
      data fdscs(  9, 8) /    .0000007066 16e0 /
      data fdscs( 10, 8) /   -.0000001543 15e0 /
      data fdscs( 11, 8) /    .0000000348 35e0 /
      data fdscs( 12, 8) /   -.0000000080 83e0 /
      data fdscs( 13, 8) /    .0000000019 20e0 /
      data fdscs( 14, 8) /   -.0000000004 65e0 /
      data fdscs( 15, 8) /    .0000000001 14e0 /
      data fdscs( 16, 8) /   -.0000000000 287e0 /
c
c series for f30m    on the interval  1.00000e+00 to  4.00000e+00
c                                        with weighted error   2.33e-11
c                                         log weighted error  10.63
c                               significant figures required  11.47
c                                    decimal places required  11.38
c
      data fdmcs(  1, 8) /   7.1719763668 6e0 /
      data fdmcs(  2, 8) /  -4.8571185185 4e0 /
      data fdmcs(  3, 8) /   2.9113197795 3e0 /
      data fdmcs(  4, 8) /  -1.5684642300 9e0 /
      data fdmcs(  5, 8) /    .7869998941 23e0 /
      data fdmcs(  6, 8) /   -.3749544936 90e0 /
      data fdmcs(  7, 8) /    .1716880079 87e0 /
      data fdmcs(  8, 8) /   -.0761760305 76e0 /
      data fdmcs(  9, 8) /    .0329421355 00e0 /
      data fdmcs( 10, 8) /   -.0139450242 81e0 /
      data fdmcs( 11, 8) /    .0057976754 88e0 /
      data fdmcs( 12, 8) /   -.0023734227 03e0 /
      data fdmcs( 13, 8) /    .0009586830 99e0 /
      data fdmcs( 14, 8) /   -.0003827164 22e0 /
      data fdmcs( 15, 8) /    .0001512084 34e0 /
      data fdmcs( 16, 8) /   -.0000591925 75e0 /
      data fdmcs( 17, 8) /    .0000229809 46e0 /
      data fdmcs( 18, 8) /   -.0000088559 00e0 /
      data fdmcs( 19, 8) /    .0000033897 35e0 /
      data fdmcs( 20, 8) /   -.0000012895 26e0 /
      data fdmcs( 21, 8) /    .0000004878 14e0 /
      data fdmcs( 22, 8) /   -.0000001835 85e0 /
      data fdmcs( 23, 8) /    .0000000687 64e0 /
      data fdmcs( 24, 8) /   -.0000000256 43e0 /
      data fdmcs( 25, 8) /    .0000000095 24e0 /
      data fdmcs( 26, 8) /   -.0000000035 23e0 /
      data fdmcs( 27, 8) /    .0000000012 99e0 /
      data fdmcs( 28, 8) /   -.0000000004 77e0 /
      data fdmcs( 29, 8) /    .0000000001 74e0 /
      data fdmcs( 30, 8) /   -.0000000000 638e0 /
      data fdmcs( 31, 8) /    .0000000000 232e0 /
c
c series for f30b    on the interval  0.00000e-01 to  2.50000e-01
c                                        with weighted error   5.48e-11
c                                         log weighted error  10.26
c                               significant figures required   9.88
c                                    decimal places required  10.85
c
      data fdbcs(  1, 8) /    .7554309639 2e0 /
      data fdbcs(  2, 8) /    .1734863060 3e0 /
      data fdbcs(  3, 8) /    .0481579481 5e0 /
      data fdbcs(  4, 8) /    .0027149874 6e0 /
      data fdbcs(  5, 8) /    .0003217541 9e0 /
      data fdbcs(  6, 8) /   -.0000071412 9e0 /
      data fdbcs(  7, 8) /   -.0000009676 6e0 /
      data fdbcs(  8, 8) /    .0000001160 1e0 /
      data fdbcs(  9, 8) /    .0000000450 4e0 /
      data fdbcs( 10, 8) /   -.0000000103 6e0 /
      data fdbcs( 11, 8) /   -.0000000025 9e0 /
      data fdbcs( 12, 8) /    .0000000014 6e0 /
      data fdbcs( 13, 8) /   -.0000000000 04e0 /
      data fdbcs( 14, 8) /   -.0000000001 8e0 /
      data fdbcs( 15, 8) /    .0000000000 54e0 /
c
c series for f35s    on the interval  0.00000e-01 to  2.71828e+00
c                                        with weighted error   7.30e-11
c                                         log weighted error  10.14
c                               significant figures required  11.18
c                                    decimal places required  10.72
c
      data fdscs(  1, 9) /  22.1653046970e0 /
      data fdscs(  2, 9) /   -.5086526539 69e0 /
      data fdscs(  3, 9) /    .0358871327 23e0 /
      data fdscs(  4, 9) /   -.0038993959 73e0 /
      data fdscs(  5, 9) /    .0005339699 07e0 /
      data fdscs(  6, 9) /   -.0000845770 08e0 /
      data fdscs(  7, 9) /    .0000148157 47e0 /
      data fdscs(  8, 9) /   -.0000027951 08e0 /
      data fdscs(  9, 9) /    .0000005582 82e0 /
      data fdscs( 10, 9) /   -.0000001166 84e0 /
      data fdscs( 11, 9) /    .0000000253 06e0 /
      data fdscs( 12, 9) /   -.0000000056 60e0 /
      data fdscs( 13, 9) /    .0000000012 99e0 /
      data fdscs( 14, 9) /   -.0000000003 05e0 /
      data fdscs( 15, 9) /    .0000000000 730e0 /
c
c series for f35m    on the interval  1.00000e+00 to  4.00000e+00
c                                        with weighted error   3.73e-11
c                                         log weighted error  10.43
c                               significant figures required  11.56
c                                    decimal places required  11.18
c
      data fdmcs(  1, 9) /  12.6701567503 6e0 /
      data fdmcs(  2, 9) /  -9.4636920164 00e0 /
      data fdmcs(  3, 9) /   6.0480654287 69e0 /
      data fdmcs(  4, 9) /  -3.4530339209 92e0 /
      data fdmcs(  5, 9) /   1.8250457226 69e0 /
      data fdmcs(  6, 9) /   -.9112870648 186e0 /
      data fdmcs(  7, 9) /    .4354953493 280e0 /
      data fdmcs(  8, 9) /   -.2009635658 884e0 /
      data fdmcs(  9, 9) /    .0901210173 526e0 /
      data fdmcs( 10, 9) /   -.0394612435 160e0 /
      data fdmcs( 11, 9) /    .0169328410 948e0 /
      data fdmcs( 12, 9) /   -.0071407017 340e0 /
      data fdmcs( 13, 9) /    .0029661522 591e0 /
      data fdmcs( 14, 9) /   -.0012158829 523e0 /
      data fdmcs( 15, 9) /    .0004926051 670e0 /
      data fdmcs( 16, 9) /   -.0001975006 123e0 /
      data fdmcs( 17, 9) /    .0000784453 353e0 /
      data fdmcs( 18, 9) /   -.0000308953 181e0 /
      data fdmcs( 19, 9) /    .0000120749 876e0 /
      data fdmcs( 20, 9) /   -.0000046864 594e0 /
      data fdmcs( 21, 9) /    .0000018072 775e0 /
      data fdmcs( 22, 9) /   -.0000006928 714e0 /
      data fdmcs( 23, 9) /    .0000002641 967e0 /
      data fdmcs( 24, 9) /   -.0000001002 365e0 /
      data fdmcs( 25, 9) /    .0000000378 535e0 /
      data fdmcs( 26, 9) /   -.0000000142 333e0 /
      data fdmcs( 27, 9) /    .0000000053 303e0 /
      data fdmcs( 28, 9) /   -.0000000019 887e0 /
      data fdmcs( 29, 9) /    .0000000007 393e0 /
      data fdmcs( 30, 9) /   -.0000000002 739e0 /
      data fdmcs( 31, 9) /    .0000000001 011e0 /
      data fdmcs( 32, 9) /   -.0000000000 372e0 /
c
c series for f35b    on the interval  0.00000e-01 to  2.50000e-01
c                                        with weighted error   5.48e-11
c                                         log weighted error  10.26
c                               significant figures required   9.91
c                                    decimal places required  10.85
c
      data fdbcs(  1, 9) /    .7665702804 4e0 /
      data fdbcs(  2, 9) /    .2216671315 0e0 /
      data fdbcs(  3, 9) /    .0657600143 6e0 /
      data fdbcs(  4, 9) /    .0058625466 3e0 /
      data fdbcs(  5, 9) /    .0006969431 7e0 /
      data fdbcs(  6, 9) /   -.0000101109 1e0 /
      data fdbcs(  7, 9) /   -.0000000660 4e0 /
      data fdbcs(  8, 9) /    .0000002536 7e0 /
      data fdbcs(  9, 9) /   -.0000000021 4e0 /
      data fdbcs( 10, 9) /   -.0000000132 6e0 /
      data fdbcs( 11, 9) /    .0000000015 0e0 /
      data fdbcs( 12, 9) /    .0000000009 5e0 /
      data fdbcs( 13, 9) /   -.0000000003 4e0 /
      data fdbcs( 14, 9) /   -.0000000000 30e0 /
      data fdbcs( 15, 9) /    .0000000000 54e0 /
c
c series for f40s    on the interval  0.00000e-01 to  2.71828e+00
c                                        with weighted error   4.91e-11
c                                         log weighted error  10.31
c                               significant figures required  11.67
c                                    decimal places required  10.90
c
      data fdscs(  1,10) /  46.3350918683 9e0 /
      data fdscs(  2,10) /   -.7807026145 628e0 /
      data fdscs(  3,10) /    .0465789224 743e0 /
      data fdscs(  4,10) /   -.0045080435 976e0 /
      data fdscs(  5,10) /    .0005646381 627e0 /
      data fdscs(  6,10) /   -.0000831331 627e0 /
      data fdscs(  7,10) /    .0000136850 757e0 /
      data fdscs(  8,10) /   -.0000024454 133e0 /
      data fdscs(  9,10) /    .0000004654 205e0 /
      data fdscs( 10,10) /   -.0000000931 320e0 /
      data fdscs( 11,10) /    .0000000194 128e0 /
      data fdscs( 12,10) /   -.0000000041 863e0 /
      data fdscs( 13,10) /    .0000000009 290e0 /
      data fdscs( 14,10) /   -.0000000002 113e0 /
      data fdscs( 15,10) /    .0000000000 491e0 /
c
c series for f40m    on the interval  1.00000e+00 to  4.00000e+00
c                                        with weighted error   6.13e-11
c                                         log weighted error  10.21
c                               significant figures required  11.66
c                                    decimal places required  10.97
c
      data fdmcs(  1,10) /  24.0980879457 2e0 /
      data fdmcs(  2,10) / -19.2973238247 9e0 /
      data fdmcs(  3,10) /  13.0411335433 2e0 /
      data fdmcs(  4,10) /  -7.8442177530 69e0 /
      data fdmcs(  5,10) /   4.3484777309 21e0 /
      data fdmcs(  6,10) /  -2.2682065310 65e0 /
      data fdmcs(  7,10) /   1.1283901506 72e0 /
      data fdmcs(  8,10) /   -.5404257776 187e0 /
      data fdmcs(  9,10) /    .2508755478 873e0 /
      data fdmcs( 10,10) /   -.1134573260 212e0 /
      data fdmcs( 11,10) /    .0501830996 299e0 /
      data fdmcs( 12,10) /   -.0217756382 802e0 /
      data fdmcs( 13,10) /    .0092927921 972e0 /
      data fdmcs( 14,10) /   -.0039080375 664e0 /
      data fdmcs( 15,10) /    .0016223028 722e0 /
      data fdmcs( 16,10) /   -.0006656886 564e0 /
      data fdmcs( 17,10) /    .0002703262 001e0 /
      data fdmcs( 18,10) /   -.0001087475 762e0 /
      data fdmcs( 19,10) /    .0000433751 765e0 /
      data fdmcs( 20,10) /   -.0000171663 866e0 /
      data fdmcs( 21,10) /    .0000067455 409e0 /
      data fdmcs( 22,10) /   -.0000026333 256e0 /
      data fdmcs( 23,10) /    .0000010217 923e0 /
      data fdmcs( 24,10) /   -.0000003942 636e0 /
      data fdmcs( 25,10) /    .0000001513 391e0 /
      data fdmcs( 26,10) /   -.0000000578 112e0 /
      data fdmcs( 27,10) /    .0000000219 841e0 /
      data fdmcs( 28,10) /   -.0000000083 247e0 /
      data fdmcs( 29,10) /    .0000000031 398e0 /
      data fdmcs( 30,10) /   -.0000000011 798e0 /
      data fdmcs( 31,10) /    .0000000004 418e0 /
      data fdmcs( 32,10) /   -.0000000001 648e0 /
      data fdmcs( 33,10) /    .0000000000 613e0 /
c
c series for f40b    on the interval  0.00000e-01 to  2.50000e-01
c                                        with weighted error   1.68e-11
c                                         log weighted error  10.77
c                               significant figures required  10.47
c                                    decimal places required  11.36
c
      data fdbcs(  1,10) /    .8056894147 6e0 /
      data fdbcs(  2,10) /    .2834447403 7e0 /
      data fdbcs(  3,10) /    .0903522185 7e0 /
      data fdbcs(  4,10) /    .0111606016 1e0 /
      data fdbcs(  5,10) /    .0014164744 6e0 /
      data fdbcs(  6,10) /    .0000100892 5e0 /
      data fdbcs(  7,10) /    .0000022449 5e0 /
      data fdbcs(  8,10) /    .0000001741 4e0 /
      data fdbcs(  9,10) /   -.0000000486 2e0 /
      data fdbcs( 10,10) /   -.0000000054 1e0 /
      data fdbcs( 11,10) /    .0000000035 1e0 /
      data fdbcs( 12,10) /   -.0000000000 82e0 /
      data fdbcs( 13,10) /   -.0000000003 1e0 /
      data fdbcs( 14,10) /    .0000000000 81e0 /
      data fdbcs( 15,10) /    .0000000000 16e0 /
c
      data nscs / 21, 20, 19, 18, 17, 17, 16, 16, 15, 15 /
      data nmcs / 21, 22, 24, 25, 26, 28, 29, 31, 32, 33 /
      data nbcs / 26, 22, 20, 18, 18, 17, 16, 15, 15, 15 /
c
      data exp1, alfsml, alfbig, alfmax / 13*0.0 /
c
      if (exp1.ne.0.0) go to 20
      exp1 = 2.0/exp(1.0)
c
      alfsml = dlog (r1mach(1)/0.8863)
      alfbig = dexp (dmin1 (-dlog(r1mach(1)),dlog(r1mach(2)))-log(2.0))
c
      alfmax(1) = r1mach(2)
      alfmax(2) = r1mach(2)
      do ndx=3,10
        xk = (ndx-2)*0.5
        alfmax(ndx) = r1mach(2)**(1.0/(xk+1.0)) * 0.99
      enddo
c
 20   ndx = (xnu+1.0)*2.0 + 0.01
      fd = 0.0
      if (alpha.lt.alfsml) return
c
      if (ndx.ne.2) go to 30
c
c calculate the fermi-dirac function for xnu = 0
c
      if (alpha.lt.10.0) fd = alnrel (exp(alpha))
      if (alpha.ge.10.0) fd = alpha + alnrel (exp(-alpha))
      return
c
 30   if (alpha.ge.1.0) go to 40
      expalf = exp (alpha)
      fd = expalf * csevl (expalf*exp1-1.0, fdscs(1,ndx), nscs(ndx))
      return
c
 40   if (alpha.ge.4.0) go to 50
      fd = alpha**(xnu+1.0) * csevl ((alpha-2.5)/1.5, fdmcs(1,ndx),
     1  nmcs(ndx))
      return
c
 50   alfi = 0.0
      if (alpha.lt.alfbig) alfi = 1.0/alpha
      fd = alpha**(xnu+1.0) * csevl ((alfi-0.125)*8.0, fdbcs(1,ndx),
     1  nbcs(ndx))
      return
c
      end
 

      function alnrel(x)
      implicit real*8 (a-h,o-z)
c
c ****  description
c
c     alnrel(x) evaluates ln(1+x) accurately in the sense of relative
c     error when x is very small.  this routine must be used to
c     maintain relative error accuracy whenever x is small and
c     accurately known.
c
c series for alnr       on the interval -3.75000d-01 to  3.75000d-01
c                                        with weighted error   1.93e-17
c                                         log weighted error  16.72
c                               significant figures required  16.44
c                                    decimal places required  17.40


      dimension alnrcs(23)
      data alnrcs( 1) /   1.0378693562 743770e0 /
      data alnrcs( 2) /   -.1336430150 4908918e0 /
      data alnrcs( 3) /    .0194082491 35520563e0 /
      data alnrcs( 4) /   -.0030107551 12753577e0 /
      data alnrcs( 5) /    .0004869461 47971548e0 /
      data alnrcs( 6) /   -.0000810548 81893175e0 /
      data alnrcs( 7) /    .0000137788 47799559e0 /
      data alnrcs( 8) /   -.0000023802 21089435e0 /
      data alnrcs( 9) /    .0000004164 04162138e0 /
      data alnrcs(10) /   -.0000000735 95828378e0 /
      data alnrcs(11) /    .0000000131 17611876e0 /
      data alnrcs(12) /   -.0000000023 54670931e0 /
      data alnrcs(13) /    .0000000004 25227732e0 /
      data alnrcs(14) /   -.0000000000 77190894e0 /
      data alnrcs(15) /    .0000000000 14075746e0 /
      data alnrcs(16) /   -.0000000000 02576907e0 /
      data alnrcs(17) /    .0000000000 00473424e0 /
      data alnrcs(18) /   -.0000000000 00087249e0 /
      data alnrcs(19) /    .0000000000 00016124e0 /
      data alnrcs(20) /   -.0000000000 00002987e0 /
      data alnrcs(21) /    .0000000000 00000554e0 /
      data alnrcs(22) /   -.0000000000 00000103e0 /
      data alnrcs(23) /    .0000000000 00000019e0 /
      data nlnrel, xmin /0, 0./
c
c **** first executable statement  alnrel
c
      if (nlnrel.ne.0) go to 10
      nlnrel = inits (alnrcs, 23, 0.1*r1mach(3))
      xmin = -1.0 + sqrt(r1mach(4))
c
 10   continue
 
c
      if (abs(x).le.0.375) alnrel = x*(1. -
     1  x*csevl (x/.375, alnrcs, nlnrel))
      if (abs(x).gt.0.375) alnrel = dlog (1.0+x)
c
      return
      end

      function inits(os,nos,eta)
      implicit real*8 (a-h,o-z)
      dimension os(nos)
c 
c .... initialize,orthogonal series,special function
c      initializes an orthogonal series so that it defines the
c     number of terms to carry in the series to meet a specified
c     error.
c
c     input arguments --
c     os     array of nos coefficients in an orthogonal series.
c     nos    number of coefficients in os.
c     eta    requested accuracy of series.


      err = 0.
      do ii=1,nos
        i = nos + 1 - ii
        err = err + abs(os(i))
        if (err.gt.eta) go to 20
      enddo

 20   continue
      inits = i
c
      return
      end


      function csevl(x,cs,n)
      implicit real*8 (a-h,o-z)
c
c ....description:
c     evaluate the n-term chebyshev series cs at x.  adapted from
c     r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973). also see fox
c     and parker, chebyshev polynomials in numerical analysis, oxford press
c     page 56.
c
c     input arguments --
c     x    value at which the series is to be evaluated.
c     cs   array of n terms of a chebyshev series.  in eval-
c          uating cs, only half the first coefficient is summed.
c     n    number of terms in array cs.
c     
       dimension cs(1)

       b1=0.
       b0=0.
       twox=2.*x
       do i=1,n
          b2=b1
          b1=b0
          ni=n+1-i
          b0=twox*b1-b2+cs(ni)
       enddo
c
       csevl = 0.5 * (b0-b2)
c
       return
      end



      function r1mach(i)
      implicit real*8 (a-h,o-z)
c      real rmach(5)

c  single-precision machine constants
c
c  machine constants for the cray 1, cray x-mp
c  smallest useable real number
c      data rmach(1) / 200034000000000000000b /
c  largest useable real number
c      data rmach(2) / 577767777777777777776b /
c  smallest eps such that 1.+eps .ne. 1.
c      data rmach(3) / 377224000000000000000b /
c  2.*rmach(3)
c      data rmach(4) / 377234000000000000000b /
c  log10(2)
c      data rmach(5) / 377774642023241175720b /


c
c for  HP 730
c

c      DATA RMACH(1) / Z'00800000' /
c      DATA RMACH(2) / Z'7F7FFFFF' /
c      DATA RMACH(3) / Z'33800000' /
c      DATA RMACH(4) / Z'34000000' /
c      DATA RMACH(5) / Z'3E9A209B' /
 

c      r1mach = rmach(i)

      r1mach = D1MACH(i)
      return
      end

      
C     ALGORITHM 745, COLLECTED ALGORITHMS FROM ACM.
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
        DOUBLE PRECISION BETA, EPS, XBIG, XMIN, XMAX, GAMMA
        PARAMETER (BETA = TWO, EPS = BETA**MACHEP, XMIN = BETA**MINEXP,
     &             XMAX = (BETA**(MAXEXP-1) -
     &       BETA**(MAXEXP+NEGEXP-1))*BETA)

        
        if (abs(ORD*2-int(ORD*2+0.1)) .lt. 1e-5 .and.
     &       ORD .ge. -0.5 .and. ORD .le. 4.0) then
           IERR = 0
           call fdm(ORD, X, fd)           
           CALL GAMMAC(ORD+1, EPS, XMAX, GAMMA, IERR)
           fd = fd/GAMMA
           IERR = 0
           return
        endif
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
