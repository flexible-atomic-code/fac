      subroutine nrrfit(iz, in, t, r)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C RR rate coefficients of M. F. Gu, ApJ ...
C iz, Nuclear charge of the recombining ion.
C in, number of electrons. 0 to iz-1
C t,  temperature in eV.
C r,  results in 10^-10 cm3 s-1.
C     if r < 0, then data not available.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      integer iz, in
      double precision t, r
      common /nrr/ ind(28), rrpar(6, 70)
      
      r = -1.0
      if (ind(iz) .lt. 0) then
         return
      endif
      if (in .ge. 10) then 
         return
      endif

      i = ind(iz) + in
      
      a = rrpar(1, i)
      b = rrpar(2, i)
      t0 = rrpar(3, i)
      t1 = rrpar(4, i)
      b1 = rrpar(5, i)
      t2 = rrpar(6, i)

      if (t2 .gt. 0) then
         b = b + b1*exp(-t2/t)
      endif

      x0 = sqrt(t/t0)
      x1 = sqrt(t/t1)
      r = x0 * (1+x0)**(1.0-b) * (1+x1)**(1.0+b)
      r = a/r
      return
      end

      
      block data nrrdata
      common /nrr/ ind(28), rrpar(6, 70)
      
      data ind/ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
     +          -1,  1, -1, 11, -1, 21, -1, 31, -1, 41,
     +          -1, -1, -1, -1, -1, 51, -1, 61 /
      data rrpar(1, 1)/ 3.959E+01/
      data rrpar(2, 1)/ 7.778E-01/
      data rrpar(3, 1)/ 3.243E-03/
      data rrpar(4, 1)/ 7.437E+03/
      data rrpar(5, 1)/ 0.000E+00/
      data rrpar(6, 1)/ 0.000E+00/
      data rrpar(1, 2)/ 2.616E+01/
      data rrpar(2, 2)/ 7.478E-01/
      data rrpar(3, 2)/ 4.678E-03/
      data rrpar(4, 2)/ 5.144E+03/
      data rrpar(5, 2)/ 0.000E+00/
      data rrpar(6, 2)/ 0.000E+00/
      data rrpar(1, 3)/ 2.082E+01/
      data rrpar(2, 3)/ 7.213E-01/
      data rrpar(3, 3)/ 4.392E-03/
      data rrpar(4, 3)/ 2.008E+03/
      data rrpar(5, 3)/ 0.000E+00/
      data rrpar(6, 3)/ 0.000E+00/
      data rrpar(1, 4)/ 2.148E+01/
      data rrpar(2, 4)/ 7.415E-01/
      data rrpar(3, 4)/ 2.805E-03/
      data rrpar(4, 4)/ 1.503E+03/
      data rrpar(5, 4)/ 0.000E+00/
      data rrpar(6, 4)/ 0.000E+00/
      data rrpar(1, 5)/ 2.038E+01/
      data rrpar(2, 5)/ 7.524E-01/
      data rrpar(3, 5)/ 2.026E-03/
      data rrpar(4, 5)/ 1.207E+03/
      data rrpar(5, 5)/ 0.000E+00/
      data rrpar(6, 5)/ 0.000E+00/
      data rrpar(1, 6)/ 1.705E+01/
      data rrpar(2, 6)/ 7.543E-01/
      data rrpar(3, 6)/ 1.681E-03/
      data rrpar(4, 6)/ 1.267E+03/
      data rrpar(5, 6)/ 0.000E+00/
      data rrpar(6, 6)/ 0.000E+00/
      data rrpar(1, 7)/ 1.165E+01/
      data rrpar(2, 7)/ 7.390E-01/
      data rrpar(3, 7)/ 1.946E-03/
      data rrpar(4, 7)/ 1.513E+03/
      data rrpar(5, 7)/ 0.000E+00/
      data rrpar(6, 7)/ 0.000E+00/
      data rrpar(1, 8)/ 5.314E+00/
      data rrpar(2, 8)/ 6.838E-01/
      data rrpar(3, 8)/ 4.443E-03/
      data rrpar(4, 8)/ 1.256E+03/
      data rrpar(5, 8)/ 5.602E-02/
      data rrpar(6, 8)/ 6.008E+01/
      data rrpar(1, 9)/ 4.878E+00/
      data rrpar(2, 9)/ 7.060E-01/
      data rrpar(3, 9)/ 1.768E-03/
      data rrpar(4, 9)/ 1.024E+03/
      data rrpar(5, 9)/ 8.149E-02/
      data rrpar(6, 9)/ 5.967E+01/
      data rrpar(1,10)/ 2.001E+00/
      data rrpar(2,10)/ 6.339E-01/
      data rrpar(3,10)/ 2.693E-03/
      data rrpar(4,10)/ 1.071E+03/
      data rrpar(5,10)/ 1.464E-01/
      data rrpar(6,10)/ 4.992E+01/
      data rrpar(1,11)/ 5.130E+01/
      data rrpar(2,11)/ 7.802E-01/
      data rrpar(3,11)/ 3.665E-03/
      data rrpar(4,11)/ 9.832E+03/
      data rrpar(5,11)/ 0.000E+00/
      data rrpar(6,11)/ 0.000E+00/
      data rrpar(1,12)/ 3.625E+01/
      data rrpar(2,12)/ 7.534E-01/
      data rrpar(3,12)/ 4.901E-03/
      data rrpar(4,12)/ 6.471E+03/
      data rrpar(5,12)/ 0.000E+00/
      data rrpar(6,12)/ 0.000E+00/
      data rrpar(1,13)/ 3.048E+01/
      data rrpar(2,13)/ 7.316E-01/
      data rrpar(3,13)/ 4.435E-03/
      data rrpar(4,13)/ 2.474E+03/
      data rrpar(5,13)/ 0.000E+00/
      data rrpar(6,13)/ 0.000E+00/
      data rrpar(1,14)/ 3.118E+01/
      data rrpar(2,14)/ 7.471E-01/
      data rrpar(3,14)/ 3.073E-03/
      data rrpar(4,14)/ 1.906E+03/
      data rrpar(5,14)/ 0.000E+00/
      data rrpar(6,14)/ 0.000E+00/
      data rrpar(1,15)/ 3.016E+01/
      data rrpar(2,15)/ 7.563E-01/
      data rrpar(3,15)/ 2.309E-03/
      data rrpar(4,15)/ 1.531E+03/
      data rrpar(5,15)/ 0.000E+00/
      data rrpar(6,15)/ 0.000E+00/
      data rrpar(1,16)/ 2.556E+01/
      data rrpar(2,16)/ 7.559E-01/
      data rrpar(3,16)/ 2.080E-03/
      data rrpar(4,16)/ 1.566E+03/
      data rrpar(5,16)/ 0.000E+00/
      data rrpar(6,16)/ 0.000E+00/
      data rrpar(1,17)/ 1.856E+01/
      data rrpar(2,17)/ 7.434E-01/
      data rrpar(3,17)/ 2.426E-03/
      data rrpar(4,17)/ 1.731E+03/
      data rrpar(5,17)/ 0.000E+00/
      data rrpar(6,17)/ 0.000E+00/
      data rrpar(1,18)/ 1.217E+01/
      data rrpar(2,18)/ 7.186E-01/
      data rrpar(3,18)/ 3.207E-03/
      data rrpar(4,18)/ 2.131E+03/
      data rrpar(5,18)/ 0.000E+00/
      data rrpar(6,18)/ 0.000E+00/
      data rrpar(1,19)/ 8.021E+00/
      data rrpar(2,19)/ 6.894E-01/
      data rrpar(3,19)/ 3.847E-03/
      data rrpar(4,19)/ 2.899E+03/
      data rrpar(5,19)/ 0.000E+00/
      data rrpar(6,19)/ 0.000E+00/
      data rrpar(1,20)/ 4.156E+00/
      data rrpar(2,20)/ 6.323E-01/
      data rrpar(3,20)/ 6.301E-03/
      data rrpar(4,20)/ 2.072E+03/
      data rrpar(5,20)/ 6.492E-02/
      data rrpar(6,20)/ 1.063E+02/
      data rrpar(1,21)/ 6.429E+01/
      data rrpar(2,21)/ 7.824E-01/
      data rrpar(3,21)/ 4.064E-03/
      data rrpar(4,21)/ 1.245E+04/
      data rrpar(5,21)/ 0.000E+00/
      data rrpar(6,21)/ 0.000E+00/
      data rrpar(1,22)/ 4.753E+01/
      data rrpar(2,22)/ 7.579E-01/
      data rrpar(3,22)/ 5.178E-03/
      data rrpar(4,22)/ 7.902E+03/
      data rrpar(5,22)/ 0.000E+00/
      data rrpar(6,22)/ 0.000E+00/
      data rrpar(1,23)/ 4.173E+01/
      data rrpar(2,23)/ 7.393E-01/
      data rrpar(3,23)/ 4.532E-03/
      data rrpar(4,23)/ 3.000E+03/
      data rrpar(5,23)/ 0.000E+00/
      data rrpar(6,23)/ 0.000E+00/
      data rrpar(1,24)/ 4.266E+01/
      data rrpar(2,24)/ 7.520E-01/
      data rrpar(3,24)/ 3.292E-03/
      data rrpar(4,24)/ 2.351E+03/
      data rrpar(5,24)/ 0.000E+00/
      data rrpar(6,24)/ 0.000E+00/
      data rrpar(1,25)/ 4.177E+01/
      data rrpar(2,25)/ 7.599E-01/
      data rrpar(3,25)/ 2.544E-03/
      data rrpar(4,25)/ 1.898E+03/
      data rrpar(5,25)/ 0.000E+00/
      data rrpar(6,25)/ 0.000E+00/
      data rrpar(1,26)/ 3.559E+01/
      data rrpar(2,26)/ 7.577E-01/
      data rrpar(3,26)/ 2.438E-03/
      data rrpar(4,26)/ 1.912E+03/
      data rrpar(5,26)/ 0.000E+00/
      data rrpar(6,26)/ 0.000E+00/
      data rrpar(1,27)/ 2.695E+01/
      data rrpar(2,27)/ 7.464E-01/
      data rrpar(3,27)/ 2.847E-03/
      data rrpar(4,27)/ 2.040E+03/
      data rrpar(5,27)/ 0.000E+00/
      data rrpar(6,27)/ 0.000E+00/
      data rrpar(1,28)/ 1.729E+01/
      data rrpar(2,28)/ 7.192E-01/
      data rrpar(3,28)/ 4.337E-03/
      data rrpar(4,28)/ 2.626E+03/
      data rrpar(5,28)/ 0.000E+00/
      data rrpar(6,28)/ 0.000E+00/
      data rrpar(1,29)/ 1.296E+01/
      data rrpar(2,29)/ 7.025E-01/
      data rrpar(3,29)/ 4.654E-03/
      data rrpar(4,29)/ 2.734E+03/
      data rrpar(5,29)/ 0.000E+00/
      data rrpar(6,29)/ 0.000E+00/
      data rrpar(1,30)/ 8.221E+00/
      data rrpar(2,30)/ 6.645E-01/
      data rrpar(3,30)/ 6.275E-03/
      data rrpar(4,30)/ 3.606E+03/
      data rrpar(5,30)/ 0.000E+00/
      data rrpar(6,30)/ 0.000E+00/
      data rrpar(1,31)/ 7.835E+01/
      data rrpar(2,31)/ 7.844E-01/
      data rrpar(3,31)/ 4.459E-03/
      data rrpar(4,31)/ 1.528E+04/
      data rrpar(5,31)/ 0.000E+00/
      data rrpar(6,31)/ 0.000E+00/
      data rrpar(1,32)/ 6.018E+01/
      data rrpar(2,32)/ 7.618E-01/
      data rrpar(3,32)/ 5.441E-03/
      data rrpar(4,32)/ 9.417E+03/
      data rrpar(5,32)/ 0.000E+00/
      data rrpar(6,32)/ 0.000E+00/
      data rrpar(1,33)/ 5.422E+01/
      data rrpar(2,33)/ 7.452E-01/
      data rrpar(3,33)/ 4.702E-03/
      data rrpar(4,33)/ 3.580E+03/
      data rrpar(5,33)/ 0.000E+00/
      data rrpar(6,33)/ 0.000E+00/
      data rrpar(1,34)/ 5.522E+01/
      data rrpar(2,34)/ 7.557E-01/
      data rrpar(3,34)/ 3.560E-03/
      data rrpar(4,34)/ 2.844E+03/
      data rrpar(5,34)/ 0.000E+00/
      data rrpar(6,34)/ 0.000E+00/
      data rrpar(1,35)/ 5.449E+01/
      data rrpar(2,35)/ 7.626E-01/
      data rrpar(3,35)/ 2.815E-03/
      data rrpar(4,35)/ 2.311E+03/
      data rrpar(5,35)/ 0.000E+00/
      data rrpar(6,35)/ 0.000E+00/
      data rrpar(1,36)/ 4.726E+01/
      data rrpar(2,36)/ 7.598E-01/
      data rrpar(3,36)/ 2.745E-03/
      data rrpar(4,36)/ 2.294E+03/
      data rrpar(5,36)/ 0.000E+00/
      data rrpar(6,36)/ 0.000E+00/
      data rrpar(1,37)/ 3.716E+01/
      data rrpar(2,37)/ 7.498E-01/
      data rrpar(3,37)/ 3.158E-03/
      data rrpar(4,37)/ 2.386E+03/
      data rrpar(5,37)/ 0.000E+00/
      data rrpar(6,37)/ 0.000E+00/
      data rrpar(1,38)/ 2.738E+01/
      data rrpar(2,38)/ 7.331E-01/
      data rrpar(3,38)/ 3.972E-03/
      data rrpar(4,38)/ 2.569E+03/
      data rrpar(5,38)/ 0.000E+00/
      data rrpar(6,38)/ 0.000E+00/
      data rrpar(1,39)/ 1.976E+01/
      data rrpar(2,39)/ 7.129E-01/
      data rrpar(3,39)/ 4.989E-03/
      data rrpar(4,39)/ 2.842E+03/
      data rrpar(5,39)/ 0.000E+00/
      data rrpar(6,39)/ 0.000E+00/
      data rrpar(1,40)/ 1.317E+01/
      data rrpar(2,40)/ 6.813E-01/
      data rrpar(3,40)/ 6.877E-03/
      data rrpar(4,40)/ 3.339E+03/
      data rrpar(5,40)/ 0.000E+00/
      data rrpar(6,40)/ 0.000E+00/
      data rrpar(1,41)/ 9.351E+01/
      data rrpar(2,41)/ 7.863E-01/
      data rrpar(3,41)/ 4.843E-03/
      data rrpar(4,41)/ 1.828E+04/
      data rrpar(5,41)/ 0.000E+00/
      data rrpar(6,41)/ 0.000E+00/
      data rrpar(1,42)/ 7.368E+01/
      data rrpar(2,42)/ 7.651E-01/
      data rrpar(3,42)/ 5.760E-03/
      data rrpar(4,42)/ 1.103E+04/
      data rrpar(5,42)/ 0.000E+00/
      data rrpar(6,42)/ 0.000E+00/
      data rrpar(1,43)/ 6.783E+01/
      data rrpar(2,43)/ 7.499E-01/
      data rrpar(3,43)/ 4.917E-03/
      data rrpar(4,43)/ 4.212E+03/
      data rrpar(5,43)/ 0.000E+00/
      data rrpar(6,43)/ 0.000E+00/
      data rrpar(1,44)/ 6.919E+01/
      data rrpar(2,44)/ 7.590E-01/
      data rrpar(3,44)/ 3.811E-03/
      data rrpar(4,44)/ 3.376E+03/
      data rrpar(5,44)/ 0.000E+00/
      data rrpar(6,44)/ 0.000E+00/
      data rrpar(1,45)/ 6.844E+01/
      data rrpar(2,45)/ 7.648E-01/
      data rrpar(3,45)/ 3.089E-03/
      data rrpar(4,45)/ 2.766E+03/
      data rrpar(5,45)/ 0.000E+00/
      data rrpar(6,45)/ 0.000E+00/
      data rrpar(1,46)/ 5.986E+01/
      data rrpar(2,46)/ 7.615E-01/
      data rrpar(3,46)/ 3.075E-03/
      data rrpar(4,46)/ 2.726E+03/
      data rrpar(5,46)/ 0.000E+00/
      data rrpar(6,46)/ 0.000E+00/
      data rrpar(1,47)/ 4.857E+01/
      data rrpar(2,47)/ 7.526E-01/
      data rrpar(3,47)/ 3.477E-03/
      data rrpar(4,47)/ 2.772E+03/
      data rrpar(5,47)/ 0.000E+00/
      data rrpar(6,47)/ 0.000E+00/
      data rrpar(1,48)/ 3.743E+01/
      data rrpar(2,48)/ 7.384E-01/
      data rrpar(3,48)/ 4.223E-03/
      data rrpar(4,48)/ 2.889E+03/
      data rrpar(5,48)/ 0.000E+00/
      data rrpar(6,48)/ 0.000E+00/
      data rrpar(1,49)/ 2.828E+01/
      data rrpar(2,49)/ 7.215E-01/
      data rrpar(3,49)/ 5.165E-03/
      data rrpar(4,49)/ 3.052E+03/
      data rrpar(5,49)/ 0.000E+00/
      data rrpar(6,49)/ 0.000E+00/
      data rrpar(1,50)/ 2.003E+01/
      data rrpar(2,50)/ 6.960E-01/
      data rrpar(3,50)/ 6.849E-03/
      data rrpar(4,50)/ 3.323E+03/
      data rrpar(5,50)/ 0.000E+00/
      data rrpar(6,50)/ 0.000E+00/
      data rrpar(1,51)/ 1.438E+02/
      data rrpar(2,51)/ 7.911E-01/
      data rrpar(3,51)/ 6.050E-03/
      data rrpar(4,51)/ 2.816E+04/
      data rrpar(5,51)/ 0.000E+00/
      data rrpar(6,51)/ 0.000E+00/
      data rrpar(1,52)/ 1.206E+02/
      data rrpar(2,52)/ 7.730E-01/
      data rrpar(3,52)/ 6.701E-03/
      data rrpar(4,52)/ 1.627E+04/
      data rrpar(5,52)/ 0.000E+00/
      data rrpar(6,52)/ 0.000E+00/
      data rrpar(1,53)/ 1.140E+02/
      data rrpar(2,53)/ 7.598E-01/
      data rrpar(3,53)/ 5.760E-03/
      data rrpar(4,53)/ 6.428E+03/
      data rrpar(5,53)/ 0.000E+00/
      data rrpar(6,53)/ 0.000E+00/
      data rrpar(1,54)/ 1.153E+02/
      data rrpar(2,54)/ 7.656E-01/
      data rrpar(3,54)/ 4.784E-03/
      data rrpar(4,54)/ 5.286E+03/
      data rrpar(5,54)/ 0.000E+00/
      data rrpar(6,54)/ 0.000E+00/
      data rrpar(1,55)/ 1.142E+02/
      data rrpar(2,55)/ 7.692E-01/
      data rrpar(3,55)/ 4.096E-03/
      data rrpar(4,55)/ 4.400E+03/
      data rrpar(5,55)/ 0.000E+00/
      data rrpar(6,55)/ 0.000E+00/
      data rrpar(1,56)/ 1.031E+02/
      data rrpar(2,56)/ 7.658E-01/
      data rrpar(3,56)/ 4.116E-03/
      data rrpar(4,56)/ 4.254E+03/
      data rrpar(5,56)/ 0.000E+00/
      data rrpar(6,56)/ 0.000E+00/
      data rrpar(1,57)/ 8.931E+01/
      data rrpar(2,57)/ 7.592E-01/
      data rrpar(3,57)/ 4.440E-03/
      data rrpar(4,57)/ 4.153E+03/
      data rrpar(5,57)/ 0.000E+00/
      data rrpar(6,57)/ 0.000E+00/
      data rrpar(1,58)/ 7.534E+01/
      data rrpar(2,58)/ 7.500E-01/
      data rrpar(3,58)/ 4.966E-03/
      data rrpar(4,58)/ 4.083E+03/
      data rrpar(5,58)/ 0.000E+00/
      data rrpar(6,58)/ 0.000E+00/
      data rrpar(1,59)/ 6.252E+01/
      data rrpar(2,59)/ 7.389E-01/
      data rrpar(3,59)/ 5.639E-03/
      data rrpar(4,59)/ 4.014E+03/
      data rrpar(5,59)/ 0.000E+00/
      data rrpar(6,59)/ 0.000E+00/
      data rrpar(1,60)/ 5.039E+01/
      data rrpar(2,60)/ 7.244E-01/
      data rrpar(3,60)/ 6.646E-03/
      data rrpar(4,60)/ 3.899E+03/
      data rrpar(5,60)/ 0.000E+00/
      data rrpar(6,60)/ 0.000E+00/
      data rrpar(1,61)/ 1.622E+02/
      data rrpar(2,61)/ 7.924E-01/
      data rrpar(3,61)/ 6.448E-03/
      data rrpar(4,61)/ 3.170E+04/
      data rrpar(5,61)/ 0.000E+00/
      data rrpar(6,61)/ 0.000E+00/
      data rrpar(1,62)/ 1.380E+02/
      data rrpar(2,62)/ 7.752E-01/
      data rrpar(3,62)/ 7.031E-03/
      data rrpar(4,62)/ 1.814E+04/
      data rrpar(5,62)/ 0.000E+00/
      data rrpar(6,62)/ 0.000E+00/
      data rrpar(1,63)/ 1.307E+02/
      data rrpar(2,63)/ 7.622E-01/
      data rrpar(3,63)/ 6.108E-03/
      data rrpar(4,63)/ 7.270E+03/
      data rrpar(5,63)/ 0.000E+00/
      data rrpar(6,63)/ 0.000E+00/
      data rrpar(1,64)/ 1.320E+02/
      data rrpar(2,64)/ 7.672E-01/
      data rrpar(3,64)/ 5.143E-03/
      data rrpar(4,64)/ 6.018E+03/
      data rrpar(5,64)/ 0.000E+00/
      data rrpar(6,64)/ 0.000E+00/
      data rrpar(1,65)/ 1.306E+02/
      data rrpar(2,65)/ 7.703E-01/
      data rrpar(3,65)/ 4.471E-03/
      data rrpar(4,65)/ 5.036E+03/
      data rrpar(5,65)/ 0.000E+00/
      data rrpar(6,65)/ 0.000E+00/
      data rrpar(1,66)/ 1.187E+02/
      data rrpar(2,66)/ 7.668E-01/
      data rrpar(3,66)/ 4.508E-03/
      data rrpar(4,66)/ 4.852E+03/
      data rrpar(5,66)/ 0.000E+00/
      data rrpar(6,66)/ 0.000E+00/
      data rrpar(1,67)/ 1.045E+02/
      data rrpar(2,67)/ 7.608E-01/
      data rrpar(3,67)/ 4.796E-03/
      data rrpar(4,67)/ 4.683E+03/
      data rrpar(5,67)/ 0.000E+00/
      data rrpar(6,67)/ 0.000E+00/
      data rrpar(1,68)/ 8.981E+01/
      data rrpar(2,68)/ 7.526E-01/
      data rrpar(3,68)/ 5.271E-03/
      data rrpar(4,68)/ 4.554E+03/
      data rrpar(5,68)/ 0.000E+00/
      data rrpar(6,68)/ 0.000E+00/
      data rrpar(1,69)/ 7.616E+01/
      data rrpar(2,69)/ 7.428E-01/
      data rrpar(3,69)/ 5.871E-03/
      data rrpar(4,69)/ 4.409E+03/
      data rrpar(5,69)/ 0.000E+00/
      data rrpar(6,69)/ 0.000E+00/
      data rrpar(1,70)/ 6.346E+01/
      data rrpar(2,70)/ 7.306E-01/
      data rrpar(3,70)/ 6.655E-03/
      data rrpar(4,70)/ 4.210E+03/
      data rrpar(5,70)/ 0.000E+00/
      data rrpar(6,70)/ 0.000E+00/

      end
