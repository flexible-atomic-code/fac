      subroutine ndrfit(iz, in, t, r)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C DR rate coefficients of M. F. Gu, ApJ ...
C iz, Nuclear charge of the recombining ion.
C in, number of electrons. 0 to iz-1
C t,  temperature in eV.
C r,  results in 10^-10 cm3 s-1.
C     if r < 0, then data not available.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      integer iz, in
      double precision t, r
      common /ndr/ ind(28), dre(16, 70), drc(16, 70)
      
      r = -1.0
      if (ind(iz) .lt. 0) then 
         return
      endif
      if (in .eq. 0 .or. in .gt. 10) then 
         return
      endif
      
      i = ind(iz) + in - 1
      
      r = 0.0
      do j = 1, 16
         if (dre(j, i) .gt. 0) then 
            r = r + drc(j, i)*exp(-dre(j, i)/t)
         endif
      enddo
       
      r = r/(t**(1.5))
      end


      block data ndrdata
      common /ndr/ ind(28), dre(16, 70), drc(16, 70)
      
      data ind/ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
     +          -1,  1, -1, 11, -1, 21, -1, 31, -1, 41,
     +          -1, -1, -1, -1, -1, 51, -1, 61 /
      data drc( 1, 1)/ 0.000E+00/
      data drc( 2, 1)/ 0.000E+00/
      data drc( 3, 1)/ 0.000E+00/
      data drc( 4, 1)/ 0.000E+00/
      data drc( 5, 1)/ 0.000E+00/
      data drc( 6, 1)/ 0.000E+00/
      data drc( 7, 1)/ 0.000E+00/
      data drc( 8, 1)/ 0.000E+00/
      data drc( 9, 1)/ 1.564E+02/
      data drc(10, 1)/ 7.123E+02/
      data drc(11, 1)/ 2.186E+02/
      data drc(12, 1)/ 0.000E+00/
      data drc(13, 1)/ 0.000E+00/
      data drc(14, 1)/ 0.000E+00/
      data drc(15, 1)/ 0.000E+00/
      data drc(16, 1)/ 0.000E+00/
      data dre( 1, 1)/ 0.000E+00/
      data dre( 2, 1)/ 0.000E+00/
      data dre( 3, 1)/ 0.000E+00/
      data dre( 4, 1)/ 0.000E+00/
      data dre( 5, 1)/ 0.000E+00/
      data dre( 6, 1)/ 0.000E+00/
      data dre( 7, 1)/ 0.000E+00/
      data dre( 8, 1)/ 0.000E+00/
      data dre( 9, 1)/ 1.038E+03/
      data dre(10, 1)/ 1.350E+03/
      data dre(11, 1)/ 1.551E+03/
      data dre(12, 1)/ 0.000E+00/
      data dre(13, 1)/ 0.000E+00/
      data dre(14, 1)/ 0.000E+00/
      data dre(15, 1)/ 0.000E+00/
      data dre(16, 1)/ 0.000E+00/
      data drc( 1, 2)/ 0.000E+00/
      data drc( 2, 2)/ 0.000E+00/
      data drc( 3, 2)/ 0.000E+00/
      data drc( 4, 2)/ 0.000E+00/
      data drc( 5, 2)/ 0.000E+00/
      data drc( 6, 2)/ 0.000E+00/
      data drc( 7, 2)/ 0.000E+00/
      data drc( 8, 2)/ 0.000E+00/
      data drc( 9, 2)/ 2.062E+02/
      data drc(10, 2)/ 8.684E+02/
      data drc(11, 2)/ 5.554E+02/
      data drc(12, 2)/ 0.000E+00/
      data drc(13, 2)/ 0.000E+00/
      data drc(14, 2)/ 0.000E+00/
      data drc(15, 2)/ 0.000E+00/
      data drc(16, 2)/ 0.000E+00/
      data dre( 1, 2)/ 0.000E+00/
      data dre( 2, 2)/ 0.000E+00/
      data dre( 3, 2)/ 0.000E+00/
      data dre( 4, 2)/ 0.000E+00/
      data dre( 5, 2)/ 0.000E+00/
      data dre( 6, 2)/ 0.000E+00/
      data dre( 7, 2)/ 0.000E+00/
      data dre( 8, 2)/ 0.000E+00/
      data dre( 9, 2)/ 9.749E+02/
      data dre(10, 2)/ 1.200E+03/
      data dre(11, 2)/ 1.353E+03/
      data dre(12, 2)/ 0.000E+00/
      data dre(13, 2)/ 0.000E+00/
      data dre(14, 2)/ 0.000E+00/
      data dre(15, 2)/ 0.000E+00/
      data dre(16, 2)/ 0.000E+00/
      data drc( 1, 3)/ 6.930E-02/
      data drc( 2, 3)/ 2.095E+00/
      data drc( 3, 3)/ 2.668E+00/
      data drc( 4, 3)/ 4.620E+00/
      data drc( 5, 3)/ 1.019E+01/
      data drc( 6, 3)/ 5.223E+01/
      data drc( 7, 3)/ 0.000E+00/
      data drc( 8, 3)/ 0.000E+00/
      data drc( 9, 3)/ 6.467E-07/
      data drc(10, 3)/ 4.092E+00/
      data drc(11, 3)/ 7.705E+00/
      data drc(12, 3)/ 4.117E+01/
      data drc(13, 3)/ 1.686E+02/
      data drc(14, 3)/ 2.178E+01/
      data drc(15, 3)/ 2.183E+02/
      data drc(16, 3)/ 4.811E+02/
      data dre( 1, 3)/ 2.064E+00/
      data dre( 2, 3)/ 2.610E+00/
      data dre( 3, 3)/ 3.157E+00/
      data dre( 4, 3)/ 6.982E+00/
      data dre( 5, 3)/ 1.270E+01/
      data dre( 6, 3)/ 1.929E+01/
      data dre( 7, 3)/ 0.000E+00/
      data dre( 8, 3)/ 0.000E+00/
      data dre( 9, 3)/ 2.835E+01/
      data dre(10, 3)/ 7.163E+01/
      data dre(11, 3)/ 8.198E+01/
      data dre(12, 3)/ 1.418E+02/
      data dre(13, 3)/ 1.901E+02/
      data dre(14, 3)/ 2.466E+02/
      data dre(15, 3)/ 1.024E+03/
      data dre(16, 3)/ 1.274E+03/
      data drc( 1, 4)/ 5.782E-01/
      data drc( 2, 4)/ 6.911E-01/
      data drc( 3, 4)/ 2.991E+00/
      data drc( 4, 4)/ 2.602E+01/
      data drc( 5, 4)/ 1.522E+02/
      data drc( 6, 4)/ 0.000E+00/
      data drc( 7, 4)/ 0.000E+00/
      data drc( 8, 4)/ 0.000E+00/
      data drc( 9, 4)/ 2.025E-10/
      data drc(10, 4)/ 1.183E+00/
      data drc(11, 4)/ 1.116E+01/
      data drc(12, 4)/ 2.688E+01/
      data drc(13, 4)/ 1.643E+02/
      data drc(14, 4)/ 3.895E+01/
      data drc(15, 4)/ 6.554E+01/
      data drc(16, 4)/ 1.470E+02/
      data dre( 1, 4)/ 3.280E-02/
      data dre( 2, 4)/ 8.783E-02/
      data dre( 3, 4)/ 3.196E+00/
      data dre( 4, 4)/ 8.997E+00/
      data dre( 5, 4)/ 3.035E+01/
      data dre( 6, 4)/ 0.000E+00/
      data dre( 7, 4)/ 0.000E+00/
      data dre( 8, 4)/ 0.000E+00/
      data dre( 9, 4)/ 3.703E+00/
      data dre(10, 4)/ 7.182E+01/
      data dre(11, 4)/ 8.723E+01/
      data dre(12, 4)/ 1.310E+02/
      data dre(13, 4)/ 1.765E+02/
      data dre(14, 4)/ 2.250E+02/
      data dre(15, 4)/ 1.039E+03/
      data dre(16, 4)/ 1.255E+03/
      data drc( 1, 5)/ 2.795E-01/
      data drc( 2, 5)/ 4.931E-01/
      data drc( 3, 5)/ 2.054E+00/
      data drc( 4, 5)/ 6.509E+00/
      data drc( 5, 5)/ 1.673E+02/
      data drc( 6, 5)/ 0.000E+00/
      data drc( 7, 5)/ 0.000E+00/
      data drc( 8, 5)/ 0.000E+00/
      data drc( 9, 5)/ 5.775E-01/
      data drc(10, 5)/ 1.834E+01/
      data drc(11, 5)/ 1.315E+02/
      data drc(12, 5)/ 2.483E+02/
      data drc(13, 5)/ 3.599E+00/
      data drc(14, 5)/ 2.008E+02/
      data drc(15, 5)/ 0.000E+00/
      data drc(16, 5)/ 0.000E+00/
      data dre( 1, 5)/ 1.812E-02/
      data dre( 2, 5)/ 1.133E-01/
      data dre( 3, 5)/ 1.298E+00/
      data dre( 4, 5)/ 5.689E+00/
      data dre( 5, 5)/ 3.156E+01/
      data dre( 6, 5)/ 0.000E+00/
      data dre( 7, 5)/ 0.000E+00/
      data dre( 8, 5)/ 0.000E+00/
      data dre( 9, 5)/ 5.895E+01/
      data dre(10, 5)/ 7.966E+01/
      data dre(11, 5)/ 1.262E+02/
      data dre(12, 5)/ 1.630E+02/
      data dre(13, 5)/ 4.575E+02/
      data dre(14, 5)/ 1.224E+03/
      data dre(15, 5)/ 0.000E+00/
      data dre(16, 5)/ 0.000E+00/
      data drc( 1, 6)/ 1.650E-02/
      data drc( 2, 6)/ 1.247E-01/
      data drc( 3, 6)/ 2.498E-01/
      data drc( 4, 6)/ 3.083E-01/
      data drc( 5, 6)/ 1.959E+00/
      data drc( 6, 6)/ 1.163E+01/
      data drc( 7, 6)/ 1.241E+02/
      data drc( 8, 6)/ 0.000E+00/
      data drc( 9, 6)/ 0.000E+00/
      data drc(10, 6)/ 1.128E-03/
      data drc(11, 6)/ 2.245E+00/
      data drc(12, 6)/ 3.211E+01/
      data drc(13, 6)/ 2.031E+02/
      data drc(14, 6)/ 8.883E+01/
      data drc(15, 6)/ 7.282E+01/
      data drc(16, 6)/ 1.071E+02/
      data dre( 1, 6)/ 7.146E-03/
      data dre( 2, 6)/ 7.240E-02/
      data dre( 3, 6)/ 1.505E-01/
      data dre( 4, 6)/ 8.442E-01/
      data dre( 5, 6)/ 3.055E+00/
      data dre( 6, 6)/ 1.444E+01/
      data dre( 7, 6)/ 3.606E+01/
      data dre( 8, 6)/ 0.000E+00/
      data dre( 9, 6)/ 0.000E+00/
      data dre(10, 6)/ 3.850E+01/
      data dre(11, 6)/ 6.307E+01/
      data dre(12, 6)/ 8.512E+01/
      data dre(13, 6)/ 1.252E+02/
      data dre(14, 6)/ 1.613E+02/
      data dre(15, 6)/ 1.166E+03/
      data dre(16, 6)/ 1.305E+03/
      data drc( 1, 7)/ 6.576E-04/
      data drc( 2, 7)/ 2.095E-01/
      data drc( 3, 7)/ 4.084E-01/
      data drc( 4, 7)/ 2.982E+00/
      data drc( 5, 7)/ 3.343E+00/
      data drc( 6, 7)/ 9.289E+01/
      data drc( 7, 7)/ 0.000E+00/
      data drc( 8, 7)/ 0.000E+00/
      data drc( 9, 7)/ 3.618E-10/
      data drc(10, 7)/ 2.125E-01/
      data drc(11, 7)/ 7.275E+00/
      data drc(12, 7)/ 4.375E+01/
      data drc(13, 7)/ 8.032E+01/
      data drc(14, 7)/ 1.881E+01/
      data drc(15, 7)/ 8.663E+01/
      data drc(16, 7)/ 0.000E+00/
      data dre( 1, 7)/ 1.171E-01/
      data dre( 2, 7)/ 7.081E-01/
      data dre( 3, 7)/ 1.471E+00/
      data dre( 4, 7)/ 3.974E+00/
      data dre( 5, 7)/ 9.465E+00/
      data dre( 6, 7)/ 2.970E+01/
      data dre( 7, 7)/ 0.000E+00/
      data dre( 8, 7)/ 0.000E+00/
      data dre( 9, 7)/ 2.939E+00/
      data dre(10, 7)/ 5.337E+01/
      data dre(11, 7)/ 7.081E+01/
      data dre(12, 7)/ 9.443E+01/
      data dre(13, 7)/ 1.192E+02/
      data dre(14, 7)/ 1.494E+02/
      data dre(15, 7)/ 1.263E+03/
      data dre(16, 7)/ 0.000E+00/
      data drc( 1, 8)/ 1.795E-02/
      data drc( 2, 8)/ 1.592E-02/
      data drc( 3, 8)/ 4.362E-02/
      data drc( 4, 8)/ 3.320E-01/
      data drc( 5, 8)/ 5.755E+01/
      data drc( 6, 8)/ 0.000E+00/
      data drc( 7, 8)/ 0.000E+00/
      data drc( 8, 8)/ 0.000E+00/
      data drc( 9, 8)/ 3.401E-10/
      data drc(10, 8)/ 7.459E-02/
      data drc(11, 8)/ 2.346E+00/
      data drc(12, 8)/ 2.842E+01/
      data drc(13, 8)/ 2.070E+01/
      data drc(14, 8)/ 2.940E+00/
      data drc(15, 8)/ 3.901E+01/
      data drc(16, 8)/ 0.000E+00/
      data dre( 1, 8)/ 1.021E-02/
      data dre( 2, 8)/ 9.071E-02/
      data dre( 3, 8)/ 6.347E-01/
      data dre( 4, 8)/ 3.139E+00/
      data dre( 5, 8)/ 3.308E+01/
      data dre( 6, 8)/ 0.000E+00/
      data dre( 7, 8)/ 0.000E+00/
      data dre( 8, 8)/ 0.000E+00/
      data dre( 9, 8)/ 2.245E+00/
      data dre(10, 8)/ 4.037E+01/
      data dre(11, 8)/ 5.699E+01/
      data dre(12, 8)/ 7.732E+01/
      data dre(13, 8)/ 9.917E+01/
      data dre(14, 8)/ 1.291E+02/
      data dre(15, 8)/ 1.256E+03/
      data dre(16, 8)/ 0.000E+00/
      data drc( 1, 9)/ 5.321E-04/
      data drc( 2, 9)/ 1.580E-03/
      data drc( 3, 9)/ 4.153E-03/
      data drc( 4, 9)/ 2.587E-01/
      data drc( 5, 9)/ 2.430E+01/
      data drc( 6, 9)/ 0.000E+00/
      data drc( 7, 9)/ 0.000E+00/
      data drc( 8, 9)/ 0.000E+00/
      data drc( 9, 9)/ 8.985E-11/
      data drc(10, 9)/ 9.842E-02/
      data drc(11, 9)/ 1.784E+00/
      data drc(12, 9)/ 1.532E+01/
      data drc(13, 9)/ 2.062E+00/
      data drc(14, 9)/ 2.427E-02/
      data drc(15, 9)/ 1.387E+01/
      data drc(16, 9)/ 0.000E+00/
      data dre( 1, 9)/ 2.038E-02/
      data dre( 2, 9)/ 6.697E-02/
      data dre( 3, 9)/ 1.980E-01/
      data dre( 4, 9)/ 1.900E+01/
      data dre( 5, 9)/ 3.765E+01/
      data dre( 6, 9)/ 0.000E+00/
      data dre( 7, 9)/ 0.000E+00/
      data dre( 8, 9)/ 0.000E+00/
      data dre( 9, 9)/ 1.957E+00/
      data dre(10, 9)/ 3.845E+01/
      data dre(11, 9)/ 5.154E+01/
      data dre(12, 9)/ 6.534E+01/
      data dre(13, 9)/ 8.955E+01/
      data dre(14, 9)/ 2.075E+02/
      data dre(15, 9)/ 1.249E+03/
      data dre(16, 9)/ 0.000E+00/
      data drc( 1,10)/ 0.000E+00/
      data drc( 2,10)/ 0.000E+00/
      data drc( 3,10)/ 0.000E+00/
      data drc( 4,10)/ 0.000E+00/
      data drc( 5,10)/ 0.000E+00/
      data drc( 6,10)/ 0.000E+00/
      data drc( 7,10)/ 0.000E+00/
      data drc( 8,10)/ 0.000E+00/
      data drc( 9,10)/ 1.702E-02/
      data drc(10,10)/ 8.035E-01/
      data drc(11,10)/ 5.588E+00/
      data drc(12,10)/ 6.642E-02/
      data drc(13,10)/ 1.545E+00/
      data drc(14,10)/ 0.000E+00/
      data drc(15,10)/ 0.000E+00/
      data drc(16,10)/ 0.000E+00/
      data dre( 1,10)/ 0.000E+00/
      data dre( 2,10)/ 0.000E+00/
      data dre( 3,10)/ 0.000E+00/
      data dre( 4,10)/ 0.000E+00/
      data dre( 5,10)/ 0.000E+00/
      data dre( 6,10)/ 0.000E+00/
      data dre( 7,10)/ 0.000E+00/
      data dre( 8,10)/ 0.000E+00/
      data dre( 9,10)/ 3.517E+01/
      data dre(10,10)/ 4.357E+01/
      data dre(11,10)/ 5.029E+01/
      data dre(12,10)/ 7.794E+01/
      data dre(13,10)/ 1.307E+03/
      data dre(14,10)/ 0.000E+00/
      data dre(15,10)/ 0.000E+00/
      data dre(16,10)/ 0.000E+00/
      data drc( 1,11)/ 0.000E+00/
      data drc( 2,11)/ 0.000E+00/
      data drc( 3,11)/ 0.000E+00/
      data drc( 4,11)/ 0.000E+00/
      data drc( 5,11)/ 0.000E+00/
      data drc( 6,11)/ 0.000E+00/
      data drc( 7,11)/ 0.000E+00/
      data drc( 8,11)/ 0.000E+00/
      data drc( 9,11)/ 2.573E+02/
      data drc(10,11)/ 9.152E+02/
      data drc(11,11)/ 2.616E+02/
      data drc(12,11)/ 0.000E+00/
      data drc(13,11)/ 0.000E+00/
      data drc(14,11)/ 0.000E+00/
      data drc(15,11)/ 0.000E+00/
      data drc(16,11)/ 0.000E+00/
      data dre( 1,11)/ 0.000E+00/
      data dre( 2,11)/ 0.000E+00/
      data dre( 3,11)/ 0.000E+00/
      data dre( 4,11)/ 0.000E+00/
      data dre( 5,11)/ 0.000E+00/
      data dre( 6,11)/ 0.000E+00/
      data dre( 7,11)/ 0.000E+00/
      data dre( 8,11)/ 0.000E+00/
      data dre( 9,11)/ 1.404E+03/
      data dre(10,11)/ 1.827E+03/
      data dre(11,11)/ 2.116E+03/
      data dre(12,11)/ 0.000E+00/
      data dre(13,11)/ 0.000E+00/
      data dre(14,11)/ 0.000E+00/
      data dre(15,11)/ 0.000E+00/
      data dre(16,11)/ 0.000E+00/
      data drc( 1,12)/ 0.000E+00/
      data drc( 2,12)/ 0.000E+00/
      data drc( 3,12)/ 0.000E+00/
      data drc( 4,12)/ 0.000E+00/
      data drc( 5,12)/ 0.000E+00/
      data drc( 6,12)/ 0.000E+00/
      data drc( 7,12)/ 0.000E+00/
      data drc( 8,12)/ 0.000E+00/
      data drc( 9,12)/ 3.818E+02/
      data drc(10,12)/ 1.347E+03/
      data drc(11,12)/ 5.637E+02/
      data drc(12,12)/ 0.000E+00/
      data drc(13,12)/ 0.000E+00/
      data drc(14,12)/ 0.000E+00/
      data drc(15,12)/ 0.000E+00/
      data drc(16,12)/ 0.000E+00/
      data dre( 1,12)/ 0.000E+00/
      data dre( 2,12)/ 0.000E+00/
      data dre( 3,12)/ 0.000E+00/
      data dre( 4,12)/ 0.000E+00/
      data dre( 5,12)/ 0.000E+00/
      data dre( 6,12)/ 0.000E+00/
      data dre( 7,12)/ 0.000E+00/
      data dre( 8,12)/ 0.000E+00/
      data dre( 9,12)/ 1.335E+03/
      data dre(10,12)/ 1.662E+03/
      data dre(11,12)/ 1.892E+03/
      data dre(12,12)/ 0.000E+00/
      data dre(13,12)/ 0.000E+00/
      data dre(14,12)/ 0.000E+00/
      data dre(15,12)/ 0.000E+00/
      data dre(16,12)/ 0.000E+00/
      data drc( 1,13)/ 2.217E+00/
      data drc( 2,13)/ 5.564E+00/
      data drc( 3,13)/ 7.610E+00/
      data drc( 4,13)/ 1.674E+01/
      data drc( 5,13)/ 7.010E+01/
      data drc( 6,13)/ 0.000E+00/
      data drc( 7,13)/ 0.000E+00/
      data drc( 8,13)/ 0.000E+00/
      data drc( 9,13)/ 5.706E-10/
      data drc(10,13)/ 8.994E+00/
      data drc(11,13)/ 1.478E+01/
      data drc(12,13)/ 7.819E+01/
      data drc(13,13)/ 3.732E+02/
      data drc(14,13)/ 5.156E+01/
      data drc(15,13)/ 4.140E+02/
      data drc(16,13)/ 8.922E+02/
      data dre( 1,13)/ 3.366E+00/
      data dre( 2,13)/ 4.388E+00/
      data dre( 3,13)/ 8.322E+00/
      data dre( 4,13)/ 1.481E+01/
      data dre( 5,13)/ 2.322E+01/
      data dre( 6,13)/ 0.000E+00/
      data dre( 7,13)/ 0.000E+00/
      data dre( 8,13)/ 0.000E+00/
      data dre( 9,13)/ 4.482E+00/
      data dre(10,13)/ 9.563E+01/
      data dre(11,13)/ 1.089E+02/
      data dre(12,13)/ 1.958E+02/
      data dre(13,13)/ 2.665E+02/
      data dre(14,13)/ 3.433E+02/
      data dre(15,13)/ 1.394E+03/
      data dre(16,13)/ 1.775E+03/
      data drc( 1,14)/ 5.687E-01/
      data drc( 2,14)/ 1.622E+00/
      data drc( 3,14)/ 3.915E+00/
      data drc( 4,14)/ 1.223E+01/
      data drc( 5,14)/ 3.492E+01/
      data drc( 6,14)/ 2.233E+02/
      data drc( 7,14)/ 0.000E+00/
      data drc( 8,14)/ 0.000E+00/
      data drc( 9,14)/ 1.317E+00/
      data drc(10,14)/ 2.692E+01/
      data drc(11,14)/ 3.365E+01/
      data drc(12,14)/ 3.536E+02/
      data drc(13,14)/ 1.930E+02/
      data drc(14,14)/ 1.262E+02/
      data drc(15,14)/ 5.204E+02/
      data drc(16,14)/ 0.000E+00/
      data dre( 1,14)/ 4.089E-02/
      data dre( 2,14)/ 9.832E-02/
      data dre( 3,14)/ 4.713E-01/
      data dre( 4,14)/ 1.699E+00/
      data dre( 5,14)/ 6.449E+00/
      data dre( 6,14)/ 3.315E+01/
      data dre( 7,14)/ 0.000E+00/
      data dre( 8,14)/ 0.000E+00/
      data dre( 9,14)/ 9.277E+01/
      data dre(10,14)/ 1.123E+02/
      data dre(11,14)/ 1.622E+02/
      data dre(12,14)/ 2.390E+02/
      data dre(13,14)/ 2.977E+02/
      data dre(14,14)/ 1.391E+03/
      data dre(15,14)/ 1.756E+03/
      data dre(16,14)/ 0.000E+00/
      data drc( 1,15)/ 3.512E-02/
      data drc( 2,15)/ 5.453E-01/
      data drc( 3,15)/ 3.179E+00/
      data drc( 4,15)/ 6.666E+00/
      data drc( 5,15)/ 1.898E+01/
      data drc( 6,15)/ 2.356E+02/
      data drc( 7,15)/ 0.000E+00/
      data drc( 8,15)/ 0.000E+00/
      data drc( 9,15)/ 5.496E-01/
      data drc(10,15)/ 2.753E+01/
      data drc(11,15)/ 5.901E+01/
      data drc(12,15)/ 6.198E+02/
      data drc(13,15)/ 2.914E+02/
      data drc(14,15)/ 5.778E+01/
      data drc(15,15)/ 4.350E+02/
      data drc(16,15)/ 0.000E+00/
      data dre( 1,15)/ 1.459E-02/
      data dre( 2,15)/ 8.405E-02/
      data dre( 3,15)/ 3.945E-01/
      data dre( 4,15)/ 1.455E+00/
      data dre( 5,15)/ 5.260E+00/
      data dre( 6,15)/ 3.747E+01/
      data dre( 7,15)/ 0.000E+00/
      data dre( 8,15)/ 0.000E+00/
      data dre( 9,15)/ 7.762E+01/
      data dre(10,15)/ 1.017E+02/
      data dre(11,15)/ 1.346E+02/
      data dre(12,15)/ 2.082E+02/
      data dre(13,15)/ 2.611E+02/
      data dre(14,15)/ 1.346E+03/
      data dre(15,15)/ 1.741E+03/
      data dre(16,15)/ 0.000E+00/
      data drc( 1,16)/ 8.829E-01/
      data drc( 2,16)/ 3.325E+00/
      data drc( 3,16)/ 7.595E+00/
      data drc( 4,16)/ 4.950E+00/
      data drc( 5,16)/ 7.758E+00/
      data drc( 6,16)/ 1.981E+02/
      data drc( 7,16)/ 0.000E+00/
      data drc( 8,16)/ 0.000E+00/
      data drc( 9,16)/ 7.591E-01/
      data drc(10,16)/ 4.272E+01/
      data drc(11,16)/ 1.506E+02/
      data drc(12,16)/ 7.341E+02/
      data drc(13,16)/ 6.236E+01/
      data drc(14,16)/ 2.883E+02/
      data drc(15,16)/ 0.000E+00/
      data drc(16,16)/ 0.000E+00/
      data dre( 1,16)/ 1.977E-02/
      data dre( 2,16)/ 6.566E-02/
      data dre( 3,16)/ 2.250E-01/
      data dre( 4,16)/ 8.196E-01/
      data dre( 5,16)/ 6.105E+00/
      data dre( 6,16)/ 4.024E+01/
      data dre( 7,16)/ 0.000E+00/
      data dre( 8,16)/ 0.000E+00/
      data dre( 9,16)/ 7.634E+01/
      data dre(10,16)/ 1.040E+02/
      data dre(11,16)/ 1.493E+02/
      data dre(12,16)/ 2.043E+02/
      data dre(13,16)/ 2.886E+02/
      data dre(14,16)/ 1.711E+03/
      data dre(15,16)/ 0.000E+00/
      data dre(16,16)/ 0.000E+00/
      data drc( 1,17)/ 5.128E-02/
      data drc( 2,17)/ 2.208E-01/
      data drc( 3,17)/ 1.960E-01/
      data drc( 4,17)/ 2.780E-01/
      data drc( 5,17)/ 1.316E+00/
      data drc( 6,17)/ 4.094E+00/
      data drc( 7,17)/ 1.041E+01/
      data drc( 8,17)/ 1.458E+02/
      data drc( 9,17)/ 3.017E-09/
      data drc(10,17)/ 9.798E-01/
      data drc(11,17)/ 3.039E+01/
      data drc(12,17)/ 1.123E+02/
      data drc(13,17)/ 4.049E+02/
      data drc(14,17)/ 8.913E+01/
      data drc(15,17)/ 1.503E+02/
      data drc(16,17)/ 0.000E+00/
      data dre( 1,17)/ 1.690E-01/
      data dre( 2,17)/ 2.308E-01/
      data dre( 3,17)/ 3.604E-01/
      data dre( 4,17)/ 1.254E+00/
      data dre( 5,17)/ 2.740E+00/
      data dre( 6,17)/ 7.361E+00/
      data dre( 7,17)/ 1.742E+01/
      data dre( 8,17)/ 3.758E+01/
      data dre( 9,17)/ 4.321E+00/
      data dre(10,17)/ 7.718E+01/
      data dre(11,17)/ 1.020E+02/
      data dre(12,17)/ 1.382E+02/
      data dre(13,17)/ 1.807E+02/
      data dre(14,17)/ 2.305E+02/
      data dre(15,17)/ 1.743E+03/
      data dre(16,17)/ 0.000E+00/
      data drc( 1,18)/ 5.305E-03/
      data drc( 2,18)/ 1.567E-02/
      data drc( 3,18)/ 9.574E-02/
      data drc( 4,18)/ 2.006E-01/
      data drc( 5,18)/ 1.583E+00/
      data drc( 6,18)/ 1.019E+02/
      data drc( 7,18)/ 0.000E+00/
      data drc( 8,18)/ 0.000E+00/
      data drc( 9,18)/ 1.099E-01/
      data drc(10,18)/ 3.876E+00/
      data drc(11,18)/ 5.307E+01/
      data drc(12,18)/ 2.061E+02/
      data drc(13,18)/ 8.457E+01/
      data drc(14,18)/ 7.473E+01/
      data drc(15,18)/ 0.000E+00/
      data drc(16,18)/ 0.000E+00/
      data dre( 1,18)/ 1.383E-02/
      data dre( 2,18)/ 8.427E-02/
      data dre( 3,18)/ 3.112E-01/
      data dre( 4,18)/ 1.490E+00/
      data dre( 5,18)/ 8.091E+00/
      data dre( 6,18)/ 4.189E+01/
      data dre( 7,18)/ 0.000E+00/
      data dre( 8,18)/ 0.000E+00/
      data dre( 9,18)/ 5.717E+01/
      data dre(10,18)/ 7.739E+01/
      data dre(11,18)/ 1.067E+02/
      data dre(12,18)/ 1.452E+02/
      data dre(13,18)/ 1.904E+02/
      data dre(14,18)/ 1.735E+03/
      data dre(15,18)/ 0.000E+00/
      data dre(16,18)/ 0.000E+00/
      data drc( 1,19)/ 1.960E-03/
      data drc( 2,19)/ 7.042E-03/
      data drc( 3,19)/ 5.463E-02/
      data drc( 4,19)/ 1.805E-01/
      data drc( 5,19)/ 3.228E-01/
      data drc( 6,19)/ 5.018E+01/
      data drc( 7,19)/ 0.000E+00/
      data drc( 8,19)/ 0.000E+00/
      data drc( 9,19)/ 5.739E-10/
      data drc(10,19)/ 2.971E-01/
      data drc(11,19)/ 7.153E+00/
      data drc(12,19)/ 6.213E+01/
      data drc(13,19)/ 6.306E+01/
      data drc(14,19)/ 1.384E+01/
      data drc(15,19)/ 3.038E+01/
      data drc(16,19)/ 0.000E+00/
      data dre( 1,19)/ 2.828E-02/
      data dre( 2,19)/ 9.223E-02/
      data dre( 3,19)/ 4.217E-01/
      data dre( 4,19)/ 1.136E+00/
      data dre( 5,19)/ 1.164E+01/
      data dre( 6,19)/ 4.781E+01/
      data dre( 7,19)/ 0.000E+00/
      data dre( 8,19)/ 0.000E+00/
      data dre( 9,19)/ 3.176E+00/
      data dre(10,19)/ 5.808E+01/
      data dre(11,19)/ 8.032E+01/
      data dre(12,19)/ 1.072E+02/
      data dre(13,19)/ 1.389E+02/
      data dre(14,19)/ 1.776E+02/
      data dre(15,19)/ 1.726E+03/
      data dre(16,19)/ 0.000E+00/
      data drc( 1,20)/ 0.000E+00/
      data drc( 2,20)/ 0.000E+00/
      data drc( 3,20)/ 0.000E+00/
      data drc( 4,20)/ 0.000E+00/
      data drc( 5,20)/ 0.000E+00/
      data drc( 6,20)/ 0.000E+00/
      data drc( 7,20)/ 0.000E+00/
      data drc( 8,20)/ 0.000E+00/
      data drc( 9,20)/ 1.070E-02/
      data drc(10,20)/ 4.512E-01/
      data drc(11,20)/ 1.471E+01/
      data drc(12,20)/ 4.667E+01/
      data drc(13,20)/ 1.378E+01/
      data drc(14,20)/ 7.862E+00/
      data drc(15,20)/ 0.000E+00/
      data drc(16,20)/ 0.000E+00/
      data dre( 1,20)/ 0.000E+00/
      data dre( 2,20)/ 0.000E+00/
      data dre( 3,20)/ 0.000E+00/
      data dre( 4,20)/ 0.000E+00/
      data dre( 5,20)/ 0.000E+00/
      data dre( 6,20)/ 0.000E+00/
      data dre( 7,20)/ 0.000E+00/
      data dre( 8,20)/ 0.000E+00/
      data dre( 9,20)/ 4.788E+01/
      data dre(10,20)/ 5.901E+01/
      data dre(11,20)/ 8.348E+01/
      data dre(12,20)/ 1.022E+02/
      data dre(13,20)/ 1.354E+02/
      data dre(14,20)/ 1.839E+03/
      data dre(15,20)/ 0.000E+00/
      data dre(16,20)/ 0.000E+00/
      data drc( 1,21)/ 0.000E+00/
      data drc( 2,21)/ 0.000E+00/
      data drc( 3,21)/ 0.000E+00/
      data drc( 4,21)/ 0.000E+00/
      data drc( 5,21)/ 0.000E+00/
      data drc( 6,21)/ 0.000E+00/
      data drc( 7,21)/ 0.000E+00/
      data drc( 8,21)/ 0.000E+00/
      data drc( 9,21)/ 3.897E+02/
      data drc(10,21)/ 1.122E+03/
      data drc(11,21)/ 2.994E+02/
      data drc(12,21)/ 0.000E+00/
      data drc(13,21)/ 0.000E+00/
      data drc(14,21)/ 0.000E+00/
      data drc(15,21)/ 0.000E+00/
      data drc(16,21)/ 0.000E+00/
      data dre( 1,21)/ 0.000E+00/
      data dre( 2,21)/ 0.000E+00/
      data dre( 3,21)/ 0.000E+00/
      data dre( 4,21)/ 0.000E+00/
      data dre( 5,21)/ 0.000E+00/
      data dre( 6,21)/ 0.000E+00/
      data dre( 7,21)/ 0.000E+00/
      data dre( 8,21)/ 0.000E+00/
      data dre( 9,21)/ 1.826E+03/
      data dre(10,21)/ 2.379E+03/
      data dre(11,21)/ 2.774E+03/
      data dre(12,21)/ 0.000E+00/
      data dre(13,21)/ 0.000E+00/
      data dre(14,21)/ 0.000E+00/
      data dre(15,21)/ 0.000E+00/
      data dre(16,21)/ 0.000E+00/
      data drc( 1,22)/ 0.000E+00/
      data drc( 2,22)/ 0.000E+00/
      data drc( 3,22)/ 0.000E+00/
      data drc( 4,22)/ 0.000E+00/
      data drc( 5,22)/ 0.000E+00/
      data drc( 6,22)/ 0.000E+00/
      data drc( 7,22)/ 0.000E+00/
      data drc( 8,22)/ 0.000E+00/
      data drc( 9,22)/ 6.266E+02/
      data drc(10,22)/ 1.796E+03/
      data drc(11,22)/ 5.979E+02/
      data drc(12,22)/ 0.000E+00/
      data drc(13,22)/ 0.000E+00/
      data drc(14,22)/ 0.000E+00/
      data drc(15,22)/ 0.000E+00/
      data drc(16,22)/ 0.000E+00/
      data dre( 1,22)/ 0.000E+00/
      data dre( 2,22)/ 0.000E+00/
      data dre( 3,22)/ 0.000E+00/
      data dre( 4,22)/ 0.000E+00/
      data dre( 5,22)/ 0.000E+00/
      data dre( 6,22)/ 0.000E+00/
      data dre( 7,22)/ 0.000E+00/
      data dre( 8,22)/ 0.000E+00/
      data dre( 9,22)/ 1.749E+03/
      data dre(10,22)/ 2.195E+03/
      data dre(11,22)/ 2.520E+03/
      data dre(12,22)/ 0.000E+00/
      data dre(13,22)/ 0.000E+00/
      data dre(14,22)/ 0.000E+00/
      data dre(15,22)/ 0.000E+00/
      data dre(16,22)/ 0.000E+00/
      data drc( 1,23)/ 1.900E-02/
      data drc( 2,23)/ 2.221E-01/
      data drc( 3,23)/ 1.613E+00/
      data drc( 4,23)/ 7.881E+00/
      data drc( 5,23)/ 6.133E+00/
      data drc( 6,23)/ 1.275E+01/
      data drc( 7,23)/ 2.714E+01/
      data drc( 8,23)/ 8.928E+01/
      data drc( 9,23)/ 3.087E-02/
      data drc(10,23)/ 2.482E+01/
      data drc(11,23)/ 2.000E+01/
      data drc(12,23)/ 1.623E+02/
      data drc(13,23)/ 6.449E+02/
      data drc(14,23)/ 5.974E+01/
      data drc(15,23)/ 6.264E+02/
      data drc(16,23)/ 1.146E+03/
      data dre( 1,23)/ 5.723E-01/
      data dre( 2,23)/ 6.415E-01/
      data dre( 3,23)/ 1.017E+00/
      data dre( 4,23)/ 1.291E+00/
      data dre( 5,23)/ 5.014E+00/
      data dre( 6,23)/ 8.605E+00/
      data dre( 7,23)/ 1.683E+01/
      data dre( 8,23)/ 2.738E+01/
      data dre( 9,23)/ 9.864E+01/
      data dre(10,23)/ 1.251E+02/
      data dre(11,23)/ 1.446E+02/
      data dre(12,23)/ 2.658E+02/
      data dre(13,23)/ 3.628E+02/
      data dre(14,23)/ 4.846E+02/
      data dre(15,23)/ 1.817E+03/
      data dre(16,23)/ 2.333E+03/
      data drc( 1,24)/ 1.327E-01/
      data drc( 2,24)/ 1.379E+00/
      data drc( 3,24)/ 7.191E+00/
      data drc( 4,24)/ 5.902E+01/
      data drc( 5,24)/ 2.909E+02/
      data drc( 6,24)/ 0.000E+00/
      data drc( 7,24)/ 0.000E+00/
      data drc( 8,24)/ 0.000E+00/
      data drc( 9,24)/ 3.588E+00/
      data drc(10,24)/ 5.808E+01/
      data drc(11,24)/ 5.970E+01/
      data drc(12,24)/ 6.341E+02/
      data drc(13,24)/ 4.036E+02/
      data drc(14,24)/ 2.239E+02/
      data drc(15,24)/ 7.586E+02/
      data drc(16,24)/ 0.000E+00/
      data dre( 1,24)/ 2.202E-02/
      data dre( 2,24)/ 3.530E-01/
      data dre( 3,24)/ 1.265E+00/
      data dre( 4,24)/ 8.368E+00/
      data dre( 5,24)/ 3.701E+01/
      data dre( 6,24)/ 0.000E+00/
      data dre( 7,24)/ 0.000E+00/
      data dre( 8,24)/ 0.000E+00/
      data dre( 9,24)/ 1.210E+02/
      data dre(10,24)/ 1.444E+02/
      data dre(11,24)/ 2.134E+02/
      data dre(12,24)/ 3.171E+02/
      data dre(13,24)/ 3.980E+02/
      data dre(14,24)/ 1.803E+03/
      data dre(15,24)/ 2.304E+03/
      data dre(16,24)/ 0.000E+00/
      data drc( 1,25)/ 7.756E-02/
      data drc( 2,25)/ 1.580E-01/
      data drc( 3,25)/ 1.047E+00/
      data drc( 4,25)/ 1.866E+00/
      data drc( 5,25)/ 3.685E+01/
      data drc( 6,25)/ 3.154E+02/
      data drc( 7,25)/ 0.000E+00/
      data drc( 8,25)/ 0.000E+00/
      data drc( 9,25)/ 1.700E+00/
      data drc(10,25)/ 6.715E+01/
      data drc(11,25)/ 9.961E+01/
      data drc(12,25)/ 1.089E+03/
      data drc(13,25)/ 7.750E+02/
      data drc(14,25)/ 7.008E+01/
      data drc(15,25)/ 6.792E+02/
      data drc(16,25)/ 0.000E+00/
      data dre( 1,25)/ 1.586E-02/
      data dre( 2,25)/ 1.284E-01/
      data dre( 3,25)/ 5.914E-01/
      data dre( 4,25)/ 1.861E+00/
      data dre( 5,25)/ 1.031E+01/
      data dre( 6,25)/ 4.449E+01/
      data dre( 7,25)/ 0.000E+00/
      data dre( 8,25)/ 0.000E+00/
      data dre( 9,25)/ 1.044E+02/
      data dre(10,25)/ 1.341E+02/
      data dre(11,25)/ 1.765E+02/
      data dre(12,25)/ 2.798E+02/
      data dre(13,25)/ 3.505E+02/
      data dre(14,25)/ 1.634E+03/
      data dre(15,25)/ 2.267E+03/
      data dre(16,25)/ 0.000E+00/
      data drc( 1,26)/ 1.056E-01/
      data drc( 2,26)/ 7.798E-01/
      data drc( 3,26)/ 1.255E+00/
      data drc( 4,26)/ 2.045E+00/
      data drc( 5,26)/ 2.124E+01/
      data drc( 6,26)/ 2.604E+02/
      data drc( 7,26)/ 0.000E+00/
      data drc( 8,26)/ 0.000E+00/
      data drc( 9,26)/ 3.013E+00/
      data drc(10,26)/ 1.224E+02/
      data drc(11,26)/ 3.143E+02/
      data drc(12,26)/ 1.504E+03/
      data drc(13,26)/ 1.218E+02/
      data drc(14,26)/ 4.176E+02/
      data drc(15,26)/ 0.000E+00/
      data drc(16,26)/ 0.000E+00/
      data dre( 1,26)/ 1.450E-02/
      data dre( 2,26)/ 8.243E-02/
      data dre( 3,26)/ 2.267E-01/
      data dre( 4,26)/ 1.182E+00/
      data dre( 5,26)/ 8.096E+00/
      data dre( 6,26)/ 4.556E+01/
      data dre( 7,26)/ 0.000E+00/
      data dre( 8,26)/ 0.000E+00/
      data dre( 9,26)/ 1.047E+02/
      data dre(10,26)/ 1.404E+02/
      data dre(11,26)/ 2.097E+02/
      data dre(12,26)/ 2.872E+02/
      data dre(13,26)/ 4.070E+02/
      data dre(14,26)/ 2.240E+03/
      data dre(15,26)/ 0.000E+00/
      data dre(16,26)/ 0.000E+00/
      data drc( 1,27)/ 1.892E-01/
      data drc( 2,27)/ 1.178E-01/
      data drc( 3,27)/ 2.889E-01/
      data drc( 4,27)/ 2.313E+00/
      data drc( 5,27)/ 1.049E+01/
      data drc( 6,27)/ 1.064E+01/
      data drc( 7,27)/ 1.760E+01/
      data drc( 8,27)/ 2.086E+02/
      data drc( 9,27)/ 7.324E-07/
      data drc(10,27)/ 3.037E+00/
      data drc(11,27)/ 7.947E+01/
      data drc(12,27)/ 2.159E+02/
      data drc(13,27)/ 1.059E+03/
      data drc(14,27)/ 2.328E+02/
      data drc(15,27)/ 2.216E+02/
      data drc(16,27)/ 0.000E+00/
      data dre( 1,27)/ 1.136E-01/
      data dre( 2,27)/ 1.701E-01/
      data dre( 3,27)/ 3.541E-01/
      data dre( 4,27)/ 1.089E+00/
      data dre( 5,27)/ 2.692E+00/
      data dre( 6,27)/ 6.056E+00/
      data dre( 7,27)/ 1.696E+01/
      data dre( 8,27)/ 4.491E+01/
      data dre( 9,27)/ 3.187E+01/
      data dre(10,27)/ 1.041E+02/
      data dre(11,27)/ 1.367E+02/
      data dre(12,27)/ 1.890E+02/
      data dre(13,27)/ 2.550E+02/
      data dre(14,27)/ 3.308E+02/
      data dre(15,27)/ 2.300E+03/
      data dre(16,27)/ 0.000E+00/
      data drc( 1,28)/ 1.093E-01/
      data drc( 2,28)/ 1.158E-01/
      data drc( 3,28)/ 1.571E+00/
      data drc( 4,28)/ 1.768E+00/
      data drc( 5,28)/ 1.955E+00/
      data drc( 6,28)/ 1.509E+02/
      data drc( 7,28)/ 0.000E+00/
      data drc( 8,28)/ 0.000E+00/
      data drc( 9,28)/ 5.200E-01/
      data drc(10,28)/ 2.075E+01/
      data drc(11,28)/ 1.621E+02/
      data drc(12,28)/ 7.323E+02/
      data drc(13,28)/ 2.613E+02/
      data drc(14,28)/ 1.171E+02/
      data drc(15,28)/ 0.000E+00/
      data drc(16,28)/ 0.000E+00/
      data dre( 1,28)/ 1.879E-02/
      data dre( 2,28)/ 4.330E-02/
      data dre( 3,28)/ 4.894E-01/
      data dre( 4,28)/ 1.042E+00/
      data dre( 5,28)/ 7.896E+00/
      data dre( 6,28)/ 5.008E+01/
      data dre( 7,28)/ 0.000E+00/
      data dre( 8,28)/ 0.000E+00/
      data dre( 9,28)/ 8.179E+01/
      data dre(10,28)/ 1.132E+02/
      data dre(11,28)/ 1.519E+02/
      data dre(12,28)/ 2.175E+02/
      data dre(13,28)/ 2.837E+02/
      data dre(14,28)/ 2.290E+03/
      data dre(15,28)/ 0.000E+00/
      data dre(16,28)/ 0.000E+00/
      data drc( 1,29)/ 1.610E-03/
      data drc( 2,29)/ 1.350E-02/
      data drc( 3,29)/ 3.504E-02/
      data drc( 4,29)/ 1.130E-01/
      data drc( 5,29)/ 1.261E+00/
      data drc( 6,29)/ 7.944E+01/
      data drc( 7,29)/ 0.000E+00/
      data drc( 8,29)/ 0.000E+00/
      data drc( 9,29)/ 7.746E-09/
      data drc(10,29)/ 1.266E+00/
      data drc(11,29)/ 4.619E+01/
      data drc(12,29)/ 2.495E+02/
      data drc(13,29)/ 3.432E+02/
      data drc(14,29)/ 2.466E+01/
      data drc(15,29)/ 3.855E-03/
      data drc(16,29)/ 5.095E+01/
      data dre( 1,29)/ 5.532E-02/
      data dre( 2,29)/ 1.052E-01/
      data dre( 3,29)/ 3.096E-01/
      data dre( 4,29)/ 8.876E-01/
      data dre( 5,29)/ 1.525E+01/
      data dre( 6,29)/ 5.814E+01/
      data dre( 7,29)/ 0.000E+00/
      data dre( 8,29)/ 0.000E+00/
      data dre( 9,29)/ 4.520E+00/
      data dre(10,29)/ 8.343E+01/
      data dre(11,29)/ 1.208E+02/
      data dre(12,29)/ 1.672E+02/
      data dre(13,29)/ 2.228E+02/
      data dre(14,29)/ 2.987E+02/
      data dre(15,29)/ 9.110E+02/
      data dre(16,29)/ 2.278E+03/
      data drc( 1,30)/ 0.000E+00/
      data drc( 2,30)/ 0.000E+00/
      data drc( 3,30)/ 0.000E+00/
      data drc( 4,30)/ 0.000E+00/
      data drc( 5,30)/ 0.000E+00/
      data drc( 6,30)/ 0.000E+00/
      data drc( 7,30)/ 0.000E+00/
      data drc( 8,30)/ 0.000E+00/
      data drc( 9,30)/ 2.377E-09/
      data drc(10,30)/ 5.554E-01/
      data drc(11,30)/ 1.608E+01/
      data drc(12,30)/ 1.030E+02/
      data drc(13,30)/ 2.034E+02/
      data drc(14,30)/ 1.064E+02/
      data drc(15,30)/ 1.959E+01/
      data drc(16,30)/ 0.000E+00/
      data dre( 1,30)/ 0.000E+00/
      data dre( 2,30)/ 0.000E+00/
      data dre( 3,30)/ 0.000E+00/
      data dre( 4,30)/ 0.000E+00/
      data dre( 5,30)/ 0.000E+00/
      data dre( 6,30)/ 0.000E+00/
      data dre( 7,30)/ 0.000E+00/
      data dre( 8,30)/ 0.000E+00/
      data dre( 9,30)/ 1.070E+01/
      data dre(10,30)/ 7.740E+01/
      data dre(11,30)/ 1.081E+02/
      data dre(12,30)/ 1.402E+02/
      data dre(13,30)/ 1.812E+02/
      data dre(14,30)/ 2.222E+02/
      data dre(15,30)/ 2.465E+03/
      data dre(16,30)/ 0.000E+00/
      data drc( 1,31)/ 0.000E+00/
      data drc( 2,31)/ 0.000E+00/
      data drc( 3,31)/ 0.000E+00/
      data drc( 4,31)/ 0.000E+00/
      data drc( 5,31)/ 0.000E+00/
      data drc( 6,31)/ 0.000E+00/
      data drc( 7,31)/ 0.000E+00/
      data drc( 8,31)/ 0.000E+00/
      data drc( 9,31)/ 5.521E+02/
      data drc(10,31)/ 1.309E+03/
      data drc(11,31)/ 3.424E+02/
      data drc(12,31)/ 0.000E+00/
      data drc(13,31)/ 0.000E+00/
      data drc(14,31)/ 0.000E+00/
      data drc(15,31)/ 0.000E+00/
      data drc(16,31)/ 0.000E+00/
      data dre( 1,31)/ 0.000E+00/
      data dre( 2,31)/ 0.000E+00/
      data dre( 3,31)/ 0.000E+00/
      data dre( 4,31)/ 0.000E+00/
      data dre( 5,31)/ 0.000E+00/
      data dre( 6,31)/ 0.000E+00/
      data dre( 7,31)/ 0.000E+00/
      data dre( 8,31)/ 0.000E+00/
      data dre( 9,31)/ 2.304E+03/
      data dre(10,31)/ 3.004E+03/
      data dre(11,31)/ 3.525E+03/
      data dre(12,31)/ 0.000E+00/
      data dre(13,31)/ 0.000E+00/
      data dre(14,31)/ 0.000E+00/
      data dre(15,31)/ 0.000E+00/
      data dre(16,31)/ 0.000E+00/
      data drc( 1,32)/ 0.000E+00/
      data drc( 2,32)/ 0.000E+00/
      data drc( 3,32)/ 0.000E+00/
      data drc( 4,32)/ 0.000E+00/
      data drc( 5,32)/ 0.000E+00/
      data drc( 6,32)/ 0.000E+00/
      data drc( 7,32)/ 0.000E+00/
      data drc( 8,32)/ 0.000E+00/
      data drc( 9,32)/ 9.355E+02/
      data drc(10,32)/ 2.178E+03/
      data drc(11,32)/ 6.700E+02/
      data drc(12,32)/ 0.000E+00/
      data drc(13,32)/ 0.000E+00/
      data drc(14,32)/ 0.000E+00/
      data drc(15,32)/ 0.000E+00/
      data drc(16,32)/ 0.000E+00/
      data dre( 1,32)/ 0.000E+00/
      data dre( 2,32)/ 0.000E+00/
      data dre( 3,32)/ 0.000E+00/
      data dre( 4,32)/ 0.000E+00/
      data dre( 5,32)/ 0.000E+00/
      data dre( 6,32)/ 0.000E+00/
      data dre( 7,32)/ 0.000E+00/
      data dre( 8,32)/ 0.000E+00/
      data dre( 9,32)/ 2.218E+03/
      data dre(10,32)/ 2.799E+03/
      data dre(11,32)/ 3.228E+03/
      data dre(12,32)/ 0.000E+00/
      data dre(13,32)/ 0.000E+00/
      data dre(14,32)/ 0.000E+00/
      data dre(15,32)/ 0.000E+00/
      data dre(16,32)/ 0.000E+00/
      data drc( 1,33)/ 1.168E-01/
      data drc( 2,33)/ 1.057E+00/
      data drc( 3,33)/ 5.268E+00/
      data drc( 4,33)/ 1.434E+01/
      data drc( 5,33)/ 1.795E+01/
      data drc( 6,33)/ 3.950E+01/
      data drc( 7,33)/ 1.160E+02/
      data drc( 8,33)/ 0.000E+00/
      data drc( 9,33)/ 1.345E-09/
      data drc(10,33)/ 3.181E+01/
      data drc(11,33)/ 4.208E+01/
      data drc(12,33)/ 2.225E+02/
      data drc(13,33)/ 9.391E+02/
      data drc(14,33)/ 1.476E+02/
      data drc(15,33)/ 8.369E+02/
      data drc(16,33)/ 1.145E+03/
      data dre( 1,33)/ 6.772E-01/
      data dre( 2,33)/ 1.012E+00/
      data dre( 3,33)/ 1.244E+00/
      data dre( 4,33)/ 4.448E+00/
      data dre( 5,33)/ 9.371E+00/
      data dre( 6,33)/ 1.845E+01/
      data dre( 7,33)/ 3.171E+01/
      data dre( 8,33)/ 0.000E+00/
      data dre( 9,33)/ 7.042E+00/
      data dre(10,33)/ 1.539E+02/
      data dre(11,33)/ 1.739E+02/
      data dre(12,33)/ 3.331E+02/
      data dre(13,33)/ 4.602E+02/
      data dre(14,33)/ 6.000E+02/
      data dre(15,33)/ 2.284E+03/
      data dre(16,33)/ 2.942E+03/
      data drc( 1,34)/ 7.848E-01/
      data drc( 2,34)/ 1.293E+00/
      data drc( 3,34)/ 4.132E+00/
      data drc( 4,34)/ 2.157E+01/
      data drc( 5,34)/ 2.099E+01/
      data drc( 6,34)/ 3.083E+01/
      data drc( 7,34)/ 9.761E+01/
      data drc( 8,34)/ 3.238E+02/
      data drc( 9,34)/ 7.151E+00/
      data drc(10,34)/ 1.044E+02/
      data drc(11,34)/ 8.845E+01/
      data drc(12,34)/ 8.897E+02/
      data drc(13,34)/ 6.782E+02/
      data drc(14,34)/ 3.252E+02/
      data drc(15,34)/ 7.570E+02/
      data drc(16,34)/ 0.000E+00/
      data dre( 1,34)/ 1.528E-01/
      data dre( 2,34)/ 2.311E-01/
      data dre( 3,34)/ 3.719E-01/
      data dre( 4,34)/ 1.311E+00/
      data dre( 5,34)/ 2.342E+00/
      data dre( 6,34)/ 7.658E+00/
      data dre( 7,34)/ 1.890E+01/
      data dre( 8,34)/ 5.101E+01/
      data dre( 9,34)/ 1.517E+02/
      data dre(10,34)/ 1.799E+02/
      data dre(11,34)/ 2.672E+02/
      data dre(12,34)/ 4.010E+02/
      data dre(13,34)/ 5.161E+02/
      data dre(14,34)/ 2.246E+03/
      data dre(15,34)/ 2.861E+03/
      data dre(16,34)/ 0.000E+00/
      data drc( 1,35)/ 4.609E-01/
      data drc( 2,35)/ 1.516E+00/
      data drc( 3,35)/ 2.946E+00/
      data drc( 4,35)/ 2.387E+01/
      data drc( 5,35)/ 4.808E+01/
      data drc( 6,35)/ 4.109E+02/
      data drc( 7,35)/ 0.000E+00/
      data drc( 8,35)/ 0.000E+00/
      data drc( 9,35)/ 1.482E+01/
      data drc(10,35)/ 1.902E+02/
      data drc(11,35)/ 7.533E+02/
      data drc(12,35)/ 2.401E+03/
      data drc(13,35)/ 4.651E+01/
      data drc(14,35)/ 7.227E+02/
      data drc(15,35)/ 0.000E+00/
      data drc(16,35)/ 0.000E+00/
      data dre( 1,35)/ 2.465E-02/
      data dre( 2,35)/ 1.180E-01/
      data dre( 3,35)/ 5.184E-01/
      data dre( 4,35)/ 2.116E+00/
      data dre( 5,35)/ 7.426E+00/
      data dre( 6,35)/ 4.758E+01/
      data dre( 7,35)/ 0.000E+00/
      data dre( 8,35)/ 0.000E+00/
      data dre( 9,35)/ 1.432E+02/
      data dre(10,35)/ 1.802E+02/
      data dre(11,35)/ 3.129E+02/
      data dre(12,35)/ 4.241E+02/
      data dre(13,35)/ 8.560E+02/
      data dre(14,35)/ 2.723E+03/
      data dre(15,35)/ 0.000E+00/
      data dre(16,35)/ 0.000E+00/
      data drc( 1,36)/ 2.286E-01/
      data drc( 2,36)/ 2.634E-01/
      data drc( 3,36)/ 3.902E+00/
      data drc( 4,36)/ 6.168E+00/
      data drc( 5,36)/ 2.564E+01/
      data drc( 6,36)/ 3.403E+02/
      data drc( 7,36)/ 0.000E+00/
      data drc( 8,36)/ 0.000E+00/
      data drc( 9,36)/ 6.298E+00/
      data drc(10,36)/ 2.302E+02/
      data drc(11,36)/ 4.504E+02/
      data drc(12,36)/ 2.722E+03/
      data drc(13,36)/ 4.537E+02/
      data drc(14,36)/ 3.532E+02/
      data drc(15,36)/ 0.000E+00/
      data drc(16,36)/ 0.000E+00/
      data dre( 1,36)/ 1.369E-02/
      data dre( 2,36)/ 9.753E-02/
      data dre( 3,36)/ 3.720E-01/
      data dre( 4,36)/ 1.697E+00/
      data dre( 5,36)/ 9.005E+00/
      data dre( 6,36)/ 5.257E+01/
      data dre( 7,36)/ 0.000E+00/
      data dre( 8,36)/ 0.000E+00/
      data dre( 9,36)/ 1.337E+02/
      data dre(10,36)/ 1.768E+02/
      data dre(11,36)/ 2.645E+02/
      data dre(12,36)/ 3.706E+02/
      data dre(13,36)/ 4.920E+02/
      data dre(14,36)/ 2.737E+03/
      data dre(15,36)/ 0.000E+00/
      data dre(16,36)/ 0.000E+00/
      data drc( 1,37)/ 4.656E-02/
      data drc( 2,37)/ 3.016E-01/
      data drc( 3,37)/ 3.969E-01/
      data drc( 4,37)/ 7.910E-01/
      data drc( 5,37)/ 3.332E+00/
      data drc( 6,37)/ 2.357E+01/
      data drc( 7,37)/ 2.831E+01/
      data drc( 8,37)/ 2.757E+02/
      data drc( 9,37)/ 3.816E+00/
      data drc(10,37)/ 1.501E+02/
      data drc(11,37)/ 3.062E+02/
      data drc(12,37)/ 2.099E+03/
      data drc(13,37)/ 7.481E+02/
      data drc(14,37)/ 2.998E+02/
      data drc(15,37)/ 0.000E+00/
      data drc(16,37)/ 0.000E+00/
      data dre( 1,37)/ 1.279E-01/
      data dre( 2,37)/ 2.053E-01/
      data dre( 3,37)/ 4.183E-01/
      data dre( 4,37)/ 9.327E-01/
      data dre( 5,37)/ 2.582E+00/
      data dre( 6,37)/ 8.247E+00/
      data dre( 7,37)/ 2.075E+01/
      data dre( 8,37)/ 5.253E+01/
      data dre( 9,37)/ 1.316E+02/
      data dre(10,37)/ 1.704E+02/
      data dre(11,37)/ 2.261E+02/
      data dre(12,37)/ 3.312E+02/
      data dre(13,37)/ 4.244E+02/
      data dre(14,37)/ 2.918E+03/
      data dre(15,37)/ 0.000E+00/
      data dre(16,37)/ 0.000E+00/
      data drc( 1,38)/ 7.478E-03/
      data drc( 2,38)/ 2.994E-02/
      data drc( 3,38)/ 1.450E-01/
      data drc( 4,38)/ 1.039E+00/
      data drc( 5,38)/ 4.018E+00/
      data drc( 6,38)/ 1.047E+01/
      data drc( 7,38)/ 2.022E+02/
      data drc( 8,38)/ 0.000E+00/
      data drc( 9,38)/ 1.288E+00/
      data drc(10,38)/ 5.744E+01/
      data drc(11,38)/ 3.256E+02/
      data drc(12,38)/ 1.564E+03/
      data drc(13,38)/ 6.938E+02/
      data drc(14,38)/ 1.649E+02/
      data drc(15,38)/ 0.000E+00/
      data drc(16,38)/ 0.000E+00/
      data dre( 1,38)/ 2.904E-02/
      data dre( 2,38)/ 6.746E-02/
      data dre( 3,38)/ 2.791E-01/
      data dre( 4,38)/ 9.292E-01/
      data dre( 5,38)/ 2.709E+00/
      data dre( 6,38)/ 9.703E+00/
      data dre( 7,38)/ 5.889E+01/
      data dre( 8,38)/ 0.000E+00/
      data dre( 9,38)/ 1.083E+02/
      data dre(10,38)/ 1.506E+02/
      data dre(11,38)/ 1.975E+02/
      data dre(12,38)/ 2.944E+02/
      data dre(13,38)/ 3.828E+02/
      data dre(14,38)/ 2.902E+03/
      data dre(15,38)/ 0.000E+00/
      data dre(16,38)/ 0.000E+00/
      data drc( 1,39)/ 1.522E-03/
      data drc( 2,39)/ 3.091E-02/
      data drc( 3,39)/ 5.125E-02/
      data drc( 4,39)/ 1.401E-01/
      data drc( 5,39)/ 5.304E-01/
      data drc( 6,39)/ 1.346E+00/
      data drc( 7,39)/ 5.391E+00/
      data drc( 8,39)/ 1.083E+02/
      data drc( 9,39)/ 6.585E-07/
      data drc(10,39)/ 2.031E+00/
      data drc(11,39)/ 7.411E+01/
      data drc(12,39)/ 2.865E+02/
      data drc(13,39)/ 1.110E+03/
      data drc(14,39)/ 3.692E+02/
      data drc(15,39)/ 7.531E+01/
      data drc(16,39)/ 0.000E+00/
      data dre( 1,39)/ 8.101E-02/
      data dre( 2,39)/ 1.506E-01/
      data dre( 3,39)/ 3.621E-01/
      data dre( 4,39)/ 8.794E-01/
      data dre( 5,39)/ 2.204E+00/
      data dre( 6,39)/ 5.522E+00/
      data dre( 7,39)/ 3.378E+01/
      data dre( 8,39)/ 7.086E+01/
      data dre( 9,39)/ 3.282E+01/
      data dre(10,39)/ 1.079E+02/
      data dre(11,39)/ 1.510E+02/
      data dre(12,39)/ 2.000E+02/
      data dre(13,39)/ 2.800E+02/
      data dre(14,39)/ 3.618E+02/
      data dre(15,39)/ 2.905E+03/
      data dre(16,39)/ 0.000E+00/
      data drc( 1,40)/ 0.000E+00/
      data drc( 2,40)/ 0.000E+00/
      data drc( 3,40)/ 0.000E+00/
      data drc( 4,40)/ 0.000E+00/
      data drc( 5,40)/ 0.000E+00/
      data drc( 6,40)/ 0.000E+00/
      data drc( 7,40)/ 0.000E+00/
      data drc( 8,40)/ 0.000E+00/
      data drc( 9,40)/ 1.695E-05/
      data drc(10,40)/ 1.782E+00/
      data drc(11,40)/ 6.525E+01/
      data drc(12,40)/ 2.708E+02/
      data drc(13,40)/ 7.781E+02/
      data drc(14,40)/ 2.841E+02/
      data drc(15,40)/ 3.775E+01/
      data drc(16,40)/ 0.000E+00/
      data dre( 1,40)/ 0.000E+00/
      data dre( 2,40)/ 0.000E+00/
      data dre( 3,40)/ 0.000E+00/
      data dre( 4,40)/ 0.000E+00/
      data dre( 5,40)/ 0.000E+00/
      data dre( 6,40)/ 0.000E+00/
      data dre( 7,40)/ 0.000E+00/
      data dre( 8,40)/ 0.000E+00/
      data dre( 9,40)/ 5.178E+01/
      data dre(10,40)/ 1.050E+02/
      data dre(11,40)/ 1.493E+02/
      data dre(12,40)/ 1.971E+02/
      data dre(13,40)/ 2.636E+02/
      data dre(14,40)/ 3.285E+02/
      data dre(15,40)/ 3.182E+03/
      data dre(16,40)/ 0.000E+00/
      data drc( 1,41)/ 0.000E+00/
      data drc( 2,41)/ 0.000E+00/
      data drc( 3,41)/ 0.000E+00/
      data drc( 4,41)/ 0.000E+00/
      data drc( 5,41)/ 0.000E+00/
      data drc( 6,41)/ 0.000E+00/
      data drc( 7,41)/ 0.000E+00/
      data drc( 8,41)/ 0.000E+00/
      data drc( 9,41)/ 7.375E+02/
      data drc(10,41)/ 1.486E+03/
      data drc(11,41)/ 3.741E+02/
      data drc(12,41)/ 0.000E+00/
      data drc(13,41)/ 0.000E+00/
      data drc(14,41)/ 0.000E+00/
      data drc(15,41)/ 0.000E+00/
      data drc(16,41)/ 0.000E+00/
      data dre( 1,41)/ 0.000E+00/
      data dre( 2,41)/ 0.000E+00/
      data dre( 3,41)/ 0.000E+00/
      data dre( 4,41)/ 0.000E+00/
      data dre( 5,41)/ 0.000E+00/
      data dre( 6,41)/ 0.000E+00/
      data dre( 7,41)/ 0.000E+00/
      data dre( 8,41)/ 0.000E+00/
      data dre( 9,41)/ 2.838E+03/
      data dre(10,41)/ 3.703E+03/
      data dre(11,41)/ 4.372E+03/
      data dre(12,41)/ 0.000E+00/
      data dre(13,41)/ 0.000E+00/
      data dre(14,41)/ 0.000E+00/
      data dre(15,41)/ 0.000E+00/
      data dre(16,41)/ 0.000E+00/
      data drc( 1,42)/ 0.000E+00/
      data drc( 2,42)/ 0.000E+00/
      data drc( 3,42)/ 0.000E+00/
      data drc( 4,42)/ 0.000E+00/
      data drc( 5,42)/ 0.000E+00/
      data drc( 6,42)/ 0.000E+00/
      data drc( 7,42)/ 0.000E+00/
      data drc( 8,42)/ 0.000E+00/
      data drc( 9,42)/ 1.281E+03/
      data drc(10,42)/ 2.456E+03/
      data drc(11,42)/ 7.871E+02/
      data drc(12,42)/ 0.000E+00/
      data drc(13,42)/ 0.000E+00/
      data drc(14,42)/ 0.000E+00/
      data drc(15,42)/ 0.000E+00/
      data drc(16,42)/ 0.000E+00/
      data dre( 1,42)/ 0.000E+00/
      data dre( 2,42)/ 0.000E+00/
      data dre( 3,42)/ 0.000E+00/
      data dre( 4,42)/ 0.000E+00/
      data dre( 5,42)/ 0.000E+00/
      data dre( 6,42)/ 0.000E+00/
      data dre( 7,42)/ 0.000E+00/
      data dre( 8,42)/ 0.000E+00/
      data dre( 9,42)/ 2.741E+03/
      data dre(10,42)/ 3.469E+03/
      data dre(11,42)/ 4.014E+03/
      data dre(12,42)/ 0.000E+00/
      data dre(13,42)/ 0.000E+00/
      data dre(14,42)/ 0.000E+00/
      data dre(15,42)/ 0.000E+00/
      data dre(16,42)/ 0.000E+00/
      data drc( 1,43)/ 3.538E-01/
      data drc( 2,43)/ 3.636E+00/
      data drc( 3,43)/ 1.646E+01/
      data drc( 4,43)/ 8.326E+00/
      data drc( 5,43)/ 2.849E+01/
      data drc( 6,43)/ 5.659E+01/
      data drc( 7,43)/ 1.464E+02/
      data drc( 8,43)/ 0.000E+00/
      data drc( 9,43)/ 1.627E-09/
      data drc(10,43)/ 5.046E+01/
      data drc(11,43)/ 6.495E+01/
      data drc(12,43)/ 3.197E+02/
      data drc(13,43)/ 1.280E+03/
      data drc(14,43)/ 2.337E+02/
      data drc(15,43)/ 1.090E+03/
      data drc(16,43)/ 1.356E+03/
      data dre( 1,43)/ 1.020E+00/
      data dre( 2,43)/ 1.413E+00/
      data dre( 3,43)/ 1.706E+00/
      data dre( 4,43)/ 3.514E+00/
      data dre( 5,43)/ 9.323E+00/
      data dre( 6,43)/ 2.013E+01/
      data dre( 7,43)/ 3.644E+01/
      data dre( 8,43)/ 0.000E+00/
      data dre( 9,43)/ 8.467E+00/
      data dre(10,43)/ 1.879E+02/
      data dre(11,43)/ 2.114E+02/
      data dre(12,43)/ 4.138E+02/
      data dre(13,43)/ 5.737E+02/
      data dre(14,43)/ 7.472E+02/
      data dre(15,43)/ 2.816E+03/
      data dre(16,43)/ 3.650E+03/
      data drc( 1,44)/ 4.393E-01/
      data drc( 2,44)/ 5.066E-01/
      data drc( 3,44)/ 2.672E+01/
      data drc( 4,44)/ 9.813E+01/
      data drc( 5,44)/ 4.562E+02/
      data drc( 6,44)/ 0.000E+00/
      data drc( 7,44)/ 0.000E+00/
      data drc( 8,44)/ 0.000E+00/
      data drc( 9,44)/ 1.181E+01/
      data drc(10,44)/ 1.690E+02/
      data drc(11,44)/ 1.086E+02/
      data drc(12,44)/ 1.237E+03/
      data drc(13,44)/ 1.205E+03/
      data drc(14,44)/ 4.495E+02/
      data drc(15,44)/ 1.050E+03/
      data drc(16,44)/ 0.000E+00/
      data dre( 1,44)/ 8.471E-02/
      data dre( 2,44)/ 5.363E-01/
      data dre( 3,44)/ 1.940E+00/
      data dre( 4,44)/ 1.029E+01/
      data dre( 5,44)/ 4.632E+01/
      data dre( 6,44)/ 0.000E+00/
      data dre( 7,44)/ 0.000E+00/
      data dre( 8,44)/ 0.000E+00/
      data dre( 9,44)/ 1.851E+02/
      data dre(10,44)/ 2.180E+02/
      data dre(11,44)/ 3.188E+02/
      data dre(12,44)/ 4.900E+02/
      data dre(13,44)/ 6.375E+02/
      data dre(14,44)/ 2.725E+03/
      data dre(15,44)/ 3.532E+03/
      data dre(16,44)/ 0.000E+00/
      data drc( 1,45)/ 3.898E-01/
      data drc( 2,45)/ 3.520E-01/
      data drc( 3,45)/ 1.249E+00/
      data drc( 4,45)/ 9.695E+00/
      data drc( 5,45)/ 1.233E+01/
      data drc( 6,45)/ 5.694E+01/
      data drc( 7,45)/ 1.266E+02/
      data drc( 8,45)/ 4.327E+02/
      data drc( 9,45)/ 7.030E+00/
      data drc(10,45)/ 2.343E+02/
      data drc(11,45)/ 2.104E+02/
      data drc(12,45)/ 2.149E+03/
      data drc(13,45)/ 2.523E+03/
      data drc(14,45)/ 1.539E+02/
      data drc(15,45)/ 8.513E+02/
      data drc(16,45)/ 0.000E+00/
      data dre( 1,45)/ 3.094E-01/
      data dre( 2,45)/ 5.315E-01/
      data dre( 3,45)/ 1.021E+00/
      data dre( 4,45)/ 2.521E+00/
      data dre( 5,45)/ 5.289E+00/
      data dre( 6,45)/ 1.222E+01/
      data dre( 7,45)/ 3.228E+01/
      data dre( 8,45)/ 6.752E+01/
      data dre( 9,45)/ 1.662E+02/
      data dre(10,45)/ 2.084E+02/
      data dre(11,45)/ 2.663E+02/
      data dre(12,45)/ 4.425E+02/
      data dre(13,45)/ 5.715E+02/
      data dre(14,45)/ 2.480E+03/
      data dre(15,45)/ 3.480E+03/
      data dre(16,45)/ 0.000E+00/
      data drc( 1,46)/ 5.545E-02/
      data drc( 2,46)/ 6.068E-01/
      data drc( 3,46)/ 1.363E+00/
      data drc( 4,46)/ 6.528E+00/
      data drc( 5,46)/ 6.888E+00/
      data drc( 6,46)/ 2.570E+01/
      data drc( 7,46)/ 7.716E+01/
      data drc( 8,46)/ 3.930E+02/
      data drc( 9,46)/ 1.569E+01/
      data drc(10,46)/ 4.173E+02/
      data drc(11,46)/ 8.374E+02/
      data drc(12,46)/ 4.092E+03/
      data drc(13,46)/ 5.437E+02/
      data drc(14,46)/ 4.946E+02/
      data drc(15,46)/ 0.000E+00/
      data drc(16,46)/ 0.000E+00/
      data dre( 1,46)/ 8.999E-02/
      data dre( 2,46)/ 1.578E-01/
      data dre( 3,46)/ 3.553E-01/
      data dre( 4,46)/ 9.498E-01/
      data dre( 5,46)/ 1.918E+00/
      data dre( 6,46)/ 6.880E+00/
      data dre( 7,46)/ 2.512E+01/
      data dre( 8,46)/ 6.881E+01/
      data dre( 9,46)/ 1.689E+02/
      data dre(10,46)/ 2.193E+02/
      data dre(11,46)/ 3.455E+02/
      data dre(12,46)/ 4.810E+02/
      data dre(13,46)/ 6.530E+02/
      data dre(14,46)/ 3.390E+03/
      data dre(15,46)/ 0.000E+00/
      data dre(16,46)/ 0.000E+00/
      data drc( 1,47)/ 1.989E+00/
      data drc( 2,47)/ 4.823E+00/
      data drc( 3,47)/ 1.569E+00/
      data drc( 4,47)/ 1.811E+01/
      data drc( 5,47)/ 4.395E+01/
      data drc( 6,47)/ 3.689E+02/
      data drc( 7,47)/ 0.000E+00/
      data drc( 8,47)/ 0.000E+00/
      data drc( 9,47)/ 1.192E+01/
      data drc(10,47)/ 3.450E+02/
      data drc(11,47)/ 5.330E+02/
      data drc(12,47)/ 3.630E+03/
      data drc(13,47)/ 1.092E+03/
      data drc(14,47)/ 3.838E+02/
      data drc(15,47)/ 0.000E+00/
      data drc(16,47)/ 0.000E+00/
      data dre( 1,47)/ 3.629E-02/
      data dre( 2,47)/ 7.742E-02/
      data dre( 3,47)/ 3.552E-01/
      data dre( 4,47)/ 2.241E+00/
      data dre( 5,47)/ 1.013E+01/
      data dre( 6,47)/ 5.460E+01/
      data dre( 7,47)/ 0.000E+00/
      data dre( 8,47)/ 0.000E+00/
      data dre( 9,47)/ 1.675E+02/
      data dre(10,47)/ 2.161E+02/
      data dre(11,47)/ 2.967E+02/
      data dre(12,47)/ 4.335E+02/
      data dre(13,47)/ 5.592E+02/
      data dre(14,47)/ 3.631E+03/
      data dre(15,47)/ 0.000E+00/
      data dre(16,47)/ 0.000E+00/
      data drc( 1,48)/ 4.445E-03/
      data drc( 2,48)/ 8.783E-02/
      data drc( 3,48)/ 1.649E+00/
      data drc( 4,48)/ 4.602E+00/
      data drc( 5,48)/ 2.725E+00/
      data drc( 6,48)/ 2.501E+00/
      data drc( 7,48)/ 2.235E+01/
      data drc( 8,48)/ 2.565E+02/
      data drc( 9,48)/ 3.155E+00/
      data drc(10,48)/ 1.528E+02/
      data drc(11,48)/ 5.445E+02/
      data drc(12,48)/ 3.125E+03/
      data drc(13,48)/ 1.200E+03/
      data drc(14,48)/ 2.195E+02/
      data drc(15,48)/ 0.000E+00/
      data drc(16,48)/ 0.000E+00/
      data dre( 1,48)/ 6.155E-02/
      data dre( 2,48)/ 1.306E-01/
      data dre( 3,48)/ 4.261E-01/
      data dre( 4,48)/ 6.894E-01/
      data dre( 5,48)/ 2.015E+00/
      data dre( 6,48)/ 6.533E+00/
      data dre( 7,48)/ 2.658E+01/
      data dre( 8,48)/ 7.151E+01/
      data dre( 9,48)/ 1.392E+02/
      data dre(10,48)/ 1.945E+02/
      data dre(11,48)/ 2.545E+02/
      data dre(12,48)/ 3.909E+02/
      data dre(13,48)/ 5.086E+02/
      data dre(14,48)/ 3.617E+03/
      data dre(15,48)/ 0.000E+00/
      data dre(16,48)/ 0.000E+00/
      data drc( 1,49)/ 4.780E-03/
      data drc( 2,49)/ 7.510E-02/
      data drc( 3,49)/ 1.180E-01/
      data drc( 4,49)/ 2.574E-01/
      data drc( 5,49)/ 5.720E-01/
      data drc( 6,49)/ 3.289E+00/
      data drc( 7,49)/ 1.107E+01/
      data drc( 8,49)/ 1.393E+02/
      data drc( 9,49)/ 2.359E+00/
      data drc(10,49)/ 1.079E+02/
      data drc(11,49)/ 5.016E+02/
      data drc(12,49)/ 2.307E+03/
      data drc(13,49)/ 1.008E+03/
      data drc(14,49)/ 1.036E+02/
      data drc(15,49)/ 0.000E+00/
      data drc(16,49)/ 0.000E+00/
      data dre( 1,49)/ 2.168E-01/
      data dre( 2,49)/ 3.170E-01/
      data dre( 3,49)/ 6.789E-01/
      data dre( 4,49)/ 1.429E+00/
      data dre( 5,49)/ 2.919E+00/
      data dre( 6,49)/ 1.793E+01/
      data dre( 7,49)/ 4.345E+01/
      data dre( 8,49)/ 8.305E+01/
      data dre( 9,49)/ 1.325E+02/
      data dre(10,49)/ 1.841E+02/
      data dre(11,49)/ 2.417E+02/
      data dre(12,49)/ 3.617E+02/
      data dre(13,49)/ 4.682E+02/
      data dre(14,49)/ 3.558E+03/
      data dre(15,49)/ 0.000E+00/
      data dre(16,49)/ 0.000E+00/
      data drc( 1,50)/ 0.000E+00/
      data drc( 2,50)/ 0.000E+00/
      data drc( 3,50)/ 0.000E+00/
      data drc( 4,50)/ 0.000E+00/
      data drc( 5,50)/ 0.000E+00/
      data drc( 6,50)/ 0.000E+00/
      data drc( 7,50)/ 0.000E+00/
      data drc( 8,50)/ 0.000E+00/
      data drc( 9,50)/ 1.231E-01/
      data drc(10,50)/ 6.652E+00/
      data drc(11,50)/ 2.173E+02/
      data drc(12,50)/ 6.010E+02/
      data drc(13,50)/ 1.979E+03/
      data drc(14,50)/ 4.831E+02/
      data drc(15,50)/ 6.177E+01/
      data drc(16,50)/ 0.000E+00/
      data dre( 1,50)/ 0.000E+00/
      data dre( 2,50)/ 0.000E+00/
      data dre( 3,50)/ 0.000E+00/
      data dre( 4,50)/ 0.000E+00/
      data dre( 5,50)/ 0.000E+00/
      data dre( 6,50)/ 0.000E+00/
      data dre( 7,50)/ 0.000E+00/
      data dre( 8,50)/ 0.000E+00/
      data dre( 9,50)/ 1.143E+02/
      data dre(10,50)/ 1.422E+02/
      data dre(11,50)/ 1.987E+02/
      data dre(12,50)/ 2.697E+02/
      data dre(13,50)/ 3.609E+02/
      data dre(14,50)/ 4.560E+02/
      data dre(15,50)/ 3.990E+03/
      data dre(16,50)/ 0.000E+00/
      data drc( 1,51)/ 0.000E+00/
      data drc( 2,51)/ 0.000E+00/
      data drc( 3,51)/ 0.000E+00/
      data drc( 4,51)/ 0.000E+00/
      data drc( 5,51)/ 0.000E+00/
      data drc( 6,51)/ 0.000E+00/
      data drc( 7,51)/ 0.000E+00/
      data drc( 8,51)/ 0.000E+00/
      data drc( 9,51)/ 1.250E+03/
      data drc(10,51)/ 1.879E+03/
      data drc(11,51)/ 4.997E+02/
      data drc(12,51)/ 0.000E+00/
      data drc(13,51)/ 0.000E+00/
      data drc(14,51)/ 0.000E+00/
      data drc(15,51)/ 0.000E+00/
      data drc(16,51)/ 0.000E+00/
      data dre( 1,51)/ 0.000E+00/
      data dre( 2,51)/ 0.000E+00/
      data dre( 3,51)/ 0.000E+00/
      data dre( 4,51)/ 0.000E+00/
      data dre( 5,51)/ 0.000E+00/
      data dre( 6,51)/ 0.000E+00/
      data dre( 7,51)/ 0.000E+00/
      data dre( 8,51)/ 0.000E+00/
      data dre( 9,51)/ 4.780E+03/
      data dre(10,51)/ 6.225E+03/
      data dre(11,51)/ 7.425E+03/
      data dre(12,51)/ 0.000E+00/
      data dre(13,51)/ 0.000E+00/
      data dre(14,51)/ 0.000E+00/
      data dre(15,51)/ 0.000E+00/
      data dre(16,51)/ 0.000E+00/
      data drc( 1,52)/ 0.000E+00/
      data drc( 2,52)/ 0.000E+00/
      data drc( 3,52)/ 0.000E+00/
      data drc( 4,52)/ 0.000E+00/
      data drc( 5,52)/ 0.000E+00/
      data drc( 6,52)/ 0.000E+00/
      data drc( 7,52)/ 0.000E+00/
      data drc( 8,52)/ 0.000E+00/
      data drc( 9,52)/ 2.172E+03/
      data drc(10,52)/ 2.591E+03/
      data drc(11,52)/ 1.592E+03/
      data drc(12,52)/ 0.000E+00/
      data drc(13,52)/ 0.000E+00/
      data drc(14,52)/ 0.000E+00/
      data drc(15,52)/ 0.000E+00/
      data drc(16,52)/ 0.000E+00/
      data dre( 1,52)/ 0.000E+00/
      data dre( 2,52)/ 0.000E+00/
      data dre( 3,52)/ 0.000E+00/
      data dre( 4,52)/ 0.000E+00/
      data dre( 5,52)/ 0.000E+00/
      data dre( 6,52)/ 0.000E+00/
      data dre( 7,52)/ 0.000E+00/
      data dre( 8,52)/ 0.000E+00/
      data dre( 9,52)/ 4.641E+03/
      data dre(10,52)/ 5.797E+03/
      data dre(11,52)/ 6.740E+03/
      data dre(12,52)/ 0.000E+00/
      data dre(13,52)/ 0.000E+00/
      data dre(14,52)/ 0.000E+00/
      data dre(15,52)/ 0.000E+00/
      data dre(16,52)/ 0.000E+00/
      data drc( 1,53)/ 1.818E+01/
      data drc( 2,53)/ 4.192E+01/
      data drc( 3,53)/ 1.699E+01/
      data drc( 4,53)/ 7.392E+01/
      data drc( 5,53)/ 1.527E+02/
      data drc( 6,53)/ 2.267E+02/
      data drc( 7,53)/ 0.000E+00/
      data drc( 8,53)/ 0.000E+00/
      data drc( 9,53)/ 2.601E-09/
      data drc(10,53)/ 1.389E+02/
      data drc(11,53)/ 1.839E+02/
      data drc(12,53)/ 7.362E+02/
      data drc(13,53)/ 2.498E+03/
      data drc(14,53)/ 5.327E+02/
      data drc(15,53)/ 1.958E+03/
      data drc(16,53)/ 2.239E+03/
      data dre( 1,53)/ 4.701E+00/
      data dre( 2,53)/ 5.360E+00/
      data dre( 3,53)/ 1.077E+01/
      data dre( 4,53)/ 1.744E+01/
      data dre( 5,53)/ 3.465E+01/
      data dre( 6,53)/ 5.778E+01/
      data dre( 7,53)/ 0.000E+00/
      data dre( 8,53)/ 0.000E+00/
      data dre( 9,53)/ 1.176E+01/
      data dre(10,53)/ 3.108E+02/
      data dre(11,53)/ 3.470E+02/
      data dre(12,53)/ 7.126E+02/
      data dre(13,53)/ 9.973E+02/
      data dre(14,53)/ 1.323E+03/
      data dre(15,53)/ 4.780E+03/
      data dre(16,53)/ 6.220E+03/
      data drc( 1,54)/ 6.664E-01/
      data drc( 2,54)/ 6.977E-01/
      data drc( 3,54)/ 1.681E+01/
      data drc( 4,54)/ 1.282E+01/
      data drc( 5,54)/ 4.379E+01/
      data drc( 6,54)/ 1.420E+02/
      data drc( 7,54)/ 2.644E+02/
      data drc( 8,54)/ 6.321E+02/
      data drc( 9,54)/ 2.946E+01/
      data drc(10,54)/ 5.203E+02/
      data drc(11,54)/ 1.905E+02/
      data drc(12,54)/ 2.580E+03/
      data drc(13,54)/ 3.116E+03/
      data drc(14,54)/ 7.440E+02/
      data drc(15,54)/ 1.972E+03/
      data drc(16,54)/ 0.000E+00/
      data dre( 1,54)/ 4.425E-01/
      data dre( 2,54)/ 6.601E-01/
      data dre( 3,54)/ 1.185E+00/
      data dre( 4,54)/ 3.512E+00/
      data dre( 5,54)/ 8.418E+00/
      data dre( 6,54)/ 1.475E+01/
      data dre( 7,54)/ 3.681E+01/
      data dre( 8,54)/ 8.463E+01/
      data dre( 9,54)/ 3.020E+02/
      data dre(10,54)/ 3.533E+02/
      data dre(11,54)/ 4.946E+02/
      data dre(12,54)/ 8.116E+02/
      data dre(13,54)/ 1.101E+03/
      data dre(14,54)/ 4.374E+03/
      data dre(15,54)/ 5.933E+03/
      data dre(16,54)/ 0.000E+00/
      data drc( 1,55)/ 1.591E+00/
      data drc( 2,55)/ 1.542E+01/
      data drc( 3,55)/ 4.739E+01/
      data drc( 4,55)/ 6.388E+00/
      data drc( 5,55)/ 2.709E+01/
      data drc( 6,55)/ 1.351E+02/
      data drc( 7,55)/ 2.677E+02/
      data drc( 8,55)/ 8.320E+02/
      data drc( 9,55)/ 1.212E-09/
      data drc(10,55)/ 8.360E+01/
      data drc(11,55)/ 1.064E+03/
      data drc(12,55)/ 1.351E+03/
      data drc(13,55)/ 5.809E+03/
      data drc(14,55)/ 2.145E+03/
      data drc(15,55)/ 1.112E+03/
      data drc(16,55)/ 1.947E+03/
      data dre( 1,55)/ 2.099E-01/
      data dre( 2,55)/ 2.933E-01/
      data dre( 3,55)/ 5.167E-01/
      data dre( 4,55)/ 1.268E+00/
      data dre( 5,55)/ 4.111E+00/
      data dre( 6,55)/ 1.192E+01/
      data dre( 7,55)/ 3.514E+01/
      data dre( 8,55)/ 9.096E+01/
      data dre( 9,55)/ 1.546E+01/
      data dre(10,55)/ 2.985E+02/
      data dre(11,55)/ 3.624E+02/
      data dre(12,55)/ 6.156E+02/
      data dre(13,55)/ 8.712E+02/
      data dre(14,55)/ 9.665E+02/
      data dre(15,55)/ 1.253E+03/
      data dre(16,55)/ 5.631E+03/
      data drc( 1,56)/ 1.238E+00/
      data drc( 2,56)/ 4.916E+00/
      data drc( 3,56)/ 7.794E+00/
      data drc( 4,56)/ 6.312E+01/
      data drc( 5,56)/ 1.029E+02/
      data drc( 6,56)/ 1.423E+02/
      data drc( 7,56)/ 7.231E+02/
      data drc( 8,56)/ 0.000E+00/
      data drc( 9,56)/ 7.380E-01/
      data drc(10,56)/ 2.629E+02/
      data drc(11,56)/ 1.283E+03/
      data drc(12,56)/ 2.977E+03/
      data drc(13,56)/ 8.459E+03/
      data drc(14,56)/ 7.561E+02/
      data drc(15,56)/ 1.039E+03/
      data drc(16,56)/ 0.000E+00/
      data dre( 1,56)/ 1.006E-01/
      data dre( 2,56)/ 2.570E-01/
      data dre( 3,56)/ 6.130E-01/
      data dre( 4,56)/ 2.380E+00/
      data dre( 5,56)/ 5.742E+00/
      data dre( 6,56)/ 2.005E+01/
      data dre( 7,56)/ 8.853E+01/
      data dre( 8,56)/ 0.000E+00/
      data dre( 9,56)/ 2.355E+02/
      data dre(10,56)/ 3.131E+02/
      data dre(11,56)/ 3.749E+02/
      data dre(12,56)/ 6.582E+02/
      data dre(13,56)/ 9.025E+02/
      data dre(14,56)/ 1.279E+03/
      data dre(15,56)/ 5.776E+03/
      data dre(16,56)/ 0.000E+00/
      data drc( 1,57)/ 6.791E-01/
      data drc( 2,57)/ 1.432E+00/
      data drc( 3,57)/ 3.880E+00/
      data drc( 4,57)/ 1.939E+01/
      data drc( 5,57)/ 3.557E+01/
      data drc( 6,57)/ 1.475E+02/
      data drc( 7,57)/ 6.744E+02/
      data drc( 8,57)/ 0.000E+00/
      data drc( 9,57)/ 5.390E+01/
      data drc(10,57)/ 1.426E+03/
      data drc(11,57)/ 1.184E+03/
      data drc(12,57)/ 8.668E+03/
      data drc(13,57)/ 4.558E+03/
      data drc(14,57)/ 4.686E+02/
      data drc(15,57)/ 0.000E+00/
      data drc(16,57)/ 0.000E+00/
      data dre( 1,57)/ 4.477E-03/
      data dre( 2,57)/ 4.785E-02/
      data dre( 3,57)/ 2.833E-01/
      data dre( 4,57)/ 1.088E+00/
      data dre( 5,57)/ 3.041E+00/
      data dre( 6,57)/ 1.498E+01/
      data dre( 7,57)/ 7.697E+01/
      data dre( 8,57)/ 0.000E+00/
      data dre( 9,57)/ 2.812E+02/
      data dre(10,57)/ 3.610E+02/
      data dre(11,57)/ 4.937E+02/
      data dre(12,57)/ 7.560E+02/
      data dre(13,57)/ 9.965E+02/
      data dre(14,57)/ 6.087E+03/
      data dre(15,57)/ 0.000E+00/
      data dre(16,57)/ 0.000E+00/
      data drc( 1,58)/ 1.145E-01/
      data drc( 2,58)/ 4.176E-01/
      data drc( 3,58)/ 2.169E-01/
      data drc( 4,58)/ 3.611E+00/
      data drc( 5,58)/ 1.170E+01/
      data drc( 6,58)/ 5.369E+01/
      data drc( 7,58)/ 4.863E+02/
      data drc( 8,58)/ 0.000E+00/
      data drc( 9,58)/ 2.262E+01/
      data drc(10,58)/ 1.011E+03/
      data drc(11,58)/ 1.514E+03/
      data drc(12,58)/ 8.610E+03/
      data drc(13,58)/ 4.861E+03/
      data drc(14,58)/ 3.034E+02/
      data drc(15,58)/ 0.000E+00/
      data drc(16,58)/ 0.000E+00/
      data dre( 1,58)/ 6.682E-03/
      data dre( 2,58)/ 6.191E-02/
      data dre( 3,58)/ 2.301E-01/
      data dre( 4,58)/ 1.291E+00/
      data dre( 5,58)/ 4.986E+00/
      data dre( 6,58)/ 2.301E+01/
      data dre( 7,58)/ 9.812E+01/
      data dre( 8,58)/ 0.000E+00/
      data dre( 9,58)/ 2.516E+02/
      data dre(10,58)/ 3.424E+02/
      data dre(11,58)/ 4.388E+02/
      data dre(12,58)/ 7.093E+02/
      data dre(13,58)/ 9.376E+02/
      data dre(14,58)/ 5.988E+03/
      data dre(15,58)/ 0.000E+00/
      data dre(16,58)/ 0.000E+00/
      data drc( 1,59)/ 1.565E-02/
      data drc( 2,59)/ 1.103E-01/
      data drc( 3,59)/ 5.213E-01/
      data drc( 4,59)/ 9.714E-01/
      data drc( 5,59)/ 2.931E+00/
      data drc( 6,59)/ 1.298E+01/
      data drc( 7,59)/ 2.662E+01/
      data drc( 8,59)/ 2.600E+02/
      data drc( 9,59)/ 1.708E+01/
      data drc(10,59)/ 8.663E+02/
      data drc(11,59)/ 1.657E+03/
      data drc(12,59)/ 8.478E+03/
      data drc(13,59)/ 3.984E+03/
      data drc(14,59)/ 1.548E+02/
      data drc(15,59)/ 0.000E+00/
      data drc(16,59)/ 0.000E+00/
      data dre( 1,59)/ 2.014E-01/
      data dre( 2,59)/ 3.854E-01/
      data dre( 3,59)/ 5.807E-01/
      data dre( 4,59)/ 2.135E+00/
      data dre( 5,59)/ 5.978E+00/
      data dre( 6,59)/ 1.769E+01/
      data dre( 7,59)/ 5.147E+01/
      data dre( 8,59)/ 1.206E+02/
      data dre( 9,59)/ 2.417E+02/
      data dre(10,59)/ 3.356E+02/
      data dre(11,59)/ 4.346E+02/
      data dre(12,59)/ 6.875E+02/
      data dre(13,59)/ 9.060E+02/
      data dre(14,59)/ 6.017E+03/
      data dre(15,59)/ 0.000E+00/
      data dre(16,59)/ 0.000E+00/
      data drc( 1,60)/ 0.000E+00/
      data drc( 2,60)/ 0.000E+00/
      data drc( 3,60)/ 0.000E+00/
      data drc( 4,60)/ 0.000E+00/
      data drc( 5,60)/ 0.000E+00/
      data drc( 6,60)/ 0.000E+00/
      data drc( 7,60)/ 0.000E+00/
      data drc( 8,60)/ 0.000E+00/
      data drc( 9,60)/ 9.881E+00/
      data drc(10,60)/ 4.677E+02/
      data drc(11,60)/ 1.404E+03/
      data drc(12,60)/ 7.589E+03/
      data drc(13,60)/ 3.804E+03/
      data drc(14,60)/ 1.590E+02/
      data drc(15,60)/ 0.000E+00/
      data drc(16,60)/ 0.000E+00/
      data dre( 1,60)/ 0.000E+00/
      data dre( 2,60)/ 0.000E+00/
      data dre( 3,60)/ 0.000E+00/
      data dre( 4,60)/ 0.000E+00/
      data dre( 5,60)/ 0.000E+00/
      data dre( 6,60)/ 0.000E+00/
      data dre( 7,60)/ 0.000E+00/
      data dre( 8,60)/ 0.000E+00/
      data dre( 9,60)/ 2.307E+02/
      data dre(10,60)/ 3.245E+02/
      data dre(11,60)/ 4.177E+02/
      data dre(12,60)/ 6.600E+02/
      data dre(13,60)/ 8.621E+02/
      data dre(14,60)/ 6.898E+03/
      data dre(15,60)/ 0.000E+00/
      data dre(16,60)/ 0.000E+00/
      data drc( 1,61)/ 0.000E+00/
      data drc( 2,61)/ 0.000E+00/
      data drc( 3,61)/ 0.000E+00/
      data drc( 4,61)/ 0.000E+00/
      data drc( 5,61)/ 0.000E+00/
      data drc( 6,61)/ 0.000E+00/
      data drc( 7,61)/ 0.000E+00/
      data drc( 8,61)/ 0.000E+00/
      data drc( 9,61)/ 1.376E+03/
      data drc(10,61)/ 1.984E+03/
      data drc(11,61)/ 5.195E+02/
      data drc(12,61)/ 0.000E+00/
      data drc(13,61)/ 0.000E+00/
      data drc(14,61)/ 0.000E+00/
      data drc(15,61)/ 0.000E+00/
      data drc(16,61)/ 0.000E+00/
      data dre( 1,61)/ 0.000E+00/
      data dre( 2,61)/ 0.000E+00/
      data dre( 3,61)/ 0.000E+00/
      data dre( 4,61)/ 0.000E+00/
      data dre( 5,61)/ 0.000E+00/
      data dre( 6,61)/ 0.000E+00/
      data dre( 7,61)/ 0.000E+00/
      data dre( 8,61)/ 0.000E+00/
      data dre( 9,61)/ 5.543E+03/
      data dre(10,61)/ 7.223E+03/
      data dre(11,61)/ 8.648E+03/
      data dre(12,61)/ 0.000E+00/
      data dre(13,61)/ 0.000E+00/
      data dre(14,61)/ 0.000E+00/
      data dre(15,61)/ 0.000E+00/
      data dre(16,61)/ 0.000E+00/
      data drc( 1,62)/ 0.000E+00/
      data drc( 2,62)/ 0.000E+00/
      data drc( 3,62)/ 0.000E+00/
      data drc( 4,62)/ 0.000E+00/
      data drc( 5,62)/ 0.000E+00/
      data drc( 6,62)/ 0.000E+00/
      data drc( 7,62)/ 0.000E+00/
      data drc( 8,62)/ 0.000E+00/
      data drc( 9,62)/ 2.470E+03/
      data drc(10,62)/ 3.116E+03/
      data drc(11,62)/ 1.301E+03/
      data drc(12,62)/ 0.000E+00/
      data drc(13,62)/ 0.000E+00/
      data drc(14,62)/ 0.000E+00/
      data drc(15,62)/ 0.000E+00/
      data drc(16,62)/ 0.000E+00/
      data dre( 1,62)/ 0.000E+00/
      data dre( 2,62)/ 0.000E+00/
      data dre( 3,62)/ 0.000E+00/
      data dre( 4,62)/ 0.000E+00/
      data dre( 5,62)/ 0.000E+00/
      data dre( 6,62)/ 0.000E+00/
      data dre( 7,62)/ 0.000E+00/
      data dre( 8,62)/ 0.000E+00/
      data dre( 9,62)/ 5.399E+03/
      data dre(10,62)/ 6.845E+03/
      data dre(11,62)/ 8.000E+03/
      data dre(12,62)/ 0.000E+00/
      data dre(13,62)/ 0.000E+00/
      data dre(14,62)/ 0.000E+00/
      data dre(15,62)/ 0.000E+00/
      data dre(16,62)/ 0.000E+00/
      data drc( 1,63)/ 1.784E+00/
      data drc( 2,63)/ 1.892E+01/
      data drc( 3,63)/ 6.022E+01/
      data drc( 4,63)/ 2.707E+01/
      data drc( 5,63)/ 9.603E+01/
      data drc( 6,63)/ 1.955E+02/
      data drc( 7,63)/ 2.647E+02/
      data drc( 8,63)/ 0.000E+00/
      data drc( 9,63)/ 0.000E+00/
      data drc(10,63)/ 1.347E+02/
      data drc(11,63)/ 2.834E+02/
      data drc(12,63)/ 4.922E+02/
      data drc(13,63)/ 2.246E+03/
      data drc(14,63)/ 1.654E+03/
      data drc(15,63)/ 1.946E+03/
      data drc(16,63)/ 2.696E+03/
      data dre( 1,63)/ 2.222E+00/
      data dre( 2,63)/ 2.621E+00/
      data dre( 3,63)/ 4.584E+00/
      data dre( 4,63)/ 1.058E+01/
      data dre( 5,63)/ 1.925E+01/
      data dre( 6,63)/ 3.921E+01/
      data dre( 7,63)/ 6.718E+01/
      data dre( 8,63)/ 0.000E+00/
      data dre( 9,63)/ 1.664E+01/
      data dre(10,63)/ 3.555E+02/
      data dre(11,63)/ 3.943E+02/
      data dre(12,63)/ 7.744E+02/
      data dre(13,63)/ 1.056E+03/
      data dre(14,63)/ 1.381E+03/
      data dre(15,63)/ 5.418E+03/
      data dre(16,63)/ 7.142E+03/
      data drc( 1,64)/ 2.383E-01/
      data drc( 2,64)/ 2.368E-01/
      data drc( 3,64)/ 2.015E+00/
      data drc( 4,64)/ 2.108E+01/
      data drc( 5,64)/ 1.629E+02/
      data drc( 6,64)/ 2.992E+02/
      data drc( 7,64)/ 8.732E+02/
      data drc( 8,64)/ 0.000E+00/
      data drc( 9,64)/ 8.802E-09/
      data drc(10,64)/ 1.822E+02/
      data drc(11,64)/ 6.337E+02/
      data drc(12,64)/ 6.578E+02/
      data drc(13,64)/ 3.661E+03/
      data drc(14,64)/ 2.714E+03/
      data drc(15,64)/ 1.047E+03/
      data drc(16,64)/ 2.039E+03/
      data dre( 1,64)/ 4.223E-02/
      data dre( 2,64)/ 1.367E-01/
      data dre( 3,64)/ 4.029E-01/
      data dre( 4,64)/ 3.193E+00/
      data dre( 5,64)/ 8.043E+00/
      data dre( 6,64)/ 2.770E+01/
      data dre( 7,64)/ 8.589E+01/
      data dre( 8,64)/ 0.000E+00/
      data dre( 9,64)/ 1.733E+01/
      data dre(10,64)/ 3.705E+02/
      data dre(11,64)/ 4.233E+02/
      data dre(12,64)/ 7.436E+02/
      data dre(13,64)/ 1.023E+03/
      data dre(14,64)/ 1.342E+03/
      data dre(15,64)/ 5.262E+03/
      data dre(16,64)/ 6.984E+03/
      data drc( 1,65)/ 1.571E-01/
      data drc( 2,65)/ 2.164E+00/
      data drc( 3,65)/ 4.367E+00/
      data drc( 4,65)/ 1.126E+01/
      data drc( 5,65)/ 7.279E+01/
      data drc( 6,65)/ 1.403E+02/
      data drc( 7,65)/ 3.204E+02/
      data drc( 8,65)/ 1.007E+03/
      data drc( 9,65)/ 3.008E-09/
      data drc(10,65)/ 1.226E+02/
      data drc(11,65)/ 1.451E+03/
      data drc(12,65)/ 1.644E+03/
      data drc(13,65)/ 6.759E+03/
      data drc(14,65)/ 2.589E+03/
      data drc(15,65)/ 1.124E+03/
      data drc(16,65)/ 2.268E+03/
      data dre( 1,65)/ 1.099E-01/
      data dre( 2,65)/ 2.807E-01/
      data dre( 3,65)/ 4.580E-01/
      data dre( 4,65)/ 1.760E+00/
      data dre( 5,65)/ 4.884E+00/
      data dre( 6,65)/ 1.308E+01/
      data dre( 7,65)/ 3.570E+01/
      data dre( 8,65)/ 9.916E+01/
      data dre( 9,65)/ 1.861E+01/
      data dre(10,65)/ 3.474E+02/
      data dre(11,65)/ 4.187E+02/
      data dre(12,65)/ 7.225E+02/
      data dre(13,65)/ 1.019E+03/
      data dre(14,65)/ 1.152E+03/
      data dre(15,65)/ 1.501E+03/
      data dre(16,65)/ 6.541E+03/
      data drc( 1,66)/ 1.587E-03/
      data drc( 2,66)/ 3.183E-01/
      data drc( 3,66)/ 9.056E+00/
      data drc( 4,66)/ 1.879E+01/
      data drc( 5,66)/ 1.490E+02/
      data drc( 6,66)/ 1.884E+02/
      data drc( 7,66)/ 8.249E+02/
      data drc( 8,66)/ 0.000E+00/
      data drc( 9,66)/ 1.399E+02/
      data drc(10,66)/ 1.771E+03/
      data drc(11,66)/ 2.003E+03/
      data drc(12,66)/ 1.023E+04/
      data drc(13,66)/ 2.315E+03/
      data drc(14,66)/ 1.245E+03/
      data drc(15,66)/ 0.000E+00/
      data drc(16,66)/ 0.000E+00/
      data dre( 1,66)/ 1.107E-03/
      data dre( 2,66)/ 3.114E-01/
      data dre( 3,66)/ 1.035E+00/
      data dre( 4,66)/ 2.266E+00/
      data dre( 5,66)/ 7.157E+00/
      data dre( 6,66)/ 2.533E+01/
      data dre( 7,66)/ 9.959E+01/
      data dre( 8,66)/ 0.000E+00/
      data dre( 9,66)/ 3.424E+02/
      data dre(10,66)/ 4.195E+02/
      data dre(11,66)/ 6.973E+02/
      data dre(12,66)/ 1.005E+03/
      data dre(13,66)/ 1.356E+03/
      data dre(14,66)/ 6.698E+03/
      data dre(15,66)/ 0.000E+00/
      data dre(16,66)/ 0.000E+00/
      data drc( 1,67)/ 2.202E+00/
      data drc( 2,67)/ 6.011E-01/
      data drc( 3,67)/ 3.740E+00/
      data drc( 4,67)/ 3.309E+01/
      data drc( 5,67)/ 1.154E+02/
      data drc( 6,67)/ 1.902E+02/
      data drc( 7,67)/ 7.533E+02/
      data drc( 8,67)/ 0.000E+00/
      data drc( 9,67)/ 6.864E+01/
      data drc(10,67)/ 1.791E+03/
      data drc(11,67)/ 1.363E+03/
      data drc(12,67)/ 1.040E+04/
      data drc(13,67)/ 6.103E+03/
      data drc(14,67)/ 8.055E+02/
      data drc(15,67)/ 0.000E+00/
      data drc(16,67)/ 0.000E+00/
      data dre( 1,67)/ 4.306E-03/
      data dre( 2,67)/ 4.661E-01/
      data dre( 3,67)/ 1.095E+00/
      data dre( 4,67)/ 4.058E+00/
      data dre( 5,67)/ 9.375E+00/
      data dre( 6,67)/ 2.890E+01/
      data dre( 7,67)/ 9.788E+01/
      data dre( 8,67)/ 0.000E+00/
      data dre( 9,67)/ 3.220E+02/
      data dre(10,67)/ 4.130E+02/
      data dre(11,67)/ 5.572E+02/
      data dre(12,67)/ 8.794E+02/
      data dre(13,67)/ 1.170E+03/
      data dre(14,67)/ 7.119E+03/
      data dre(15,67)/ 0.000E+00/
      data dre(16,67)/ 0.000E+00/
      data drc( 1,68)/ 6.848E-03/
      data drc( 2,68)/ 1.859E-01/
      data drc( 3,68)/ 7.646E-01/
      data drc( 4,68)/ 4.585E+00/
      data drc( 5,68)/ 1.697E+01/
      data drc( 6,68)/ 9.391E+01/
      data drc( 7,68)/ 5.812E+02/
      data drc( 8,68)/ 0.000E+00/
      data drc( 9,68)/ 3.934E+01/
      data drc(10,68)/ 1.695E+03/
      data drc(11,68)/ 1.759E+03/
      data drc(12,68)/ 1.081E+04/
      data drc(13,68)/ 6.040E+03/
      data drc(14,68)/ 4.913E+02/
      data drc(15,68)/ 0.000E+00/
      data drc(16,68)/ 0.000E+00/
      data dre( 1,68)/ 2.488E-02/
      data dre( 2,68)/ 1.864E-01/
      data dre( 3,68)/ 5.843E-01/
      data dre( 4,68)/ 1.830E+00/
      data dre( 5,68)/ 6.573E+00/
      data dre( 6,68)/ 2.113E+01/
      data dre( 7,68)/ 1.083E+02/
      data dre( 8,68)/ 0.000E+00/
      data dre( 9,68)/ 2.967E+02/
      data dre(10,68)/ 4.011E+02/
      data dre(11,68)/ 5.265E+02/
      data dre(12,68)/ 8.393E+02/
      data dre(13,68)/ 1.116E+03/
      data dre(14,68)/ 7.041E+03/
      data dre(15,68)/ 0.000E+00/
      data dre(16,68)/ 0.000E+00/
      data drc( 1,69)/ 5.393E-01/
      data drc( 2,69)/ 1.191E+00/
      data drc( 3,69)/ 6.070E+00/
      data drc( 4,69)/ 3.602E+01/
      data drc( 5,69)/ 3.244E+02/
      data drc( 6,69)/ 0.000E+00/
      data drc( 7,69)/ 0.000E+00/
      data drc( 8,69)/ 0.000E+00/
      data drc( 9,69)/ 2.697E+01/
      data drc(10,69)/ 1.376E+03/
      data drc(11,69)/ 2.121E+03/
      data drc(12,69)/ 1.061E+04/
      data drc(13,69)/ 5.344E+03/
      data drc(14,69)/ 2.341E+02/
      data drc(15,69)/ 0.000E+00/
      data drc(16,69)/ 0.000E+00/
      data dre( 1,69)/ 9.250E-02/
      data dre( 2,69)/ 7.467E-01/
      data dre( 3,69)/ 4.174E+00/
      data dre( 4,69)/ 1.268E+01/
      data dre( 5,69)/ 1.189E+02/
      data dre( 6,69)/ 0.000E+00/
      data dre( 7,69)/ 0.000E+00/
      data dre( 8,69)/ 0.000E+00/
      data dre( 9,69)/ 2.837E+02/
      data dre(10,69)/ 3.926E+02/
      data dre(11,69)/ 5.066E+02/
      data dre(12,69)/ 8.109E+02/
      data dre(13,69)/ 1.078E+03/
      data dre(14,69)/ 7.034E+03/
      data dre(15,69)/ 0.000E+00/
      data dre(16,69)/ 0.000E+00/
      data drc( 1,70)/ 0.000E+00/
      data drc( 2,70)/ 0.000E+00/
      data drc( 3,70)/ 0.000E+00/
      data drc( 4,70)/ 0.000E+00/
      data drc( 5,70)/ 0.000E+00/
      data drc( 6,70)/ 0.000E+00/
      data drc( 7,70)/ 0.000E+00/
      data drc( 8,70)/ 0.000E+00/
      data drc( 9,70)/ 1.452E+01/
      data drc(10,70)/ 7.861E+02/
      data drc(11,70)/ 1.889E+03/
      data drc(12,70)/ 9.862E+03/
      data drc(13,70)/ 5.120E+03/
      data drc(14,70)/ 1.945E+02/
      data drc(15,70)/ 0.000E+00/
      data drc(16,70)/ 0.000E+00/
      data dre( 1,70)/ 0.000E+00/
      data dre( 2,70)/ 0.000E+00/
      data dre( 3,70)/ 0.000E+00/
      data dre( 4,70)/ 0.000E+00/
      data dre( 5,70)/ 0.000E+00/
      data dre( 6,70)/ 0.000E+00/
      data dre( 7,70)/ 0.000E+00/
      data dre( 8,70)/ 0.000E+00/
      data dre( 9,70)/ 2.700E+02/
      data dre(10,70)/ 3.810E+02/
      data dre(11,70)/ 4.861E+02/
      data dre(12,70)/ 7.828E+02/
      data dre(13,70)/ 1.033E+03/
      data dre(14,70)/ 8.024E+03/
      data dre(15,70)/ 0.000E+00/
      data dre(16,70)/ 0.000E+00/
      
      end
