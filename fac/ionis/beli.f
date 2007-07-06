      subroutine ebeli(z, n, e)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     return the ionization threshold of the ion
c     input:
c     z, nuclear charge of the element.
c     n, number of electrons of the ion.
c     output:
c     e, the ionization threshold in eV.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer z, n
      double precision e
      double precision error(406), cdi(7,406), cea(8,406)
      common /idata/ error, cdi, cea

      integer k

      k = 1 + (z*(z-1)/2) + (z-n)
      e = cdi(1, k)
      return
      end
      
      subroutine cbeli(z, n, ene, ea, dir, e)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculate the ionization cross sections.
c     Reference:
c     K. L. Bell et al. J. Phys. Chem. Ref. Data 12, 891 (1983)
c     input:
c     z, nuclear charge of the element.
c     n, number of electrons of the ion.
c     ene, energy where cross sections to be calcuated (eV)
c     output:
c     ea, excitation autoionization contribution (10^-20 cm2).
c     dir, direct ionization cross section (10^-20 cm2).
c     e, relative error estimate of the results in %
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer z, n
      double precision ene, dir, ea, e
      double precision error(406), cdi(7,406), cea(8,406)
      common /idata/ error, cdi, cea

      integer k, i
      double precision ei0, ei1, a0, a1, ea0
      double precision power, power1, x, x2, xs

      k = 1 + (z*(z-1)/2) + (z-n)
      ei0 = cdi(1, k)
      a0 = cdi(2, k)
      ea0 = cea(1, k)
      ei1 = cea(2, k)
      a1 = cea(3, k)

      e = error(k)

      if (1+ei0 .eq. 1) then
         ea = -1.0
         dir = -1.0
         return
      endif

      dir = 0.0
      if (1+ei0 .ne. 1 .and. ene .gt. ei0) then
         x = ei0/ene
         x2 = 1.0/x
         xs = a0*dlog(x2)
         power1 = 1.0 - x
         power = power1
         do i = 3, 7
            if (cdi(i, k) .eq. 0) goto 10
            xs = xs + cdi(i, k)*power
            power = power*power1
         enddo
 10      dir = 1E7*xs/(ene*ei0)
      endif

      ea = 0.0
      if (1+ea0 .ne. 1 .and. ene .gt. ea0) then
         x = ei1/ene
         x2 = 1.0/x
         xs = a1*dlog(x2)
         power1 = 1.0 - x
         power = power1
         do i = 4, 8
            if (cea(i, k) .eq. 0) goto 20
            xs = xs + cea(i, k)*power
            power = power * power1
         enddo
 20      ea = 1E7*xs/(ene*ei1) - dir
      endif

      return
      end

      subroutine rbeli(z, n, t, ea, dir)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculate the ionization rate coeff..
c     Reference:
c     K. L. Bell et al. J. Phys. Chem. Ref. Data 12, 891 (1983)
c     input:
c     z, nuclear charge of the element.
c     n, number of electrons of the ion.
c     t, temperature where rate coeff. to be calcuated (eV)
c     output:
c     ea, excitation autoionization contribution (10^-10 cm3/s).
c     dir, direct ionization rate coeff. (10^-10 cm3/s).
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer z, n
      double precision t, ea, dir
      double precision error(406), cdi(7,406), cea(8,406)
      common /idata/ error, cdi, cea
      integer ng
      parameter(ng=15)
      double precision xg(ng), wg(ng)
      data xg /.933078120172818E-01,
     +         .492691740301883E+00,
     +         .121559541207095E+01,
     +         .226994952620374E+01,
     +         .366762272175144E+01,
     +         .542533662741355E+01,
     +         .756591622661307E+01,
     +         .101202285680191E+02,
     +         .131302824821757E+02,
     +         .166544077083300E+02,
     +         .207764788994488E+02,
     +         .256238942267288E+02,
     +         .314075191697539E+02,
     +         .385306833064860E+02,
     +         .480260855726858E+02/
      data wg /.218234885940086E+00,
     +         .342210177922884E+00,
     +         .263027577941681E+00,
     +         .126425818105931E+00,
     +         .402068649210010E-01,
     +         .856387780361184E-02,
     +         .121243614721425E-02,
     +         .111674392344251E-03,
     +         .645992676202287E-05,
     +         .222631690709627E-06,
     +         .422743038497938E-08,
     +         .392189726704110E-10,
     +         .145651526407313E-12,
     +         .148302705111330E-15,
     +         .160059490621113E-19/
      integer k, i
      double precision ei0, a, b, c, x0, ene
      double precision mc
c     mc = 1.12837967*c*sqrt(2.0/mc^2)
      parameter(mc = 0.0066923847825)

      k = 1 + (z*(z-1)/2) + (z-n)
      ei0 = cdi(1, k)
      if (1+ei0 .eq. 1) then
         ea = -1.0
         dir = -1.0
         return
      endif

      ea = 0.0
      dir = 0.0

      x0 = ei0/t
      do i = 1, ng
         ene = t*xg(i) + ei0
         call cbeli(z, n, ene, a, b, c)
         ea = ea + wg(i)*a*(x0+xg(i))
         dir = dir + wg(i)*b*(x0+xg(i))
      enddo

      a = sqrt(t)
      ea = a * mc * ea * exp(-x0)
      dir = a * mc * dir * exp(-x0)
      
      return
      end

      block data belid
      double precision error(406), cdi(7,406), cea(8,406)
      common /idata/ error, cdi, cea
C     Z= 1 [+ 0]
      data error(  1) /  7.0/
      data cdi(1,  1) / 1.36000E+01/
      data cdi(2,  1) / 1.85000E-01/
      data cdi(3,  1) /-1.90000E-02/
      data cdi(4,  1) / 1.23000E-01/
      data cdi(5,  1) /-1.90000E-01/
      data cdi(6,  1) / 9.53000E-01/
C     Z= 2 [+ 0]
      data error(  2) /  5.0/
      data cdi(1,  2) / 2.46000E+01/
      data cdi(2,  2) / 5.72000E-01/
      data cdi(3,  2) /-3.44000E-01/
      data cdi(4,  2) /-5.23000E-01/
      data cdi(5,  2) / 3.44500E+00/
      data cdi(6,  2) /-6.82100E+00/
      data cdi(7,  2) / 5.57800E+00/
C     Z= 2 [+ 1]
      data error(  3) / 10.0/
      data cdi(1,  3) / 5.44000E+01/
      data cdi(2,  3) / 1.85000E-01/
      data cdi(3,  3) / 8.90000E-02/
      data cdi(4,  3) / 1.31000E-01/
      data cdi(5,  3) / 3.88000E-01/
      data cdi(6,  3) /-1.09100E+00/
      data cdi(7,  3) / 1.35400E+00/
C     Z= 3 [+ 0]
      data error(  4) / 10.0/
      data cdi(1,  4) / 5.40000E+00/
      data cdi(2,  4) / 8.50000E-02/
      data cdi(3,  4) /-4.00000E-03/
      data cdi(4,  4) / 7.57000E-01/
      data cdi(5,  4) /-1.78000E-01/
C     Z= 3 [+ 1]
      data error(  5) / 12.0/
      data cdi(1,  5) / 7.56000E+01/
      data cdi(2,  5) / 7.22000E-01/
      data cdi(3,  5) /-1.49000E-01/
      data cdi(4,  5) /-1.30100E+00/
      data cdi(5,  5) / 1.94400E+00/
C     Z= 3 [+ 2]
      data error(  6) / 10.0/
      data cdi(1,  6) / 1.22400E+02/
      data cdi(2,  6) / 4.00000E-01/
C     Z= 4 [+ 0]
      data error(  7) / 20.0/
      data cdi(1,  7) / 9.30000E+00/
      data cdi(2,  7) / 9.24000E-01/
      data cdi(3,  7) /-7.70000E-01/
      data cdi(4,  7) / 3.62000E-01/
C     Z= 4 [+ 1]
      data error(  8) / 20.0/
      data cdi(1,  8) / 1.82000E+01/
      data cdi(2,  8) / 2.69000E-01/
      data cdi(3,  8) / 3.89000E-01/
      data cdi(4,  8) /-1.83600E+00/
      data cdi(5,  8) / 3.93900E+00/
      data cdi(6,  8) /-2.27500E+00/
      data cea(1,  8) / 9.28200E+01/
      data cea(2,  8) / 1.82000E+01/
      data cea(3,  8) / 7.53000E-01/
      data cea(4,  8) /-5.82000E-01/
      data cea(5,  8) / 6.43000E-01/
      data cea(6,  8) /-9.66000E-01/
C     Z= 4 [+ 2]
      data error(  9) / 20.0/
      data cdi(1,  9) / 1.53900E+02/
      data cdi(2,  9) / 7.96000E-01/
      data cdi(3,  9) /-5.00000E-01/
      data cdi(4,  9) / 8.84000E-01/
C     Z= 4 [+ 3]
      data error( 10) / 20.0/
      data cdi(1, 10) / 2.17700E+02/
      data cdi(2, 10) / 4.00000E-01/
C     Z= 5 [+ 0]
      data error( 11) / 20.0/
      data cdi(1, 11) / 8.30000E+00/
      data cdi(2, 11) / 1.10600E+00/
      data cdi(3, 11) /-1.06900E+00/
      data cdi(4, 11) /-8.80000E-02/
C     Z= 5 [+ 1]
      data error( 12) / 20.0/
      data cdi(1, 12) / 2.51000E+01/
      data cdi(2, 12) / 9.07000E-01/
      data cdi(3, 12) /-4.77000E-01/
      data cdi(4, 12) / 1.97000E-01/
C     Z= 5 [+ 2]
      data error( 13) / 20.0/
      data cdi(1, 13) / 3.79000E+01/
      data cdi(2, 13) / 3.93000E-01/
      data cdi(3, 13) /-8.20000E-02/
      data cdi(4, 13) /-3.03000E-01/
      data cdi(5, 13) / 2.63000E-01/
C     Z= 5 [+ 3]
      data error( 14) / 20.0/
      data cdi(1, 14) / 2.59400E+02/
      data cdi(2, 14) / 7.96000E-01/
      data cdi(3, 14) /-5.00000E-01/
      data cdi(4, 14) / 8.84000E-01/
C     Z= 5 [+ 4]
      data error( 15) / 20.0/
      data cdi(1, 15) / 3.40200E+02/
      data cdi(2, 15) / 4.00000E-01/
C     Z= 6 [+ 0]
      data error( 16) /  5.0/
      data cdi(1, 16) / 1.13000E+01/
      data cdi(2, 16) / 2.11400E+00/
      data cdi(3, 16) /-1.96500E+00/
      data cdi(4, 16) /-6.08000E-01/
C     Z= 6 [+ 1]
      data error( 17) / 10.0/
      data cdi(1, 17) / 2.44000E+01/
      data cdi(2, 17) / 1.08200E+00/
      data cdi(3, 17) /-1.61000E-01/
      data cdi(4, 17) /-8.56000E-01/
      data cdi(5, 17) / 9.06000E-01/
C     Z= 6 [+ 2]
      data error( 18) / 10.0/
      data cdi(1, 18) / 4.79000E+01/
      data cdi(2, 18) / 7.15000E-01/
      data cdi(3, 18) /-4.10000E-02/
      data cdi(4, 18) / 1.75000E-01/
C     Z= 6 [+ 3]
      data error( 19) / 20.0/
      data cdi(1, 19) / 6.45000E+01/
      data cdi(2, 19) / 4.50000E-01/
      data cdi(3, 19) /-3.18000E-01/
      data cdi(4, 19) / 1.02600E+00/
      data cdi(5, 19) /-2.85900E+00/
      data cdi(6, 19) / 1.99500E+00/
C     Z= 6 [+ 4]
      data error( 20) / 20.0/
      data cdi(1, 20) / 3.92100E+02/
      data cdi(2, 20) / 7.96000E-01/
      data cdi(3, 20) /-5.00000E-01/
      data cdi(4, 20) / 8.84000E-01/
C     Z= 6 [+ 5]
      data error( 21) / 10.0/
      data cdi(1, 21) / 4.90000E+02/
      data cdi(2, 21) / 4.00000E-01/
C     Z= 7 [+ 0]
      data error( 22) /  5.0/
      data cdi(1, 22) / 1.45000E+01/
      data cdi(2, 22) / 2.26500E+00/
      data cdi(3, 22) /-1.71000E+00/
      data cdi(4, 22) /-2.32200E+00/
      data cdi(5, 22) / 1.73200E+00/
C     Z= 7 [+ 1]
      data error( 23) / 10.0/
      data cdi(1, 23) / 2.96000E+01/
      data cdi(2, 23) / 1.07600E+00/
      data cdi(3, 23) /-8.29000E-01/
      data cdi(4, 23) / 8.72000E-01/
      data cdi(5, 23) /-1.62000E-01/
      data cdi(6, 23) / 1.53300E+00/
C     Z= 7 [+ 2]
      data error( 24) / 10.0/
      data cdi(1, 24) / 4.75000E+01/
      data cdi(2, 24) / 5.00000E-01/
      data cdi(3, 24) / 2.23000E-01/
      data cdi(4, 24) / 2.20700E+00/
      data cdi(5, 24) /-4.15500E+00/
      data cdi(6, 24) / 3.76900E+00/
C     Z= 7 [+ 3]
      data error( 25) / 10.0/
      data cdi(1, 25) / 7.75000E+01/
      data cdi(2, 25) / 8.13000E-01/
      data cdi(3, 25) /-7.00000E-03/
      data cdi(4, 25) /-4.60000E-02/
C     Z= 7 [+ 4]
      data error( 26) / 10.0/
      data cdi(1, 26) / 9.79000E+01/
      data cdi(2, 26) / 2.18000E-01/
      data cdi(3, 26) / 2.38000E-01/
      data cdi(4, 26) /-2.20000E-01/
      data cdi(5, 26) /-4.46000E-01/
      data cdi(6, 26) / 2.52300E+00/
      data cdi(7, 26) /-1.90200E+00/
      data cea(1, 26) / 4.11180E+02/
      data cea(2, 26) / 9.79000E+01/
      data cea(3, 26) / 8.37000E-01/
      data cea(4, 26) /-2.14000E-01/
      data cea(5, 26) /-2.53800E+00/
      data cea(6, 26) / 7.48800E+00/
      data cea(7, 26) /-1.10060E+01/
      data cea(8, 26) / 5.52300E+00/
C     Z= 7 [+ 5]
      data error( 27) / 20.0/
      data cdi(1, 27) / 5.52100E+02/
      data cdi(2, 27) / 7.96000E-01/
      data cdi(3, 27) /-5.00000E-01/
      data cdi(4, 27) / 8.84000E-01/
C     Z= 7 [+ 6]
      data error( 28) / 10.0/
      data cdi(1, 28) / 6.67000E+02/
      data cdi(2, 28) / 4.00000E-01/
C     Z= 8 [+ 0]
      data error( 29) /  5.0/
      data cdi(1, 29) / 1.36000E+01/
      data cdi(2, 29) / 2.45500E+00/
      data cdi(3, 29) /-2.18100E+00/
      data cdi(4, 29) /-1.57000E+00/
C     Z= 8 [+ 1]
      data error( 30) / 10.0/
      data cdi(1, 30) / 3.51000E+01/
      data cdi(2, 30) / 1.52600E+00/
      data cdi(3, 30) /-5.93000E-01/
      data cdi(4, 30) /-3.99000E-01/
      data cdi(5, 30) /-5.83000E-01/
      data cdi(6, 30) / 3.23500E+00/
C     Z= 8 [+ 2]
      data error( 31) /  7.0/
      data cdi(1, 31) / 5.49000E+01/
      data cdi(2, 31) / 1.06600E+00/
      data cdi(3, 31) / 4.42000E-01/
      data cdi(4, 31) / 4.75000E-01/
      data cdi(5, 31) /-2.96100E+00/
      data cdi(6, 31) / 4.47000E+00/
C     Z= 8 [+ 3]
      data error( 32) / 10.0/
      data cdi(1, 32) / 7.74000E+01/
      data cdi(2, 32) / 1.04500E+00/
      data cdi(3, 32) /-6.52000E-01/
      data cdi(4, 32) / 1.29900E+00/
C     Z= 8 [+ 4]
      data error( 33) / 10.0/
      data cdi(1, 33) / 1.13900E+02/
      data cdi(2, 33) / 7.27000E-01/
      data cdi(3, 33) / 9.10000E-02/
      data cdi(4, 33) / 2.20000E-02/
C     Z= 8 [+ 5]
      data error( 34) / 20.0/
      data cdi(1, 34) / 1.38100E+02/
      data cdi(2, 34) / 3.36000E-01/
      data cdi(3, 34) / 8.00000E-02/
      data cdi(4, 34) / 1.43000E-01/
      data cdi(5, 34) /-7.31000E-01/
      data cdi(6, 34) / 1.33600E+00/
      data cdi(7, 34) /-7.85000E-01/
      data cea(1, 34) / 4.41920E+02/
      data cea(2, 34) / 3.20000E+02/
      data cea(3, 34) / 8.01000E-01/
      data cea(4, 34) / 3.86500E+00/
      data cea(5, 34) /-4.60500E+00/
      data cea(6, 34) / 1.54300E+00/
C     Z= 8 [+ 6]
      data error( 35) / 20.0/
      data cdi(1, 35) / 7.39300E+02/
      data cdi(2, 35) / 7.96000E-01/
      data cdi(3, 35) /-5.00000E-01/
      data cdi(4, 35) / 8.84000E-01/
C     Z= 8 [+ 7]
      data error( 36) / 10.0/
      data cdi(1, 36) / 8.71400E+02/
      data cdi(2, 36) / 4.00000E-01/
C     Z= 9 [+ 0]
      data error( 37) / 20.0/
      data cdi(1, 37) / 1.74000E+01/
      data cdi(2, 37) / 2.79000E+00/
      data cdi(3, 37) / 4.69000E-01/
      data cdi(4, 37) /-1.29000E+01/
      data cdi(5, 37) / 2.62600E+01/
      data cdi(6, 37) /-1.34300E+01/
      data cea(1, 37) / 7.98700E+01/
      data cea(2, 37) / 3.00000E+01/
      data cea(3, 37) / 3.92500E+00/
      data cea(4, 37) /-9.47000E-01/
      data cea(5, 37) /-5.68800E+00/
      data cea(6, 37) / 4.91100E+00/
      data cea(7, 37) /-8.30000E-02/
      data cea(8, 37) /-1.60000E-01/
C     Z= 9 [+ 1]
      data error( 38) / 60.0/
      data cdi(1, 38) / 3.50000E+01/
      data cdi(2, 38) / 2.01900E+00/
      data cdi(3, 38) /-1.32000E+00/
      data cdi(4, 38) / 1.70000E+00/
C     Z= 9 [+ 2]
      data error( 39) / 10.0/
      data cdi(1, 39) / 6.27000E+01/
      data cdi(2, 39) / 2.04200E+00/
      data cdi(3, 39) /-5.86000E-01/
      data cdi(4, 39) /-6.19000E-01/
      data cdi(5, 39) / 2.07200E+00/
C     Z= 9 [+ 3]
      data error( 40) / 40.0/
      data cdi(1, 40) / 8.71000E+01/
      data cdi(2, 40) / 1.06600E+00/
      data cdi(3, 40) / 4.42000E-01/
      data cdi(4, 40) / 4.75000E-01/
      data cdi(5, 40) /-2.96100E+00/
      data cdi(6, 40) / 4.47000E+00/
C     Z= 9 [+ 4]
      data error( 41) / 60.0/
      data cdi(1, 41) / 1.14200E+02/
      data cdi(2, 41) / 1.20000E+00/
      data cdi(3, 41) /-6.52000E-01/
      data cdi(4, 41) / 1.29900E+00/
C     Z= 9 [+ 5]
      data error( 42) / 20.0/
      data cdi(1, 42) / 1.57200E+02/
      data cdi(2, 42) / 6.99000E-01/
      data cdi(3, 42) / 2.20000E-02/
      data cdi(4, 42) / 2.26000E-01/
      data cdi(5, 42) /-1.79000E-01/
C     Z= 9 [+ 6]
      data error( 43) / 80.0/
      data cdi(1, 43) / 1.85200E+02/
      data cdi(2, 43) / 3.36000E-01/
      data cdi(3, 43) / 8.00000E-02/
      data cdi(4, 43) / 1.43000E-01/
      data cdi(5, 43) /-7.31000E-01/
      data cdi(6, 43) / 1.33600E+00/
      data cdi(7, 43) /-7.85000E-01/
C     Z= 9 [+ 7]
      data error( 44) / 50.0/
      data cdi(1, 44) / 9.53900E+02/
      data cdi(2, 44) / 7.96000E-01/
      data cdi(3, 44) /-5.00000E-01/
      data cdi(4, 44) / 8.84000E-01/
C     Z= 9 [+ 8]
      data error( 45) / 50.0/
      data cdi(1, 45) / 1.10310E+03/
      data cdi(2, 45) / 4.00000E-01/
C     Z=10 [+ 0]
      data error( 46) / 10.0/
      data cdi(1, 46) / 2.16000E+01/
      data cdi(2, 46) / 2.19200E+00/
      data cdi(3, 46) /-4.47000E-01/
      data cdi(4, 46) /-7.00600E+00/
      data cdi(5, 46) / 5.92700E+00/
C     Z=10 [+ 1]
      data error( 47) / 12.0/
      data cdi(1, 47) / 4.11000E+01/
      data cdi(2, 47) / 2.70500E+00/
      data cdi(3, 47) /-2.94600E+00/
      data cdi(4, 47) / 4.86200E+00/
      data cdi(5, 47) /-1.50700E+01/
      data cdi(6, 47) / 1.77800E+01/
      data cdi(7, 47) /-5.71600E+00/
C     Z=10 [+ 2]
      data error( 48) / 15.0/
      data cdi(1, 48) / 6.35000E+01/
      data cdi(2, 48) / 3.70100E+00/
      data cdi(3, 48) /-1.12800E+00/
      data cdi(4, 48) /-6.34400E+00/
      data cdi(5, 48) / 4.84200E+00/
C     Z=10 [+ 3]
      data error( 49) / 10.0/
      data cdi(1, 49) / 9.25000E+01/
      data cdi(2, 49) / 7.85000E-01/
      data cdi(3, 49) / 1.70900E+00/
      data cdi(4, 49) /-1.08500E+01/
      data cdi(5, 49) / 4.15000E+01/
      data cdi(6, 49) /-5.80500E+01/
      data cdi(7, 49) / 3.07200E+01/
C     Z=10 [+ 4]
      data error( 50) / 40.0/
      data cdi(1, 50) / 1.26200E+02/
      data cdi(2, 50) / 1.06600E+00/
      data cdi(3, 50) / 4.42000E-01/
      data cdi(4, 50) / 4.75000E-01/
      data cdi(5, 50) /-2.96100E+00/
      data cdi(6, 50) / 4.47000E+00/
C     Z=10 [+ 5]
      data error( 51) / 60.0/
      data cdi(1, 51) / 1.57900E+02/
      data cdi(2, 51) / 1.04500E+00/
      data cdi(3, 51) /-6.52000E-01/
      data cdi(4, 51) / 1.29900E+00/
C     Z=10 [+ 6]
      data error( 52) / 20.0/
      data cdi(1, 52) / 2.07500E+02/
      data cdi(2, 52) / 7.37000E-01/
      data cdi(3, 52) / 3.80000E-02/
      data cdi(4, 52) / 2.73000E-01/
      data cdi(5, 52) /-3.84000E-01/
      data cdi(6, 52) / 1.33000E-01/
C     Z=10 [+ 7]
      data error( 53) / 40.0/
      data cdi(1, 53) / 2.39100E+02/
      data cdi(2, 53) / 4.20000E-02/
      data cdi(3, 53) /-1.12000E-01/
      data cdi(4, 53) / 2.32500E+00/
      data cdi(5, 53) /-1.74300E+00/
      data cdi(6, 53) / 1.50000E-02/
      data cea(1, 53) / 8.96620E+02/
      data cea(2, 53) / 2.39100E+02/
      data cea(3, 53) / 6.47000E-01/
      data cea(4, 53) /-7.45000E-01/
      data cea(5, 53) / 2.76800E+00/
      data cea(6, 53) /-4.05900E+00/
      data cea(7, 53) / 1.52000E+00/
C     Z=10 [+ 8]
      data error( 54) / 15.0/
      data cdi(1, 54) / 1.19600E+03/
      data cdi(2, 54) / 8.84000E-01/
      data cdi(3, 54) /-3.48000E-01/
      data cdi(4, 54) / 8.72000E-01/
      data cdi(5, 54) /-8.44300E+00/
      data cdi(6, 54) / 1.71600E+01/
      data cdi(7, 54) /-9.17900E+00/
C     Z=10 [+ 9]
      data error( 55) / 20.0/
      data cdi(1, 55) / 1.36060E+03/
      data cdi(2, 55) / 4.34000E-01/
      data cdi(3, 55) /-6.90000E-02/
      data cdi(4, 55) / 9.30000E-02/
      data cdi(5, 55) /-3.90000E-02/
C     Z=11 [+ 0]
      data error( 56) / 20.0/
      data cdi(1, 56) / 5.10000E+00/
      data cdi(2, 56) / 7.96000E-01/
      data cdi(3, 56) /-7.70000E-01/
      data cdi(4, 56) / 2.32300E+00/
      data cdi(5, 56) /-3.92900E+00/
      data cdi(6, 56) / 1.56200E+00/
C     Z=11 [+ 1]
      data error( 57) / 10.0/
      data cdi(1, 57) / 4.73000E+01/
      data cdi(2, 57) / 3.76300E+00/
      data cdi(3, 57) /-1.31700E+00/
      data cdi(4, 57) /-1.54400E+01/
      data cdi(5, 57) / 2.73700E+01/
      data cdi(6, 57) /-1.49200E+01/
C     Z=11 [+ 2]
      data error( 58) / 30.0/
      data cdi(1, 58) / 7.16000E+01/
      data cdi(2, 58) / 2.79400E+00/
      data cdi(3, 58) / 4.69000E-01/
      data cdi(4, 58) /-1.29400E+01/
      data cdi(5, 58) / 2.62600E+01/
      data cdi(6, 58) /-1.34300E+01/
C     Z=11 [+ 3]
      data error( 59) / 60.0/
      data cdi(1, 59) / 9.89000E+01/
      data cdi(2, 59) / 2.01900E+00/
      data cdi(3, 59) /-1.32000E+00/
      data cdi(4, 59) / 1.70000E+00/
C     Z=11 [+ 4]
      data error( 60) / 40.0/
      data cdi(1, 60) / 1.38200E+02/
      data cdi(2, 60) / 7.85000E-01/
      data cdi(3, 60) / 1.70900E+00/
      data cdi(4, 60) /-1.08500E+01/
      data cdi(5, 60) / 4.15000E+01/
      data cdi(6, 60) /-5.80500E+01/
      data cdi(7, 60) / 3.07200E+01/
C     Z=11 [+ 5]
      data error( 61) / 40.0/
      data cdi(1, 61) / 1.72100E+02/
      data cdi(2, 61) / 1.06600E+00/
      data cdi(3, 61) / 4.42000E-01/
      data cdi(4, 61) / 4.75000E-01/
      data cdi(5, 61) /-2.96100E+00/
      data cdi(6, 61) / 4.47000E+00/
C     Z=11 [+ 6]
      data error( 62) / 60.0/
      data cdi(1, 62) / 2.08500E+02/
      data cdi(2, 62) / 1.20000E+00/
      data cdi(3, 62) /-6.52000E-01/
      data cdi(4, 62) / 1.29900E+00/
C     Z=11 [+ 7]
      data error( 63) / 40.0/
      data cdi(1, 63) / 2.64200E+02/
      data cdi(2, 63) / 7.27000E-01/
      data cdi(3, 63) / 9.10000E-02/
      data cdi(4, 63) / 2.20000E-02/
C     Z=11 [+ 8]
      data error( 64) / 80.0/
      data cdi(1, 64) / 2.99900E+02/
      data cdi(2, 64) / 3.36000E-01/
      data cdi(3, 64) / 8.00000E-02/
      data cdi(4, 64) / 1.43000E-01/
      data cdi(5, 64) /-7.31000E-01/
      data cdi(6, 64) / 1.33600E+00/
      data cdi(7, 64) /-7.85000E-01/
C     Z=11 [+ 9]
      data error( 65) / 30.0/
      data cdi(1, 65) / 1.46510E+03/
      data cdi(2, 65) / 8.40000E-01/
      data cdi(3, 65) /-3.10000E-02/
      data cdi(4, 65) /-7.13000E-01/
      data cdi(5, 65) / 2.44800E+00/
      data cdi(6, 65) /-2.40100E+00/
      data cdi(7, 65) / 7.15000E-01/
C     Z=11 [+10]
      data error( 66) / 50.0/
      data cdi(1, 66) / 1.64870E+03/
      data cdi(2, 66) / 4.00000E-01/
C     Z=12 [+ 0]
      data error( 67) / 15.0/
      data cdi(1, 67) / 7.60000E+00/
      data cdi(2, 67) / 4.84000E-01/
      data cdi(3, 67) / 1.75000E+00/
      data cdi(4, 67) /-1.56200E+00/
      data cdi(5, 67) / 3.78700E+00/
C     Z=12 [+ 1]
      data error( 68) / 15.0/
      data cdi(1, 68) / 1.52000E+01/
      data cdi(2, 68) / 9.14000E-01/
      data cdi(3, 68) /-1.88000E-01/
      data cdi(4, 68) /-4.38300E+00/
      data cdi(5, 68) / 1.38300E+01/
      data cdi(6, 68) /-2.05700E+01/
      data cdi(7, 68) / 1.00100E+01/
C     Z=12 [+ 2]
      data error( 69) / 10.0/
      data cdi(1, 69) / 7.80000E+01/
      data cdi(2, 69) / 4.21500E+00/
      data cdi(3, 69) /-1.29500E+00/
      data cdi(4, 69) /-1.66600E+01/
      data cdi(5, 69) / 3.34300E+01/
      data cdi(6, 69) /-1.95700E+01/
C     Z=12 [+ 3]
      data error( 70) / 70.0/
      data cdi(1, 70) / 1.09200E+02/
      data cdi(2, 70) / 2.79000E+00/
      data cdi(3, 70) / 4.69000E-01/
      data cdi(4, 70) /-1.29400E+01/
      data cdi(5, 70) / 2.62600E+01/
      data cdi(6, 70) /-1.34300E+01/
C     Z=12 [+ 4]
      data error( 71) / 60.0/
      data cdi(1, 71) / 1.41300E+02/
      data cdi(2, 71) / 2.01900E+00/
      data cdi(3, 71) /-1.32000E+00/
      data cdi(4, 71) / 1.70000E+00/
C     Z=12 [+ 5]
      data error( 72) / 40.0/
      data cdi(1, 72) / 1.86500E+02/
      data cdi(2, 72) / 7.85000E-01/
      data cdi(3, 72) / 1.70900E+00/
      data cdi(4, 72) /-1.08500E+01/
      data cdi(5, 72) / 4.15000E+01/
      data cdi(6, 72) /-5.80500E+01/
      data cdi(7, 72) / 3.07200E+01/
C     Z=12 [+ 6]
      data error( 73) / 40.0/
      data cdi(1, 73) / 2.24900E+02/
      data cdi(2, 73) / 1.06600E+00/
      data cdi(3, 73) / 4.42000E-01/
      data cdi(4, 73) / 4.75000E-01/
      data cdi(5, 73) /-2.96100E+00/
      data cdi(6, 73) / 4.47000E+00/
C     Z=12 [+ 7]
      data error( 74) / 60.0/
      data cdi(1, 74) / 2.65900E+02/
      data cdi(2, 74) / 1.20000E+00/
      data cdi(3, 74) /-6.52000E-01/
      data cdi(4, 74) / 1.29900E+00/
C     Z=12 [+ 8]
      data error( 75) / 40.0/
      data cdi(1, 75) / 3.28000E+02/
      data cdi(2, 75) / 7.27000E-01/
      data cdi(3, 75) / 9.10000E-02/
      data cdi(4, 75) / 2.20000E-02/
C     Z=12 [+ 9]
      data error( 76) / 40.0/
      data cdi(1, 76) / 3.67700E+02/
      data cdi(2, 76) / 3.61000E-01/
      data cdi(3, 76) / 7.90000E-02/
      data cdi(4, 76) /-1.12000E-01/
      data cdi(5, 76) / 5.00000E-02/
C     Z=12 [+10]
      data error( 77) / 30.0/
      data cdi(1, 77) / 1.76180E+03/
      data cdi(2, 77) / 9.30000E-01/
      data cdi(3, 77) /-2.34000E-01/
      data cdi(4, 77) /-5.13000E-01/
      data cdi(5, 77) /-4.62200E+00/
      data cdi(6, 77) / 1.24000E+01/
      data cdi(7, 77) /-7.04500E+00/
C     Z=12 [+11]
      data error( 78) / 50.0/
      data cdi(1, 78) / 1.96260E+03/
      data cdi(2, 78) / 4.00000E-01/
C     Z=13 [+ 0]
      data error( 79) / 60.0/
      data cdi(1, 79) / 6.00000E+00/
      data cdi(2, 79) / 1.17000E+00/
      data cdi(3, 79) / 8.43000E-01/
      data cdi(4, 79) /-2.87700E+00/
      data cdi(5, 79) / 1.95800E+00/
C     Z=13 [+ 1]
      data error( 80) / 10.0/
      data cdi(1, 80) / 1.88000E+01/
      data cdi(2, 80) / 4.96000E-01/
      data cdi(3, 80) /-9.00000E-02/
      data cdi(4, 80) / 4.72400E+00/
      data cdi(5, 80) /-1.38100E+01/
      data cdi(6, 80) / 1.50000E+01/
      data cdi(7, 80) /-5.73900E+00/
C     Z=13 [+ 2]
      data error( 81) / 10.0/
      data cdi(1, 81) / 2.80000E+01/
      data cdi(2, 81) / 9.20000E-02/
      data cdi(3, 81) / 3.30000E-02/
      data cdi(4, 81) / 2.63000E-01/
      data cdi(5, 81) / 7.60000E-02/
      data cdi(6, 81) /-2.18000E-01/
      data cea(1, 81) / 8.00800E+01/
      data cea(2, 81) / 4.50000E+01/
      data cea(3, 81) / 1.99700E+00/
      data cea(4, 81) /-1.25900E+00/
      data cea(5, 81) /-7.50000E-02/
      data cea(6, 81) /-4.60000E+00/
      data cea(7, 81) / 4.01400E+00/
C     Z=13 [+ 3]
      data error( 82) / 40.0/
      data cdi(1, 82) / 1.20000E+02/
      data cdi(2, 82) / 3.71000E+00/
      data cdi(3, 82) /-2.36000E+00/
      data cdi(4, 82) / 2.35000E+00/
      data cdi(5, 82) /-6.60000E-01/
C     Z=13 [+ 4]
      data error( 83) / 70.0/
      data cdi(1, 83) / 1.53700E+02/
      data cdi(2, 83) / 2.79400E+00/
      data cdi(3, 83) / 4.69000E-01/
      data cdi(4, 83) /-1.29400E+01/
      data cdi(5, 83) / 2.62600E+01/
      data cdi(6, 83) /-1.34300E+01/
C     Z=13 [+ 5]
      data error( 84) / 40.0/
      data cdi(1, 84) / 1.90500E+02/
      data cdi(2, 84) / 2.41900E+00/
      data cdi(3, 84) /-1.32000E+00/
      data cdi(4, 84) / 1.70000E+00/
C     Z=13 [+ 6]
      data error( 85) / 40.0/
      data cdi(1, 85) / 2.41400E+02/
      data cdi(2, 85) / 7.85000E-01/
      data cdi(3, 85) / 1.70900E+00/
      data cdi(4, 85) /-1.08500E+01/
      data cdi(5, 85) / 4.15000E+01/
      data cdi(6, 85) /-5.80500E+01/
      data cdi(7, 85) / 3.07200E+01/
C     Z=13 [+ 7]
      data error( 86) / 30.0/
      data cdi(1, 86) / 2.84600E+02/
      data cdi(2, 86) / 2.00300E+00/
      data cdi(3, 86) /-4.05000E-01/
      data cdi(4, 86) /-6.51900E+00/
      data cdi(5, 86) / 1.69800E+01/
      data cdi(6, 86) /-1.17600E+01/
C     Z=13 [+ 8]
      data error( 87) / 30.0/
      data cdi(1, 87) / 3.90200E+02/
      data cdi(2, 87) / 7.80000E-01/
      data cdi(3, 87) / 9.07000E-01/
      data cdi(4, 87) /-3.01000E-01/
      data cdi(5, 87) /-2.33700E+00/
      data cdi(6, 87) / 9.42000E+00/
      data cdi(7, 87) /-6.11900E+00/
C     Z=13 [+ 9]
      data error( 88) / 20.0/
      data cdi(1, 88) / 3.98000E+02/
      data cdi(2, 88) / 6.01000E-01/
      data cdi(3, 88) /-5.49000E-01/
      data cdi(4, 88) / 1.09200E+00/
      data cdi(5, 88) / 9.10000E-02/
C     Z=13 [+10]
      data error( 89) / 40.0/
      data cdi(1, 89) / 4.42100E+02/
      data cdi(2, 89) / 7.18000E-01/
      data cdi(3, 89) /-4.32000E-01/
      data cdi(4, 89) /-3.04000E-01/
      data cdi(5, 89) / 1.81000E-01/
C     Z=13 [+11]
      data error( 90) / 30.0/
      data cdi(1, 90) / 2.08600E+03/
      data cdi(2, 90) / 1.01200E+00/
      data cdi(3, 90) /-6.78000E-01/
      data cdi(4, 90) /-3.60000E-01/
      data cdi(5, 90) / 2.83000E+00/
      data cdi(6, 90) /-1.92700E+00/
C     Z=13 [+12]
      data error( 91) / 50.0/
      data cdi(1, 91) / 2.30410E+03/
      data cdi(2, 91) / 4.00000E-01/
C     Z=14 [+ 0]
      data error( 92) / 60.0/
      data cdi(1, 92) / 8.20000E+00/
      data cdi(2, 92) / 1.57300E+00/
      data cdi(3, 92) / 7.22000E-01/
      data cdi(4, 92) /-2.68700E+00/
      data cdi(5, 92) / 1.85600E+00/
C     Z=14 [+ 1]
      data error( 93) / 60.0/
      data cdi(1, 93) / 1.63000E+01/
      data cdi(2, 93) / 1.17000E+00/
      data cdi(3, 93) / 8.43000E-01/
      data cdi(4, 93) /-2.87700E+00/
      data cdi(5, 93) / 1.95800E+00/
C     Z=14 [+ 2]
      data error( 94) /100.0/
      data cdi(1, 94) / 3.35000E+01/
      data cdi(2, 94) / 1.08800E+00/
      data cdi(3, 94) / 2.03000E-01/
      data cdi(4, 94) /-2.43000E-01/
      data cdi(5, 94) / 9.90000E-02/
C     Z=14 [+ 3]
      data error( 95) / 15.0/
      data cdi(1, 95) / 5.40000E+01/
      data cdi(2, 95) / 2.46000E-01/
      data cdi(3, 95) / 3.98000E-01/
      data cdi(4, 95) /-8.26000E-01/
      data cdi(5, 95) / 2.14700E+00/
      data cdi(6, 95) /-1.59700E+00/
      data cea(1, 95) / 1.09080E+02/
      data cea(2, 95) / 9.00000E+01/
      data cea(3, 95) / 3.55700E+00/
      data cea(4, 95) /-2.20000E-01/
      data cea(5, 95) /-5.64400E+00/
      data cea(6, 95) / 5.85600E+00/
      data cea(7, 95) /-1.89000E+00/
C     Z=14 [+ 4]
      data error( 96) / 50.0/
      data cdi(1, 96) / 1.66800E+02/
      data cdi(2, 96) / 4.21500E+00/
      data cdi(3, 96) /-1.29500E+00/
      data cdi(4, 96) /-1.66600E+01/
      data cdi(5, 96) / 3.34300E+01/
      data cdi(6, 96) /-1.95700E+01/
C     Z=14 [+ 5]
      data error( 97) / 70.0/
      data cdi(1, 97) / 2.05100E+02/
      data cdi(2, 97) / 2.79400E+00/
      data cdi(3, 97) / 4.69000E-01/
      data cdi(4, 97) /-1.29400E+01/
      data cdi(5, 97) / 2.62600E+01/
      data cdi(6, 97) /-1.34300E+01/
C     Z=14 [+ 6]
      data error( 98) / 60.0/
      data cdi(1, 98) / 2.46500E+02/
      data cdi(2, 98) / 2.01900E+00/
      data cdi(3, 98) /-1.32000E+00/
      data cdi(4, 98) / 1.70000E+00/
C     Z=14 [+ 7]
      data error( 99) / 40.0/
      data cdi(1, 99) / 3.03200E+02/
      data cdi(2, 99) / 7.85000E-01/
      data cdi(3, 99) / 1.70900E+00/
      data cdi(4, 99) /-1.08500E+01/
      data cdi(5, 99) / 4.15000E+01/
      data cdi(6, 99) /-5.80500E+01/
      data cdi(7, 99) / 3.07200E+01/
C     Z=14 [+ 8]
      data error(100) / 40.0/
      data cdi(1,100) / 3.51100E+02/
      data cdi(2,100) / 1.06600E+00/
      data cdi(3,100) / 4.42000E-01/
      data cdi(4,100) / 4.75000E-01/
      data cdi(5,100) /-2.96100E+00/
      data cdi(6,100) / 4.47000E+00/
C     Z=14 [+ 9]
      data error(101) / 60.0/
      data cdi(1,101) / 4.01400E+02/
      data cdi(2,101) / 1.20000E+00/
      data cdi(3,101) /-6.52000E-01/
      data cdi(4,101) / 1.29900E+00/
C     Z=14 [+10]
      data error(102) / 40.0/
      data cdi(1,102) / 4.76100E+02/
      data cdi(2,102) / 7.27000E-01/
      data cdi(3,102) / 9.10000E-02/
      data cdi(4,102) / 2.20000E-02/
C     Z=14 [+11]
      data error(103) / 80.0/
      data cdi(1,103) / 5.23500E+02/
      data cdi(2,103) / 3.36000E-01/
      data cdi(3,103) / 8.00000E-02/
      data cdi(4,103) / 1.43000E-01/
      data cdi(5,103) /-7.31000E-01/
      data cdi(6,103) / 1.33600E+00/
      data cdi(7,103) /-7.85000E-01/
C     Z=14 [+12]
      data error(104) / 50.0/
      data cdi(1,104) / 2.43770E+03/
      data cdi(2,104) / 7.96000E-01/
      data cdi(3,104) /-5.00000E-01/
      data cdi(4,104) / 8.84000E-01/
C     Z=14 [+13]
      data error(105) / 50.0/
      data cdi(1,105) / 2.67310E+03/
      data cdi(2,105) / 4.00000E-01/
C     Z=15 [+ 0]
      data error(106) / 50.0/
      data cdi(1,106) / 1.05000E+01/
      data cdi(2,106) / 1.70400E+00/
      data cdi(3,106) / 1.51800E+00/
      data cdi(4,106) /-2.98200E+00/
      data cdi(5,106) / 1.77400E+00/
C     Z=15 [+ 1]
      data error(107) / 60.0/
      data cdi(1,107) / 1.97000E+01/
      data cdi(2,107) / 1.57300E+00/
      data cdi(3,107) / 7.22000E-01/
      data cdi(4,107) /-2.68700E+00/
      data cdi(5,107) / 1.85600E+00/
C     Z=15 [+ 2]
      data error(108) / 60.0/
      data cdi(1,108) / 3.02000E+01/
      data cdi(2,108) / 1.17000E+00/
      data cdi(3,108) / 8.43000E-01/
      data cdi(4,108) /-2.87700E+00/
      data cdi(5,108) / 1.95800E+00/
C     Z=15 [+ 3]
      data error(109) /100.0/
      data cdi(1,109) / 5.14000E+01/
      data cdi(2,109) / 1.08800E+00/
      data cdi(3,109) / 2.03000E-01/
      data cdi(4,109) /-2.43000E-01/
      data cdi(5,109) / 9.90000E-02/
C     Z=15 [+ 4]
      data error(110) /200.0/
      data cdi(1,110) / 6.50000E+01/
      data cdi(2,110) / 9.07000E-01/
      data cdi(3,110) / 1.05000E-01/
      data cdi(4,110) / 1.47000E-01/
      data cdi(5,110) /-7.30000E-02/
      data cdi(6,110) /-1.00000E-03/
C     Z=15 [+ 5]
      data error(111) / 50.0/
      data cdi(1,111) / 2.20400E+02/
      data cdi(2,111) / 3.50000E+00/
      data cdi(3,111) /-2.36000E+00/
      data cdi(4,111) / 2.35000E+00/
      data cdi(5,111) /-6.60000E-01/
C     Z=15 [+ 6]
      data error(112) / 70.0/
      data cdi(1,112) / 2.63200E+02/
      data cdi(2,112) / 2.79400E+00/
      data cdi(3,112) / 4.69000E-01/
      data cdi(4,112) /-1.29400E+01/
      data cdi(5,112) / 2.62600E+01/
      data cdi(6,112) /-1.34300E+01/
C     Z=15 [+ 7]
      data error(113) / 60.0/
      data cdi(1,113) / 3.09400E+02/
      data cdi(2,113) / 2.01900E+00/
      data cdi(3,113) /-1.32000E+00/
      data cdi(4,113) / 1.70000E+00/
C     Z=15 [+ 8]
      data error(114) / 40.0/
      data cdi(1,114) / 3.71700E+02/
      data cdi(2,114) / 7.85000E-01/
      data cdi(3,114) / 1.70900E+00/
      data cdi(4,114) /-1.08500E+01/
      data cdi(5,114) / 4.15000E+01/
      data cdi(6,114) /-5.80500E+01/
      data cdi(7,114) / 3.07200E+01/
C     Z=15 [+ 9]
      data error(115) / 40.0/
      data cdi(1,115) / 4.24500E+02/
      data cdi(2,115) / 1.06600E+00/
      data cdi(3,115) / 4.42000E-01/
      data cdi(4,115) / 4.75000E-01/
      data cdi(5,115) /-2.96100E+00/
      data cdi(6,115) / 4.47000E+00/
C     Z=15 [+10]
      data error(116) / 60.0/
      data cdi(1,116) / 4.79600E+02/
      data cdi(2,116) / 1.20000E+00/
      data cdi(3,116) /-6.52000E-01/
      data cdi(4,116) / 1.29900E+00/
C     Z=15 [+11]
      data error(117) / 40.0/
      data cdi(1,117) / 5.60400E+02/
      data cdi(2,117) / 6.68000E-01/
      data cdi(3,117) / 1.97000E-01/
      data cdi(4,117) / 1.20000E-01/
      data cdi(5,117) /-1.78000E-01/
C     Z=15 [+12]
      data error(118) / 80.0/
      data cdi(1,118) / 6.11800E+02/
      data cdi(2,118) / 3.36000E-01/
      data cdi(3,118) / 8.00000E-02/
      data cdi(4,118) / 1.43000E-01/
      data cdi(5,118) /-7.31000E-01/
      data cdi(6,118) / 1.33600E+00/
      data cdi(7,118) /-7.85000E-01/
C     Z=15 [+13]
      data error(119) / 50.0/
      data cdi(1,119) / 2.81690E+03/
      data cdi(2,119) / 7.96000E-01/
      data cdi(3,119) /-5.00000E-01/
      data cdi(4,119) / 8.84000E-01/
C     Z=15 [+14]
      data error(120) / 50.0/
      data cdi(1,120) / 3.06980E+03/
      data cdi(2,120) / 4.00000E-01/
C     Z=16 [+ 0]
      data error(121) / 20.0/
      data cdi(1,121) / 1.04000E+01/
      data cdi(2,121) / 3.15000E+00/
      data cdi(3,121) /-2.35000E+00/
      data cdi(4,121) /-2.03200E+00/
C     Z=16 [+ 1]
      data error(122) / 50.0/
      data cdi(1,122) / 2.33000E+01/
      data cdi(2,122) / 1.70400E+00/
      data cdi(3,122) / 1.51800E+00/
      data cdi(4,122) /-2.98200E+00/
      data cdi(5,122) / 1.77400E+00/
C     Z=16 [+ 2]
      data error(123) / 60.0/
      data cdi(1,123) / 3.48000E+01/
      data cdi(2,123) / 1.66600E+00/
      data cdi(3,123) / 2.99000E-01/
      data cdi(4,123) /-2.14800E+00/
      data cdi(5,123) / 2.01400E+00/
C     Z=16 [+ 3]
      data error(124) / 60.0/
      data cdi(1,124) / 4.73000E+01/
      data cdi(2,124) / 1.17000E+00/
      data cdi(3,124) / 8.43000E-01/
      data cdi(4,124) /-2.87700E+00/
      data cdi(5,124) / 1.95800E+00/
C     Z=16 [+ 4]
      data error(125) / 10.0/
      data cdi(1,125) / 7.27000E+01/
      data cdi(2,125) / 5.98000E-01/
      data cdi(3,125) / 1.98000E-01/
      data cdi(4,125) /-1.76000E-01/
      data cdi(5,125) / 5.10000E-02/
      data cdi(6,125) / 1.10000E-02/
      data cdi(7,125) / 1.00000E-02/
      data cea(1,125) / 1.49760E+02/
      data cea(2,125) / 1.20000E+02/
      data cea(3,125) / 2.14700E+00/
      data cea(4,125) / 1.54200E+00/
      data cea(5,125) /-1.50800E+00/
      data cea(6,125) / 8.74000E-01/
      data cea(7,125) /-2.20000E-02/
      data cea(8,125) /-5.00000E-02/
C     Z=16 [+ 5]
      data error(126) /200.0/
      data cdi(1,126) / 8.81000E+01/
      data cdi(2,126) / 9.07000E-01/
      data cdi(3,126) / 1.05000E-01/
      data cdi(4,126) / 1.47000E-01/
      data cdi(5,126) /-7.30000E-02/
      data cdi(6,126) /-1.00000E-03/
C     Z=16 [+ 6]
      data error(127) / 50.0/
      data cdi(1,127) / 2.80900E+02/
      data cdi(2,127) / 3.50000E+00/
      data cdi(3,127) /-2.36000E+00/
      data cdi(4,127) / 2.35000E+00/
      data cdi(5,127) /-6.60000E-01/
C     Z=16 [+ 7]
      data error(128) / 70.0/
      data cdi(1,128) / 3.28200E+02/
      data cdi(2,128) / 2.79400E+00/
      data cdi(3,128) / 4.69000E-01/
      data cdi(4,128) /-1.29400E+01/
      data cdi(5,128) / 2.62600E+01/
      data cdi(6,128) /-1.34300E+01/
C     Z=16 [+ 8]
      data error(129) / 60.0/
      data cdi(1,129) / 3.79100E+02/
      data cdi(2,129) / 2.01900E+00/
      data cdi(3,129) /-1.32000E+00/
      data cdi(4,129) / 1.70000E+00/
C     Z=16 [+ 9]
      data error(130) / 40.0/
      data cdi(1,130) / 4.47100E+02/
      data cdi(2,130) / 7.85000E-01/
      data cdi(3,130) / 1.70900E+00/
      data cdi(4,130) /-1.08500E+01/
      data cdi(5,130) / 4.15000E+01/
      data cdi(6,130) /-5.80500E+01/
      data cdi(7,130) / 3.07200E+01/
C     Z=16 [+10]
      data error(131) / 40.0/
      data cdi(1,131) / 5.04800E+02/
      data cdi(2,131) / 1.06600E+00/
      data cdi(3,131) / 4.42000E-01/
      data cdi(4,131) / 4.75000E-01/
      data cdi(5,131) /-2.96100E+00/
      data cdi(6,131) / 4.47000E+00/
C     Z=16 [+11]
      data error(132) / 60.0/
      data cdi(1,132) / 5.64700E+02/
      data cdi(2,132) / 1.20000E+00/
      data cdi(3,132) /-6.52000E-01/
      data cdi(4,132) / 1.29900E+00/
C     Z=16 [+12]
      data error(133) / 40.0/
      data cdi(1,133) / 6.51600E+02/
      data cdi(2,133) / 6.68000E-01/
      data cdi(3,133) / 1.97000E-01/
      data cdi(4,133) / 1.20000E-01/
      data cdi(5,133) /-1.78000E-01/
C     Z=16 [+13]
      data error(134) / 80.0/
      data cdi(1,134) / 7.07100E+02/
      data cdi(2,134) / 3.36000E-01/
      data cdi(3,134) / 8.00000E-02/
      data cdi(4,134) / 1.43000E-01/
      data cdi(5,134) /-7.31000E-01/
      data cdi(6,134) / 1.33600E+00/
      data cdi(7,134) /-7.85000E-01/
C     Z=16 [+14]
      data error(135) / 50.0/
      data cdi(1,135) / 3.22380E+03/
      data cdi(2,135) / 7.96000E-01/
      data cdi(3,135) /-5.00000E-01/
      data cdi(4,135) / 8.84000E-01/
C     Z=16 [+15]
      data error(136) / 50.0/
      data cdi(1,136) / 3.49310E+03/
      data cdi(2,136) / 4.00000E-01/
C     Z=17 [+ 0]
      data error(137) / 15.0/
      data cdi(1,137) / 1.30000E+01/
      data cdi(2,137) / 2.24100E+00/
      data cdi(3,137) /-1.94800E+00/
      data cdi(4,137) / 1.18900E+00/
      data cdi(5,137) /-2.30000E-02/
      data cdi(6,137) / 5.00000E-03/
      data cea(1,137) / 5.01800E+01/
      data cea(2,137) / 2.00000E+01/
      data cea(3,137) / 6.10200E+00/
      data cea(4,137) / 2.19500E+00/
      data cea(5,137) /-1.56160E+01/
      data cea(6,137) / 1.03100E+01/
      data cea(7,137) /-5.88000E-01/
C     Z=17 [+ 1]
      data error(138) / 50.0/
      data cdi(1,138) / 2.38000E+01/
      data cdi(2,138) / 2.08600E+00/
      data cdi(3,138) / 1.07700E+00/
      data cdi(4,138) /-2.17200E+00/
      data cdi(5,138) / 8.09000E-01/
C     Z=17 [+ 2]
      data error(139) / 10.0/
      data cdi(1,139) / 3.96000E+01/
      data cdi(2,139) / 1.70400E+00/
      data cdi(3,139) / 1.51800E+00/
      data cdi(4,139) /-2.98200E+00/
      data cdi(5,139) / 1.77400E+00/
C     Z=17 [+ 3]
      data error(140) / 60.0/
      data cdi(1,140) / 5.35000E+01/
      data cdi(2,140) / 1.66600E+00/
      data cdi(3,140) / 2.99000E-01/
      data cdi(4,140) /-2.14800E+00/
      data cdi(5,140) / 2.01400E+00/
C     Z=17 [+ 4]
      data error(141) / 60.0/
      data cdi(1,141) / 6.78000E+01/
      data cdi(2,141) / 1.17000E+00/
      data cdi(3,141) / 8.43000E-01/
      data cdi(4,141) /-2.87700E+00/
      data cdi(5,141) / 1.95800E+00/
C     Z=17 [+ 5]
      data error(142) / 10.0/
      data cdi(1,142) / 9.70000E+01/
      data cdi(2,142) / 9.82000E-01/
      data cdi(3,142) /-6.33000E-01/
      data cdi(4,142) / 3.73000E-01/
      data cdi(5,142) / 5.50000E-02/
      data cdi(6,142) / 1.30000E-02/
      data cdi(7,142) / 1.20000E-02/
      data cea(1,142) / 1.99820E+02/
      data cea(2,142) / 1.70000E+02/
      data cea(3,142) / 3.57200E+00/
      data cea(4,142) / 4.30100E+00/
      data cea(5,142) /-1.07090E+01/
      data cea(6,142) / 5.99300E+00/
      data cea(7,142) / 4.60000E-02/
      data cea(8,142) / 5.60000E-02/
C     Z=17 [+ 6]
      data error(143) /200.0/
      data cdi(1,143) / 1.14200E+02/
      data cdi(2,143) / 9.07000E-01/
      data cdi(3,143) / 1.05000E-01/
      data cdi(4,143) / 1.47000E-01/
      data cdi(5,143) /-7.30000E-02/
      data cdi(6,143) /-1.00000E-03/
C     Z=17 [+ 7]
      data error(144) / 50.0/
      data cdi(1,144) / 3.48300E+02/
      data cdi(2,144) / 3.50000E+00/
      data cdi(3,144) /-2.36000E+00/
      data cdi(4,144) / 2.35000E+00/
      data cdi(5,144) /-6.60000E-01/
C     Z=17 [+ 8]
      data error(145) / 70.0/
      data cdi(1,145) / 4.00000E+02/
      data cdi(2,145) / 2.79400E+00/
      data cdi(3,145) / 4.69000E-01/
      data cdi(4,145) /-1.29400E+01/
      data cdi(5,145) / 2.62600E+01/
      data cdi(6,145) /-1.34300E+01/
C     Z=17 [+ 9]
      data error(146) / 60.0/
      data cdi(1,146) / 4.55600E+02/
      data cdi(2,146) / 2.01900E+00/
      data cdi(3,146) /-1.32000E+00/
      data cdi(4,146) / 1.70000E+00/
C     Z=17 [+10]
      data error(147) / 40.0/
      data cdi(1,147) / 5.29300E+02/
      data cdi(2,147) / 7.85000E-01/
      data cdi(3,147) / 1.70900E+00/
      data cdi(4,147) /-1.08500E+01/
      data cdi(5,147) / 4.15000E+01/
      data cdi(6,147) /-5.80500E+01/
      data cdi(7,147) / 3.07200E+01/
C     Z=17 [+11]
      data error(148) / 40.0/
      data cdi(1,148) / 5.92000E+02/
      data cdi(2,148) / 1.06600E+00/
      data cdi(3,148) / 4.42000E-01/
      data cdi(4,148) / 4.75000E-01/
      data cdi(5,148) /-2.96100E+00/
      data cdi(6,148) / 4.47000E+00/
C     Z=17 [+12]
      data error(149) / 60.0/
      data cdi(1,149) / 6.59700E+02/
      data cdi(2,149) / 1.20000E+00/
      data cdi(3,149) /-6.52000E-01/
      data cdi(4,149) / 1.29900E+00/
C     Z=17 [+13]
      data error(150) / 40.0/
      data cdi(1,150) / 7.49700E+02/
      data cdi(2,150) / 6.68000E-01/
      data cdi(3,150) / 1.97000E-01/
      data cdi(4,150) / 1.20000E-01/
      data cdi(5,150) /-1.78000E-01/
C     Z=17 [+14]
      data error(151) / 80.0/
      data cdi(1,151) / 8.09400E+02/
      data cdi(2,151) / 3.36000E-01/
      data cdi(3,151) / 8.00000E-02/
      data cdi(4,151) / 1.43000E-01/
      data cdi(5,151) /-7.31000E-01/
      data cdi(6,151) / 1.33600E+00/
      data cdi(7,151) /-7.85000E-01/
C     Z=17 [+15]
      data error(152) / 50.0/
      data cdi(1,152) / 3.65840E+03/
      data cdi(2,152) / 7.96000E-01/
      data cdi(3,152) /-5.00000E-01/
      data cdi(4,152) / 8.84000E-01/
C     Z=17 [+16]
      data error(153) / 50.0/
      data cdi(1,153) / 3.94620E+03/
      data cdi(2,153) / 4.00000E-01/
C     Z=18 [+ 0]
      data error(154) / 15.0/
      data cdi(1,154) / 1.58000E+01/
      data cdi(2,154) / 2.53200E+00/
      data cdi(3,154) /-2.67200E+00/
      data cdi(4,154) / 2.54300E+00/
      data cdi(5,154) /-7.69000E-01/
      data cdi(6,154) / 8.00000E-03/
      data cdi(7,154) / 6.00000E-03/
      data cea(1,154) / 6.19500E+01/
      data cea(2,154) / 1.58000E+01/
      data cea(3,154) / 4.33700E+00/
      data cea(4,154) / 3.09200E+00/
      data cea(5,154) /-2.12530E+01/
      data cea(6,154) / 1.46260E+01/
      data cea(7,154) / 1.80000E-02/
      data cea(8,154) / 3.10000E-02/
C     Z=18 [+ 1]
      data error(155) / 15.0/
      data cdi(1,155) / 2.74000E+01/
      data cdi(2,155) / 2.89600E+00/
      data cdi(3,155) / 7.77000E-01/
      data cdi(4,155) /-4.44700E+00/
      data cdi(5,155) / 2.86700E+00/
      data cdi(6,155) / 7.20000E-02/
C     Z=18 [+ 2]
      data error(156) / 15.0/
      data cdi(1,156) / 4.07000E+01/
      data cdi(2,156) / 2.08600E+00/
      data cdi(3,156) / 1.07700E+00/
      data cdi(4,156) /-2.17200E+00/
      data cdi(5,156) / 8.09000E-01/
C     Z=18 [+ 3]
      data error(157) / 10.0/
      data cdi(1,157) / 5.23000E+01/
      data cdi(2,157) / 1.18600E+00/
      data cdi(3,157) /-1.18000E+00/
      data cdi(4,157) / 1.10500E+01/
      data cdi(5,157) /-3.07900E+01/
      data cdi(6,157) / 3.66200E+01/
      data cdi(7,157) /-1.54200E+01/
C     Z=18 [+ 4]
      data error(158) / 15.0/
      data cdi(1,158) / 7.50000E+01/
      data cdi(2,158) / 1.57400E+00/
      data cdi(3,158) / 7.22000E-01/
      data cdi(4,158) /-2.68700E+00/
      data cdi(5,158) / 1.85600E+00/
      data cea(1,158) / 2.19750E+02/
      data cea(2,158) / 1.26000E+02/
      data cea(3,158) / 2.79800E+00/
      data cea(4,158) / 4.11400E+00/
      data cea(5,158) /-3.10300E+00/
      data cea(6,158) / 4.38000E-01/
C     Z=18 [+ 5]
      data error(159) / 10.0/
      data cdi(1,159) / 9.10000E+01/
      data cdi(2,159) / 1.17000E+00/
      data cdi(3,159) / 8.43000E-01/
      data cdi(4,159) /-2.87700E+00/
      data cdi(5,159) / 1.95800E+00/
      data cea(1,159) / 2.30230E+02/
      data cea(2,159) / 2.00000E+02/
      data cea(3,159) / 3.77100E+00/
      data cea(4,159) / 1.61630E+01/
      data cea(5,159) /-3.49520E+01/
      data cea(6,159) / 2.08530E+01/
C     Z=18 [+ 6]
      data error(160) / 10.0/
      data cdi(1,160) / 1.24300E+02/
      data cdi(2,160) / 9.68000E-01/
      data cdi(3,160) /-3.06000E-01/
      data cdi(4,160) / 2.23000E-01/
      data cdi(5,160) /-5.00000E-03/
      data cdi(6,160) / 1.80000E-02/
      data cdi(7,160) / 1.70000E-02/
      data cea(1,160) / 2.39320E+02/
      data cea(2,160) / 2.20000E+02/
      data cea(3,160) / 3.73900E+00/
      data cea(4,160) / 6.81700E+00/
      data cea(5,160) /-1.36650E+01/
      data cea(6,160) / 7.35700E+00/
      data cea(7,160) /-1.17000E-01/
      data cea(8,160) /-1.39000E-01/
C     Z=18 [+ 7]
      data error(161) /200.0/
      data cdi(1,161) / 1.43500E+02/
      data cdi(2,161) / 9.07000E-01/
      data cdi(3,161) / 1.05000E-01/
      data cdi(4,161) / 1.47000E-01/
      data cdi(5,161) /-7.30000E-02/
      data cdi(6,161) /-1.00000E-03/
C     Z=18 [+ 8]
      data error(162) / 40.0/
      data cdi(1,162) / 4.22500E+02/
      data cdi(2,162) / 3.51200E+00/
      data cdi(3,162) /-3.57000E-01/
      data cdi(4,162) / 1.15800E+00/
      data cdi(5,162) /-7.76000E-01/
C     Z=18 [+ 9]
      data error(163) / 70.0/
      data cdi(1,163) / 4.78700E+02/
      data cdi(2,163) / 2.79400E+00/
      data cdi(3,163) / 4.69000E-01/
      data cdi(4,163) /-1.29400E+01/
      data cdi(5,163) / 2.62600E+01/
      data cdi(6,163) /-1.34300E+01/
C     Z=18 [+10]
      data error(164) / 60.0/
      data cdi(1,164) / 5.39000E+02/
      data cdi(2,164) / 2.01900E+00/
      data cdi(3,164) /-1.32000E+00/
      data cdi(4,164) / 1.70000E+00/
C     Z=18 [+11]
      data error(165) / 40.0/
      data cdi(1,165) / 6.18200E+02/
      data cdi(2,165) / 7.85000E-01/
      data cdi(3,165) / 1.70900E+00/
      data cdi(4,165) /-1.08500E+01/
      data cdi(5,165) / 4.15000E+01/
      data cdi(6,165) /-5.80500E+01/
      data cdi(7,165) / 3.07200E+01/
C     Z=18 [+12]
      data error(166) / 40.0/
      data cdi(1,166) / 6.86100E+02/
      data cdi(2,166) / 1.06600E+00/
      data cdi(3,166) / 4.42000E-01/
      data cdi(4,166) / 4.75000E-01/
      data cdi(5,166) /-2.96100E+00/
      data cdi(6,166) / 4.47000E+00/
C     Z=18 [+13]
      data error(167) / 60.0/
      data cdi(1,167) / 7.55700E+02/
      data cdi(2,167) / 1.20000E+00/
      data cdi(3,167) /-6.52000E-01/
      data cdi(4,167) / 1.29900E+00/
C     Z=18 [+14]
      data error(168) / 20.0/
      data cdi(1,168) / 8.54800E+02/
      data cdi(2,168) / 6.68000E-01/
      data cdi(3,168) / 1.97000E-01/
      data cdi(4,168) / 1.20000E-01/
      data cdi(5,168) /-1.78000E-01/
C     Z=18 [+15]
      data error(169) / 80.0/
      data cdi(1,169) / 9.18000E+02/
      data cdi(2,169) / 3.36000E-01/
      data cdi(3,169) / 8.00000E-02/
      data cdi(4,169) / 1.43000E-01/
      data cdi(5,169) /-7.31000E-01/
      data cdi(6,169) / 1.33600E+00/
      data cdi(7,169) /-7.85000E-01/
C     Z=18 [+16]
      data error(170) / 50.0/
      data cdi(1,170) / 4.12080E+03/
      data cdi(2,170) / 7.96000E-01/
      data cdi(3,170) /-5.00000E-01/
      data cdi(4,170) / 8.84000E-01/
C     Z=18 [+17]
      data error(171) / 50.0/
      data cdi(1,171) / 4.42610E+03/
      data cdi(2,171) / 4.00000E-01/
C     Z=19 [+ 0]
      data error(172) / 20.0/
      data cdi(1,172) / 4.30000E+00/
      data cdi(2,172) / 3.73000E-01/
      data cdi(3,172) /-3.49000E-01/
      data cdi(4,172) / 2.71200E+00/
      data cdi(5,172) /-5.25300E+00/
      data cdi(6,172) / 2.87100E+00/
      data cea(1,172) / 2.02800E+01/
      data cea(2,172) / 1.50000E+01/
      data cea(3,172) / 3.94900E+00/
      data cea(4,172) / 7.06100E+00/
      data cea(5,172) /-1.63400E+01/
      data cea(6,172) / 8.11800E+00/
      data cea(7,172) / 1.25700E+00/
C     Z=19 [+ 1]
      data error(173) / 10.0/
      data cdi(1,173) / 3.18000E+01/
      data cdi(2,173) / 2.71400E+00/
      data cdi(3,173) /-4.04000E-01/
      data cdi(4,173) / 2.45600E+00/
      data cdi(5,173) /-1.85500E+00/
      data cdi(6,173) / 2.39000E-01/
C     Z=19 [+ 2]
      data error(174) / 30.0/
      data cdi(1,174) / 4.58000E+01/
      data cdi(2,174) / 1.66200E+00/
      data cdi(3,174) /-4.00000E-02/
      data cdi(4,174) /-2.77000E-01/
      data cdi(5,174) / 2.88000E-01/
C     Z=19 [+ 3]
      data error(175) / 50.0/
      data cdi(1,175) / 6.09000E+01/
      data cdi(2,175) / 2.08600E+00/
      data cdi(3,175) / 1.07700E+00/
      data cdi(4,175) /-2.17200E+00/
      data cdi(5,175) / 8.09000E-01/
C     Z=19 [+ 4]
      data error(176) / 60.0/
      data cdi(1,176) / 8.27000E+01/
      data cdi(2,176) / 1.75500E+00/
      data cdi(3,176) /-2.70000E-01/
      data cdi(4,176) / 7.72000E-01/
      data cdi(5,176) /-3.47000E-01/
C     Z=19 [+ 5]
      data error(177) / 60.0/
      data cdi(1,177) / 1.00000E+02/
      data cdi(2,177) / 1.24900E+00/
      data cdi(3,177) /-1.74000E-01/
      data cdi(4,177) / 1.03200E+00/
      data cdi(5,177) /-5.75000E-01/
C     Z=19 [+ 6]
      data error(178) / 60.0/
      data cdi(1,178) / 1.17600E+02/
      data cdi(2,178) / 9.60000E-01/
      data cdi(3,178) /-2.49000E-01/
      data cdi(4,178) / 6.87000E-01/
      data cdi(5,178) /-3.26000E-01/
C     Z=19 [+ 7]
      data error(179) /100.0/
      data cdi(1,179) / 1.54900E+02/
      data cdi(2,179) / 1.00800E+00/
      data cdi(3,179) / 2.03000E-01/
      data cdi(4,179) /-2.43000E-01/
      data cdi(5,179) / 9.90000E-02/
C     Z=19 [+ 8]
      data error(180) /200.0/
      data cdi(1,180) / 1.75800E+02/
      data cdi(2,180) / 9.07000E-01/
      data cdi(3,180) / 1.05000E-01/
      data cdi(4,180) / 1.47000E-01/
      data cdi(5,180) /-7.30000E-02/
      data cdi(6,180) /-1.00000E-03/
C     Z=19 [+ 9]
      data error(181) / 50.0/
      data cdi(1,181) / 5.03400E+02/
      data cdi(2,181) / 3.51200E+00/
      data cdi(3,181) /-3.57000E-01/
      data cdi(4,181) / 1.15800E+00/
      data cdi(5,181) /-7.76000E-01/
C     Z=19 [+10]
      data error(182) / 70.0/
      data cdi(1,182) / 5.64100E+02/
      data cdi(2,182) / 2.79400E+00/
      data cdi(3,182) / 4.69000E-01/
      data cdi(4,182) /-1.29400E+01/
      data cdi(5,182) / 2.62600E+01/
      data cdi(6,182) /-1.34300E+01/
C     Z=19 [+11]
      data error(183) / 60.0/
      data cdi(1,183) / 6.29100E+02/
      data cdi(2,183) / 2.01900E+00/
      data cdi(3,183) /-1.32000E+00/
      data cdi(4,183) / 1.70000E+00/
C     Z=19 [+12]
      data error(184) / 40.0/
      data cdi(1,184) / 7.14000E+02/
      data cdi(2,184) / 7.85000E-01/
      data cdi(3,184) / 1.70900E+00/
      data cdi(4,184) /-1.08500E+01/
      data cdi(5,184) / 4.15000E+01/
      data cdi(6,184) /-5.80500E+01/
      data cdi(7,184) / 3.07200E+01/
C     Z=19 [+13]
      data error(185) / 40.0/
      data cdi(1,185) / 7.87100E+02/
      data cdi(2,185) / 1.06600E+00/
      data cdi(3,185) / 4.42000E-01/
      data cdi(4,185) / 4.75000E-01/
      data cdi(5,185) /-2.96100E+00/
      data cdi(6,185) / 4.47000E+00/
C     Z=19 [+14]
      data error(186) / 60.0/
      data cdi(1,186) / 8.61800E+02/
      data cdi(2,186) / 1.20000E+00/
      data cdi(3,186) /-6.52000E-01/
      data cdi(4,186) / 1.29900E+00/
C     Z=19 [+15]
      data error(187) / 40.0/
      data cdi(1,187) / 9.68000E+02/
      data cdi(2,187) / 6.68000E-01/
      data cdi(3,187) / 1.97000E-01/
      data cdi(4,187) / 1.20000E-01/
      data cdi(5,187) /-1.78000E-01/
C     Z=19 [+16]
      data error(188) / 80.0/
      data cdi(1,188) / 1.03400E+03/
      data cdi(2,188) / 3.36000E-01/
      data cdi(3,188) / 8.00000E-02/
      data cdi(4,188) / 1.43000E-01/
      data cdi(5,188) /-7.31000E-01/
      data cdi(6,188) / 1.33600E+00/
      data cdi(7,188) /-7.85000E-01/
C     Z=19 [+17]
      data error(189) / 50.0/
      data cdi(1,189) / 4.61100E+03/
      data cdi(2,189) / 7.96000E-01/
      data cdi(3,189) /-5.00000E-01/
      data cdi(4,189) / 8.84000E-01/
C     Z=19 [+18]
      data error(190) / 50.0/
      data cdi(1,190) / 4.93390E+03/
      data cdi(2,190) / 4.00000E-01/
C     Z=20 [+ 0]
      data error(191) / 40.0/
      data cdi(1,191) / 6.10000E+00/
      data cdi(2,191) / 4.82000E-01/
      data cdi(3,191) /-4.30000E-02/
      data cdi(4,191) / 2.01200E+00/
      data cdi(5,191) /-1.31000E+00/
C     Z=20 [+ 1]
      data error(192) / 15.0/
      data cdi(1,192) / 1.19000E+01/
      data cdi(2,192) / 8.80000E-02/
      data cdi(3,192) / 4.00000E-02/
      data cdi(4,192) / 5.49000E-01/
      data cdi(5,192) /-4.97000E-01/
      data cea(1,192) / 2.22200E+01/
      data cea(2,192) / 1.40000E+01/
      data cea(3,192) / 8.68000E-01/
      data cea(4,192) / 9.46000E-01/
      data cea(5,192) /-2.23300E+00/
      data cea(6,192) /-3.04000E+00/
      data cea(7,192) / 6.35100E+00/
C     Z=20 [+ 2]
      data error(193) / 50.0/
      data cdi(1,193) / 5.09000E+01/
      data cdi(2,193) / 2.85700E+00/
      data cdi(3,193) /-6.57000E-01/
      data cdi(4,193) / 1.21400E+00/
      data cdi(5,193) /-6.89000E-01/
C     Z=20 [+ 3]
      data error(194) / 50.0/
      data cdi(1,194) / 6.71000E+01/
      data cdi(2,194) / 2.53400E+00/
      data cdi(3,194) / 1.37700E+00/
      data cdi(4,194) /-5.72100E+00/
      data cdi(5,194) / 7.13200E+00/
      data cdi(6,194) /-2.73100E+00/
C     Z=20 [+ 4]
      data error(195) / 50.0/
      data cdi(1,195) / 8.44000E+01/
      data cdi(2,195) / 1.03500E+00/
      data cdi(3,195) / 1.84600E+00/
      data cdi(4,195) / 6.57500E+00/
      data cdi(5,195) /-3.53400E+01/
      data cdi(6,195) / 5.35800E+01/
      data cdi(7,195) /-2.48400E+01/
C     Z=20 [+ 5]
      data error(196) / 60.0/
      data cdi(1,196) / 1.08800E+02/
      data cdi(2,196) / 1.64700E+00/
      data cdi(3,196) /-2.70000E-01/
      data cdi(4,196) / 7.72000E-01/
      data cdi(5,196) /-3.47000E-01/
C     Z=20 [+ 6]
      data error(197) / 60.0/
      data cdi(1,197) / 1.27700E+02/
      data cdi(2,197) / 1.14900E+00/
      data cdi(3,197) /-1.74000E-01/
      data cdi(4,197) / 1.03200E+00/
      data cdi(5,197) /-5.75000E-01/
C     Z=20 [+ 7]
      data error(198) / 60.0/
      data cdi(1,198) / 1.47200E+02/
      data cdi(2,198) / 9.60000E-01/
      data cdi(3,198) /-2.49000E-01/
      data cdi(4,198) / 6.87000E-01/
      data cdi(5,198) /-3.26000E-01/
C     Z=20 [+ 8]
      data error(199) /100.0/
      data cdi(1,199) / 1.88500E+02/
      data cdi(2,199) / 1.00800E+00/
      data cdi(3,199) / 2.03000E-01/
      data cdi(4,199) /-2.43000E-01/
      data cdi(5,199) / 9.90000E-02/
C     Z=20 [+ 9]
      data error(200) /200.0/
      data cdi(1,200) / 2.11300E+02/
      data cdi(2,200) / 9.07000E-01/
      data cdi(3,200) / 1.05000E-01/
      data cdi(4,200) / 1.47000E-01/
      data cdi(5,200) /-7.30000E-02/
      data cdi(6,200) /-1.00000E-03/
C     Z=20 [+10]
      data error(201) / 50.0/
      data cdi(1,201) / 5.91300E+02/
      data cdi(2,201) / 3.50000E+00/
      data cdi(3,201) /-2.36000E+00/
      data cdi(4,201) / 2.35000E+00/
      data cdi(5,201) /-6.60000E-01/
C     Z=20 [+11]
      data error(202) / 70.0/
      data cdi(1,202) / 6.56400E+02/
      data cdi(2,202) / 2.79400E+00/
      data cdi(3,202) / 4.69000E-01/
      data cdi(4,202) /-1.29400E+01/
      data cdi(5,202) / 2.62600E+01/
      data cdi(6,202) /-1.34300E+01/
C     Z=20 [+12]
      data error(203) / 60.0/
      data cdi(1,203) / 7.26000E+02/
      data cdi(2,203) / 2.01900E+00/
      data cdi(3,203) /-1.32000E+00/
      data cdi(4,203) / 1.70000E+00/
C     Z=20 [+13]
      data error(204) / 40.0/
      data cdi(1,204) / 8.16600E+02/
      data cdi(2,204) / 7.85000E-01/
      data cdi(3,204) / 1.70900E+00/
      data cdi(4,204) /-1.08500E+01/
      data cdi(5,204) / 4.15000E+01/
      data cdi(6,204) /-5.80500E+01/
      data cdi(7,204) / 3.07200E+01/
C     Z=20 [+14]
      data error(205) / 40.0/
      data cdi(1,205) / 8.95100E+02/
      data cdi(2,205) / 1.06600E+00/
      data cdi(3,205) / 4.42000E-01/
      data cdi(4,205) / 4.75000E-01/
      data cdi(5,205) /-2.96100E+00/
      data cdi(6,205) / 4.47000E+00/
C     Z=20 [+15]
      data error(206) / 60.0/
      data cdi(1,206) / 9.74000E+02/
      data cdi(2,206) / 1.20000E+00/
      data cdi(3,206) /-6.52000E-01/
      data cdi(4,206) / 1.29900E+00/
C     Z=20 [+16]
      data error(207) / 40.0/
      data cdi(1,207) / 1.08700E+03/
      data cdi(2,207) / 6.68000E-01/
      data cdi(3,207) / 1.97000E-01/
      data cdi(4,207) / 1.20000E-01/
      data cdi(5,207) /-1.78000E-01/
C     Z=20 [+17]
      data error(208) / 80.0/
      data cdi(1,208) / 1.15700E+03/
      data cdi(2,208) / 3.36000E-01/
      data cdi(3,208) / 8.00000E-02/
      data cdi(4,208) / 1.43000E-01/
      data cdi(5,208) /-7.31000E-01/
      data cdi(6,208) / 1.33600E+00/
      data cdi(7,208) /-7.85000E-01/
C     Z=20 [+18]
      data error(209) / 50.0/
      data cdi(1,209) / 5.12900E+03/
      data cdi(2,209) / 7.96000E-01/
      data cdi(3,209) /-5.00000E-01/
      data cdi(4,209) / 8.84000E-01/
C     Z=20 [+19]
      data error(210) / 50.0/
      data cdi(1,210) / 5.46970E+03/
      data cdi(2,210) / 4.00000E-01/
C     Z=21 [+ 0]
      data error(211) /200.0/
      data cdi(1,211) / 6.50000E+00/
      data cdi(2,211) / 2.88800E+00/
      data cdi(3,211) / 7.53000E-01/
      data cdi(4,211) /-2.12400E+01/
      data cdi(5,211) / 6.62900E+01/
      data cdi(6,211) /-7.52300E+01/
      data cdi(7,211) / 2.90800E+01/
C     Z=21 [+ 1]
      data error(212) /200.0/
      data cdi(1,212) / 1.28000E+01/
      data cdi(2,212) / 2.87600E+00/
      data cdi(3,212) /-2.32100E+00/
      data cdi(4,212) / 5.91000E+00/
      data cdi(5,212) /-1.03060E+01/
      data cdi(6,212) /-1.33300E+00/
      data cdi(7,212) / 7.73800E+00/
C     Z=21 [+ 2]
      data error(213) /300.0/
      data cdi(1,213) / 2.48000E+01/
      data cdi(2,213) / 2.36900E+00/
      data cdi(3,213) / 5.31000E-01/
      data cdi(4,213) /-1.32300E+00/
      data cdi(5,213) /-7.77800E+00/
      data cdi(6,213) / 1.49200E+01/
      data cdi(7,213) /-6.09200E+00/
C     Z=21 [+ 3]
      data error(214) / 30.0/
      data cdi(1,214) / 7.35000E+01/
      data cdi(2,214) / 2.62200E+00/
      data cdi(3,214) /-1.25600E+00/
      data cdi(4,214) / 2.21500E+00/
      data cdi(5,214) /-1.05000E+00/
C     Z=21 [+ 4]
      data error(215) / 30.0/
      data cdi(1,215) / 9.17000E+01/
      data cdi(2,215) / 2.40200E+00/
      data cdi(3,215) /-1.02800E+00/
      data cdi(4,215) / 1.49500E+00/
      data cdi(5,215) /-4.20000E-01/
C     Z=21 [+ 5]
      data error(216) / 30.0/
      data cdi(1,216) / 1.11100E+02/
      data cdi(2,216) / 1.84200E+00/
      data cdi(3,216) /-1.23000E-01/
      data cdi(4,216) / 9.08000E-01/
      data cdi(5,216) /-4.39000E-01/
C     Z=21 [+ 6]
      data error(217) / 30.0/
      data cdi(1,217) / 1.38000E+02/
      data cdi(2,217) / 1.75500E+00/
      data cdi(3,217) /-2.70000E-01/
      data cdi(4,217) / 7.72000E-01/
      data cdi(5,217) /-3.47000E-01/
C     Z=21 [+ 7]
      data error(218) / 30.0/
      data cdi(1,218) / 1.58700E+02/
      data cdi(2,218) / 1.24900E+00/
      data cdi(3,218) /-1.74000E-01/
      data cdi(4,218) / 1.03200E+00/
      data cdi(5,218) /-5.75000E-01/
C     Z=21 [+ 8]
      data error(219) / 30.0/
      data cdi(1,219) / 1.80000E+02/
      data cdi(2,219) / 9.60000E-01/
      data cdi(3,219) /-2.49000E-01/
      data cdi(4,219) / 6.87000E-01/
      data cdi(5,219) /-3.26000E-01/
C     Z=21 [+ 9]
      data error(220) / 30.0/
      data cdi(1,220) / 2.25300E+02/
      data cdi(2,220) / 1.00800E+00/
      data cdi(3,220) / 2.03000E-01/
      data cdi(4,220) /-2.43000E-01/
      data cdi(5,220) / 9.90000E-02/
C     Z=21 [+10]
      data error(221) /200.0/
      data cdi(1,221) / 2.49800E+02/
      data cdi(2,221) / 9.07000E-01/
      data cdi(3,221) / 1.05000E-01/
      data cdi(4,221) / 1.47000E-01/
      data cdi(5,221) /-7.30000E-02/
      data cdi(6,221) /-1.00000E-03/
C     Z=21 [+11]
      data error(222) / 50.0/
      data cdi(1,222) / 6.85900E+02/
      data cdi(2,222) / 3.44100E+00/
      data cdi(3,222) / 1.06600E+00/
      data cdi(4,222) /-2.77400E+00/
      data cdi(5,222) / 1.48900E+00/
C     Z=21 [+12]
      data error(223) / 70.0/
      data cdi(1,223) / 7.55500E+02/
      data cdi(2,223) / 2.49200E+00/
      data cdi(3,223) / 1.04300E+00/
      data cdi(4,223) /-5.28000E-01/
      data cdi(5,223) /-1.13000E-01/
C     Z=21 [+13]
      data error(224) / 60.0/
      data cdi(1,224) / 8.29800E+02/
      data cdi(2,224) / 1.91200E+00/
      data cdi(3,224) / 1.07900E+00/
      data cdi(4,224) /-3.50000E-02/
      data cdi(5,224) /-4.47000E-01/
C     Z=21 [+14]
      data error(225) / 60.0/
      data cdi(1,225) / 9.26000E+02/
      data cdi(2,225) / 1.63300E+00/
      data cdi(3,225) / 6.45000E-01/
      data cdi(4,225) /-1.01000E-01/
      data cdi(5,225) /-2.39000E-01/
C     Z=21 [+15]
      data error(226) / 60.0/
      data cdi(1,226) / 1.00900E+03/
      data cdi(2,226) / 1.09100E+00/
      data cdi(3,226) / 2.76000E-01/
      data cdi(4,226) / 1.05200E+00/
      data cdi(5,226) /-8.92000E-01/
C     Z=21 [+16]
      data cdi(1,227) / 1.09400E+03/
      data cdi(2,227) / 1.89000E-01/
      data cdi(3,227) / 3.27000E-01/
      data cdi(4,227) / 1.67900E+00/
      data cdi(5,227) /-1.16000E+00/
C     Z=21 [+17]
      data error(228) / 60.0/
      data cdi(1,228) / 1.21300E+03/
      data cdi(2,228) / 4.14000E-01/
      data cdi(3,228) / 9.70000E-02/
      data cdi(4,228) / 1.37000E-01/
      data cdi(5,228) /-1.44000E-01/
C     Z=21 [+18]
      data cdi(1,229) / 1.28800E+03/
      data cdi(2,229) / 2.89000E-01/
      data cdi(3,229) / 1.10000E-01/
      data cdi(4,229) /-4.41000E-01/
      data cdi(5,229) / 2.69000E-01/
C     Z=21 [+19]
      data error(230) / 70.0/
      data cdi(1,230) / 5.67530E+03/
      data cdi(2,230) / 5.25000E-01/
      data cdi(3,230) /-1.88000E-01/
      data cdi(4,230) / 4.89000E-01/
      data cdi(5,230) /-2.92000E-01/
C     Z=21 [+20]
      data error(231) / 50.0/
      data cdi(1,231) / 6.03380E+03/
      data cdi(2,231) / 2.60000E-01/
      data cdi(3,231) /-2.20000E-02/
      data cdi(4,231) / 5.10000E-02/
      data cdi(5,231) /-2.30000E-02/
C     Z=22 [+ 0]
      data error(232) / 40.0/
      data cdi(1,232) / 6.80000E+00/
      data cdi(2,232) / 8.27000E-01/
      data cdi(3,232) / 1.27000E+00/
      data cdi(4,232) /-8.80600E+00/
      data cdi(5,232) / 2.45200E+01/
      data cdi(6,232) /-2.67600E+01/
      data cdi(7,232) / 1.00400E+01/
C     Z=22 [+ 1]
      data error(233) / 10.0/
      data cdi(1,233) / 1.36000E+01/
      data cdi(2,233) / 6.18000E-01/
      data cdi(3,233) / 4.80000E-02/
      data cdi(4,233) / 1.85000E-01/
      data cdi(5,233) /-7.40000E-01/
      data cdi(6,233) / 5.22000E-01/
      data cea(1,233) / 3.25000E+01/
      data cea(2,233) / 2.80000E+01/
      data cea(3,233) / 4.02600E+00/
      data cea(4,233) / 8.92000E+00/
      data cea(5,233) /-3.00700E+01/
      data cea(6,233) / 3.17300E+01/
      data cea(7,233) /-1.79700E+01/
      data cea(8,233) / 5.94500E+00/
C     Z=22 [+ 2]
      data error(234) /100.0/
      data cdi(1,234) / 2.75000E+01/
      data cdi(2,234) / 2.87600E+00/
      data cdi(3,234) /-2.32100E+00/
      data cdi(4,234) / 5.91000E+00/
      data cdi(5,234) /-1.03100E+01/
      data cdi(6,234) /-1.33300E+00/
      data cdi(7,234) / 7.73800E+00/
C     Z=22 [+ 3]
      data error(235) / 10.0/
      data cdi(1,235) / 4.31000E+01/
      data cdi(2,235) / 3.57000E+00/
      data cdi(3,235) / 5.31000E-01/
      data cdi(4,235) /-1.32300E+00/
      data cdi(5,235) /-7.77800E+00/
      data cdi(6,235) / 1.49200E+01/
      data cdi(7,235) /-6.09200E+00/
C     Z=22 [+ 4]
      data error(236) /100.0/
      data cdi(1,236) / 9.92000E+01/
      data cdi(2,236) / 2.85700E+00/
      data cdi(3,236) / 2.09000E-01/
      data cdi(4,236) / 4.88000E-01/
      data cdi(5,236) /-3.79000E-01/
C     Z=22 [+ 5]
      data error(237) /100.0/
      data cdi(1,237) / 1.19400E+02/
      data cdi(2,237) / 2.45700E+00/
      data cdi(3,237) / 2.89000E-01/
      data cdi(4,237) / 4.84000E-01/
      data cdi(5,237) /-4.13000E-01/
C     Z=22 [+ 6]
      data error(238) /100.0/
      data cdi(1,238) / 1.40800E+02/
      data cdi(2,238) / 2.03400E+00/
      data cdi(3,238) / 4.10000E-01/
      data cdi(4,238) / 2.53000E-01/
      data cdi(5,238) /-3.21100E-01/
C     Z=22 [+ 7]
      data error(239) / 50.0/
      data cdi(1,239) / 1.68500E+02/
      data cdi(2,239) / 3.15900E+00/
      data cdi(3,239) /-8.78000E-01/
      data cdi(4,239) /-3.01800E+00/
      data cdi(5,239) / 3.19600E+00/
      data cdi(6,239) /-2.44000E+00/
      data cdi(7,239) /-9.60000E+00/
      data cea(1,239) / 3.56160E+02/
      data cea(2,239) / 3.05600E+02/
      data cea(3,239) / 2.80000E-02/
      data cea(4,239) / 2.82280E+01/
      data cea(5,239) /-4.82170E+01/
      data cea(6,239) / 3.24560E+01/
      data cea(7,239) / 4.38000E-01/
      data cea(8,239) / 6.47000E-01/
C     Z=22 [+ 8]
      data error(240) /100.0/
      data cdi(1,240) / 1.93200E+02/
      data cdi(2,240) / 1.93600E+00/
      data cdi(3,240) / 5.42000E-01/
      data cdi(4,240) / 5.95000E-01/
      data cdi(5,240) /-5.18000E-01/
C     Z=22 [+ 9]
      data error(241) /100.0/
      data cdi(1,241) / 2.15900E+02/
      data cdi(2,241) / 7.72000E-01/
      data cdi(3,241) / 4.72000E-01/
      data cdi(4,241) / 6.95000E-01/
      data cdi(5,241) /-8.25000E-01/
      data cdi(6,241) /-1.80000E-02/
      data cea(1,241) / 3.78400E+02/
      data cea(2,241) / 3.63300E+02/
      data cea(3,241) / 2.50100E+00/
      data cea(4,241) / 2.83120E+01/
      data cea(5,241) /-5.54650E+01/
      data cea(6,241) / 2.74410E+01/
      data cea(7,241) / 1.65100E+00/
      data cea(8,241) / 2.23200E+00/
C     Z=22 [+10]
      data error(242) /200.0/
      data cdi(1,242) / 2.65200E+02/
      data cdi(2,242) / 2.13300E+00/
      data cdi(3,242) / 1.91000E-01/
      data cdi(4,242) /-1.50000E-01/
      data cdi(5,242) / 3.20000E-02/
C     Z=22 [+11]
      data error(243) /200.0/
      data cdi(1,243) / 2.91500E+02/
      data cdi(2,243) / 1.50700E+00/
      data cdi(3,243) / 1.05000E-01/
      data cdi(4,243) / 1.47000E-01/
      data cdi(5,243) /-7.30000E-02/
      data cdi(6,243) /-1.00000E-03/
C     Z=22 [+12]
      data error(244) / 50.0/
      data cdi(1,244) / 7.87300E+02/
      data cdi(2,244) / 3.44100E+00/
      data cdi(3,244) / 1.06600E+00/
      data cdi(4,244) /-2.77400E+00/
      data cdi(5,244) / 1.48900E+00/
C     Z=22 [+13]
      data error(245) / 70.0/
      data cdi(1,245) / 8.61300E+02/
      data cdi(2,245) / 2.49200E+00/
      data cdi(3,245) / 1.04300E+00/
      data cdi(4,245) /-5.28000E-01/
      data cdi(5,245) /-1.13000E-01/
C     Z=22 [+14]
      data error(246) / 60.0/
      data cdi(1,246) / 9.40400E+02/
      data cdi(2,246) / 1.91200E+00/
      data cdi(3,246) / 1.07900E+00/
      data cdi(4,246) /-3.50000E-02/
      data cdi(5,246) /-4.47000E-01/
C     Z=22 [+15]
      data error(247) / 60.0/
      data cdi(1,247) / 1.04200E+03/
      data cdi(2,247) / 1.63300E+00/
      data cdi(3,247) / 6.45000E-01/
      data cdi(4,247) /-1.01000E-01/
      data cdi(5,247) /-2.39000E-01/
C     Z=22 [+16]
      data error(248) / 60.0/
      data cdi(1,248) / 1.13100E+03/
      data cdi(2,248) / 1.09100E+00/
      data cdi(3,248) / 2.76000E-01/
      data cdi(4,248) / 1.05200E+00/
      data cdi(5,248) /-8.92000E-01/
C     Z=22 [+17]
      data error(249) /100.0/
      data cdi(1,249) / 1.22000E+03/
      data cdi(2,249) / 1.89000E-01/
      data cdi(3,249) / 3.27000E-01/
      data cdi(4,249) / 1.67900E+00/
      data cdi(5,249) /-1.16000E+00/
C     Z=22 [+18]
      data error(250) / 60.0/
      data cdi(1,250) / 1.34200E+03/
      data cdi(2,250) / 4.14000E-01/
      data cdi(3,250) / 9.70000E-02/
      data cdi(4,250) / 1.37000E-01/
      data cdi(5,250) /-1.44000E-01/
C     Z=22 [+19]
      data error(251) /100.0/
      data cdi(1,251) / 1.42500E+03/
      data cdi(2,251) / 2.89000E-01/
      data cdi(3,251) / 1.10000E-01/
      data cdi(4,251) /-4.41000E-01/
      data cdi(5,251) / 2.69000E-01/
C     Z=22 [+20]
      data error(252) / 70.0/
      data cdi(1,252) / 6.24900E+03/
      data cdi(2,252) / 5.25000E-01/
      data cdi(3,252) /-1.88000E-01/
      data cdi(4,252) / 4.89000E-01/
      data cdi(5,252) /-2.92000E-01/
C     Z=22 [+21]
      data error(253) / 50.0/
      data cdi(1,253) / 6.62560E+03/
      data cdi(2,253) / 2.60000E-01/
      data cdi(3,253) /-2.20000E-02/
      data cdi(4,253) / 5.10000E-02/
      data cdi(5,253) /-2.30000E-02/
C     Z=23 [+ 0]
      data error(254) /100.0/
      data cdi(1,254) / 6.70000E+00/
      data cdi(2,254) / 2.25800E+00/
      data cdi(3,254) /-8.43000E-01/
      data cdi(4,254) / 6.22000E-01/
      data cdi(5,254) / 1.00000E-02/
      data cea(1,254) / 1.19400E+01/
      data cea(2,254) / 9.80000E+00/
      data cea(3,254) / 6.11300E+00/
      data cea(4,254) /-1.14800E+00/
      data cea(5,254) / 1.58300E+00/
      data cea(6,254) /-3.40000E-02/
C     Z=23 [+ 1]
      data error(255) /100.0/
      data cdi(1,255) / 1.46000E+01/
      data cdi(2,255) / 2.65900E+00/
      data cdi(3,255) /-1.06800E+00/
      data cdi(4,255) /-3.96000E-01/
      data cdi(5,255) / 1.03600E+00/
      data cea(1,255) / 2.43600E+01/
      data cea(2,255) / 1.95000E+01/
      data cea(3,255) / 6.50300E+00/
      data cea(4,255) /-2.78000E-01/
      data cea(5,255) /-5.66300E+00/
      data cea(6,255) / 3.86000E+00/
C     Z=23 [+ 2]
      data error(256) /200.0/
      data cdi(1,256) / 2.93000E+01/
      data cdi(2,256) / 2.88800E+00/
      data cdi(3,256) / 7.53000E-01/
      data cdi(4,256) /-2.12400E+01/
      data cdi(5,256) / 6.62900E+01/
      data cdi(6,256) /-7.52300E+01/
      data cdi(7,256) / 2.90800E+01/
C     Z=23 [+ 3]
      data error(257) /200.0/
      data cdi(1,257) / 4.67000E+01/
      data cdi(2,257) / 2.87600E+00/
      data cdi(3,257) /-2.32100E+00/
      data cdi(4,257) / 5.91000E+00/
      data cdi(5,257) /-1.03060E+01/
      data cdi(6,257) /-1.33300E+00/
      data cdi(7,257) / 7.73800E+00/
C     Z=23 [+ 4]
      data error(258) /300.0/
      data cdi(1,258) / 6.52000E+01/
      data cdi(2,258) / 2.36900E+00/
      data cdi(3,258) / 5.31000E-01/
      data cdi(4,258) /-1.32300E+00/
      data cdi(5,258) /-7.77800E+00/
      data cdi(6,258) / 1.49200E+01/
      data cdi(7,258) /-6.09200E+00/
C     Z=23 [+ 5]
      data error(259) / 30.0/
      data cdi(1,259) / 1.28100E+02/
      data cdi(2,259) / 2.96300E+00/
      data cdi(3,259) /-3.03000E-01/
      data cdi(4,259) /-5.19000E-01/
      data cdi(5,259) / 7.67000E-01/
C     Z=23 [+ 6]
      data error(260) /100.0/
      data cdi(1,260) / 1.50200E+02/
      data cdi(2,260) / 2.45700E+00/
      data cdi(3,260) / 2.89000E-01/
      data cdi(4,260) / 4.84000E-01/
      data cdi(5,260) /-4.13000E-01/
C     Z=23 [+ 7]
      data error(261) /100.0/
      data cdi(1,261) / 1.73700E+02/
      data cdi(2,261) / 2.03400E+00/
      data cdi(3,261) / 4.09000E-01/
      data cdi(4,261) / 2.53000E-01/
      data cdi(5,261) / 3.21000E-01/
C     Z=23 [+ 8]
      data error(262) /100.0/
      data cdi(1,262) / 2.05800E+02/
      data cdi(2,262) / 3.15900E+00/
      data cdi(3,262) /-8.78000E-01/
      data cdi(4,262) /-3.01800E+00/
      data cdi(5,262) / 3.19600E+00/
      data cdi(6,262) /-2.44000E+00/
      data cdi(7,262) /-9.60000E+00/
      data cea(1,262) / 4.34600E+02/
      data cea(2,262) / 3.73000E+02/
      data cea(3,262) / 2.80000E-02/
      data cea(4,262) / 2.82280E+01/
      data cea(5,262) /-4.82170E+01/
      data cea(6,262) / 3.24560E+01/
      data cea(7,262) / 4.38000E-01/
      data cea(8,262) / 6.47000E-01/
C     Z=23 [+ 9]
      data error(263) /100.0/
      data cdi(1,263) / 2.30500E+02/
      data cdi(2,263) / 1.93600E+00/
      data cdi(3,263) / 5.42000E-01/
      data cdi(4,263) / 5.95000E-01/
      data cdi(5,263) /-5.18000E-01/
C     Z=23 [+10]
      data error(264) /100.0/
      data cdi(1,264) / 2.55000E+02/
      data cdi(2,264) / 7.72000E-01/
      data cdi(3,264) / 4.72000E-01/
      data cdi(4,264) / 6.95000E-01/
      data cdi(5,264) /-8.25000E-01/
      data cdi(6,264) /-1.80000E-02/
      data cea(1,264) / 4.48800E+02/
      data cea(2,264) / 4.28400E+02/
      data cea(3,264) / 2.50100E+00/
      data cea(4,264) / 2.83120E+01/
      data cea(5,264) /-5.54650E+01/
      data cea(6,264) / 2.74410E+01/
      data cea(7,264) / 1.65100E+00/
      data cea(8,264) / 2.23200E+00/
C     Z=23 [+11]
      data error(265) /200.0/
      data cdi(1,265) / 3.08300E+02/
      data cdi(2,265) / 2.13300E+00/
      data cdi(3,265) / 1.91000E-01/
      data cdi(4,265) /-1.50000E-01/
      data cdi(5,265) /-7.30000E-02/
      data cdi(6,265) / 3.20000E-02/
C     Z=23 [+12]
      data error(266) /200.0/
      data cdi(1,266) / 3.36300E+02/
      data cdi(2,266) / 1.50700E+00/
      data cdi(3,266) / 1.05000E-01/
      data cdi(4,266) / 1.47000E-01/
      data cdi(5,266) /-7.30000E-02/
      data cdi(6,266) /-1.00000E-03/
C     Z=23 [+13]
      data error(267) / 50.0/
      data cdi(1,267) / 8.95600E+02/
      data cdi(2,267) / 3.44100E+00/
      data cdi(3,267) / 1.06600E+00/
      data cdi(4,267) /-2.77400E+00/
      data cdi(5,267) / 1.48900E+00/
C     Z=23 [+14]
      data error(268) / 70.0/
      data cdi(1,268) / 9.74000E+02/
      data cdi(2,268) / 2.49200E+00/
      data cdi(3,268) / 1.04300E+00/
      data cdi(4,268) /-5.28000E-01/
      data cdi(5,268) /-1.13000E-01/
C     Z=23 [+15]
C     Z=23 [+16]
C     Z=23 [+17]
C     Z=23 [+18]
C     Z=23 [+19]
C     Z=23 [+20]
C     Z=23 [+21]
C     Z=23 [+22]
C     Z=24 [+ 0]
      data error(277) / 40.0/
      data cdi(1,277) / 6.80000E+00/
      data cdi(2,277) / 1.02900E+00/
      data cdi(3,277) /-4.29000E-01/
      data cdi(4,277) / 5.20000E-02/
      data cdi(5,277) / 7.40000E-02/
C     Z=24 [+ 1]
      data error(278) / 20.0/
      data cdi(1,278) / 1.50000E+01/
      data cdi(2,278) / 7.25000E-01/
      data cdi(3,278) /-3.44000E-01/
      data cdi(4,278) / 1.42000E-01/
      data cdi(5,278) / 6.70000E-02/
      data cea(1,278) / 4.00500E+01/
      data cea(2,278) / 3.50000E+01/
      data cea(3,278) / 5.37000E+00/
      data cea(4,278) / 5.85600E+00/
      data cea(5,278) /-1.94050E+01/
      data cea(6,278) / 8.23300E+00/
      data cea(7,278) / 6.07000E-01/
      data cea(8,278) / 2.76800E+00/
C     Z=24 [+ 2]
      data error(279) /100.0/
      data cdi(1,279) / 3.10000E+01/
      data cdi(2,279) / 2.65900E+00/
      data cdi(3,279) /-1.06800E+00/
      data cdi(4,279) /-3.96000E-01/
      data cdi(5,279) / 1.03600E+00/
C     Z=24 [+ 3]
      data error(280) /200.0/
      data cdi(1,280) / 4.91000E+01/
      data cdi(2,280) / 2.88800E+00/
      data cdi(3,280) / 7.53000E-01/
      data cdi(4,280) /-2.12400E+01/
      data cdi(5,280) / 6.62900E+01/
      data cdi(6,280) /-7.52300E+01/
      data cdi(7,280) / 2.90800E+01/
C     Z=24 [+ 4]
      data error(281) /200.0/
      data cdi(1,281) / 6.93000E+01/
      data cdi(2,281) / 2.87600E+00/
      data cdi(3,281) /-2.32100E+00/
      data cdi(4,281) / 5.91000E+00/
      data cdi(5,281) /-1.03060E+01/
      data cdi(6,281) /-1.33300E+00/
      data cdi(7,281) / 7.73800E+00/
C     Z=24 [+ 5]
      data error(282) /300.0/
      data cdi(1,282) / 9.06000E+01/
      data cdi(2,282) / 2.36900E+00/
      data cdi(3,282) / 5.31000E-01/
      data cdi(4,282) /-1.32300E+00/
      data cdi(5,282) /-7.77800E+00/
      data cdi(6,282) / 1.49200E+01/
      data cdi(7,282) /-6.09200E+00/
C     Z=24 [+ 6]
      data error(283) /100.0/
      data cdi(1,283) / 1.61100E+02/
      data cdi(2,283) / 2.85700E+00/
      data cdi(3,283) / 2.09000E-01/
      data cdi(4,283) / 4.88000E-01/
      data cdi(5,283) /-3.79000E-01/
C     Z=24 [+ 7]
      data error(284) /100.0/
      data cdi(1,284) / 1.84700E+02/
      data cdi(2,284) / 2.45700E+00/
      data cdi(3,284) / 2.89000E-01/
      data cdi(4,284) / 4.84000E-01/
      data cdi(5,284) /-4.13000E-01/
C     Z=24 [+ 8]
      data error(285) /100.0/
      data cdi(1,285) / 2.09300E+02/
      data cdi(2,285) / 2.03400E+00/
      data cdi(3,285) / 4.09000E-01/
      data cdi(4,285) / 2.53000E-01/
      data cdi(5,285) /-3.21000E-01/
C     Z=24 [+ 9]
      data error(286) /100.0/
      data cdi(1,286) / 2.44400E+02/
      data cdi(2,286) / 3.15900E+00/
      data cdi(3,286) /-8.78000E-01/
      data cdi(4,286) /-3.01800E+00/
      data cdi(5,286) / 3.19600E+00/
      data cdi(6,286) /-2.44000E+00/
      data cdi(7,286) /-9.60000E+00/
      data cea(1,286) / 5.17280E+02/
      data cea(2,286) / 4.42400E+02/
      data cea(3,286) / 2.80000E-02/
      data cea(4,286) / 2.82280E+01/
      data cea(5,286) /-4.82170E+01/
      data cea(6,286) / 3.24560E+01/
      data cea(7,286) / 4.38000E-01/
      data cea(8,286) / 6.47000E-01/
C     Z=24 [+10]
      data error(287) /100.0/
      data cdi(1,287) / 2.70800E+02/
      data cdi(2,287) / 1.93600E+00/
      data cdi(3,287) / 5.42000E-01/
      data cdi(4,287) / 5.95000E-01/
      data cdi(5,287) /-5.18000E-01/
C     Z=24 [+11]
      data error(288) /100.0/
      data cdi(1,288) / 2.98000E+02/
      data cdi(2,288) / 7.72000E-01/
      data cdi(3,288) / 4.72000E-01/
      data cdi(4,288) / 6.95000E-01/
      data cdi(5,288) /-8.25000E-01/
      data cdi(6,288) /-1.80000E-02/
      data cea(1,288) / 5.24480E+02/
      data cea(2,288) / 5.00600E+02/
      data cea(3,288) / 2.50100E+00/
      data cea(4,288) / 2.83120E+01/
      data cea(5,288) /-5.54650E+01/
      data cea(6,288) / 2.74410E+01/
      data cea(7,288) / 1.65100E+00/
      data cea(8,288) / 2.23200E+00/
C     Z=24 [+12]
      data error(289) /200.0/
      data cdi(1,289) / 3.55000E+02/
      data cdi(2,289) / 2.13300E+00/
      data cdi(3,289) / 1.91000E-01/
      data cdi(4,289) /-1.50000E-01/
      data cdi(5,289) / 3.20000E-02/
C     Z=24 [+13]
      data error(290) /200.0/
      data cdi(1,290) / 3.84300E+02/
      data cdi(2,290) / 1.50700E+00/
      data cdi(3,290) / 1.05000E-01/
      data cdi(4,290) / 1.47000E-01/
      data cdi(5,290) /-7.30000E-02/
      data cdi(6,290) /-1.00000E-03/
C     Z=24 [+14]
      data error(291) / 50.0/
      data cdi(1,291) / 1.01060E+03/
      data cdi(2,291) / 3.44100E+00/
      data cdi(3,291) / 1.06600E+00/
      data cdi(4,291) /-2.77400E+00/
      data cdi(5,291) / 1.48900E+00/
C     Z=24 [+15]
      data error(292) / 70.0/
      data cdi(1,292) / 1.09700E+03/
      data cdi(2,292) / 2.49200E+00/
      data cdi(3,292) / 1.04300E+00/
      data cdi(4,292) /-5.28000E-01/
      data cdi(5,292) /-1.13000E-01/
C     Z=24 [+16]
      data error(293) / 60.0/
      data cdi(1,293) / 1.18500E+03/
      data cdi(2,293) / 1.91200E+00/
      data cdi(3,293) / 1.07900E+00/
      data cdi(4,293) /-3.50000E-02/
      data cdi(5,293) /-4.47000E-01/
C     Z=24 [+17]
      data error(294) / 60.0/
      data cdi(1,294) / 1.29900E+03/
      data cdi(2,294) / 1.63300E+00/
      data cdi(3,294) / 6.45000E-01/
      data cdi(4,294) /-1.01000E-01/
      data cdi(5,294) /-2.39000E-01/
C     Z=24 [+18]
      data error(295) / 60.0/
      data cdi(1,295) / 1.39600E+03/
      data cdi(2,295) / 1.09100E+00/
      data cdi(3,295) / 2.76000E-01/
      data cdi(4,295) / 1.05200E+00/
      data cdi(5,295) /-8.92000E-01/
C     Z=24 [+19]
      data error(296) /100.0/
      data cdi(1,296) / 1.49600E+03/
      data cdi(2,296) / 1.89000E-01/
      data cdi(3,296) / 3.27000E-01/
      data cdi(4,296) / 1.67900E+00/
      data cdi(5,296) /-1.16000E+00/
C     Z=24 [+20]
      data error(297) / 60.0/
      data cdi(1,297) / 1.63400E+03/
      data cdi(2,297) / 4.14000E-01/
      data cdi(3,297) / 9.70000E-02/
      data cdi(4,297) / 1.37000E-01/
      data cdi(5,297) /-1.44000E-01/
C     Z=24 [+21]
      data error(298) /100.0/
      data cdi(1,298) / 1.72140E+03/
      data cdi(2,298) / 2.89000E-01/
      data cdi(3,298) / 1.10000E-01/
      data cdi(4,298) /-4.41000E-01/
      data cdi(5,298) / 2.69000E-01/
C     Z=24 [+22]
      data error(299) / 70.0/
      data cdi(1,299) / 7.48200E+03/
      data cdi(2,299) / 5.25000E-01/
      data cdi(3,299) /-1.88000E-01/
      data cdi(4,299) / 4.89000E-01/
      data cdi(5,299) /-2.92000E-01/
C     Z=24 [+23]
C     Z=25 [+ 0]
      data error(301) /200.0/
      data cdi(1,301) / 7.40000E+00/
      data cdi(2,301) / 2.12400E+00/
      data cdi(3,301) /-5.30000E-01/
      data cdi(4,301) /-9.61700E+00/
      data cdi(5,301) / 2.91700E+01/
      data cdi(6,301) /-4.08100E+01/
      data cdi(7,301) / 1.92900E+01/
C     Z=25 [+ 1]
      data error(302) /200.0/
      data cdi(1,302) / 1.56000E+01/
      data cdi(2,302) / 4.42000E+00/
      data cdi(3,302) /-4.81000E+00/
      data cdi(4,302) / 5.69200E+00/
      data cdi(5,302) /-1.29770E+01/
      data cdi(6,302) / 6.69800E+00/
C     Z=25 [+ 2]
      data error(303) /100.0/
      data cdi(1,303) / 3.37000E+01/
      data cdi(2,303) / 2.25800E+00/
      data cdi(3,303) /-8.43000E-01/
      data cdi(4,303) / 6.22000E-01/
      data cdi(5,303) / 1.00000E-02/
      data cea(1,303) / 6.53400E+01/
      data cea(2,303) / 4.92000E+01/
      data cea(3,303) / 6.11300E+00/
      data cea(4,303) /-1.14800E+00/
      data cea(5,303) / 1.58300E+00/
      data cea(6,303) /-3.40000E-02/
C     Z=25 [+ 3]
      data error(304) /100.0/
      data cdi(1,304) / 5.12000E+01/
      data cdi(2,304) / 2.65900E+00/
      data cdi(3,304) /-1.06800E+00/
      data cdi(4,304) /-3.96000E-01/
      data cdi(5,304) / 1.03600E+00/
      data cea(1,304) / 8.87400E+01/
      data cea(2,304) / 6.81000E+01/
      data cea(3,304) / 6.50300E+00/
      data cea(4,304) /-2.78000E-01/
      data cea(5,304) /-5.66300E+00/
      data cea(6,304) / 3.86000E+00/
C     Z=25 [+ 4]
      data error(305) /200.0/
      data cdi(1,305) / 7.24000E+01/
      data cdi(2,305) / 2.88800E+00/
      data cdi(3,305) / 7.53000E-01/
      data cdi(4,305) /-2.12400E+01/
      data cdi(5,305) / 6.62900E+01/
      data cdi(6,305) /-7.52300E+01/
      data cdi(7,305) / 2.90800E+01/
C     Z=25 [+ 5]
      data error(306) /100.0/
      data cdi(1,306) / 9.50000E+01/
      data cdi(2,306) / 2.87600E+00/
      data cdi(3,306) /-2.32100E+00/
      data cdi(4,306) / 5.91000E+00/
      data cdi(5,306) /-1.03060E+01/
      data cdi(6,306) /-1.33300E+00/
      data cdi(7,306) / 7.73800E+00/
C     Z=25 [+ 6]
      data error(307) /300.0/
      data cdi(1,307) / 1.19300E+02/
      data cdi(2,307) / 2.36900E+00/
      data cdi(3,307) / 5.31000E-01/
      data cdi(4,307) /-1.32300E+00/
      data cdi(5,307) /-7.77800E+00/
      data cdi(6,307) / 1.49200E+01/
      data cdi(7,307) /-6.09200E+00/
C     Z=25 [+ 7]
      data error(308) /100.0/
      data cdi(1,308) / 1.96500E+02/
      data cdi(2,308) / 2.85600E+00/
      data cdi(3,308) / 2.09000E-01/
      data cdi(4,308) / 4.88000E-01/
      data cdi(5,308) /-3.79000E-01/
C     Z=25 [+ 8]
      data error(309) /100.0/
      data cdi(1,309) / 2.21800E+02/
      data cdi(2,309) / 2.45700E+00/
      data cdi(3,309) / 2.89000E-01/
      data cdi(4,309) / 4.84000E-01/
      data cdi(5,309) /-4.13000E-01/
C     Z=25 [+ 9]
      data error(310) /100.0/
      data cdi(1,310) / 2.48300E+02/
      data cdi(2,310) / 2.03400E+00/
      data cdi(3,310) / 4.09000E-01/
      data cdi(4,310) / 2.53000E-01/
      data cdi(5,310) /-3.21000E-01/
C     Z=25 [+10]
      data error(311) /100.0/
      data cdi(1,311) / 2.86000E+02/
      data cdi(2,311) / 3.15900E+00/
      data cdi(3,311) /-8.78000E-01/
      data cdi(4,311) /-3.01800E+00/
      data cdi(5,311) / 3.19600E+00/
      data cdi(6,311) /-2.44000E+00/
      data cdi(7,311) /-9.60000E+00/
      data cea(1,311) / 6.06320E+02/
      data cea(2,311) / 5.17600E+02/
      data cea(3,311) / 2.80000E-02/
      data cea(4,311) / 2.82280E+01/
      data cea(5,311) /-4.82170E+01/
      data cea(6,311) / 3.24560E+01/
      data cea(7,311) / 4.38000E-01/
      data cea(8,311) / 6.47000E-01/
C     Z=25 [+11]
      data error(312) /100.0/
      data cdi(1,312) / 3.14400E+02/
      data cdi(2,312) / 1.93600E+00/
      data cdi(3,312) / 5.42000E-01/
      data cdi(4,312) / 5.95000E-01/
      data cdi(5,312) /-5.18000E-01/
C     Z=25 [+12]
      data error(313) /100.0/
      data cdi(1,313) / 3.43600E+02/
      data cdi(2,313) / 7.72000E-01/
      data cdi(3,313) / 4.72000E-01/
      data cdi(4,313) / 6.95000E-01/
      data cdi(5,313) /-8.25000E-01/
      data cdi(6,313) /-1.80000E-02/
      data cea(1,313) / 6.03680E+02/
      data cea(2,313) / 5.77300E+02/
      data cea(3,313) / 2.50100E+00/
      data cea(4,313) / 2.83120E+01/
      data cea(5,313) /-5.54650E+01/
      data cea(6,313) / 2.74410E+01/
      data cea(7,313) / 1.65100E+00/
      data cea(8,313) / 2.23200E+00/
C     Z=25 [+13]
      data error(314) /200.0/
      data cdi(1,314) / 4.04000E+02/
      data cdi(2,314) / 2.13300E+00/
      data cdi(3,314) / 1.91000E-01/
      data cdi(4,314) /-1.50000E-01/
      data cdi(5,314) / 3.20000E-02/
C     Z=25 [+14]
      data error(315) / 70.0/
      data cdi(1,315) / 1.24400E+03/
      data cdi(2,315) / 2.49200E+00/
      data cdi(3,315) / 1.04300E+00/
      data cdi(4,315) /-5.28000E-01/
      data cdi(5,315) /-1.13000E-01/
      data cdi(6,315) /-1.00000E-03/
C     Z=25 [+15]
      data error(316) / 60.0/
      data cdi(1,316) / 1.31700E+03/
      data cdi(2,316) / 1.91200E+00/
      data cdi(3,316) / 1.07900E+00/
      data cdi(4,316) /-3.50000E-02/
      data cdi(5,316) /-4.47000E-01/
C     Z=25 [+16]
      data error(317) / 60.0/
      data cdi(1,317) / 1.43700E+03/
      data cdi(2,317) / 1.63300E+00/
      data cdi(3,317) / 6.45000E-01/
      data cdi(4,317) /-1.01000E-01/
      data cdi(5,317) /-2.39000E-01/
C     Z=25 [+17]
      data error(318) / 60.0/
      data cdi(1,318) / 1.53900E+03/
      data cdi(2,318) / 1.09100E+00/
      data cdi(3,318) / 2.76000E-01/
      data cdi(4,318) / 1.05200E+00/
      data cdi(5,318) /-8.92000E-01/
C     Z=25 [+18]
      data error(319) /100.0/
      data cdi(1,319) / 1.64400E+03/
      data cdi(2,319) / 1.89000E-01/
      data cdi(3,319) / 3.27000E-01/
      data cdi(4,319) / 1.67900E+00/
      data cdi(5,319) /-1.16000E+00/
C     Z=25 [+19]
      data error(320) / 60.0/
      data cdi(1,320) / 1.78800E+03/
      data cdi(2,320) / 4.14000E-01/
      data cdi(3,320) / 9.70000E-02/
      data cdi(4,320) / 1.37000E-01/
      data cdi(5,320) /-1.44000E-01/
C     Z=25 [+20]
      data error(321) /100.0/
      data cdi(1,321) / 1.88000E+03/
      data cdi(2,321) / 2.89000E-01/
      data cdi(3,321) / 1.10000E-01/
      data cdi(4,321) /-4.41000E-01/
      data cdi(5,321) / 2.69000E-01/
C     Z=25 [+21]
      data error(322) / 70.0/
      data cdi(1,322) / 8.14140E+03/
      data cdi(2,322) / 5.25000E-01/
      data cdi(3,322) /-1.88000E-01/
      data cdi(4,322) / 4.89000E-01/
      data cdi(5,322) /-2.92000E-01/
C     Z=25 [+22]
C     Z=25 [+23]
C     Z=25 [+24]
C     Z=26 [+ 0]
      data error(326) / 40.0/
      data cdi(1,326) / 7.90000E+00/
      data cdi(2,326) / 1.14200E+00/
      data cdi(3,326) /-9.20000E-01/
      data cdi(4,326) / 1.78200E+00/
      data cdi(5,326) /-1.69400E+00/
C     Z=26 [+ 1]
      data error(327) / 10.0/
      data cdi(1,327) / 1.62000E+01/
      data cdi(2,327) / 2.12400E+00/
      data cdi(3,327) /-5.30000E-01/
      data cdi(4,327) /-9.61700E+00/
      data cdi(5,327) / 2.91700E+01/
      data cdi(6,327) /-4.08100E+01/
      data cdi(7,327) / 1.92900E+01/
C     Z=26 [+ 2]
      data error(328) / 10.0/
      data cdi(1,328) / 3.06000E+01/
      data cdi(2,328) / 4.42000E+00/
      data cdi(3,328) /-2.96800E+00/
      data cdi(4,328) /-3.71200E+00/
      data cdi(5,328) / 1.21900E+00/
      data cdi(6,328) / 4.30000E-02/
C     Z=26 [+ 3]
      data error(329) / 40.0/
      data cdi(1,329) / 5.48000E+01/
      data cdi(2,329) / 2.25800E+00/
      data cdi(3,329) /-8.43000E-01/
      data cdi(4,329) / 6.22000E-01/
      data cdi(5,329) / 1.00000E-02/
      data cea(1,329) / 1.13400E+02/
      data cea(2,329) / 8.00000E+01/
      data cea(3,329) / 6.11300E+00/
      data cea(4,329) /-1.14800E+00/
      data cea(5,329) / 1.58300E+00/
      data cea(6,329) /-3.40000E-02/
C     Z=26 [+ 4]
      data error(330) / 40.0/
      data cdi(1,330) / 7.50000E+01/
      data cdi(2,330) / 2.65900E+00/
      data cdi(3,330) /-1.06800E+00/
      data cdi(4,330) /-3.96000E-01/
      data cdi(5,330) / 1.03600E+00/
      data cea(1,330) / 1.34250E+02/
      data cea(2,330) / 1.00000E+02/
      data cea(3,330) / 6.50300E+00/
      data cea(4,330) /-2.78000E-01/
      data cea(5,330) /-5.66300E+00/
      data cea(6,330) / 3.86000E+00/
C     Z=26 [+ 5]
      data error(331) / 10.0/
      data cdi(1,331) / 9.90000E+01/
      data cdi(2,331) / 4.39700E+00/
      data cdi(3,331) / 9.80000E-02/
      data cdi(4,331) /-2.68100E+00/
      data cdi(5,331) / 3.04800E+00/
C     Z=26 [+ 6]
      data error(332) / 10.0/
      data cdi(1,332) / 1.25000E+02/
      data cdi(2,332) / 4.61200E+00/
      data cdi(3,332) / 8.19000E-01/
      data cdi(4,332) /-5.12600E+00/
      data cdi(5,332) / 4.02900E+00/
C     Z=26 [+ 7]
      data error(333) /300.0/
      data cdi(1,333) / 1.51100E+02/
      data cdi(2,333) / 3.57000E+00/
      data cdi(3,333) / 5.31000E-01/
      data cdi(4,333) /-1.32300E+00/
      data cdi(5,333) /-7.77800E+00/
      data cdi(6,333) / 1.49210E+01/
      data cdi(7,333) /-6.09200E+00/
C     Z=26 [+ 8]
      data error(334) / 30.0/
      data cdi(1,334) / 2.33600E+02/
      data cdi(2,334) / 2.85600E+00/
      data cdi(3,334) / 2.09000E-01/
      data cdi(4,334) / 4.88000E-01/
      data cdi(5,334) /-3.79000E-01/
C     Z=26 [+ 9]
      data error(335) / 30.0/
      data cdi(1,335) / 2.62100E+02/
      data cdi(2,335) / 2.45700E+00/
      data cdi(3,335) / 2.89000E-01/
      data cdi(4,335) / 4.84000E-01/
      data cdi(5,335) /-4.13000E-01/
C     Z=26 [+10]
      data error(336) / 30.0/
      data cdi(1,336) / 2.90300E+02/
      data cdi(2,336) / 2.03400E+00/
      data cdi(3,336) / 4.10000E-01/
      data cdi(4,336) / 2.53000E-01/
      data cdi(5,336) /-3.21000E-01/
C     Z=26 [+11]
      data error(337) / 20.0/
      data cdi(1,337) / 3.30800E+02/
      data cdi(2,337) / 3.15900E+00/
      data cdi(3,337) /-8.78000E-01/
      data cdi(4,337) /-3.01800E+00/
      data cdi(5,337) / 3.19600E+00/
      data cdi(6,337) /-2.44000E+00/
      data cdi(7,337) /-9.60000E+00/
      data cea(1,337) / 6.99600E+02/
      data cea(2,337) / 6.00000E+02/
      data cea(3,337) / 2.80000E-02/
      data cea(4,337) / 2.82280E+01/
      data cea(5,337) /-4.82170E+01/
      data cea(6,337) / 3.24560E+01/
      data cea(7,337) / 4.38000E-01/
      data cea(8,337) / 6.47000E-01/
C     Z=26 [+12]
      data error(338) / 30.0/
      data cdi(1,338) / 3.61000E+02/
      data cdi(2,338) / 1.93600E+00/
      data cdi(3,338) / 5.42000E-01/
      data cdi(4,338) / 5.95000E-01/
      data cdi(5,338) /-5.18000E-01/
C     Z=26 [+13]
      data error(339) / 25.0/
      data cdi(1,339) / 3.92200E+02/
      data cdi(2,339) / 7.72000E-01/
      data cdi(3,339) / 4.72000E-01/
      data cdi(4,339) / 6.95000E-01/
      data cdi(5,339) /-8.25000E-01/
      data cdi(6,339) /-1.80000E-02/
      data cea(1,339) / 6.89920E+02/
      data cea(2,339) / 6.60000E+02/
      data cea(3,339) / 2.50100E+00/
      data cea(4,339) / 2.83120E+01/
      data cea(5,339) /-5.54650E+01/
      data cea(6,339) / 2.74410E+01/
      data cea(7,339) / 1.65100E+00/
      data cea(8,339) / 2.33200E+00/
C     Z=26 [+14]
      data error(340) / 30.0/
      data cdi(1,340) / 4.57000E+02/
      data cdi(2,340) / 2.13300E+00/
      data cdi(3,340) / 1.91000E-01/
      data cdi(4,340) /-1.50000E-01/
      data cdi(5,340) / 3.20000E-02/
C     Z=26 [+15]
      data error(341) / 25.0/
      data cdi(1,341) / 4.89000E+02/
      data cdi(2,341) / 8.91000E-01/
      data cdi(3,341) /-6.71000E-01/
      data cdi(4,341) / 6.00000E-01/
      data cdi(5,341) /-7.60000E-02/
      data cea(1,341) / 7.28610E+02/
      data cea(2,341) / 7.20000E+02/
      data cea(3,341) / 3.78000E+00/
      data cea(4,341) / 1.70260E+01/
      data cea(5,341) /-5.86440E+01/
      data cea(6,341) / 4.45470E+01/
      data cea(7,341) / 5.70000E-02/
      data cea(8,341) / 5.70000E-02/
C     Z=26 [+16]
      data error(342) / 40.0/
      data cdi(1,342) / 1.26600E+03/
      data cdi(2,342) / 3.44100E+00/
      data cdi(3,342) / 1.06600E+00/
      data cdi(4,342) /-2.77400E+00/
      data cdi(5,342) / 1.48900E+00/
C     Z=26 [+17]
      data error(343) / 30.0/
      data cdi(1,343) / 1.36600E+03/
      data cdi(2,343) / 2.49200E+00/
      data cdi(3,343) / 1.04300E+00/
      data cdi(4,343) /-5.28000E-01/
      data cdi(5,343) /-1.13000E-01/
C     Z=26 [+18]
      data error(344) / 40.0/
      data cdi(1,344) / 1.46300E+03/
      data cdi(2,344) / 1.91200E+00/
      data cdi(3,344) / 1.07900E+00/
      data cdi(4,344) /-3.50000E-02/
      data cdi(5,344) /-4.47000E-01/
C     Z=26 [+19]
      data error(345) / 30.0/
      data cdi(1,345) / 1.58300E+03/
      data cdi(2,345) / 1.63300E+00/
      data cdi(3,345) / 6.45000E-01/
      data cdi(4,345) /-1.01000E-01/
      data cdi(5,345) /-2.39000E-01/
C     Z=26 [+20]
      data error(346) / 30.0/
      data cdi(1,346) / 1.67800E+03/
      data cdi(2,346) / 1.09100E+00/
      data cdi(3,346) / 2.76000E-01/
      data cdi(4,346) / 1.05200E+00/
      data cdi(5,346) /-8.92000E-01/
C     Z=26 [+21]
      data error(347) / 30.0/
      data cdi(1,347) / 1.78900E+03/
      data cdi(2,347) / 1.89000E-01/
      data cdi(3,347) / 3.27000E-01/
      data cdi(4,347) / 1.67900E+00/
      data cdi(5,347) /-1.16000E+00/
C     Z=26 [+22]
      data error(348) / 30.0/
      data cdi(1,348) / 1.95000E+03/
      data cdi(2,348) / 4.14000E-01/
      data cdi(3,348) / 9.70000E-02/
      data cdi(4,348) / 1.37000E-01/
      data cdi(5,348) /-1.44000E-01/
C     Z=26 [+23]
      data error(349) / 40.0/
      data cdi(1,349) / 2.04600E+03/
      data cdi(2,349) / 2.89000E-01/
      data cdi(3,349) / 1.10000E-01/
      data cdi(4,349) /-4.41000E-01/
      data cdi(5,349) / 2.69000E-01/
C     Z=26 [+24]
      data error(350) / 40.0/
      data cdi(1,350) / 8.82900E+03/
      data cdi(2,350) / 5.25000E-01/
      data cdi(3,350) /-1.88000E-01/
      data cdi(4,350) / 4.89000E-01/
      data cdi(5,350) /-2.92000E-01/
C     Z=26 [+25]
      data error(351) / 30.0/
      data cdi(1,351) / 9.27700E+03/
      data cdi(2,351) / 2.60000E-01/
      data cdi(3,351) /-2.20000E-02/
      data cdi(4,351) / 5.10000E-02/
      data cdi(5,351) /-2.30000E-02/
C     Z=27 [+ 0]
      data error(352) /200.0/
      data cdi(1,352) / 7.90000E+00/
      data cdi(2,352) / 2.67900E+00/
      data cdi(3,352) /-2.86500E+00/
      data cdi(4,352) / 6.28900E+00/
      data cdi(5,352) /-1.59470E+01/
      data cdi(6,352) / 9.54700E+00/
C     Z=27 [+ 1]
      data error(353) /200.0/
      data cdi(1,353) / 1.71000E+01/
      data cdi(2,353) / 3.45300E+00/
      data cdi(3,353) /-2.00800E+00/
      data cdi(4,353) /-1.62100E+00/
      data cdi(5,353) / 8.89000E-01/
C     Z=27 [+ 2]
      data error(354) /200.0/
      data cdi(1,354) / 3.35000E+01/
      data cdi(2,354) / 3.49700E+00/
      data cdi(3,354) /-1.69800E+00/
      data cdi(4,354) /-1.72300E+00/
      data cdi(5,354) / 2.32900E+00/
      data cea(1,354) / 8.18400E+01/
      data cea(2,354) / 4.20000E+01/
      data cea(3,354) / 6.04400E+00/
      data cea(4,354) / 2.13100E+00/
      data cea(5,354) /-1.83700E+01/
      data cea(6,354) / 1.54590E+01/
C     Z=27 [+ 3]
      data error(355) /200.0/
      data cdi(1,355) / 5.13000E+01/
      data cdi(2,355) / 4.42000E+00/
      data cdi(3,355) /-4.81000E+00/
      data cdi(4,355) / 5.69200E+00/
      data cdi(5,355) /-1.29770E+01/
      data cdi(6,355) / 6.69800E+00/
C     Z=27 [+ 4]
      data error(356) /100.0/
      data cdi(1,356) / 7.95000E+01/
      data cdi(2,356) / 2.25800E+00/
      data cdi(3,356) /-8.43000E-01/
      data cdi(4,356) / 6.22000E-01/
      data cdi(5,356) / 1.00000E-02/
      data cea(1,356) / 1.58000E+02/
      data cea(2,356) / 1.19300E+02/
      data cea(3,356) / 6.61300E+00/
      data cea(4,356) /-1.14800E+00/
      data cea(5,356) / 1.58300E+00/
      data cea(6,356) /-3.40000E-02/
C     Z=27 [+ 5]
      data error(357) /100.0/
      data cdi(1,357) / 1.02000E+02/
      data cdi(2,357) / 2.65900E+00/
      data cdi(3,357) /-1.06800E+00/
      data cdi(4,357) /-3.96000E-01/
      data cdi(5,357) / 1.03600E+00/
      data cea(1,357) / 1.83600E+02/
      data cea(2,357) / 1.27500E+02/
      data cea(3,357) / 6.50300E+00/
      data cea(4,357) /-2.78000E-01/
      data cea(5,357) /-5.66300E+00/
      data cea(6,357) / 3.86000E+00/
C     Z=27 [+ 6]
      data error(358) /200.0/
      data cdi(1,358) / 1.29000E+02/
      data cdi(2,358) / 2.88800E+00/
      data cdi(3,358) / 7.53000E-01/
      data cdi(4,358) /-2.12400E+01/
      data cdi(5,358) / 6.62900E+01/
      data cdi(6,358) /-7.52300E+01/
      data cdi(7,358) / 2.90800E+01/
C     Z=27 [+ 7]
      data error(359) /200.0/
      data cdi(1,359) / 1.57000E+02/
      data cdi(2,359) / 2.87600E+00/
      data cdi(3,359) /-2.32100E+00/
      data cdi(4,359) / 5.91000E+00/
      data cdi(5,359) /-1.03060E+01/
      data cdi(6,359) /-1.33300E+00/
      data cdi(7,359) / 7.73800E+00/
C     Z=27 [+ 8]
      data error(360) /300.0/
      data cdi(1,360) / 1.86100E+02/
      data cdi(2,360) / 2.36900E+00/
      data cdi(3,360) / 5.31000E-01/
      data cdi(4,360) /-1.32300E+00/
      data cdi(5,360) /-7.77800E+00/
      data cdi(6,360) / 1.49200E+01/
      data cdi(7,360) /-6.09200E+00/
C     Z=27 [+ 9]
      data error(361) /100.0/
      data cdi(1,361) / 2.76000E+02/
      data cdi(2,361) / 2.85600E+00/
      data cdi(3,361) / 2.09000E-01/
      data cdi(4,361) / 4.88000E-01/
      data cdi(5,361) /-3.79000E-01/
C     Z=27 [+10]
      data error(362) /100.0/
      data cdi(1,362) / 3.05000E+02/
      data cdi(2,362) / 2.45700E+00/
      data cdi(3,362) / 2.89000E-01/
      data cdi(4,362) / 4.84000E-01/
      data cdi(5,362) /-4.13000E-01/
C     Z=27 [+11]
      data error(363) /100.0/
      data cdi(1,363) / 3.36000E+02/
      data cdi(2,363) / 2.03400E+00/
      data cdi(3,363) / 4.09000E-01/
      data cdi(4,363) / 2.53000E-01/
      data cdi(5,363) /-3.21000E-01/
C     Z=27 [+12]
      data error(364) /100.0/
      data cdi(1,364) / 3.79000E+02/
      data cdi(2,364) / 3.15900E+00/
      data cdi(3,364) /-8.78000E-01/
      data cdi(4,364) /-3.01800E+00/
      data cdi(5,364) / 3.19600E+00/
      data cdi(6,364) /-2.44000E+00/
      data cdi(7,364) /-9.60000E+00/
      data cea(1,364) / 8.03480E+02/
      data cea(2,364) / 6.85900E+02/
      data cea(3,364) / 2.80000E-02/
      data cea(4,364) / 2.82280E+01/
      data cea(5,364) /-4.82170E+01/
      data cea(6,364) / 3.24560E+01/
      data cea(7,364) / 4.38000E-01/
      data cea(8,364) / 6.47000E-01/
C     Z=27 [+13]
      data error(365) /100.0/
      data cdi(1,365) / 4.11000E+02/
      data cdi(2,365) / 1.93600E+00/
      data cdi(3,365) / 5.42000E-01/
      data cdi(4,365) / 5.95000E-01/
      data cdi(5,365) /-5.18000E-01/
C     Z=27 [+14]
      data error(366) /100.0/
      data cdi(1,366) / 4.44000E+02/
      data cdi(2,366) / 7.72000E-01/
      data cdi(3,366) / 4.72000E-01/
      data cdi(4,366) / 6.95000E-01/
      data cdi(5,366) /-8.25000E-01/
      data cdi(6,366) /-1.80000E-02/
      data cea(1,366) / 7.81440E+02/
      data cea(2,366) / 7.45900E+02/
      data cea(3,366) / 2.50100E+00/
      data cea(4,366) / 2.83120E+01/
      data cea(5,366) /-5.54650E+01/
      data cea(6,366) / 2.74410E+01/
      data cea(7,366) / 1.65100E+00/
      data cea(8,366) / 2.23200E+00/
C     Z=27 [+15]
      data error(367) /200.0/
      data cdi(1,367) / 5.12000E+02/
      data cdi(2,367) / 2.13300E+00/
      data cdi(3,367) / 1.91000E-01/
      data cdi(4,367) /-1.50000E-01/
      data cdi(5,367) / 3.20000E-02/
C     Z=27 [+16]
      data error(368) /200.0/
      data cdi(1,368) / 5.46800E+02/
      data cdi(2,368) / 1.50700E+00/
      data cdi(3,368) /-1.05000E-01/
      data cdi(4,368) / 1.47000E-01/
      data cdi(5,368) /-7.30000E-02/
      data cdi(6,368) /-1.00000E-03/
C     Z=27 [+17]
      data error(369) / 50.0/
      data cdi(1,369) / 1.40300E+03/
      data cdi(2,369) / 3.44100E+00/
      data cdi(3,369) / 1.06600E+00/
      data cdi(4,369) /-2.77400E+00/
      data cdi(5,369) / 1.48900E+00/
C     Z=27 [+18]
C     Z=27 [+19]
C     Z=27 [+20]
C     Z=27 [+21]
C     Z=27 [+22]
C     Z=27 [+23]
C     Z=27 [+24]
C     Z=27 [+25]
C     Z=27 [+26]
C     Z=28 [+ 0]
      data error(379) /200.0/
      data cdi(1,379) / 7.60000E+00/
      data cdi(2,379) / 8.32000E-01/
      data cdi(3,379) /-7.66000E-01/
      data cdi(4,379) / 1.25800E+00/
      data cdi(5,379) /-6.79000E-01/
C     Z=28 [+ 1]
      data error(380) /200.0/
      data cdi(1,380) / 1.82000E+01/
      data cdi(2,380) / 2.67900E+00/
      data cdi(3,380) /-2.86500E+00/
      data cdi(4,380) / 6.28900E+00/
      data cdi(5,380) /-1.59470E+01/
      data cdi(6,380) / 9.54700E+00/
C     Z=28 [+ 2]
      data error(381) / 40.0/
      data cdi(1,381) / 3.62000E+01/
      data cdi(2,381) / 3.45300E+00/
      data cdi(3,381) /-2.00800E+00/
      data cdi(4,381) /-1.62100E+00/
      data cdi(5,381) / 8.89000E-01/
C     Z=28 [+ 3]
      data error(382) / 10.0/
      data cdi(1,382) / 5.49000E+01/
      data cdi(2,382) / 3.49700E+00/
      data cdi(3,382) /-1.69800E+00/
      data cdi(4,382) /-1.72300E+00/
      data cdi(5,382) / 2.32900E+00/
      data cea(1,382) / 1.22580E+02/
      data cea(2,382) / 7.00000E+01/
      data cea(3,382) / 6.04400E+00/
      data cea(4,382) / 2.13100E+00/
      data cea(5,382) /-1.83700E+01/
      data cea(6,382) / 1.54590E+01/
C     Z=28 [+ 4]
      data error(383) /200.0/
      data cdi(1,383) / 7.61000E+01/
      data cdi(2,383) / 4.42000E+00/
      data cdi(3,383) /-4.81000E+00/
      data cdi(4,383) / 5.69200E+00/
      data cdi(5,383) /-1.29770E+01/
      data cdi(6,383) / 6.69800E+00/
C     Z=28 [+ 5]
      data error(384) /100.0/
      data cdi(1,384) / 1.08000E+02/
      data cdi(2,384) / 2.25800E+00/
      data cdi(3,384) /-8.43000E-01/
      data cdi(4,384) / 6.22000E-01/
      data cdi(5,384) / 1.00000E-02/
      data cea(1,384) / 2.26800E+02/
      data cea(2,384) / 1.60600E+02/
      data cea(3,384) / 6.11300E+00/
      data cea(4,384) /-1.14800E+00/
      data cea(5,384) / 1.58300E+00/
      data cea(6,384) /-3.40000E-02/
C     Z=28 [+ 6]
      data error(385) /100.0/
      data cdi(1,385) / 1.33000E+02/
      data cdi(2,385) / 2.65900E+00/
      data cdi(3,385) /-1.06800E+00/
      data cdi(4,385) /-3.96000E-01/
      data cdi(5,385) / 1.03600E+00/
      data cea(1,385) / 2.39400E+02/
      data cea(2,385) / 1.66000E+02/
      data cea(3,385) / 6.50300E+00/
      data cea(4,385) /-2.78000E-01/
      data cea(5,385) /-5.66300E+00/
      data cea(6,385) / 3.86000E+00/
C     Z=28 [+ 7]
      data error(386) /200.0/
      data cdi(1,386) / 1.62000E+02/
      data cdi(2,386) / 2.88800E+00/
      data cdi(3,386) / 7.53000E-01/
      data cdi(4,386) /-2.12400E+01/
      data cdi(5,386) / 6.62900E+01/
      data cdi(6,386) /-7.52300E+01/
      data cdi(7,386) / 2.90800E+01/
C     Z=28 [+ 8]
      data error(387) /200.0/
      data cdi(1,387) / 1.93000E+02/
      data cdi(2,387) / 2.87600E+00/
      data cdi(3,387) /-2.32100E+00/
      data cdi(4,387) / 5.91000E+00/
      data cdi(5,387) /-1.03060E+01/
      data cdi(6,387) /-1.33300E+00/
      data cdi(7,387) / 7.73800E+00/
C     Z=28 [+ 9]
      data error(388) /300.0/
      data cdi(1,388) / 2.24600E+02/
      data cdi(2,388) / 2.36900E+00/
      data cdi(3,388) / 5.31000E-01/
      data cdi(4,388) /-1.32300E+00/
      data cdi(5,388) /-7.77800E+00/
      data cdi(6,388) / 1.49200E+01/
      data cdi(7,388) /-6.09200E+00/
C     Z=28 [+10]
      data error(389) /100.0/
      data cdi(1,389) / 3.21000E+02/
      data cdi(2,389) / 2.85600E+00/
      data cdi(3,389) / 2.09000E-01/
      data cdi(4,389) / 4.88000E-01/
      data cdi(5,389) /-3.79000E-01/
C     Z=28 [+11]
      data error(390) /100.0/
      data cdi(1,390) / 3.52000E+02/
      data cdi(2,390) / 2.45700E+00/
      data cdi(3,390) / 2.89000E-01/
      data cdi(4,390) / 4.84000E-01/
      data cdi(5,390) /-4.13000E-01/
C     Z=28 [+12]
      data error(391) /100.0/
      data cdi(1,391) / 3.84000E+02/
      data cdi(2,391) / 2.03400E+00/
      data cdi(3,391) / 4.10000E-01/
      data cdi(4,391) / 2.53000E-01/
      data cdi(5,391) /-3.21000E-01/
C     Z=28 [+13]
      data error(392) /100.0/
      data cdi(1,392) / 4.30000E+02/
      data cdi(2,392) / 3.15900E+00/
      data cdi(3,392) /-8.78000E-01/
      data cdi(4,392) /-3.01800E+00/
      data cdi(5,392) / 3.19600E+00/
      data cdi(6,392) /-2.44000E+00/
      data cdi(7,392) /-9.60000E+00/
      data cea(1,392) / 9.11600E+02/
      data cea(2,392) / 7.78300E+02/
      data cea(3,392) / 2.80000E-02/
      data cea(4,392) / 2.82280E+01/
      data cea(5,392) /-4.82170E+01/
      data cea(6,392) / 3.24560E+01/
      data cea(7,392) / 4.38000E-01/
      data cea(8,392) / 6.47000E-01/
C     Z=28 [+14]
      data error(393) /100.0/
      data cdi(1,393) / 4.64000E+02/
      data cdi(2,393) / 1.93600E+00/
      data cdi(3,393) / 5.42000E-01/
      data cdi(4,393) / 5.95000E-01/
      data cdi(5,393) /-5.18000E-01/
C     Z=28 [+15]
      data error(394) /100.0/
      data cdi(1,394) / 4.99000E+02/
      data cdi(2,394) / 7.72000E-01/
      data cdi(3,394) / 4.72000E-01/
      data cdi(4,394) / 6.95000E-01/
      data cdi(5,394) /-8.25000E-01/
      data cdi(6,394) /-1.80000E-02/
      data cea(1,394) / 8.78240E+02/
      data cea(2,394) / 8.38300E+02/
      data cea(3,394) / 2.51000E+00/
      data cea(4,394) / 2.83120E+01/
      data cea(5,394) /-5.54650E+01/
      data cea(6,394) / 2.74410E+01/
      data cea(7,394) / 1.65100E+00/
      data cea(8,394) / 2.23200E+00/
C     Z=28 [+16]
      data error(395) /200.0/
      data cdi(1,395) / 5.71000E+02/
      data cdi(2,395) / 2.13300E+00/
      data cdi(3,395) / 1.91000E-01/
      data cdi(4,395) /-1.50000E-01/
      data cdi(5,395) / 3.20000E-02/
C     Z=28 [+17]
      data error(396) /200.0/
      data cdi(1,396) / 6.07200E+02/
      data cdi(2,396) / 1.50700E+00/
      data cdi(3,396) / 1.05000E-01/
      data cdi(4,396) / 1.47000E-01/
      data cdi(5,396) /-7.30000E-02/
      data cdi(6,396) /-1.00000E-03/
C     Z=28 [+18]
      data error(397) / 50.0/
      data cdi(1,397) / 1.54700E+03/
      data cdi(2,397) / 3.44100E+00/
      data cdi(3,397) / 1.06600E+00/
      data cdi(4,397) /-2.77400E+00/
      data cdi(5,397) / 1.48900E+00/
C     Z=28 [+19]
      data error(398) / 70.0/
      data cdi(1,398) / 1.63900E+03/
      data cdi(2,398) / 2.49200E+00/
      data cdi(3,398) / 1.04300E+00/
      data cdi(4,398) /-5.28000E-01/
      data cdi(5,398) /-1.13000E-01/
C     Z=28 [+20]
      data error(399) / 60.0/
      data cdi(1,399) / 1.74700E+03/
      data cdi(2,399) / 1.91200E+00/
      data cdi(3,399) / 1.07900E+00/
      data cdi(4,399) /-3.50000E-02/
      data cdi(5,399) /-4.47000E-01/
C     Z=28 [+21]
      data error(400) / 60.0/
      data cdi(1,400) / 1.88200E+03/
      data cdi(2,400) / 1.63300E+00/
      data cdi(3,400) / 6.45000E-01/
      data cdi(4,400) /-1.01000E-01/
      data cdi(5,400) /-2.39000E-01/
C     Z=28 [+22]
      data error(401) / 60.0/
      data cdi(1,401) / 2.00300E+03/
      data cdi(2,401) / 1.09100E+00/
      data cdi(3,401) / 2.76000E-01/
      data cdi(4,401) / 1.05200E+00/
      data cdi(5,401) /-8.92000E-01/
C     Z=28 [+23]
      data error(402) /100.0/
      data cdi(1,402) / 2.12300E+03/
      data cdi(2,402) / 1.89000E-01/
      data cdi(3,402) / 3.27000E-01/
      data cdi(4,402) / 1.67900E+00/
      data cdi(5,402) /-1.16000E+00/
C     Z=28 [+24]
      data error(403) / 60.0/
      data cdi(1,403) / 2.29900E+03/
      data cdi(2,403) / 4.14000E-01/
      data cdi(3,403) / 9.70000E-02/
      data cdi(4,403) / 1.37000E-01/
      data cdi(5,403) /-1.44000E-01/
C     Z=28 [+25]
      data error(404) / 10.0/
      data cdi(1,404) / 2.39800E+03/
      data cdi(2,404) / 2.89000E-01/
      data cdi(3,404) / 1.10000E-01/
      data cdi(4,404) /-4.41000E-01/
      data cdi(5,404) / 2.69000E-01/
C     Z=28 [+26]
C     Z=28 [+27]
      end
