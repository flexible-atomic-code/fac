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

      if (ei0 .eq. 0) then
         ea = -1.0
         dir = -1.0
         return
      endif

      dir = 0.0
      if (ei0 .gt. 0 .and. ene .gt. ei0) then
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
      if (ea0 .gt. 0 .and. ene .gt. ea0) then
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
      double precision ei0, ei1, a, b, c, x0, ene
      double precision mc
c     mc = 1.12837967*c*sqrt(2.0/mc^2)
      parameter(mc = 0.0066923847825)

      k = 1 + (z*(z-1)/2) + (z-n)
      ei0 = cdi(1, k)
      if (ei0 .eq. 0) then
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
