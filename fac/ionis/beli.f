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
c     dir, direct ionization cross section (10^-20 cm2).
c     ea, excitation autoionization contribution (10^-20 cm2).
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

      dir = 0.0
      if (ene .gt. ei0) then
         x = ei0/ene
         x2 = 1.0/x
         xs = a0*dlog(x2)
         power1 = 1.0 - x
         power = power1
         do i = 3, 7
            if (cdi(i, k) .gt. 0) then
               xs = xs + cdi(i, k)*power
               power = power*power1
            endif
         enddo
         dir = 1E7*xs/(ene*ei0)
      endif

      ea = 0.0
      if (ea0 .gt. 0 .and. ene .gt. ea0) then
         x = ei1/ene
         x2 = 1.0/x
         xs = a1*dlog(x2)
         power1 = 1.0 - x
         power = power1
         do i = 4, 8
            if (cea(i, k) .gt. 0) then
               xs = xs + cea(i, k)*power
               power = power * power1
            endif
         enddo
         ea = 1E7*xs/(ene*ei1) - dir
      endif

      return
      end
