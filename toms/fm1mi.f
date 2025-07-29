      double precision function fm1mi(f)
      implicit real*8 (a-h, o-z)
      dimension a1(9),b1(9),a2(9),b2(9)
      data an,m1,k1,m2,k2/-0.5d0,5,6,6,6/
      data a1(1)/-1.570044577033d4/
      data a1(2)/1.001958278442d4/
      data a1(3)/-2.805343454951d3/
      data a1(4)/4.121170498099d2/
      data a1(5)/-3.174780572961d1/
      data a1(6)/1.0d0/
      data b1(1)/-2.782831558471d4/
      data b1(2)/2.886114034012d4/
      data b1(3)/-1.274243093149d4/
      data b1(4)/3.063252215963d3/
      data b1(5)/-4.225615045074d2/
      data b1(6)/3.168918168284d1/
      data b1(7)/-1.008561571363d0/
      data a2(1)/2.206779160034d-8/
      data a2(2)/-1.437701234283d-6/
      data a2(3)/6.103116850636d-5/
      data a2(4)/-1.169411057416d-3/
      data a2(5)/1.814141021608d-2/
      data a2(6)/-9.588603457639d-2/
      data a2(7)/1.0d0/
      data b2(1)/8.827116613576d-8/
      data b2(2)/-5.750804196059d-6/
      data b2(3)/2.429627688357d-4/
      data b2(4)/-4.601959491394d-3/
      data b2(5)/6.932122275959d-2/
      data b2(6)/-3.217372489776d-1/
      data b2(7)/3.124344749296d0/

      if (f .lt. 4d0) then
         rn = f+a1(m1)
         do i=m1-1,1,-1
            rn = rn*f+a1(i)
         enddo
         den = b1(k1+1)
         do i = k1,1,-1
            den = den*f+b1(i)
         enddo
         fm1mi = dlog(f*rn/den)
      else
         ff = 1d0/f**(1d0/(1d0+an))
         rn = ff+a2(m2)
         do i=m2-1,1,-1
            rn = rn*ff+a2(i)
         enddo
         den = b2(k2+1)
         do i=k2,1,-1
            den = den*ff+b2(i)
         enddo
         fm1mi = rn/(den*ff)
      endif
      return
      end
      
