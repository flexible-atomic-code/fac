      double precision function fm1mi(f)
      implicit real*8 (a-h, o-z)
      dimension a1(9),b1(9),a2(9),b2(9)
      data an,m1,k1,m2,k2/-0.5d0,3,3,2,2/
      data (a1(i),i=1,4)/7.8516685d2,-1.4034065d2,1.3257418d1,1d0/
      data (b1(i),i=1,4)/1.3917278d3,-8.0463066d2,
     +     .5854806d2,-1.0640712d1/
      data (a2(i),i=1,3)/8.9742174d-3,-1.0604768d-1,1.0d0/
      data (b2(i),i=1,3)/3.5898124d-2,-4.2520975d-1,3.6612154d0/

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
      
