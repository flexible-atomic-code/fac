      double precision function fm1m(x)
      implicit real*8 (a-h, o-z)
      dimension a1(12), b1(12), a2(12), b2(12)
      data an, m1, k1, m2, k2/-0.5d0, 2, 3, 2, 2/
      data (a1(i),i=1,3)/2.31456d+1,1.37820d+1,1.00000d+0/
      data (b1(i),i=1,4)/1.30586d+1,1.70048d+1,5.07527d+0,2.36620d-1/
      data (a2(i),i=1,3)/1.53602d-2,1.46815d-1,1.00000d0/
      data (b2(i),i=1,3)/7.68015e-3,7.63700d-2,5.70485d-1/

      if (x .lt. 2d0) then
         xx=dexp(x)
         rn=xx+a1(m1)
         do i=m1-1,1,-1
            rn=rn*xx+a1(i)
         enddo
         den=b1(k1+1)
         do i=k1,1,-1
            den = den*xx+b1(i)
         enddo
         fm1m = xx*rn/den
      else
         xx=1d0/x**2
         rn=xx+a2(m2)
         do i=m2-1,1,-1
            rn = rn*xx+a2(i)
         enddo
         den=b2(k2+1)
         do i=k2,1,-1
            den=den*xx+b2(i)
         enddo
         fm1m = (x**(an+1d0))*rn/den
      endif
      return
      end
      
