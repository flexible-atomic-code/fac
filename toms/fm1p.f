      double precision function fm1p(x)
      implicit real*8 (a-h, o-z)
      dimension a1(12), b1(12), a2(12), b2(12)
      data an, m1, k1, m2, k2/0.5d0, 2, 3, 2, 2/
      data (a1(i),i=1,3)/2.18168d+1,1.31693d+1,1.00000d+0/
      data (b1(i),i=1,4)/2.46180d+1,2.35546d+1,4.76290d+1,1.34481d-1/
      data (a2(i),i=1,3)/4.73011d-2,5.48433d-1,1.00000d0/
      data (b2(i),i=1,3)/7.09478e-2,7.37041d-1,3.82065d-1/

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
         fm1p = xx*rn/den
      else
         xx=1d0/x**2
         rn=xx+a2(m2)
         do i=m2-1,1,-1
            rn=rn*xx+a2(i)
         enddo
         den=b2(k2+1)
         do i=k2,1,-1
            den=den*xx+b2(i)
         enddo
         fm1p = (x**(an+1d0))*rn/den
      endif
      return
      end
      
