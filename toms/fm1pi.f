      double precision function fm1pi(f)
      implicit real*8 (a-h, o-z)
      dimension a1(9),b1(9),a2(9),b2(9)
      data an,m1,k1,m2,k2/0.5d0,2,2,2,2/
      data (a1(i),i=1,3)/4.4593646d1,1.1288764d1,1.0d0/
      data (b1(i),i=1,3)/3.9519346d1,-5.7517464d0,2.6594291d-1/
      data (a2(i),i=1,3)/3.4873722d1,-2.6922515d1,1.0d0/
      data (b2(i),i=1,3)/2.6612832d1,-2.0452930d1,1.1808945d1/

      if (f .lt. 4d0) then
         rn = f+a1(m1)
         do i=m1-1,1,-1
            rn = rn*f+a1(i)
         enddo
         den = b1(k1+1)
         do i = k1,1,-1
            den = den*f+b1(i)
         enddo
         fm1pi = dlog(f*rn/den)
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
         fm1pi = rn/(den*ff)
      endif
      return
      end
      
