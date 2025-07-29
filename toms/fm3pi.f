      double precision function fm3pi(f)
      implicit real*8 (a-h, o-z)
      dimension a1(9),b1(9),a2(9),b2(9)
      data an,m1,k1,m2,k2/1.5d0,3,4,6,5/
      data a1(1)/1.715627994191d2/
      data a1(2)/1.125926232897d2/
      data a1(3)/2.056296753055d1/
      data a1(4)/1.0d0/
      data b1(1)/2.280653583157d2/
      data b1(2)/1.193456203021d2/
      data b1(3)/1.167743113540d1/
      data b1(4)/-3.226808804038d-1/
      data b1(5)/3.519268762788d-3/
      data a2(1)/-6.321828169799d-3/
      data a2(2)/-2.183147266896d-2/
      data a2(3)/-1.057562799320d-1/
      data a2(4)/-4.657944387545d-1/
      data a2(5)/-5.951932864088d-1/
      data a2(6)/3.684471177100d-1/
      data a2(7)/1.0d0/
      data b2(1)/-4.381942605018d-3/
      data b2(2)/-1.513236504100d-2/
      data b2(3)/-7.850001283886d-2/
      data b2(4)/-3.407561772613d-1/
      data b2(5)/-5.074812565486d-1/
      data b2(6)/-1.387107009074d-1/

      if (f .lt. 4d0) then
         rn = f+a1(m1)
         do i=m1-1,1,-1
            rn = rn*f+a1(i)
         enddo
         den = b1(k1+1)
         do i = k1,1,-1
            den = den*f+b1(i)
         enddo
         fm3pi = dlog(f*rn/den)
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
         fm3pi = rn/(den*ff)
      endif
      return
      end
      
