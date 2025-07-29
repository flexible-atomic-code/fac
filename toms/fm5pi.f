      double precision function fm5pi(f)
      implicit real*8 (a-h, o-z)
      dimension a1(9),b1(9),a2(9),b2(9)
      data an,m1,k1,m2,k2/2.5d0,2,3,6,6/
      data a1(1)/2.138969250409d2/
      data a1(2)/3.539903493971d1/
      data a1(3)/1.0d0/
      data b1(1)/7.108545512710d2/
      data b1(2)/9.873746988121d1/
      data b1(3)/1.067755522895d0/
      data b1(4)/-1.182798726503d-2/
      data a2(1)/-3.312041011227d-2/
      data a2(2)/1.315763372315d-1/
      data a2(3)/-4.820942898296d-1/
      data a2(4)/5.099038074944d-1/
      data a2(5)/5.495613498630d-1/
      data a2(6)/-1.498867562255d0/
      data a2(7)/1.0d0/
      data b2(1)/-2.315515517515d-2/
      data b2(2)/9.198776585252d-2/
      data b2(3)/-3.835879295548d-1/
      data b2(4)/5.415026856351d-1/
      data b2(5)/-3.847241692193d-1/
      data b2(6)/3.739781456585d-2/
      data b2(7)/-3.008504449098d-2/

      if (f .lt. 4d0) then
         rn = f+a1(m1)
         do i=m1-1,1,-1
            rn = rn*f+a1(i)
         enddo
         den = b1(k1+1)
         do i = k1,1,-1
            den = den*f+b1(i)
         enddo
         fm5pi = dlog(f*rn/den)
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
         fm5pi = rn/(den*ff)
      endif
      return
      end
      
