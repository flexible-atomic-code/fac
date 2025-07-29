      double precision function fm1pi(f)
      implicit real*8 (a-h, o-z)
      dimension a1(9),b1(9),a2(9),b2(9)
      data an,m1,k1,m2,k2/0.5d0,4,3,6,5/
      data a1(1)/1.999266880833d4/
      data a1(2)/5.702479099336d3/
      data a1(3)/6.610132843877d2/
      data a1(4)/3.818838129486d1/
      data a1(5)/1.0d0/
      data b1(1)/1.771804140488d4/
      data b1(2)/-2.014785161019d3/
      data b1(3)/9.130355392717d1/
      data b1(4)/-1.670718177489d0/
      data a2(1)/-1.277060388085d-2/
      data a2(2)/7.187946804945d-2/
      data a2(3)/-4.262314235106d-1/
      data a2(4)/4.997559426872d-1/
      data a2(5)/-1.285579118012d0/
      data a2(6)/-3.930805454272d-1/
      data a2(7)/1.0d0/
      data b2(1)/-9.745794806288d-3/
      data b2(2)/5.485432756838d-2/
      data b2(3)/-3.299466243260d-1/
      data b2(4)/4.077841975923d-1/
      data b2(5)/-1.145531476975d0/
      data b2(6)/-6.067091689181d-2/

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
      
