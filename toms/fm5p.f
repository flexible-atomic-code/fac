      double precision function fm5p(x)
      implicit real*8 (a-h, o-z)
      dimension a1(12), b1(12), a2(12), b2(12)
      data an, m1, k1, m2, k2/2.5d0, 6, 7, 10, 9/
      data a1(1)/6.61606300631656d4/
      data a1(2)/1.20132462801652d5/
      data a1(3)/7.67255995316812d4/
      data a1(4)/2.10427138842443d4/
      data a1(5)/2.44325236813275d3/
      data a1(6)/1.02589947781696d2/
      data a1(7)/1.0d0/
      data b1(1)/1.99078071053871d4/
      data b1(2)/3.79076097261066d4/
      data b1(3)/2.60117136841197d4/
      data b1(4)/7.97584657659364d3/
      data b1(5)/1.10886130159658d3/
      data b1(6)/6.35483623268093d1/
      data b1(7)/1.16951072617142d0/
      data b1(8)/3.31482978240026d-3/
      data a2(1)/8.42667076131315d-12/
      data a2(2)/2.31618876821567d-9/
      data a2(3)/3.54323824923987d-7/
      data a2(4)/2.77981736000034d-5/
      data a2(5)/1.14008027400645d-3/
      data a2(6)/2.32779790773633d-2/
      data a2(7)/2.39564845938301d-1/
      data a2(8)/1.24415366126179d0/
      data a2(9)/3.18831203950106d0/
      data a2(10)/3.42040216997894d0/
      data a2(11)/1.0d0/
      data b2(1)/2.94933476646033d-11/
      data b2(2)/7.68215783076936d-9/
      data b2(3)/1.12818616415947d-6/
      data b2(4)/8.09451165406274d-5/
      data b2(5)/2.81111224925648d-3/
      data b2(6)/3.99937801931919d-2/
      data b2(7)/2.27132567866839d-1/
      data b2(8)/5.31886045222680d-1/
      data b2(9)/3.70866321410385d-1/
      data b2(10)/2.27326643192516d-2/

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
         fm5p = xx*rn/den
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
         fm5p = (x**(an+1d0))*rn/den
      endif
      return
      end
      
