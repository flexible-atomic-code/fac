      double precision function fm3p(x)
      implicit real*8 (a-h, o-z)
      dimension a1(12), b1(12), a2(12), b2(12)
      data an, m1, k1, m2, k2/1.5d0, 6, 7, 9, 10/
      data a1(1)/4.32326386604283d4/
      data a1(2)/8.55472308218786d4/
      data a1(3)/5.95275291210962d4/
      data a1(4)/1.77294861572005d4/
      data a1(5)/2.21876607796460d3/
      data a1(6)/9.90562948053193d1/
      data a1(7)/1.0d0/
      data b1(1)/3.25218725353467d4/
      data b1(2)/7.01022511904373d4/
      data b1(3)/5.50859144223638d4/
      data b1(4)/1.95942074576400d4/
      data b1(5)/3.20803912586318d3/
      data b1(6)/2.20853967067789d2/
      data b1(7)/5.05580641737527d0/
      data b1(8)/1.99507945223266d-2/
      data a2(1)/2.80452693148553d-13/
      data a2(2)/8.60096863656367d-11/
      data a2(3)/1.62974620742993d-8/
      data a2(4)/1.63598843752050d-6/
      data a2(5)/9.12915407846722d-5/
      data a2(6)/2.62988766922117d-3/
      data a2(7)/3.85682997219346d-2/
      data a2(8)/2.78383256609605d-1/
      data a2(9)/9.02250179334496d-1/
      data a2(10)/1.0d0/
      data b2(1)/7.01131732871184d-13/
      data b2(2)/2.10699282897576d-10/
      data b2(3)/3.94452010378723d-8/
      data b2(4)/3.84703231868724d-6/
      data b2(5)/2.04569943213216d-4/
      data b2(6)/5.31999109566385d-3/
      data b2(7)/6.39899717779153d-2/
      data b2(8)/3.14236143831882d-1/
      data b2(9)/4.70252591891375d-1/
      data b2(10)/-2.15540156936373d-2/
      data b2(11)/2.34829436438087d-3/

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
         fm3p = xx*rn/den
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
         fm3p = (x**(an+1d0))*rn/den
      endif
      return
      end
      
