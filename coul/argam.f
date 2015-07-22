      function argam(x, y)
      real*8 x, y, argam, PI, TPI
      complex*16 clogam, logam, r          
      PARAMETER (PI=3.1415926535897932D0,TPI=PI+PI)  

      r = logam(1D-16)
      r = dcmplx(x, y)
      r = clogam(r)
      
      argam = dimag(r) 

      if (argam .gt. 0) then
         argam = dmod(argam, TPI)
      else
         argam = -dmod(-argam, TPI)
      endif

      return
      end
      
