      function dlogam(x)
      double precision x, dlogam
      complex*16 cx, clogam 

      cx = dcmplx(x, 0.0D0)
      cx = clogam(cx)
      dlogam = dble(cx)
      return
      end
      
