      function dlogam(x)
      double precision x, dlogam
      complex*16 cx, clgam 

      cx = dcmplx(x, 0.0D0)
      cx = clgam(cx)
      dlogam = dble(cx)
      return
      end
      
