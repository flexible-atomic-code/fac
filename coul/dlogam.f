      function dlogam(x)
      double precision x, dlogam
      complex*16 cx, r, clogam, logam

c      r = logam(2D-16)
      cx = dcmplx(x, 0.0D0)
      r = clogam(cx)
      dlogam = dble(r)
      return
      end
