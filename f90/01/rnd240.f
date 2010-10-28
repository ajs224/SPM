      real*8 function rnd240(u0, u1, u2)
      integer*4 u0, u1, u2
      integer*4 c0, c1, c2
      integer*4 m0, m1, m2
      real*8 x0, x1, x2

      data m0/11973/, m1/2800/, m2/2842/

      data x0/9.094947017729282379150390625D-13/,
     &x1/1.490116119384765625D-8/,
     &x2/2.44140625D-4/

      c0 = m0*u0
      c1 = m0*u1 + m1*u0
      c2 = m0*u2 + m1*u1 + m2*u0

      u0 = c0 - ishft(ishft(c0, -14), 14)
      n = c1 + ishft(c0, -14) 
      u1 = n - ishft(ishft(n, -14), 14)
      n = c2 + ishft(n, -14) 
      u2 = n - ishft(ishft(n, -12), 12)

      rnd240 = u0*x0 + u1*x1 + u2*x2

      end function