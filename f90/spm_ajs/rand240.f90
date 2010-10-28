REAL(KIND=8) FUNCTION rnd240(u0, u1, u2)

	IMPLICIT NONE

	INTEGER(KIND=4)::u0, u1, u2
	INTEGER(KIND=4)::c0, c1, c2
	INTEGER(KIND=4)::m0, m1, m2
	
	REAL(KIND=8)::x0, x1, x2

	INTEGER::n

	m0=11973
	m1=2800
	m2=2842

	x0=9.094947017729282379150390625D-13
	x1=1.490116119384765625D-8
	x2=2.44140625D-4

	c0 = m0*u0
	c1 = m0*u1 + m1*u0
	c2 = m0*u2 + m1*u1 + m2*u0

	u0 = c0 - ISHFT(ISHFT(c0, -14), 14)
	n = c1 + ISHFT(c0, -14) 
	u1 = n - ISHFT(ISHFT(n, -14), 14)
	n = c2 + ISHFT(n, -14) 
	u2 = n - ISHFT(ISHFT(n, -12), 12)

	rnd240 = u0*x0 + u1*x1 + u2*x2

END FUNCTION