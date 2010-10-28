SUBROUTINE replacement(y1, x)

	USE global_data

	IMPLICIT NONE

	REAL, INTENT(INOUT)::y1(Np)
	REAL(KIND=8), INTENT(IN)::x
	INTEGER::i

	i=1.0+Np*rnd240(n0,n1,n2)
	y1(i)=x

END SUBROUTINE replacement