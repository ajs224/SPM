SUBROUTINE collision(y, x, Phi_present, ymin)

	USE global_data

	IMPLICIT NONE

	REAL, INTENT(IN)::y(Np), Phi_present, ymin
	REAL(KIND=8), INTENT(INOUT)::x
	REAL::f_coll_part_denom, f_coll_part_numer
	INTEGER::i
	LOGICAL::chosen_i

	chosen_i=.FALSE.

	! Denominator of collision partner distribution
	f_coll_part_denom=c1*Phi_present*(x**a)*(ymin**(a-1.0))

	! Choose a particle 
	DO WHILE(chosen_i .EQV. .FALSE.)
		i=1.0+Np*rnd240(n0,n1,n2)
		! Numerator of collision partner distribution
		f_coll_part_numer=c1*Phi_present*(x**a)*(y(i)**(a-1.0))
		IF (rnd240(n0,n1,n2)<f_coll_part_numer/f_coll_part_denom) chosen_i=.TRUE.
	ENDDO

	x=x+y(i)

END SUBROUTINE collision
