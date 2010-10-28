SUBROUTINE compute_rates(y, Phi_present, x, coll_rate, out_rate, ymin, tau)

	USE global_data

	IMPLICIT NONE
	
	REAL, INTENT(IN)::y(Np), Phi_present
	REAL(KIND=8), INTENT(IN)::x
	REAL, INTENT(OUT)::coll_rate, out_rate, ymin, tau
	INTEGER::i

	! Compute collision rate
	coll_rate=0.0
	ymin=1.0e8
	DO i=1, Np
		ymin=MIN(ymin, y(i))
		coll_rate=coll_rate+c1*Phi_present*(x**a)*(y(i)**(a-1.0))/Np
	ENDDO

	! Compute outflow rate
	out_rate=c2*x

	! Compute exponentially distributed waiting time
	tau=-LOG(rnd240(n0,n1,n2))/(coll_rate+out_rate)

END SUBROUTINE compute_rates