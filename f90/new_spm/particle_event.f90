SUBROUTINE particle_event(y, y1, Phi_present, write_to_file)

	USE global_data
	IMPLICIT NONE

	REAL, INTENT(INOUT)::y(np), y1(np), Phi_present
	LOGICAL, INTENT(IN)::write_to_file
	REAL::Phi_future
	REAL(KIND=8)::x
	INTEGER:: i, j, nc
	REAL::ymin
	REAL::coll_rate, out_rate
	REAL::f_coll_part_numer, f_coll_part_denom, tau, t
	LOGICAL::collision_occurred
	
	collision_occurred=.TRUE.

	! Initialise time and iteration number
	t=0.0
	nc=0   

	! Initialise total future mass density of field particles
	Phi_future=0.0

	! Generate a test particle x
	x=rnd240(n0,n1,n2) 

	DO WHILE (collision_occurred)
		!IF(.NOT. collision_occurred) EXIT
		! Compute the rates and the waiting time tau
		CALL compute_rates(y, Phi_present, x, coll_rate, out_rate, ymin, tau)

		! Wait time tau
		t=t+tau

		! recalculate future field particles mass	
		Phi_future=Phi_future+Q_in*tau

		IF (rnd240(n0,n1,n2)<eps2*Q_in*tau/Phi_present) THEN
			! Particle replacement
			CALL replacement(y1, x)

		ENDIF
	
		collision_occurred=.FALSE.	
		
		IF (rnd240(n0,n1,n2)<coll_rate/(coll_rate+out_rate)) THEN
			! A collision occurs
			CALL collision(y, x, Phi_present, ymin)
			nc=nc+1
			collision_occurred=.TRUE.
		ENDIF
	ENDDO


	IF (write_to_file.EQV..TRUE.) WRITE(UNIT=9, FMT=*) x, Phi_present, nc, t 	

	PRINT *, 1.0*nc, Phi_present

	! Underrelax the field particles mass density
	Phi_present=eps1*Phi_future+(1.0-eps1)*Phi_present

	!
	DO i=1,Np
		y(i)=y1(i)
	ENDDO

END SUBROUTINE particle_event
