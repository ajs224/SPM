! SPM implementation for the simulation of a zero dimensional reactor.
! This implements the version of the SPM found in Sasha's SPM Paper.
! Alastair J. Smith (ajs224@cam.ac.uk)
! 20/06/04

PROGRAM spm_main

	USE global_data

	IMPLICIT NONE

	REAL::y(Np), y1(Np) ! Field particle y
	REAL(4)::ta(2), xt1, xt2
	REAL(4), EXTERNAL :: dtime
	INTEGER::k, i
	REAL::Phi_present ! Present total mass density of field particles

	! Get seeds for the random number generator
	OPEN(UNIT=5, FILE="rand_seeds", STATUS="old")
	READ(UNIT=5,FMT=*) n0, n1, n2
	CLOSE(UNIT=5)

	! Open file for data output
	OPEN(UNIT=9, FILE=output_file, STATUS="unknown")

	! Generate an initial state
	! y(1),..., y(Np)
	Phi_present=1.0
	DO i=1,Np
		y(i)=1.0
		y1(i)=1.0
	ENDDO
	
	DO k=1,12*Np
		CALL particle_event(y, y1, Phi_present, .TRUE.) 
	ENDDO
	
	xt1=dtime(ta)

	DO k=1, Nr
		CALL particle_event(y, y1, Phi_present, .TRUE.)     
	ENDDO

	CLOSE(UNIT=9)

	OPEN(UNIT=5, FILE="rand_seeds", STATUS="unknown")
	WRITE(UNIT=5,FMT=*) n0,n1,n2
	CLOSE(UNIT=5)

	! Print how long the algorithm has taken
	xt2=dtime(ta)
	PRINT *, 'Time elapsed: ', xt2

      
END PROGRAM spm_main