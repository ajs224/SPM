MODULE global_data

	IMPLICIT NONE
	SAVE

	!Constants
	INTEGER, PARAMETER::Np=1 ! Number of particles
	INTEGER, PARAMETER::Nr=20000 !20000 ! Reps.
	REAL, PARAMETER::a=1.0/3.0 ! Power occurring in kernel for stirred reactor, K(x, y)=(xy)^a 
	REAL, PARAMETER::c1=1.0 ! Free parameter
	REAL, PARAMETER::c2=0.1 ! Residence time factor (theta is proportional to mass density) 
	REAL, PARAMETER::eps1=2.0e-2, eps2=0.5 ! Relaxation times
	REAL, PARAMETER::Q_in=6.0 ! Inflow rate Q_in (could be read from the command line or from a file)
	CHARACTER(LEN=18), PARAMETER::output_file="spm_out_one_6" ! Output filename (could be read from the command line)

	! Uniform random number generator
	REAL(KIND=8), EXTERNAL::rnd240

	! Seeds for the random number generator 
	INTEGER(KIND=4)::n0, n1, n2

END MODULE global_data