MODULE global_data

  USE nrtype

  IMPLICIT NONE
  SAVE

  !Constants
  INTEGER, PARAMETER::N=100 ! Number of field particles
  INTEGER, PARAMETER::x_max=100
  REAL(DP), PARAMETER::c1=1.0 ! 
  REAL(DP), PARAMETER::Q_in=1.0 ! Inflow rate Q_in (could be read from the command line or from a file)
  CHARACTER(LEN=18), PARAMETER::fp_file="fp_dist_Qin1" ! Output filename (could be read from the command line)
  CHARACTER(LEN=18), PARAMETER::mass_dens_file="mass_dens_Qin1" ! Output filename (could be read from the command line)
  ! Uniform random number generator
  REAL(KIND=8), EXTERNAL::rnd240
  
  ! Seeds for the random number generator 
  INTEGER(KIND=4)::n0, n1, n2

END MODULE global_data
