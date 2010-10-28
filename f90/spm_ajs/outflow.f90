SUBROUTINE outflow(x, t, M, xmin, xmax, write_to_file)

	USE global_data

	IMPLICIT NONE

	REAL, INTENT(INOUT)::x(Np)
	REAL, INTENT(IN)::xmin, xmax
	LOGICAL, INTENT(IN)::write_to_file
	
	INTEGER::i,j
	LOGICAL::chosen_i,chosen_j

	chosen_i=.FALSE.
	chosen_j=.FALSE. 



END SUBROUTINE outflow