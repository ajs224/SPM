FUNCTION theta(x)

  USE global_data
  
  USE nrtype
  
  IMPLICIT NONE

  REAL(DP),INTENT(IN)::x
  REAL(DP)::theta
 
  theta=one/(c1*x)

END FUNCTION theta
