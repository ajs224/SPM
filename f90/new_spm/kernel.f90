FUNCTION K(x,y)

  USE global_data
  
  USE nrtype
  
  IMPLICIT NONE

  REAL(DP), INTENT(IN)::x,y
  REAL(DP), PARAMETER::a=1.0/3.0 
  REAL(DP):: K

  K=(x**a)*(y**a)

END FUNCTION K
