
#include "cctk.h"
#include "cctk_Arguments.h"

subroutine MagScalar_zero_rhs( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS

  rhs_Ex    = 0
  rhs_Ey    = 0
  rhs_Ez    = 0

  rhs_Ax    = 0
  rhs_Ay    = 0
  rhs_Az    = 0

  rhs_Aphi  = 0

  rhs_Zeta  = 0

end subroutine MagScalar_zero_rhs
