#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

subroutine MagScalar_InitSymBound( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT ierr

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "MagScalarEvolve::rhs_Ex" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "MagScalarEvolve::rhs_Ey" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "MagScalarEvolve::rhs_Ez" )

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "MagScalarEvolve::rhs_Ax" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "MagScalarEvolve::rhs_Ay" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "MagScalarEvolve::rhs_Az" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "MagScalarEvolve::rhs_Zeta" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "MagScalarEvolve::rhs_Aphi" )

end subroutine MagScalar_InitSymBound
