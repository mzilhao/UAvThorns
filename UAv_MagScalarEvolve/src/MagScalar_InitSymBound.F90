#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

subroutine MagScalar_InitSymBound( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT ierr

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "UAv_MagScalarEvolve::rhs_Ex" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "UAv_MagScalarEvolve::rhs_Ey" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "UAv_MagScalarEvolve::rhs_Ez" )

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "UAv_MagScalarEvolve::rhs_Ax" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "UAv_MagScalarEvolve::rhs_Ay" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "UAv_MagScalarEvolve::rhs_Az" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "UAv_MagScalarEvolve::rhs_Zeta" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "UAv_MagScalarEvolve::rhs_Aphi" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "UAv_MagScalarEvolve::rhs_phi1" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "UAv_MagScalarEvolve::rhs_phi2" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "UAv_MagScalarEvolve::rhs_Kphi1" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "UAv_MagScalarEvolve::rhs_Kphi2" )

end subroutine MagScalar_InitSymBound
