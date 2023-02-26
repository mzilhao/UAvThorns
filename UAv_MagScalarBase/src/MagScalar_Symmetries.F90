
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

subroutine MagScalar_symmetries( CCTK_ARGUMENTS )
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT ierr

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "MagScalarBase::Ex" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "MagScalarBase::Ey" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "MagScalarBase::Ez" )

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "MagScalarBase::Ax" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "MagScalarBase::Ay" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "MagScalarBase::Az" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "MagScalarBase::Zeta" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "MagScalarBase::Aphi" )

end subroutine MagScalar_symmetries
