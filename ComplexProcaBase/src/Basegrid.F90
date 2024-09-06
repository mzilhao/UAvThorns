! Basegrid.F90 : Register symmetries of the grid functions
!
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"


subroutine proca_symmetries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT ierr

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "ProcaBase::E1x" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "ProcaBase::E1y" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "ProcaBase::E1z" )

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "ProcaBase::A1x" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "ProcaBase::A1y" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "ProcaBase::A1z" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ProcaBase::Zeta1" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ProcaBase::Aphi1" )

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "ProcaBase::E2x" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "ProcaBase::E2y" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "ProcaBase::E2z" )

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "ProcaBase::A2x" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "ProcaBase::A2y" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "ProcaBase::A2z" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ProcaBase::Zeta2" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ProcaBase::Aphi2" )

end subroutine proca_symmetries
