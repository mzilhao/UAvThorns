#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

subroutine Proca_InitSymBound( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT ierr

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "ProcaEvolve::rhs_E1x" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "ProcaEvolve::rhs_E1y" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "ProcaEvolve::rhs_E1z" )

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "ProcaEvolve::rhs_A1x" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "ProcaEvolve::rhs_A1y" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "ProcaEvolve::rhs_A1z" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ProcaEvolve::rhs_Zeta1" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ProcaEvolve::rhs_Aphi1" )

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "ProcaEvolve::rhs_E2x" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "ProcaEvolve::rhs_E2y" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "ProcaEvolve::rhs_E2z" )

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "ProcaEvolve::rhs_A2x" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "ProcaEvolve::rhs_A2y" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "ProcaEvolve::rhs_A2z" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ProcaEvolve::rhs_Zeta2" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ProcaEvolve::rhs_Aphi2" )

end subroutine Proca_InitSymBound
