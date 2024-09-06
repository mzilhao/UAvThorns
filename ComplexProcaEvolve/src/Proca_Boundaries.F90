#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine Proca_Boundaries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT ierr

  CCTK_INT, parameter :: one = 1

  ! The outgoing (radiative) boundary conditions are being handled from calls to
  ! the NewRad infraestructure. Here we register all BCs as 'none', which
  ! enforces all the symmetry BCs.

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,      &
       "ProcaBase::E1i", "none")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for ProcaBase::E1i!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,      &
       "ProcaBase::A1i", "none")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for ProcaBase::A1i!")

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, one, -one,        &
       "ProcaBase::Aphi1", "none")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for ProcaBase::Aphi1!")

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, one, -one,        &
       "ProcaBase::Zeta1", "none")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for ProcaBase::Zeta1!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,      &
       "ProcaBase::E2i", "none")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for ProcaBase::E2i!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,      &
       "ProcaBase::A2i", "none")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for ProcaBase::A2i!")

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, one, -one,        &
       "ProcaBase::Aphi2", "none")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for ProcaBase::Aphi2!")

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, one, -one,        &
       "ProcaBase::Zeta2", "none")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for ProcaBase::Zeta2!")

end subroutine Proca_Boundaries
