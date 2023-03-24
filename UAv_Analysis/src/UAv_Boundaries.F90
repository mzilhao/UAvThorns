#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine UAv_Analysis_Boundaries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT ierr
  CCTK_INT, parameter :: one = 1
  CCTK_INT, parameter :: bndsize = 3

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, bndsize, -one, &
       "UAv_Analysis::densities", "flat")
  if (ierr < 0)                                                           &
       call CCTK_ERROR("Failed to register BC for UAv_Analysis::densities!")

end subroutine UAv_Analysis_Boundaries
