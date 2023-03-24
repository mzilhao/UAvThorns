
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#include "SpaceMask.h"

subroutine UAv_Analysis_gfs( CCTK_ARGUMENTS )
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL alph, beta(3), Tab(4,4)
  CCTK_REAL gd(3,3), gu(3,3), detgd

  CCTK_REAL aux, S

  CCTK_INT  i, j, k, m, n

  CCTK_INT type_bits, state_outside

  type_bits     = -1
  state_outside = -1


  if (excise_horizon /= 0) then

     call SpaceMask_GetTypeBits(type_bits, "mask")
     call SpaceMask_GetStateBits(state_outside, "mask", "outside")

     if (type_bits < 0) then
        call CCTK_WARN(0, "Thorn AHFinderDirect not activated, but excise_horizon requires it.")
     end if

     if (state_outside < 0) then
        call CCTK_WARN(0, "Error in obtaining AHFinderDirect GetStateBits")
     end if

  end if


  dE = 0

  do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
  do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
  do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)

    !--------------Get local variables ----------
    gd(1,1) = gxx(i,j,k)
    gd(1,2) = gxy(i,j,k)
    gd(1,3) = gxz(i,j,k)
    gd(2,2) = gyy(i,j,k)
    gd(2,3) = gyz(i,j,k)
    gd(3,3) = gzz(i,j,k)
    gd(2,1) = gd(1,2)
    gd(3,1) = gd(1,3)
    gd(3,2) = gd(2,3)

    alph    = alp(i,j,k)

    beta(1) = betax(i,j,k)
    beta(2) = betay(i,j,k)
    beta(3) = betaz(i,j,k)

    ! stress-energy tensor variables
    Tab = 0
    if (stress_energy_state /= 0) then
       Tab(4,4) = eTtt(i,j,k)
       Tab(4,1) = eTtx(i,j,k)
       Tab(4,2) = eTty(i,j,k)
       Tab(4,3) = eTtz(i,j,k)
       Tab(1,1) = eTxx(i,j,k)
       Tab(1,2) = eTxy(i,j,k)
       Tab(1,3) = eTxz(i,j,k)
       Tab(2,2) = eTyy(i,j,k)
       Tab(2,3) = eTyz(i,j,k)
       Tab(3,3) = eTzz(i,j,k)
       Tab(1,4) = Tab(4,1)
       Tab(2,4) = Tab(4,2)
       Tab(3,4) = Tab(4,3)
       Tab(2,1) = Tab(1,2)
       Tab(3,1) = Tab(1,3)
       Tab(3,2) = Tab(2,3)
    end if
    !--------------------------------------------


    !-------------- Invert metric ---------------
    detgd =       gd(1,1) * gd(2,2) * gd(3,3)                                &
            + 2 * gd(1,2) * gd(1,3) * gd(2,3)                                &
            -     gd(1,1) * gd(2,3) ** 2                                     &
            -     gd(2,2) * gd(1,3) ** 2                                     &
            -     gd(3,3) * gd(1,2) ** 2
    gu(1,1) = (gd(2,2) * gd(3,3) - gd(2,3) ** 2     ) / detgd
    gu(2,2) = (gd(1,1) * gd(3,3) - gd(1,3) ** 2     ) / detgd
    gu(3,3) = (gd(1,1) * gd(2,2) - gd(1,2) ** 2     ) / detgd
    gu(1,2) = (gd(1,3) * gd(2,3) - gd(1,2) * gd(3,3)) / detgd
    gu(1,3) = (gd(1,2) * gd(2,3) - gd(1,3) * gd(2,2)) / detgd
    gu(2,3) = (gd(1,3) * gd(1,2) - gd(2,3) * gd(1,1)) / detgd
    gu(2,1) = gu(1,2)
    gu(3,1) = gu(1,3)
    gu(3,2) = gu(2,3)
    !--------------------------------------------

    aux = 0
    do m = 1, 3
      do n = 1, 3
        aux = aux + beta(m) * beta(n) * Tab(m,n)
      end do
    end do
    aux = (Tab(4,4) - aux) / alph

    S = 0
    do m = 1, 3
       do n = 1, 3
          S = S + gu(m,n) * Tab(m,n)
       end do
    end do


    ! checking if outside the horizon, if asking for it to be excised
    if (SpaceMask_CheckStateBitsF90(space_mask, i, j, k, type_bits, state_outside) .or. &
         excise_horizon == 0) then

       ! dE = (alpha h^ij T_ij + T_tt / alpha - beta^i beta^j T_ij / alpha) sqrt(detgd)

       dE(i,j,k) = (alph * S + aux) * sqrt(detgd)

    end if

  end do
  end do
  end do
end subroutine UAv_Analysis_gfs


subroutine UAv_Analysis_IntegrateVol( CCTK_ARGUMENTS )
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT ierr
  CCTK_INT reduction_handle, varid

  CCTK_REAL E_int

  call CCTK_ReductionHandle(reduction_handle, 'sum')
  if( reduction_handle < 0 ) then
     call CCTK_WARN(0, 'Could not obtain a handle for sum reduction')
  end if


  ! get index to the integration array
  call CCTK_VarIndex(varid, 'KillingQuantities::dE')
  if( varid < 0 ) then
     call CCTK_WARN(0, 'Could not get index to grid array dE')
  end if

  ! do a sum over all processors
  call CCTK_Reduce(ierr, cctkGH, -1, reduction_handle, 1, CCTK_VARIABLE_REAL, &
       E_int, 1, varid)

  ! the multiplication with the volume element needs to be done here
  total_energy = E_int * cctk_delta_space(1) * cctk_delta_space(2) * cctk_delta_space(3)

  ! write(*,*) 'total_energy = ', total_energy

  ! note that the integrated value just obtained does *not* take into account
  ! any grid symmetries that may be present. so one needs to remember to
  ! multiply by the appropriate factor (ie, 8 in the case of using just an
  ! octant, etc) a posteriori.

end subroutine UAv_Analysis_IntegrateVol
