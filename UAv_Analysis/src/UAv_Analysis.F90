
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

  CCTK_REAL aux, S, rho

  CCTK_REAL x1, y1, z1

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


  dE_gf_volume   = 0
  dIxx_gf_volume = 0
  dIxy_gf_volume = 0
  dIxz_gf_volume = 0
  dIyy_gf_volume = 0
  dIyz_gf_volume = 0
  dIzz_gf_volume = 0
  density_rho    = 0

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

    x1      = x(i,j,k)
    y1      = y(i,j,k)
    z1      = z(i,j,k)

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

    ! Eulerian energy density
    rho = Tab(4,4)
    do m = 1, 3
       rho = rho - 2 * beta(m) * Tab(m,4)
       do n = 1, 3
          rho = rho + beta(m) * beta(n) * Tab(m,n)
       end do
    end do
    rho = rho / ( alph * alph )

    density_rho(i,j,k) = rho


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

       ! dE_gf_volume = (alpha h^ij T_ij + T_tt / alpha - beta^i beta^j T_ij / alpha) sqrt(detgd)

       dE_gf_volume(i,j,k)   = (alph * S + aux) * sqrt(detgd)

       ! dI_ij = rho * x^i x^j * alpha * sqrt(detgd)
       dIxx_gf_volume(i,j,k) = alph * rho * x1 * x1 * sqrt(detgd)
       dIxy_gf_volume(i,j,k) = alph * rho * x1 * y1 * sqrt(detgd)
       dIxz_gf_volume(i,j,k) = alph * rho * x1 * z1 * sqrt(detgd)
       dIyy_gf_volume(i,j,k) = alph * rho * y1 * y1 * sqrt(detgd)
       dIyz_gf_volume(i,j,k) = alph * rho * y1 * z1 * sqrt(detgd)
       dIzz_gf_volume(i,j,k) = alph * rho * z1 * z1 * sqrt(detgd)

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
  CCTK_REAL Ixx_int, Ixy_int, Ixz_int, Iyy_int, Iyz_int, Izz_int

  call CCTK_ReductionHandle(reduction_handle, 'sum')
  if (reduction_handle < 0) then
     call CCTK_WARN(0, 'Could not obtain a handle for sum reduction')
  end if

  ! energy

  ! get index to the integration array
  call CCTK_VarIndex(varid, 'UAv_Analysis::dE_gf_volume')
  if (varid < 0) then
     call CCTK_WARN(0, 'Could not get index to grid array dE_gf_volume')
  end if

  ! do a sum over all processors
  call CCTK_Reduce(ierr, cctkGH, -1, reduction_handle, 1, CCTK_VARIABLE_REAL, &
       E_int, 1, varid)
  if (ierr < 0) then
     call CCTK_WARN(0, 'Error while reducing dE_gf_volume')
  end if

  ! TODO: is there a way of doing all components at once?

  ! Ixx
  call CCTK_VarIndex(varid, 'UAv_Analysis::dIxx_gf_volume')
  if (varid < 0) then
     call CCTK_WARN(0, 'Could not get index to grid array dIxx_gf_volume')
  end if
  call CCTK_Reduce(ierr, cctkGH, -1, reduction_handle, 1, CCTK_VARIABLE_REAL, &
       Ixx_int, 1, varid)
  if (ierr < 0) then
     call CCTK_WARN(0, 'Error while reducing dIxx_gf_volume')
  end if

  ! Ixy
  call CCTK_VarIndex(varid, 'UAv_Analysis::dIxy_gf_volume')
  if (varid < 0) then
     call CCTK_WARN(0, 'Could not get index to grid array dIxy_gf_volume')
  end if
  call CCTK_Reduce(ierr, cctkGH, -1, reduction_handle, 1, CCTK_VARIABLE_REAL, &
       Ixy_int, 1, varid)
  if (ierr < 0) then
     call CCTK_WARN(0, 'Error while reducing dIxy_gf_volume')
  end if

  ! Ixz
  call CCTK_VarIndex(varid, 'UAv_Analysis::dIxz_gf_volume')
  if (varid < 0) then
     call CCTK_WARN(0, 'Could not get index to grid array dIxz_gf_volume')
  end if
  call CCTK_Reduce(ierr, cctkGH, -1, reduction_handle, 1, CCTK_VARIABLE_REAL, &
       Ixz_int, 1, varid)
  if (ierr < 0) then
     call CCTK_WARN(0, 'Error while reducing dIxz_gf_volume')
  end if

  ! Iyy
  call CCTK_VarIndex(varid, 'UAv_Analysis::dIyy_gf_volume')
  if (varid < 0) then
     call CCTK_WARN(0, 'Could not get index to grid array dIyy_gf_volume')
  end if
  call CCTK_Reduce(ierr, cctkGH, -1, reduction_handle, 1, CCTK_VARIABLE_REAL, &
       Iyy_int, 1, varid)
  if (ierr < 0) then
     call CCTK_WARN(0, 'Error while reducing dIyy_gf_volume')
  end if

  ! Iyz
  call CCTK_VarIndex(varid, 'UAv_Analysis::dIyz_gf_volume')
  if (varid < 0) then
     call CCTK_WARN(0, 'Could not get index to grid array dIyz_gf_volume')
  end if
  call CCTK_Reduce(ierr, cctkGH, -1, reduction_handle, 1, CCTK_VARIABLE_REAL, &
       Iyz_int, 1, varid)
  if (ierr < 0) then
     call CCTK_WARN(0, 'Error while reducing dIyz_gf_volume')
  end if

  ! Izz
  call CCTK_VarIndex(varid, 'UAv_Analysis::dIzz_gf_volume')
  if (varid < 0) then
     call CCTK_WARN(0, 'Could not get index to grid array dIzz_gf_volume')
  end if
  call CCTK_Reduce(ierr, cctkGH, -1, reduction_handle, 1, CCTK_VARIABLE_REAL, &
       Izz_int, 1, varid)
  if (ierr < 0) then
     call CCTK_WARN(0, 'Error while reducing dIzz_gf_volume')
  end if


  ! the multiplication with the volume element needs to be done here
  total_energy = E_int * cctk_delta_space(1) * cctk_delta_space(2) * cctk_delta_space(3)

  Ixx = Ixx_int * cctk_delta_space(1) * cctk_delta_space(2) * cctk_delta_space(3)
  Ixy = Ixy_int * cctk_delta_space(1) * cctk_delta_space(2) * cctk_delta_space(3)
  Ixz = Ixz_int * cctk_delta_space(1) * cctk_delta_space(2) * cctk_delta_space(3)
  Iyy = Iyy_int * cctk_delta_space(1) * cctk_delta_space(2) * cctk_delta_space(3)
  Iyz = Iyz_int * cctk_delta_space(1) * cctk_delta_space(2) * cctk_delta_space(3)
  Izz = Izz_int * cctk_delta_space(1) * cctk_delta_space(2) * cctk_delta_space(3)

  ! write(*,*) 'total_energy = ', total_energy

  ! note that the integrated values just obtained do *not* take into account
  ! any grid symmetries that may be present. so one needs to remember to
  ! multiply by the appropriate factor (ie, 8 in the case of using just an
  ! octant, etc) a posteriori.

end subroutine UAv_Analysis_IntegrateVol
