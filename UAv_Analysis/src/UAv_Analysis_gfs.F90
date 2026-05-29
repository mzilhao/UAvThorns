
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
  ! Auxiliaries to avoid repeating computations
  CCTK_REAL sqrt_detgd, vol3, vol4, rho_vol4

  CCTK_REAL S, rho
  CCTK_REAL mom(3)
  
  ! names x0, y0, z0 used as members of the thorn
  ! They are set in dedicated functions, to be called before this routine
  CCTK_REAL x1, y1, z1

  CCTK_INT  i, j, k, m, n

  CCTK_INT type_bits, state_outside

  logical docalc
  
  ! Volume element to be used with multipatch for integration variables
  CCTK_REAL, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: volume_form
  pointer (volume_form_ptr, volume_form)
  CCTK_REAL dV
  ! Cartesian volume element used without multipatch
  CCTK_REAL dV_cart

  
  type_bits     = -1
  state_outside = -1
  
  if (do_analysis_every .le. 0) then
    return
  end if
  
  if (MOD(cctk_iteration, do_analysis_every) .ne. 0 ) then
    return
  endif

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


  if (use_volume_form > 0) then
    call CCTK_VarDataPtr(volume_form_ptr, cctkGH, 0, "Coordinates::volume_form")
  end if
  ! If not using multipatch, we multiply by the coarse Cartesian volume element
  ! (mesh refinement is tackled by sum reduction; lower case cctk_delta_space is the base level spacing)
  dV_cart = cctk_delta_space(1) * cctk_delta_space(2) * cctk_delta_space(3)


!   write(*,*) 'Checking origin coordinates for the analysis in UAv_Analysis'
!   write(*,*) 'x0 = ', x0
!   write(*,*) 'y0 = ', y0
!   write(*,*) 'z0 = ', z0

  dE_gf_volume    = 0
  dJx_gf_volume   = 0
  dJy_gf_volume   = 0
  dJz_gf_volume   = 0
if (early_CoM == 0) then ! They were already computed if early_CoM == 1
  drho_gf_volume  = 0
  dCoMx_gf_volume = 0
  dCoMy_gf_volume = 0
  dCoMz_gf_volume = 0
end if
  dpx_gf_volume   = 0
  dpy_gf_volume   = 0
  dpz_gf_volume   = 0
  dIxx_gf_volume  = 0
  dIxy_gf_volume  = 0
  dIxz_gf_volume  = 0
  dIyy_gf_volume  = 0
  dIyz_gf_volume  = 0
  dIzz_gf_volume  = 0
  if (compute_density_rho == 1 .and. early_CoM == 0) then ! It was already computed if early_CoM == 1
    density_rho    = 0
  end if
  if (compute_density_p == 1) then
    density_px     = 0
    density_py     = 0
    density_pz     = 0
  end if

  ! Note that these loops will also exclude at least one layer of points on the physical boundary,
  ! (see https://lists.einsteintoolkit.org/pipermail/users/2024-September/009465.html )
  ! but as of now we choose to live with that
  do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
  do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
  do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)

    ! checking if outside the horizon, if asking for it to be excised
    docalc = .true.
    if (excise_horizon /= 0) then
      if (.not. SpaceMask_CheckStateBitsF90(space_mask, i, j, k, type_bits, state_outside)) then
        docalc = .false.
      end if
    end if

    ! if inside the horizon, no need to compute the rest (continue with the next
    ! iteration of the do loop)
    if (.not. docalc) cycle

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

    x1      = x(i,j,k) - x0
    y1      = y(i,j,k) - y0
    z1      = z(i,j,k) - z0

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

    sqrt_detgd = sqrt(detgd)
    !--------------------------------------------

    ! With multipatch we need to multiply by the volume element stored in the corresponding variable
    ! Else, just use the standard Cartesian one
    if (use_volume_form > 0) then
      dV = volume_form(i,j,k)
    else
      dV = dV_cart
    end if

    vol3 = sqrt_detgd * dV
    vol4 = alph * vol3

    ! Eulerian energy density

    ! If early_CoM, rho was already computed, so we can just get the value.
    ! We could do something more elaborate from the integrand GF, but we would
    ! need to be careful with divisions by sym_factor_drho and vol4. They
    ! probably don't vanish in general, but for legibility, we keep it simple.
    
    if (early_CoM == 1 .and. compute_density_rho == 1) then 
      rho = density_rho(i,j,k)
    
    else ! We need to compute it in general
      rho = Tab(4,4)
      do m = 1, 3
        rho = rho - 2 * beta(m) * Tab(m,4)
        do n = 1, 3
          rho = rho + beta(m) * beta(n) * Tab(m,n)
        end do
      end do
      rho = rho / ( alph * alph )
      
      if (compute_density_rho == 1) then
        density_rho(i,j,k) = rho
      end if
    end if ! if early_CoM
    rho_vol4 = rho * vol4
    
    ! momentum density
    do n = 1, 3
      mom(n) = Tab(4,n)
      do m = 1, 3
        mom(n) = mom(n) - beta(m) * Tab(m,n)
      end do
      mom(n) = - mom(n) / alph
    end do
   
    if (compute_density_p == 1) then
      density_px(i,j,k) = mom(1)
      density_py(i,j,k) = mom(2)
      density_pz(i,j,k) = mom(3)
    end if
   

    S = 0
    do m = 1, 3
      do n = 1, 3
        S = S + gu(m,n) * Tab(m,n)
      end do
    end do


    ! Symmetry: we multiply by the factors here, so that other thorns can use these GFs safely

    ! dE = (alpha h^ij T_ij + T_tt / alpha - beta^i beta^j T_ij / alpha) sqrt(detgd)
    !              = (alpha * (rho + S) - 2 p_i beta^i) sqrt(detgd)

    dE_gf_volume(i,j,k)   = (alph * (rho + S) - 2 * sum(beta * mom)) * vol3 * sym_factor_dE

    ! dJz = (-y p_x + x p_y) sqrt(detgd)        + permutations
    dJz_gf_volume(i,j,k)  = (-y1 * mom(1) + x1 * mom(2)) * vol3 * sym_factor_dJz
    dJx_gf_volume(i,j,k)  = (-z1 * mom(2) + y1 * mom(3)) * vol3 * sym_factor_dJx
    dJy_gf_volume(i,j,k)  = (-x1 * mom(3) + z1 * mom(1)) * vol3 * sym_factor_dJy
    
  if (early_CoM == 0) then
    ! drho = rho * alpha * sqrt(detgd)
    drho_gf_volume(i,j,k) = rho_vol4 * sym_factor_drho

    ! dCoM^i = rho * x^i * alpha * sqrt(detgd)
    ! Division by integral of density in IntegrateVol
    ! We don't use x1 here (and add x0 back in IntegrateVol), so that this GF can be used in other thorns directly
    dCoMx_gf_volume(i,j,k) = x(i,j,k) * rho_vol4 * sym_factor_dCoMx
    dCoMy_gf_volume(i,j,k) = y(i,j,k) * rho_vol4 * sym_factor_dCoMy
    dCoMz_gf_volume(i,j,k) = z(i,j,k) * rho_vol4 * sym_factor_dCoMz
  end if

    ! dp^i = p^i * alpha * sqrt(detgd)
    dpx_gf_volume(i,j,k) = mom(1) * vol4 * sym_factor_dpx
    dpy_gf_volume(i,j,k) = mom(2) * vol4 * sym_factor_dpy
    dpz_gf_volume(i,j,k) = mom(3) * vol4 * sym_factor_dpz

    ! dI_ij = rho * x^i x^j * alpha * sqrt(detgd)
    dIxx_gf_volume(i,j,k) = x1 * x1 * rho_vol4 * sym_factor_dIxx
    dIxy_gf_volume(i,j,k) = x1 * y1 * rho_vol4 * sym_factor_dIxy
    dIxz_gf_volume(i,j,k) = x1 * z1 * rho_vol4 * sym_factor_dIxz
    dIyy_gf_volume(i,j,k) = y1 * y1 * rho_vol4 * sym_factor_dIyy
    dIyz_gf_volume(i,j,k) = y1 * z1 * rho_vol4 * sym_factor_dIyz
    dIzz_gf_volume(i,j,k) = z1 * z1 * rho_vol4 * sym_factor_dIzz

  end do
  end do
  end do
end subroutine UAv_Analysis_gfs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Special function when we need to compute the center of mass earlier than tracking.
! We only compute the GFs needed for the center of mass.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine UAv_Analysis_early_CoM_gfs( CCTK_ARGUMENTS )
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS



  CCTK_REAL alph, beta(3), Tab(4,4)
  CCTK_REAL gd(3,3), detgd
  CCTK_REAL rho_alph, rho_vol4
  
  CCTK_INT  i, j, k, m, n

  CCTK_INT type_bits, state_outside

  logical docalc
  
  ! Volume element to be used with multipatch for integration variables
  CCTK_REAL, dimension(cctk_lsh(1),cctk_lsh(2),cctk_lsh(3)) :: volume_form
  pointer (volume_form_ptr, volume_form)
  CCTK_REAL dV
  ! Cartesian volume element used without multipatch
  CCTK_REAL dV_cart

  if (early_CoM == 0) then
    return
  end if
  
  if (do_analysis_every .le. 0) then
    return
  end if

  if (MOD(cctk_iteration, do_analysis_every) .ne. 0 ) then
    return
  endif

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


  if (use_volume_form > 0) then
    call CCTK_VarDataPtr(volume_form_ptr, cctkGH, 0, "Coordinates::volume_form")
  end if
  ! If not using multipatch, we multiply by the coarse Cartesian volume element
  ! (mesh refinement is tackled by sum reduction; lower case cctk_delta_space is the base level spacing)
  dV_cart = cctk_delta_space(1) * cctk_delta_space(2) * cctk_delta_space(3)


  drho_gf_volume  = 0
  dCoMx_gf_volume = 0
  dCoMy_gf_volume = 0
  dCoMz_gf_volume = 0
  if (compute_density_rho == 1) then
    density_rho    = 0
  end if

  ! Note that these loops will also exclude at least one layer of points on the physical boundary,
  ! (see https://lists.einsteintoolkit.org/pipermail/users/2024-September/009465.html )
  ! but as of now we choose to live with that
  do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
  do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
  do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)

    ! checking if outside the horizon, if asking for it to be excised
    docalc = .true.
    if (excise_horizon /= 0) then
      if (.not. SpaceMask_CheckStateBitsF90(space_mask, i, j, k, type_bits, state_outside)) then
        docalc = .false.
      end if
    end if

    ! if inside the horizon, no need to compute the rest (continue with the next
    ! iteration of the do loop)
    if (.not. docalc) cycle

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


    !------------ Metric determinant ------------
    detgd =       gd(1,1) * gd(2,2) * gd(3,3)                                &
            + 2 * gd(1,2) * gd(1,3) * gd(2,3)                                &
            -     gd(1,1) * gd(2,3) ** 2                                     &
            -     gd(2,2) * gd(1,3) ** 2                                     &
            -     gd(3,3) * gd(1,2) ** 2
    !--------------------------------------------

    ! With multipatch we need to multiply by the volume element stored in the corresponding variable
    ! Else, just use the standard Cartesian one
    if (use_volume_form > 0) then
      dV = volume_form(i,j,k)
    else
      dV = dV_cart
    end if


    ! Eulerian energy density
    ! We compute rho*alpha here
    rho_alph = Tab(4,4)
    do m = 1, 3
      rho_alph = rho_alph - 2 * beta(m) * Tab(m,4)
      do n = 1, 3
        rho_alph = rho_alph + beta(m) * beta(n) * Tab(m,n)
      end do
    end do
    rho_alph = rho_alph / alph
    rho_vol4 = rho_alph * sqrt(detgd) * dV

    if (compute_density_rho == 1) then
      density_rho(i,j,k) = rho_alph / alph
    end if
    
    ! Symmetry: we multiply by the factors here, so that other thorns can use these GFs safely

    ! drho = rho * alpha * sqrt(detgd)
    drho_gf_volume(i,j,k) = rho_vol4 * sym_factor_drho

    ! dCoM^i = rho * x^i * alpha * sqrt(detgd)
    ! Division by integral of density in UAv_Analysis_early_CoM_reduce
    ! We don't use x1 here (and add x0 back in UAv_Analysis_early_CoM_reduce) for consistency
    dCoMx_gf_volume(i,j,k) = x(i,j,k) * rho_vol4 * sym_factor_dCoMx
    dCoMy_gf_volume(i,j,k) = y(i,j,k) * rho_vol4 * sym_factor_dCoMy
    dCoMz_gf_volume(i,j,k) = z(i,j,k) * rho_vol4 * sym_factor_dCoMz

  end do
  end do
  end do

end subroutine UAv_Analysis_early_CoM_gfs
